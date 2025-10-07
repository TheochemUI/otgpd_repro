/*
 * GPR.cpp
 *
 *  Created on: 23 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "AtomicDimer.h"

#include <cmath>
#include <cstddef>
#include <limits>

#include "../managers/io/FileManager.h"
#include "Enums.h"
#include "auxiliary/Distance.h"
#include "dimer/Dimer.h"
#include "ml/GaussianProcessRegression.h"
#include "structures/Structures.h"

#ifdef WITH_HDF5

#define H5_USE_EIGEN
#include <highfive/H5Easy.hpp>

#ifdef H5_USE_EIGEN
#include <Eigen/Eigen>
#endif

#endif

namespace atmd {

bool is_subsystem(const gpr::AtomsConfiguration& atoms_config)
{
    return atoms_config.atoms_froz_active.positions.getSize() > 0 ||
           atoms_config.atoms_froz_inactive.positions.getSize() > 0;
}

AtomicDimer::AtomicDimer()
{
    clear();

    // Assume cubic cell
    cell_dimensions[0] = {15.3455999999999992, 0., 0.};
    cell_dimensions[1] = {0., 21.7020000000000017, 0.};
    cell_dimensions[2] = {0., 0., 100.0000000000000000};

    gpr_model = new gpr::GaussianProcessRegression();

    num_of_gen_potential_calls = 0;
    early_params._threshold = 0.4;
    early_params._use_adaptive_threshold = true;
    early_params._adaptive_A = 1.3;
    early_params._adaptive_floor = 0.25;
    early_params._dist_metric = DistanceMetricType::EMD;
}

AtomicDimer::~AtomicDimer()
{
    if (gpr_model != nullptr) {
        delete gpr_model;
        gpr_model = nullptr;
    }

    clear();
}

void AtomicDimer::clear()
{
    atoms_config.clear();
    cell_dimensions[0] = {};
    cell_dimensions[1] = {};
    cell_dimensions[2] = {};
    orient_init.clear();

    middle_point_init.clear();
    image1.clear();
    all_obs.clear();

    E_zero_level.clear();
    stop_citeria_dimer.clear();
    stop_citeria_gpr.clear();

    dimer_sep = 0.;

    max_iter_init_rot.setZero();
    max_iter_new_pairs.setZero();
    max_iter_relax_rot.setZero();

    ratio_at_limit = 0;
    actdist_fro = 0.;
    num_bigiter_initparam = 0;

    curvature_final = 0.;

    gpr_model = nullptr;
}

bool AtomicDimer::isAtomsConfigurationCorrect()
{
    bool res = true;

    if (atoms_config.atoms_froz_active.type.isEmpty()) {
        res = false;
    }
    if (atoms_config.atoms_mov.positions.getSize() ==
        atoms_config.positions.getSize()) {
        res = true;
    }

    return res;
}

void AtomicDimer::updateLocation(const gpr::Coord orient,
                                 const gpr::Coord& middle_point,
                                 gpr::Coord& new_point)
{
    new_point = middle_point + orient * dimer_sep;
}

bool AtomicDimer::isFinalConvergenceReached(
    const gpr::Observation& middle_point, const double stop_criteria)
{
    // TODO: check if the implementation is fine. In matlab it is done as
    // `max(abs(G_R))`
    // TODO(rg): Take an input parameter and pick either norm, max_force_on_atom
    // etc.
    double max_value = middle_point.G.norm();

    if (max_value < stop_criteria) return true;
    log_man << "Current Max is " << max_value << "\n";

    return false;
}

bool AtomicDimer::isRotationalConvergenceReached(const gpr::Coord& orient,
                                                 gpr::Observation& middle_point)
{
    dimer::Dimer dimer(dimer_sep, is_subsystem(atoms_config));
    gpr::Coord F_rot;
    gpr::Coord G01_loc;
    bool is_convergence_rached = false;

    G01_loc = middle_point.G;
    G01_loc.append(image1.G);

    dimer.rotateForce(G01_loc, orient, F_rot);

    double omega = dimer.estimateRotationalAngle(orient, G01_loc, F_rot);

    if (omega <= stop_citeria_dimer.angle_rotation)
        is_convergence_rached = true;
    else
        is_convergence_rached = false;

    return is_convergence_rached;
}

bool AtomicDimer::rotateDimer(const double stop_criteria, gpr::Coord& orient,
                              gpr::LBFGSInfo& rot_info,
                              gpr::Observation& middle_point,
                              gpr::Observation& approx_obs)
{
    dimer::Dimer dimer(dimer_sep, is_subsystem(atoms_config));
    gpr::Observation observed{
        middle_point};  // observed coordinates, energy and gradient of
                        // the middle point and the image 1
                        // For starters it can just be the approx_obs
    gpr::Observation reduced_obs;
    gpr::Coord orient_old;
    double curv_dummy = 0;  // dummy variable
    bool terminate_loop = false;

    // if necessary, rotate the dimer and re-calculate approximated energy and
    // gradient at image 1
    orient_old = orient;

    dimer.rotate(middle_point.R, orient_old, approx_obs.G, POTENTIAL_GP,
                 stop_criteria, false, all_obs, gpr_model, rot_info, orient,
                 curv_dummy, observed);

    if (!observed.R.isEmpty()) {
        // TODO: refactor, this looks like updateLocation() method
        for (gpr::Index_t n = 0; n < approx_obs.R.getNumCols(); ++n) {
            approx_obs.R(1, n) =
                middle_point.R(0, n) + orient(0, n) * dimer_sep;
        }

        // TODO: optimize!
        // Strictly speaking, this copy is unnecessary. We can modify
        // calculatePotential() method to consider only slice from the
        // input data sets
        reduced_obs.R.append(approx_obs.R, 1);

        // TODO: check, do we really need G and E?
        reduced_obs.G.append(approx_obs.G, 1);
        reduced_obs.E.resize(1, 1);
        reduced_obs.E(0, 0) = approx_obs.E(0, 1);

        gpr_model->calculatePotential(reduced_obs);

        for (gpr::Index_t j = 0; j < approx_obs.R.getNumCols(); ++j) {
            approx_obs.R(1, j) = reduced_obs.R(0, j);
            approx_obs.G(1, j) = reduced_obs.G(0, j);
        }
        approx_obs.E(0, 1) = reduced_obs.E(0, 0);

        terminate_loop =
            checkRotationalConvergence(orient, orient_old, stop_criteria);
    } else {
        terminate_loop = true;
        log_man << "Terminated due to R_obs is empty\n";
    }

    return terminate_loop;
}

bool AtomicDimer::checkRotationalConvergence(const gpr::Coord& orient,
                                             const gpr::Coord& orient_old,
                                             const double stop_criteria)
{
    if (acos(orient.dot(orient_old)) < stop_criteria) {
        return true;
    }
    return false;
}

void AtomicDimer::translateDimer(const gpr::Observation& middle_point,
                                 const gpr::Coord& orient,
                                 const gpr::Observation& approx_obs,
                                 const double curvature,
                                 gpr::LBFGSInfo& trans_info, gpr::Coord& R_new)
{
    dimer::Dimer dimer(dimer_sep, is_subsystem(atoms_config));
    gpr::Coord G01_approx_negative;

    G01_approx_negative.append(approx_obs.G, 0);
    G01_approx_negative *= -1.;

    dimer.translate(middle_point.R, orient, G01_approx_negative, curvature,
                    transition_param, trans_info, R_new);
}

bool AtomicDimer::updateActiveFrozenAtoms(const gpr::Coord& R_new,
                                          gpr::LBFGSInfo& rot_info,
                                          gpr::LBFGSInfo& trans_info)
{
    aux::ProblemSetUp setup;
    bool atoms_were_updated = false;

    if (debug_level == DEBUG_L1 || debug_level == DEBUG_L2 ||
        debug_level == DEBUG_L3) {
        log_man << "Frozen:: "
                << atoms_config.atoms_froz_active.positions.getSize()
                << " active, "
                << atoms_config.atoms_froz_inactive.positions.getSize()
                << " inactive\n"
                << "Moving:: " << atoms_config.atoms_mov.positions.getSize()
                << " Total:: " << atoms_config.positions.getSize() << "\n";
    }

    // if (std::isfinite(actdist_fro)) {
    //     atoms_were_updated = setup.activateFrozenAtoms(R_new, actdist_fro,
    //                                                    atoms_config);
    //     if (atoms_were_updated) {
    //         gpr::Field<double> dummy;
    //         gpr::Field<double> dist;

    //         dummy.resize(1, 1);
    //         dummy[0] = 1.;

    //         log_man << "More frozen atoms activated. Now "
    //                 << atoms_config.atoms_froz_active.positions.getNj() / 3
    //                 << " active and "
    //                 << atoms_config.atoms_froz_inactive.positions.getNj() / 3
    //                 << " inactive frozen atoms.\n";

    //         gpr_model->getSexpAtCovarianceFunction()->setConfInfo(atoms_config);
    //         gpr_model->setHyperparameters(all_obs, atoms_config, true, false,
    //         false);

    //         rot_info.clearOld();
    //         rot_info.cgiter_rot = 0;

    //         trans_info.clear();
    //     }
    // }

    return atoms_were_updated;
}

bool AtomicDimer::limitTranslation(const gpr::Coord& R,
                                   gpr::LBFGSInfo& trans_info,
                                   gpr::Coord& R_new)
{
    aux::Distance distance;
    gpr::Field<double> R_diff;
    gpr::Field<double> dist;
    gpr::Field<double> step_length_atomwise;
    double tmp = 0.;
    double step_coeff = 0.;
    gpr::Field<double> steplength_atomwise_limit;
    bool step_length_is_larger_than_limit = false;

    R_diff = R_new - R;

    // atom-wise step lengths
    step_length_atomwise.resize(R.getNumRows(), R.getNumCols() / 3);
    for (gpr::Index_t i = 0; i < R_new.getNumRows(); ++i) {
        for (gpr::Index_t j = 0; j < R_new.getNumCols(); j += 3) {
            double x = R_new(i, j) - R(i, j);
            double y = R_new(i, j + 1) - R(i, j + 1);
            double z = R_new(i, j + 2) - R(i, j + 2);
            step_length_atomwise[j / 3] += sqrt(x * x + y * y + z * z);
        }
    }

    distance.mindist_interatomic(R, atoms_config, dist);
    steplength_atomwise_limit.resize(dist.getSize());
    tmp = 0.5 * (1. - ratio_at_limit);
    for (gpr::Index_t n = 0; n < steplength_atomwise_limit.getSize(); ++n) {
        steplength_atomwise_limit[n] = tmp * dist[n];
    }

    for (gpr::Index_t n = 0; n < steplength_atomwise_limit.getSize(); ++n) {
        if (step_length_atomwise[n] > 0.99 * steplength_atomwise_limit[n]) {
            step_length_is_larger_than_limit = true;
            break;
        }
    }

    if (step_length_is_larger_than_limit) {
        step_coeff = std::numeric_limits<double>::max();
        for (gpr::Index_t n = 0; n < steplength_atomwise_limit.getSize(); ++n) {
            tmp = steplength_atomwise_limit[n] / step_length_atomwise[n];
            if (tmp < step_coeff) {
                step_coeff = tmp;
            }
        }
        step_coeff *= 0.99;

        R_new = R + (R_new - R) * step_coeff;

        trans_info.clear();
    }

    return step_length_is_larger_than_limit;
}

bool AtomicDimer::isInterAtomicDistanceTooBig(const gpr::Coord& R_new,
                                              const gpr::Coord& R_all,
                                              gpr::Index_t& num_es1)
{
    bool res = false;
    aux::Distance distance;
    gpr::Field<double> dist;
    double __dist{0};
    double thresh{0};
    if (early_params._use_adaptive_threshold) {
        // thresh = std::max(early_params._adaptive_floor,
        //                   early_params._adaptive_A / sqrt(num_atoms));

        // This block implements an adaptive threshold that grows with the
        // number of data points but is capped by a physically motivated,
        // size-dependent upper limit. This prevents the step from becoming too
        // large, even with many data points, and enforces a stricter learning
        // curve at the beginning of the search.

        // Assumption: Get the current number of data points in the history.
        const int num_data_points = static_cast<int>(all_obs.E.getSize());

        // Unpack parameters for clarity.
        const double t_min = 0.15;
        // early_params._adaptive_t_min;  // e.g., 0.2
        const double delta_t_explore = 0.35;  // 0.35
        // early_params._adaptive_delta_t_explore;           // e.g., 0.4
        const int n_half = 50;  // early_params._adaptive_n_half;  // e.g., 20.0
        // used to be 50

        // Calculate the rate constant k from N_half.
        const double k = std::log(2.0) / n_half;
        const double earned_thresh =
            t_min + delta_t_explore * (1.0 - std::exp(-k * num_data_points));

        const int num_atoms =
            static_cast<int>(atoms_config.atoms_mov.type.getSize());
        const double adaptive_A = early_params._adaptive_A;  // e.g., 1.3
        const double adaptive_floor =
            early_params._adaptive_floor;  // e.g., 0.2
        const double physical_limit =
            std::max(adaptive_floor, adaptive_A / std::sqrt(num_atoms));
        // const double physical_limit = 10000;

        // Part 3: The final threshold is the minimum of the two.
        // The step can be as large as the data allows, but never larger than
        // the physical limit.
        thresh = std::min(earned_thresh, physical_limit);

    } else {
        thresh = early_params._threshold;
    }
    // Now get the actual distance
    switch (early_params._dist_metric) {
        case DistanceMetricType::MAX_1D_LOG: {
            distance.dist_max1Dlog(R_new, R_all, atoms_config, dist);
            // NOTE(rg): This restores MATLAB/JCTC settings
            thresh = fabs(log(ratio_at_limit));
            break;
        }
        case DistanceMetricType::EMD: {
            distance.dist_emd(R_new, R_all, atoms_config, dist);
            break;
        }
        case DistanceMetricType::RMSD: {
            distance.dist_rmsd(R_new, R_all, atoms_config, dist);
            break;
        }
        default:
            throw("Shouldn't be here");
    }
    __dist = dist.getMinElt();

    if (__dist > thresh) {
        ++num_es1;
        res = true;
    }

    log_man << "NOTE: got " << __dist << " compared to " << thresh << ".\n";
    return res;
}

gpr::Index_t AtomicDimer::getNumEvaluations()
{
    // size(E_all,1)-N_obs_init
    return (gpr::Index_t)(all_obs.E.getNumRows() - E_all_init.getNumRows());
}

void AtomicDimer::checkIfMaxIterReachedAndPrintWarning(
    const gpr::Index_t iter_counter, const gpr::Index_t iter_limit,
    const std::string& message)
{
    if (iter_counter == iter_limit - 1) {
        log_man << "WARNING: Maximum number of iterations " << iter_counter
                << " was reached during the " << message << ".\n";
    }
}

void AtomicDimer::defineStartPath(const gpr::Coord& R_latest_conv,
                                  const gpr::Coord& orient_latest_conv,
                                  const gpr::Coord& orient_init_gpr,
                                  gpr::Coord& R_previous,
                                  gpr::Coord& orient_previous,
                                  gpr::Observation& middle_point,
                                  gpr::Coord& orient)
{
    if (assume_many_iterations || R_previous.isEmpty()) {
        if (!R_latest_conv.isEmpty()) {
            middle_point.R = R_latest_conv;
            orient = orient_latest_conv;
            log_man << "Started relaxation phase from the latest converged "
                       "dimer.\n";
        } else {
            middle_point.R = middle_point_init.R;
            orient = orient_init_gpr;
            log_man << "Started relaxation phase from the initial location.\n";
        }
    } else {
        middle_point.R = R_previous;
        orient = orient_previous;
        log_man << "Started relaxation phase where the previous one stopped.\n";
        R_previous.clear();
        orient_previous.clear();
    }
}

double AtomicDimer::performRotationAndTranslationWithGPR(
    const double stop_dimer_gp, gpr::Observation& middle_point,
    gpr::Coord& orient, gpr::Coord& R_latest_conv, gpr::Coord& R_previous,
    gpr::Coord& orient_latest_conv, gpr::Coord& orient_previous,
    gpr::Coord& R_new, gpr::Index_t& num_esmax, gpr::Index_t& num_es1,
    const size_t& outer_iter_idx)
{
    //    gpr::LBFGSInfo rot_info;               // Information on rotation
    //    gpr::LBFGSInfo trans_info;           // Information on translation
    gpr::LBFGSInfo rot_info;
    gpr::LBFGSInfo trans_info;
    gpr::Observation approx_obs;
    gpr::Coord temp;
    gpr::Index_t inner_iter;
    double curvature = 0.;

    rot_info.clear();
    rot_info.num_cg_iter = middle_point.R.getNumCols();
    rot_info.num_lbfgs_iter = middle_point.R.getNumCols();

    trans_info.clear();
    trans_info.num_cg_iter = middle_point.R.getNumCols();
    trans_info.num_lbfgs_iter = middle_point.R.getNumCols();

    for (inner_iter = 0; inner_iter < max_iter_new_pairs.inner; ++inner_iter) {
        log_man << " Inner relaxation iteration: " << inner_iter << "\n";

        // Calculate approximated energy and gradient at the middle point
        // and image 1 of the dimer
        approx_obs.R = middle_point.R;
        updateLocation(orient, middle_point.R, temp);
        approx_obs.R.append(temp);
        gpr_model->calculatePotential(approx_obs);

        // Stop the relaxation phase if converged
        // TODO: move evaluation to something similar to
        // isFinal0ConvergenceReached
        if (approx_obs.G.getMaxAbsEltInRow(0) < stop_dimer_gp) {
            R_latest_conv = middle_point.R;
            orient_latest_conv = orient;
            log_man << "Stopped relaxation phase: converged after "
                    << inner_iter + 1 << " inner iterations.\n";
            break;
        }

        //        E_R_gp.appendScalarAsSlice(approx_obs.E(0, 0));
        //        maxF_R_gp.push_back(middle_point.G.getMaxAbsElt());

        // Stop the relaxation phase if maximum number of inner iterations
        // reached
        if (inner_iter == max_iter_new_pairs.inner - 1) {
            if (!assume_many_iterations) {
                R_previous = middle_point.R;
                orient_previous = orient;
            }
            log_man << "Stopped relaxation phase: maximum number of inner "
                       "iterations "
                    << inner_iter + 1 << " reached.\n";
            ++num_esmax;
            break;
        }

        // Use LBFGS to rotate and translate the dimer
        curvature = performRotationAndTranslationPure(
            middle_point, approx_obs, orient, R_new, rot_info, trans_info);

        // Check if new active frozen atoms and update 'atoms_config' and
        // 'atoms_config'
        if (updateActiveFrozenAtoms(R_new, rot_info, trans_info)) {
            log_man << "updateActiveFrozenAtoms\n";
            gpr_model->optimize(all_obs);
        }

        // Limit the move if any atom-wise step length is larger than
        // 99 % of 0.5*(1-'ratio_at_limit') times the minimum inter-atomic
        // distance
        if (limitTranslation(middle_point.R, trans_info, R_new)) {
            log_man << "Warning: the step length of inner iteration "
                    << inner_iter + 1 << " limited.\n";
        }

        // STOPPING CRITERION FOR INTER-ATOMIC DISTANCES
        // reject the step and stop the relaxation phase if the following
        // does not hold:
        // there is an observed data point so that all inter-atomic distances
        // of the current image are more than 'ratio_at_limit' (by default 2/3)
        // but less than 1/'ratio_at_limit' (3/2) times the corresponding
        // inter-atomic distance of the observed data point, i.e.,
        // |log(r_im/r_nearobs)| < |log(ratio_at_limit)| ( = |log(2/3)| = 0.4055
        // )
        if (inner_iter > 0) {
            if (isInterAtomicDistanceTooBig(R_new, all_obs.R, num_es1)) {
                log_man
                    << "Stopped the relaxation phase after " << inner_iter + 1
                    << " inner iterations: inter-atomic distance changes too"
                       " much compared to 'nearest' observed data point.\n";
                break;
            }
        }

        // Otherwise accept the step and continue the relaxation
        middle_point.R = R_new;

        if (debug_level == DEBUG_L3) writeLatticeData(middle_point);
#ifdef WITH_HDF5
        writeHDF5(
            middle_point,
            "outer_loop/" + std::to_string(outer_iter_idx) + "/inner_loop",
            inner_iter);
#endif
    }
    if (debug_level == DEBUG_L2) writeLatticeData(middle_point);
#ifdef WITH_HDF5
    // TODO(rg): Use the debug_levels
    writeHDF5(middle_point, "outer_loop", outer_iter_idx);
#endif

    return curvature;
}

#ifdef WITH_HDF5
// TODO(rg): Take the openFlags as a parameter
void AtomicDimer::writeHDF5(const gpr::Observation& middle_point,
                            const std::string& addl_group, const size_t& idx)
{
    H5Easy::File file("gpr_optim_out.h5", H5Easy::File::OpenOrCreate);
    H5Easy::dump(file,
                 "/" + addl_group + "/" + std::to_string(idx) + "/positions",
                 middle_point.R.extractEigenVector());
    H5Easy::dump(file,
                 "/" + addl_group + "/" + std::to_string(idx) + "/gradients",
                 middle_point.G.extractEigenVector());
    H5Easy::dump(file, "/" + addl_group + "/" + std::to_string(idx) + "/energy",
                 middle_point.E.extractEigenVector());
}
#endif

void AtomicDimer::writeLatticeData(gpr::Observation& middle_point)
{
    // Print the lattice (forget about the last point if they are falling
    // beyond the cell)
    gpr::Observation tmp_obs;
    gpr::io::FileManager file_man;
    std::string full_file_name_E;
    std::string full_file_name_R;
    std::string full_file_name_G;
    gpr::vector3_reg orig_coord = middle_point.R.at(0);
    gpr::vector3_reg current_coord;
    gpr::Index_t Jmax;
    gpr::Index_t Kmax;

    orig_coord = middle_point.R.at(0);
    orig_coord.y -= debug_output_info.offset;
    orig_coord.z -= debug_output_info.offset;
    current_coord = orig_coord;

    // We consider only a small square around the middle point
    // FIXME: make it generic (for arbitrary number of moving atoms)
    Jmax = debug_output_info.offset * 2. / debug_output_info.dy + 1;
    Kmax = debug_output_info.offset * 2. / debug_output_info.dz + 1;

    // make sure that the folder exists
    // FIXME: won't work on Windows system
    int error =
        system((std::string("mkdir -p ") + debug_output_info.out_dir).c_str());
    gpr::assertMsg(error == 0, "Error! Can't create folder '" +
                                   debug_output_info.out_dir +
                                   "'. Command exited with the code " +
                                   std::to_string(error));

    log_man << "...Writing output files: " << file_counter << "\n";

    gpr::Field<double> loc_E;
    gpr::Coord loc_R, loc_G;
    gpr::Coord loc_grid;

    loc_grid.resize(Jmax, 3 * Kmax);
    loc_E.resize(Jmax, Kmax);
    loc_R.resize(Jmax * Kmax, middle_point.R.getNumCols());
    loc_G.resize(Jmax * Kmax, middle_point.G.getNumCols());

    // Assuming only 1 moving atom!
    gpr::Index_t counter = 0;
    for (gpr::Index_t j = 0; j < Jmax; ++j) {
        current_coord.y = orig_coord.y + j * debug_output_info.dy;
        current_coord.z = orig_coord.z;
        for (gpr::Index_t k = 0; k < Kmax; ++k) {
            tmp_obs = middle_point;
            current_coord.z = orig_coord.z + k * debug_output_info.dz;
            tmp_obs.R.set(0, 0, current_coord);

            gpr_model->calculatePotential(tmp_obs);

            loc_R.assignToRow(counter, tmp_obs.R.getInternalVector());
            loc_G.assignToRow(counter, tmp_obs.G.getInternalVector());
            loc_E[counter] = tmp_obs.E[0];
            loc_grid.set(j, k, current_coord);
            ++counter;
        }
    }

    full_file_name_E =
        debug_output_info.out_dir + "/" + debug_output_info.file_name_E +
        std::to_string(file_counter) + "." + debug_output_info.file_extension;

    full_file_name_R =
        debug_output_info.out_dir + "/" + debug_output_info.file_name_R +
        std::to_string(file_counter) + "." + debug_output_info.file_extension;

    full_file_name_G =
        debug_output_info.out_dir + "/" + debug_output_info.file_name_G +
        std::to_string(file_counter) + "." + debug_output_info.file_extension;

    file_man.writeXYZFile(full_file_name_R, loc_grid, std::ios::app);
    //    file_man.writeXYZFile(full_file_name_R, loc_R, std::ios::app);
    file_man.writeXYZFile(full_file_name_G, loc_G, std::ios::app);
    file_man.writePlainFile(full_file_name_E, loc_E, std::ios::app);

    ++file_counter;
}

double AtomicDimer::performRotationAndTranslationPure(
    gpr::Observation& middle_point, gpr::Observation& approx_obs,
    gpr::Coord& orient, gpr::Coord& R_new, gpr::LBFGSInfo& rot_info,
    gpr::LBFGSInfo& trans_info)
{
    double curvature = 0.;

    rot_info.deltaOrient.clear();
    rot_info.deltaF.clear();

    // If necessary, rotate the dimer and re-calculate approximated
    // energy and gradient at image 1
    for (gpr::Index_t inner_iter_rot = 0;
         inner_iter_rot < max_iter_relax_rot.inner; ++inner_iter_rot) {
        log_man << "  Inner rotational iteration: " << inner_iter_rot << "\n";
        bool terminate_loop =
            rotateDimer(stop_citeria_gpr.angle_rotation, orient, rot_info,
                        middle_point, approx_obs);

        if (terminate_loop) {
            break;
        }
    }

    // Translate the dimer
    curvature = calculateCurvature(orient, approx_obs);
    translateDimer(middle_point, orient, approx_obs, curvature, trans_info,
                   R_new);

    return curvature;
}

void AtomicDimer::pruneObservedData()
{
    if (use_prune == true) {
        // TODO: Make this generic, use parameters
        if (all_obs.E.getSize() > start_prune_at * 2) {
            log_man << "Reached twice pruning size at " << start_prune_at * 2
                    << "\n";
            prune_threshold = prune_threshold / 10;
            log_man << "Reduced prune threshold to " << prune_threshold << "\n";
            start_prune_at = start_prune_at * 2;
            log_man << "Updated pruning start point " << start_prune_at << "\n";
        }
        if (all_obs.E.getSize() > start_prune_at) {
            log_man << "Dropping values (max " << nprune_vals << ")\n";
            std::vector<double> allMax;
            allMax.resize(all_obs.E.getSize());
            for (gpr::Index_t i = 0; i < all_obs.E.getSize(); ++i) {
                allMax[i] = all_obs.G.getMaxAbsEltInRow(i);
            }
            for (gpr::Index_t i = 0; i < nprune_vals; ++i) {
                auto k = std::max_element(allMax.begin(), allMax.end() - 1);
                if (all_obs.G.getMaxAbsEltInRow((gpr::Index_t)std::distance(
                        allMax.begin(), k)) > prune_threshold) {
                    log_man << "Size is " << all_obs.E.getSize() << "\n";
                    log_man
                        << all_obs.G.getMaxAbsEltInRow(
                               (gpr::Index_t)std::distance(allMax.begin(), k))
                        << " at " << std::distance(allMax.begin(), k)
                        << " is bad; dropping\n";
                    all_obs.E.deleteRow(
                        (gpr::Index_t)std::distance(allMax.begin(), k));
                    all_obs.R.deleteRow(
                        (gpr::Index_t)std::distance(allMax.begin(), k));
                    all_obs.G.deleteRow(
                        (gpr::Index_t)std::distance(allMax.begin(), k));
                    allMax.erase(k);
                } else {
                    log_man << all_obs.G.getMaxAbsEltInRow(i) << " for " << i
                            << " is good; keeping\n";
                    log_man << "Size is now " << all_obs.E.getSize() << "\n";
                }
            }
        }
    }
}

} /* namespace atmd */
