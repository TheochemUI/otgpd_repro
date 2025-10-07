/*
 * AtomicDimer.inl
 *
 *  Created on: 18 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_ATOMICDIMER_INL_
#define GPR_ATOMICDIMER_INL_

#ifdef WITH_HDF5

#define H5_USE_EIGEN
#include <highfive/H5Easy.hpp>

#ifdef H5_USE_EIGEN
#include <Eigen/Eigen>
#endif

#endif

namespace atmd {

inline double AtomicDimer::calculateCurvature(
    const gpr::Coord& orient, const gpr::Observation& approx_obs)
{
    double curvature = 0.;

    for (gpr::Index_t j = 0; j < approx_obs.G.getNumCols(); ++j) {
        curvature += (approx_obs.G(1, j) - approx_obs.G(0, j)) * orient(0, j);
    }
    curvature /= dimer_sep;

    return curvature;
}

inline Eigen::VectorXd AtomicDimer::getFinalForceAtMidPoint()
{
    Eigen::VectorXd force(middle_point_final.G.getNumCols());

    for (gpr::Index_t n = 0; n < force.rows(); ++n)
        force(n) = middle_point_final.G(0, n);

    return force;
}

inline Eigen::VectorXd AtomicDimer::getFinalCoordOfMidPoint()
{
    Eigen::VectorXd coord(middle_point_final.R.getNumCols());

    for (gpr::Index_t n = 0; n < coord.rows(); ++n)
        coord(n) = middle_point_final.R(0, n);

    return coord;
}

inline gpr::Index_t AtomicDimer::getTotalForceCalls()
{
    return num_of_gen_potential_calls;
}

inline gpr::Index_t AtomicDimer::getTotalGPRForceCalls()
{
    return num_of_gpr_potential_calls;
}

inline gpr::Index_t AtomicDimer::getIterations()
{
    return totalIterations;
}

inline double AtomicDimer::getFinalEnergy()
{
    return middle_point_final.E(0, 0);
}

inline double AtomicDimer::getFinalCurvature()
{
    return curvature_final;
}

inline gpr::Coord* AtomicDimer::getFinalOrientation()
{
    return &orient_final;
}

inline gpr::Coord& AtomicDimer::getMidPointCoordForAllObservationPointsRef()
{
    return all_obs.R;
}

/**
 * @brief Return energies for all observation points.
 */
inline gpr::Field<double>& AtomicDimer::getEnergyForAllObservationPointsRef()
{
    return all_obs.E;
}

/**
 * @brief Return gradients for all observation points.
 */
inline gpr::Coord& AtomicDimer::getGradientForAllObservationPointsRef()
{
    return all_obs.G;
}

/**
 * @brief Get GPR model.
 */
inline gpr::GaussianProcessRegression* AtomicDimer::getGPRModel()
{
    return gpr_model;
}

inline void AtomicDimer::setAtomsConfiguration(
    const gpr::AtomsConfiguration& _atoms_config)
{
    atoms_config = _atoms_config;
}

inline void AtomicDimer::assignFinalResults(
    const gpr::Observation& middle_point,
    const std::pair<double, gpr::Coord>& eigen_data)
{
    middle_point_final.R = middle_point.R;
    middle_point_final.E = middle_point.E;
    middle_point_final.G = middle_point.G;
    curvature_final = eigen_data.first;
    orient_final = eigen_data.second;
}

template <typename T>
void AtomicDimer::callGeneralPotentialFromEON(
    const gpr::vector3_reg (&cell_dimensions)[3], gpr::Observation& mid_point,
    T& potential)
{
    long N;
    const double* R;
    const int* atomicNrs;
    double* F;
    double* U;
    const double box[9] = {
        cell_dimensions[0].x, cell_dimensions[0].y, cell_dimensions[0].z,
        cell_dimensions[1].x, cell_dimensions[1].y, cell_dimensions[1].z,
        cell_dimensions[2].x, cell_dimensions[2].y, cell_dimensions[2].z};

    // FIXME: check if the size is correct!!!
    mid_point.E.resize(mid_point.R.getNumRows(), 1);
    mid_point.G.resize(mid_point.R.getNumRows(), mid_point.R.getNumCols());

    // Total number of moving atoms
    N = mid_point.R.getNumPoints();

    // Position at middle point
    R = mid_point.R.getInternalVector().data();

    // Force at middle point
    F = mid_point.G.getInternalVector().data();

    // Energy at middle point
    U = mid_point.E.getInternalVector().data();

    // FIXME: Maybe ditch atomicNrs altogether?
    // Atom types
    atomicNrs = reinterpret_cast<const int*>(
        atoms_config.is_frozen.getInternalVector().data());

    // Call for the potential calculations
    potential.force(N, R, atomicNrs, F, U, nullptr /*variance*/, box);
}

template <typename Pot>
void AtomicDimer::evaluateAccurateEnergyAndForceAtMidPoint(
    gpr::Observation& middle_point, Pot& general_potential)
{
    aux::ProblemSetUp problem_setup;

    if (middle_point_init.E.isEmpty()) {
        // Calculate energy and gradient at the middle point of the dimer
        problem_setup.calculateGeneralPotential(atoms_config, cell_dimensions,
                                                middle_point, general_potential,
                                                num_of_gen_potential_calls);

        // OR call for
        // callGeneralPotentialFromEON(...);

        // Set zero level of biased potential to the energy of the middle point
        // of the initial dimer
        E_zero_level = middle_point.E;
        middle_point.E.setZero();
        problem_setup.cutOffEnergy(E_zero_level, all_obs.E);

        // Assemble coordinates, energy and gradients of all Observation& points
        all_obs.append(middle_point);
    } else {
        E_zero_level = middle_point_init.E;

        problem_setup.calculateGeneralPotential(atoms_config, cell_dimensions,
                                                middle_point, general_potential,
                                                num_of_gen_potential_calls);

        // OR call for
        // callGeneralPotentialFromEON(...);

        middle_point.E = middle_point_init.E;
        middle_point.G = middle_point_init.G;
        problem_setup.cutOffEnergy(E_zero_level, middle_point.E);

        all_obs.R = middle_point_init.R;
        all_obs.E = middle_point_init.E;
        all_obs.G = middle_point_init.G;
        problem_setup.cutOffEnergy(E_zero_level, all_obs.E);
    }
}

template <typename Pot>
void AtomicDimer::evaluateAccurateEnergyAndForce(gpr::Observation& obs,
                                                 Pot& general_potential)
{
    aux::ProblemSetUp problem_setup;

    problem_setup.calculateGeneralPotential(atoms_config, cell_dimensions, obs,
                                            general_potential,
                                            num_of_gen_potential_calls);

    // OR call for
    // callGeneralPotentialFromEON(...);

    problem_setup.cutOffEnergy(E_zero_level, obs.E);

    all_obs.append(obs);
}

template <typename Pot>
void AtomicDimer::performInitialRotations(gpr::Observation& middle_point,
                                          gpr::Coord& orient, Pot& potential)
{
    gpr::Coord temp;
    gpr::Index_t outer_iter;
    gpr::Coord orient_old;
    gpr::LBFGSInfo rot_info;  // Information on rotation
    gpr::Observation approx_obs;
    double stop_citeria_gpr_init;

    stop_citeria_gpr_init =
        std::min(0.01, stop_citeria_dimer.angle_rotation * 0.1);

#ifdef WITH_HDF5
    // TODO(rg): This can use the writeHDF5 function as well, its the same
    // structure, maybe pass the mode via parameter
    H5Easy::File file("gpr_optim_out.h5", H5Easy::File::Overwrite);
#endif
    for (outer_iter = 0; outer_iter < max_iter_init_rot.outer; ++outer_iter) {
        log_man << "Outer rotational iteration: " << outer_iter << "\n";

        // 5.d) Check rotational convergence using accurate forces F0 and F1
        // (not
        //    invoked in MATLAB)
        if (isRotationalConvergenceReached(orient, middle_point)) {
            log_man << "Rotated the initial dimer in " << outer_iter + 1
                    << " outer"
                       " iterations (total number of image evaluations: "
                    << getNumEvaluations() << ").\n";
            break;
        }

        // 5.a) Update the GP model based on the energy and force evaluations
        gpr_model->setHyperparameters(all_obs, atoms_config);
        // NOTE(rg): no need to optimize hyperparameters during rotation
        gpr_model->optimize(all_obs);
        // gpr_model->updateModelWithFullData(all_obs);

        // Define the initial dimer orientation for the inner iteration loop
        orient_old = orient;
        orient = orient_init;

        // Calculate approximated energy and gradient at the middle point
        // and image 1 of the dimer
        approx_obs.R = middle_point.R;
        updateLocation(orient, middle_point.R, temp);
        approx_obs.R.append(temp);
        gpr_model->calculatePotential(approx_obs);
#ifdef WITH_HDF5
        H5Easy::dump(
            file,
            "/initial_rotations/" + std::to_string(outer_iter) + "/positions",
            approx_obs.R.extractEigenVector());
        H5Easy::dump(
            file,
            "/initial_rotations/" + std::to_string(outer_iter) + "/gradients",
            approx_obs.G.extractEigenVector());
        H5Easy::dump(
            file,
            "/initial_rotations/" + std::to_string(outer_iter) + "/energy",
            approx_obs.E.extractEigenVector());
#endif

        rot_info.clear();
        rot_info.num_cg_iter = middle_point.R.getNumCols();
        rot_info.num_lbfgs_iter = middle_point.R.getNumCols();

        // 5.b) Rotate the dimer until rotational convergence using the GP
        // approximation of the energy gradient.
        gpr::Index_t inner_iter;
        for (inner_iter = 0; inner_iter < max_iter_init_rot.inner;
             ++inner_iter) {
            log_man << " Inner rotational iteration: " << inner_iter << "\n";
            bool terminate_loop =
                rotateDimer(stop_citeria_gpr_init, orient, rot_info,
                            middle_point, approx_obs);
            if (terminate_loop) break;
        }
        checkIfMaxIterReachedAndPrintWarning(
            inner_iter, max_iter_init_rot.inner, "dimer rotation");

        if (inner_iter == 0) {
            log_man << "Rotated the initial dimer in " << outer_iter + 1
                    << " outer iterations (total number of image evaluations: "
                    << getNumEvaluations() << ").\n";
            log_man << "WARNING: Dimer orientation converged on the GP surface "
                       "(stop criteria = "
                    << stop_citeria_gpr_init
                    << "), but not on the true PES (stop criteria = "
                    << stop_citeria_dimer.angle_rotation << ").\n";
            break;
        } else {
            if (outer_iter > 0) {
                if (checkRotationalConvergence(
                        orient, orient_old,
                        stop_citeria_dimer.angle_rotation)) {
                    log_man << "Rotated the initial dimer in " << outer_iter + 1
                            << " outer iterations (total number of image "
                               "evaluations: "
                            << getNumEvaluations() << ".\n";
                    break;
                }
            }
            // 5.c) Evaluate accurate energy E1 and force F1 at R1.
            updateLocation(orient, middle_point.R, image1.R);
            evaluateAccurateEnergyAndForce(image1, potential);
        }
    }
    checkIfMaxIterReachedAndPrintWarning(outer_iter, max_iter_init_rot.outer,
                                         "initial rotations");
}

template <typename Pot>
std::pair<double, gpr::Coord> AtomicDimer::performGPRIterations(
    gpr::Observation& middle_point, gpr::Coord& orient, Pot& potential)
{
    gpr::Coord orient_init_gpr;  // Initial orientation before applying GPR???
    gpr::Coord R_previous;
    gpr::Coord R_latest_conv;
    gpr::Coord orient_previous;
    gpr::Coord orient_latest_conv;
    gpr::Coord R_new;
    gpr::Observation approx_obs;  // Approximated observation at image1
    double stop_dimer_gp = 0.;
    gpr::Index_t outer_iter;
    gpr::Index_t num_esmax = 0;
    gpr::Index_t num_es1 = 0;  // FIXME:: this variable is reduntant
    double curvature = 0.;

    orient_init_gpr = orient;

    for (outer_iter = 0; outer_iter < max_iter_new_pairs.outer; ++outer_iter) {
        log_man << "Outer relaxation iteration: " << outer_iter << "\n";

        pruneObservedData();

        // 6.a) Update the GP model based on the energy and force evaluations.
        bool update_sexpat_param = false;
        if (outer_iter < 1 ||
            atoms_config.n_pt > gpr_model->getSexpAtCovarianceFunction()
                                    ->getLengthScaleRef()
                                    .getNumCols()) {
            update_sexpat_param = true;
        }
        gpr_model->setHyperparameters(all_obs, atoms_config,
                                      update_sexpat_param);
        gpr_model->optimize(all_obs);

        // Define the convergence threshold for the relaxation phase
        if (divisor_stop_criteria_gpr > 0) {
            // If this option is set on, the GP convergence threshold is
            // 1/'divisor_stop_criteria_gpr' of the smallest accurate 'maxF_R'
            // obtained so far, but not less than 1/10 of the final threshold
            stop_dimer_gp = std::max(
                middle_point.G.getMaxAbsElt() / divisor_stop_criteria_gpr,
                stop_citeria_dimer.force * 0.1);
        } else {
            // otherwise the GP convergence threshold is always 1/10 of the
            // final threshold
            stop_dimer_gp = stop_citeria_dimer.force * 0.1;
        }

        // Define the initial dimer for the relaxation phase
        middle_point.R = middle_point_init.R;
        orient = orient_init_gpr;

        // Define the start path for the relaxation phase
        defineStartPath(R_latest_conv, orient_latest_conv, orient_init_gpr,
                        R_previous, orient_previous, middle_point, orient);

        // 6.b) Rotate and translate the dimer until early stopping or
        // convergence
        //      using the GP approximation of the energy gradient.
        curvature = performRotationAndTranslationWithGPR(
            stop_dimer_gp, middle_point, orient, R_latest_conv, R_previous,
            orient_latest_conv, orient_previous, R_new, num_esmax, num_es1,
            outer_iter);

        // 6.c) Evaluate accurate energy E0 and force F0 at R0.
        evaluateAccurateEnergyAndForce(middle_point, potential);

        // 6.d) Check final convergence using accurate force F0.
        if (isFinalConvergenceReached(middle_point, stop_citeria_dimer.force)) {
            log_man << "Final convergence obtained after " << outer_iter + 1
                    << " relaxation phases (total number of image evaluations: "
                    << getNumEvaluations() << ").\n";
            break;
        }
    }

    if (outer_iter == max_iter_new_pairs.outer) {
        log_man << "Maximum number of outer iterations ("
                << max_iter_new_pairs.outer << ") reached.\n";
    }
    totalIterations = outer_iter + 1;

    return {curvature, orient};
}

template <typename Pot>
void AtomicDimer::execute(Pot& general_potential)
{
    gpr::Observation middle_point;  // Observation at the middle point:
                                    // coordinates, energy and gradient

    gpr::Field<double>
        E_gp;           // Vector gathering approximated energy of the
                        // middle point of the dimer for each inner iteration
    gpr::Coord orient;  // Unit vector along the direction of the dimer
    gpr::Coord orient_old;  // Unit vector along the direction of the dimer
                            // form the previous iteration
    std::chrono::duration<double> elp_time;
    auto start = std::chrono::steady_clock::now();
    std::pair<double, gpr::Coord> eigen_data;

    if (!isAtomsConfigurationCorrect()) {
        log_man
            << "Error! The structure of atoms configuration is incorrect!\n";
        return;
    }

    orient_init.normalizeRow(0);
    middle_point.R = middle_point_init.R;
    orient = orient_init;

    // 1) Evaluate accurate energy E0 and force F0 at the middle point R0
    evaluateAccurateEnergyAndForceAtMidPoint(middle_point, general_potential);

    // 2) Check final convergence using accurate force F0
    if (isFinalConvergenceReached(middle_point, stop_citeria_dimer.force)) {
        log_man << "Final convergence obtained in the beginning ("
                << getNumEvaluations() << " image evaluations).\n";
        eigen_data.first = 0;
        eigen_data.second = orient_init;
        assignFinalResults(middle_point, eigen_data);
        return;
    }

    if (max_iter_init_rot.outer > 0) {
        // 3) Evaluate energy E1 and force F1 at R1
        updateLocation(orient, middle_point.R, image1.R);
        evaluateAccurateEnergyAndForce(image1, general_potential);

        // 5) Repeat initial rotations until rotational convergence
        performInitialRotations(middle_point, orient, general_potential);
    }

    if (all_obs.R.getNumRows() < 2) {
        updateLocation(orient, middle_point.R, image1.R);
        evaluateAccurateEnergyAndForce(image1, general_potential);
        log_man << "Evaluated image 1 of the dimer for the initial GP model.\n";
    }

    // 6) Repeat GPR iterations until final convergence
    eigen_data = performGPRIterations(middle_point, orient, general_potential);

    assignFinalResults(middle_point, eigen_data);

    if (debug_level == DEBUG_L1) {
        log_man << "\n";
        log_man << "Call statistics:\n";
        log_man << " General potential: " << getTotalForceCalls() << "\n";
        log_man << " GPR     potential: " << getTotalGPRForceCalls() << "\n";
        log_man << "\n";
    }

    elp_time = std::chrono::steady_clock::now() - start;
    log_man << "Elapsed time: " << elp_time.count() << "s\n";

    middle_point.R.print();
    middle_point.E.print();
    middle_point.G.print();
}

} /* namespace atmd */

#endif /* GPR_ATOMICDIMER_INL_ */
