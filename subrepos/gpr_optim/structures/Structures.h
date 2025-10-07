/*
 * Structures.h
 *
 *  Created on: 17 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef STRUCTURES_STRUCTURES_H_
#define STRUCTURES_STRUCTURES_H_

#include <Eigen/Dense>
#include <string>

#include "../data_types/Coord.h"
#include "../data_types/Field.h"
#include "../gpr/Enums.h"
#include "../managers/io/LogManager.h"

namespace gpr {

/**
 * @brief Structure of a string key and template value.
 */
template <typename T>
struct KeyValuePair {
    std::string key;
    T value;
};

/**
 * @brief A simple pair with identical datatypes.
 */
template <typename T>
struct Pair {
    T first;
    T second;
};

/**
 * @brief A simple triplet with identical datatypes.
 */
template <typename T>
struct Triplet {
    T first;
    T second;
    T third;
};

/**
 * @brief Structure of directional derivatives.
 */
template <typename T>
struct Derivatives {
    T D1;
    T D2;
    T D12;
};

/**
 * @brief Structure of pair of indices.
 */
struct Indices2D {
    uint32_t i;
    uint32_t j;
};

/**
 * @brief Structure of limits for iterational loops.
 */
struct IterationsLimits {
    Index_t outer = 0;
    Index_t inner = 0;

    void setZero()
    {
        outer = inner = 0;
    }
};

/**
 * @brief Structure of GPR properties.
 */
struct GPRSetup {
    double sigma2;
    double jitter_sigma2;
    double tol_fun;
    double tol_X;
    uint8_t optimization_alg;

    GPRSetup()
    {
        clear();
    }

    void clear()
    {
        sigma2 = 0.;
        jitter_sigma2 = 0.;
        tol_fun = 0.;
        tol_X = 0.;
        optimization_alg = 0;
    }
};

/**
 * @brief Settings for the optimization algorithm.
 */
struct OptimizationAlgorithmSettings {
    bool check_derivative;  // TODO(rg): Unused?
    uint8_t report_level;
    Index_t max_iter;
    double tolerance_func;
    double tolerance_sol;
    // SCG parameters
    double lambda_limit;
    double lambda;
    // ADAM-specific parameters
    double learning_rate;
    double learning_rate_decay;
    double beta1;
    double beta2;
    double epsilon;
    double weight_decay;
    bool amsgrad;
    OptimizationAlgorithmSettings()
    {
        setDefault();
    }

    ~OptimizationAlgorithmSettings()
    {
        setDefault();
    }

    void setDefault()
    {
        check_derivative = false;
        report_level = 1;
        max_iter = 400;
        tolerance_func = 1e-4;
        tolerance_sol = 1e-4;
        lambda_limit = 1e20;
        lambda = 10.;
        // ADAM-specific parameters
        // NOTE(rg): Generally a high learning rate here is fine, just need to
        // reach a low value, fast
        learning_rate = 1e-1;  // from PyTorch
        learning_rate_decay = 0.99;
        beta1 = 0.9;
        beta2 = 0.999;
        epsilon = 1.0e-8;
        weight_decay = 0.1;
        amsgrad = true;
    }
};

/**
 * @brief Structure of energy and gradient stored in Eigen-compatible form.
 */
struct EnergyAndGradient {
    double* energy;
    Eigen::VectorXd* gradient;
};

/**
 * @brief Structure of the input parameters.
 * All input parameters are pairs of string names and values.
 */
struct InputParameters {
    KeyValuePair<int> i_dist = {"i_dist", 0};
    KeyValuePair<int> i_run = {"i_run", 0};

    KeyValuePair<double[4]> dist_sp = {"dist_sp", {0., 0., 0., 0.}};

    /**
     * Activation distance for moving+frozen atom pairs (inf if all active).
     */
    KeyValuePair<int> actdist_fro = {"actdist_fro", 5};

    /**
     * Coordinates of the initial data points (N_obs x D).
     */
    KeyValuePair<double*> R_all_init = {"R_all_init", nullptr};

    /**
     * Energies at the initial data points (N_obs x 1).
     */
    KeyValuePair<double*> E_all_init = {"E_all_init", nullptr};

    /**
     * Gradients at the initial data points (N_obs x D).
     */
    KeyValuePair<double*> G_all_init = {"G_all_init", nullptr};

    /**
     * Coordinates of the middle point of the initial dimer (1 x D).
     */
    KeyValuePair<double[3]> R_init = {"R_init", {0., 0., 0.}};

    /**
     * Energy at the middle point of the initial dimer (if not empty, 'R_init'
     * should be included in 'R_all_init').
     */
    KeyValuePair<double> E_init = {"E_init", 0.};

    /**
     * Gradient at the middle point of the initial dimer (1 x D) (if not empty,
     * 'R_init' should be included in 'R_all_init').
     */
    KeyValuePair<double[3]> G_init = {"G_init", {0., 0., 0.}};

    /**
     * Unit vector along the direction of the initial dimer (1 x D).
     */
    KeyValuePair<double[3]> orient_init = {"orient_init", {0., 0., 0.}};

    /**
     * Dimer separation (distance from the middle point of the dimer to the two
     * images).
     */
    KeyValuePair<double> dimer_sep = {"dimer_sep", 0.};

    /**
     * A function defining the rotation step.
     */
    KeyValuePair<std::string> method_rot = {"method_rot", "none"};

    /**
     * A function defining the translation step.
     */
    KeyValuePair<std::string> method_trans = {"method_trans", "none"};

    /**
     * Parameters of the translation method (shape depends on 'method_trans').
     */
    KeyValuePair<double[2]> param_trans = {"param_trans", {0., 0.}};

    /**
     * Indicator if image 1 of the dimer is evaluted (1) or not (0) after each
     * relaxation phase in addition to the middle point of the dimer.
     */
    KeyValuePair<int> eval_image1 = {"eval_image1", 0};

    /**
     * Final convergence threshold for 'maxF_R', which is the maximum component
     * of the force acting on the middle point of the dimer (i.e., the algorithm
     * is stopped when all components of the accurate force are below
     * 'T_dimer').
     */
    KeyValuePair<double> T_dimer = {"T_dimer", 0.};

    /**
     * Indicator if the initial rotations are performed without GP (1) or with
     * GP (0).
     */
    KeyValuePair<int> initrot_nogp = {"initrot_nogp", 0};

    /**
     * Convergence threshold for rotation angle in the initial rotations (the
     * dimer is not rotated when the estimated rotation angle is less than
     * this).
     */
    KeyValuePair<double> T_anglerot_init = {"T_anglerot_init", 0.};

    /**
     * Maximum number of initial rotations (0 if initial rotations skipped).
     */
    KeyValuePair<int> num_iter_initrot = {"num_iter_initrot", 0};

    /**
     * Convergence threshold for rotation angle during a relaxation phase.
     */
    KeyValuePair<double> T_anglerot_gp = {"T_anglerot_gp", 0.};

    /**
     * Maximum number of rotation iterations per translation during a relaxation
     * phase.
     */
    KeyValuePair<int> num_iter_rot_gp = {"num_iter_rot_gp", 0};

    /**
     * If this option is set on (> 0), the convergence threshold for a
     * relaxation phase is 1/'divisor_T_dimer_gp' of the smallest accurate
     * 'maxF_R' obtained so far, but not less than 1/10 of 'T_dimer' (otherwise
     * the GP convergence threshold is always 1/10 of 'T_dimer').
     */
    KeyValuePair<int> divisor_T_dimer_gp = {"divisor_T_dimer_gp", 0};

    /**
     * Maximum displacement of the middle point of the dimer from the nearest
     * observed data point (the relaxation phase is stopped if 'disp_max' is
     * reached).
     */
    KeyValuePair<double> disp_max = {"disp_max", 0.};

    /**
     * Limit for the ratio (< 1) of inter-atomic distances between image and its
     * "nearest" observed data point (the relaxation phase is stopped if
     * 'ratio_at_limit' is reached for any image).
     */
    KeyValuePair<double> ratio_at_limit = {"ratio_at_limit", 0.};

    /**
     * Number of outer iterations started from the initial location 'R_init'
     * (after that, each relaxation phase is started from the latest converged
     * dimer).
     */
    KeyValuePair<int> num_bigiter_initloc = {"num_bigiter_initloc", 0};

    /**
     * Number of outer iterations where the hyperparameter optimization is
     * started from values initialized based on the range of the current data
     * (after that, the optimization is started from the values of the previous
     * round).
     */
    KeyValuePair<int> num_bigiter_initparam = {"num_bigiter_initparam", 0};

    /**
     * Maximum number of outer iterations (new pairs of observations).
     */
    KeyValuePair<int> num_bigiter = {"num_bigiter", 0};

    /**
     * Maximum number of inner iterations (steps during a relaxation phase).
     */
    KeyValuePair<int> num_iter = {"num_iter", 0};

    /**
     * indicator if 'num_iter' is assumed to be much larger than required for
     * dimer convergence on accurate energy surface (if not, the next relaxation
     * phase is continued from the current path if 'num_iter' is reached).
     */
    KeyValuePair<int> islarge_num_iter = {"islarge_num_iter", 0};

    /**
     * Path to the data file required to continue from a cancelled run (empty if
     * started normally from the beginning).
     */
    KeyValuePair<std::string> load_file = {"load_file", "none"};

    /**
     * Path to the data file where data is saved (empty if not saved).
     */
    KeyValuePair<std::string> save_file = {"save_file", "none"};

    /* * --- FOR PRUNE * */
    /*
     * Activation
     */
    KeyValuePair<bool> use_prune = {"use_prune", false};

    /**
     * Number of values to keep before prune begins
     */
    KeyValuePair<int> start_prune_at = {"start_prune_at", 10};

    /**
     * Number of values drop
     * Note that no values will be dropped if they are smaller than the
     * threshold
     */
    KeyValuePair<int> nprune_vals = {"nprune_vals", 3};

    /**
     * Threshold below which candidates will be dropped
     */
    KeyValuePair<double> prune_threshold = {"prune_threshold", 0.5};

    /* * --- FOR DEBUGGING * */
    /**
     * Name of the output folder
     */
    KeyValuePair<std::string> debug_output_dir = {"debug_output_dir", "output"};

    /**
     * Name of the file with positions
     */
    KeyValuePair<std::string> debug_output_file_R = {"debug_output_file_R",
                                                     "position"};

    /**
     * Name of the file with energies
     */
    KeyValuePair<std::string> debug_output_file_E = {"debug_output_file_E",
                                                     "energy"};

    /**
     * Name of the file with gradients
     */
    KeyValuePair<std::string> debug_output_file_G = {"debug_output_file_G",
                                                     "gradient"};

    /**
     * Extension of the output files
     */
    KeyValuePair<std::string> debug_output_file_extension = {
        "debug_output_file_extension", "dat"};

    /**
     * Offset for a sub-box (used to constract a smaller 2D cell around the
     * middle pioint)
     */
    KeyValuePair<double> debug_offset_from_mid_point = {
        "debug_offset_from_mid_point", 3.};

    /**
     * Grid step in Y direction
     */
    KeyValuePair<double> debug_dy = {"debug_dy", 0.1};

    /**
     * Grid step in Z direction
     */
    KeyValuePair<double> debug_dz = {"debug_dz", 0.1};

    /**
     * Debug level
     */
    KeyValuePair<int> debug_level = {"debug_level", 0};
    /* * --- FOR DEBUGGING * */

    /**
     * Accurate potential and gradient function.
     */
    KeyValuePair<std::string> pot_general = {"pot_general", "none"};

    /**
     * Structure array including information about the configurations necessary
     * for the GP model
     *  - conf_info.conf_fro: coordinates of active frozen atoms (N_fro x 3)
     *  - conf_info.atomtype_mov: atomtype indices for moving atoms (1 x N_mov)
     *  - conf_info.atomtype_fro: pairtype indices for active frozen atoms (1 x
     * N_fro)
     *  - conf_info.pairtype: pairtype indices for pairs of atomtypes (n_at x
     * n_at)
     *  - conf_info.n_pt: number of active pairtypes
     */
    KeyValuePair<std::string*> conf_info = {"conf_info", nullptr};

    /**
     * Structure array including information about inactive frozen atoms
     *  - conf_info_inactive.conf_ifro: coordinates of inactive frozen atoms
     *                                  (N_ifro x 3)
     *  - conf_info_inactive.atomtype_ifro: atomtype indices for inactive frozen
     *                                      atoms (1 x N_ifro)
     */
    KeyValuePair<std::string*> conf_info_inactive = {"conf_info_inactive",
                                                     nullptr};

    /**
     * Dimensions of the computational domain.
     */
    KeyValuePair<double[9]> cell_dimensions = {
        "cell_dimensions", {0., 0., 0., 0., 0., 0., 0., 0., 0.}};

    /*
     * GPR parameters.
     */
    /**
     * Magnitude of a covariance function?
     */
    KeyValuePair<double> gp_sigma2 = {"gp_sigma2", 0.};
    KeyValuePair<double> jitter_sigma2 = {"jitter_sigma2", 0.};
    KeyValuePair<std::string> optimization_alg = {"optimization_alg",
                                                  "SCG_opt"};

    /**
     * Noise variance?
     */
    KeyValuePair<double> sigma2 = {"sigma2", 0.};

    // (calculated on a fly):
    /**
     * SexpAt covariance.
     */
    KeyValuePair<double> magnSigma2 = {"magnSigma2", 0.};

    /**
     * Constant covariance.
     */
    KeyValuePair<double> constSigma2 = {"constSigma2", 0.};

    /**
     * Prior parameters
     */
    KeyValuePair<double> prior_mu = {"prior_mu", 0.};
    KeyValuePair<double> prior_nu = {"prior_nu", 0.};
    KeyValuePair<double> prior_s2 = {"prior_s2", 0.};

    /**
     * Early stopping parameters
     */
    KeyValuePair<DistanceMetricType> es_dist_metric = {"es_dist_metric",
                                                       DistanceMetricType::EMD};
    KeyValuePair<double> es_threshold = {"es_threshold", 1.2};

    /**
     * Farthest point sampling parameters
     */
    KeyValuePair<DistanceMetricType> fps_metric = {"fps_metric",
                                                       DistanceMetricType::EMD};
    KeyValuePair<int> fps_history = {"fps_history", 5};

    /**
     * Dimer details
     */
    // TODO(rg):: Need to actually log these below
    KeyValuePair<double> rotation_removal_projection_threshold = {
        "rotation_removal_projection_threshold",
        std::numeric_limits<double>::infinity()};
    KeyValuePair<std::string> rot_opt = {"rot_opt", "lbfgs"};

    /**
     * Optimization algorithm
     */
    KeyValuePair<std::string> check_derivative = {"check_derivative", "false"};
    KeyValuePair<int> report_level = {"report_level", 1};
    KeyValuePair<int> max_iter = {"max_iter", 400};
    KeyValuePair<double> tolerance_func = {"tolerance_func", 1e-4};
    KeyValuePair<double> tolerance_sol = {"tolerance_sol", 1e-4};
    KeyValuePair<double> lambda_limit = {"lambda_limit", 1e20};
    KeyValuePair<double> lambda = {"lambda", 10.};
    // TODO(rg):: Need to actually log these below
    KeyValuePair<double> learning_rate = {"learning_rate", 0.8};
    KeyValuePair<double> learning_rate_decay = {"learning_rate_decay", 0.999};
    KeyValuePair<double> beta1 = {"beta1", 0.9};
    KeyValuePair<double> beta2 = {"beta2", 0.99};
    KeyValuePair<double> epsilon = {"epsilon", 1.0e-8};
    KeyValuePair<double> weight_decay = {"weight_decay", 0.0};
    KeyValuePair<bool> amsgrad = {"amsgrad", true};

    /**
     * @brief Print the content of the structure.
     */
    void print()
    {
        io::LogManager log_man;

        log_man << i_dist.key << " : " << i_dist.value << "\n";
        log_man << i_run.key << " : " << i_run.value << "\n";
        log_man << dist_sp.key << " : " << dist_sp.value[0] << " "
                << dist_sp.value[1] << " " << dist_sp.value[2] << " "
                << dist_sp.value[3] << "\n";
        log_man << actdist_fro.key << " : " << actdist_fro.value << "\n";
        log_man << dimer_sep.key << " : " << dimer_sep.value << "\n";
        log_man << method_rot.key << " : " << method_rot.value << "\n";
        log_man << method_trans.key << " : " << method_trans.value << "\n";
        log_man << param_trans.key << " : " << param_trans.value[0] << " "
                << param_trans.value[1] << "\n";
        log_man << eval_image1.key << " : " << eval_image1.value << "\n";
        log_man << T_dimer.key << " : " << T_dimer.value << "\n";
        log_man << initrot_nogp.key << " : " << initrot_nogp.value << "\n";
        log_man << T_anglerot_init.key << " : " << T_anglerot_init.value
                << "\n";
        log_man << num_iter_initrot.key << " : " << num_iter_initrot.value
                << "\n";
        log_man << T_anglerot_gp.key << " : " << T_anglerot_gp.value << "\n";
        log_man << num_iter_rot_gp.key << " : " << num_iter_rot_gp.value
                << "\n";
        log_man << divisor_T_dimer_gp.key << " : " << divisor_T_dimer_gp.value
                << "\n";
        log_man << disp_max.key << " : " << disp_max.value << "\n";
        log_man << ratio_at_limit.key << " : " << ratio_at_limit.value << "\n";
        log_man << num_bigiter.key << " : " << num_bigiter.value << "\n";
        log_man << num_iter.key << " : " << num_iter.value << "\n";
        log_man << islarge_num_iter.key << " : " << islarge_num_iter.value
                << "\n";
        log_man << load_file.key << " : " << load_file.value << "\n";
        log_man << save_file.key << " : " << save_file.value << "\n";

        log_man << debug_output_dir.key << " : " << debug_output_dir.value
                << "\n";
        log_man << debug_output_file_R.key << " : " << debug_output_file_R.value
                << "\n";
        log_man << debug_output_file_E.key << " : " << debug_output_file_E.value
                << "\n";
        log_man << debug_output_file_G.key << " : " << debug_output_file_G.value
                << "\n";
        log_man << debug_output_file_extension.key << " : "
                << debug_output_file_extension.value << "\n";
        log_man << debug_offset_from_mid_point.key << " : "
                << debug_offset_from_mid_point.value << "\n";
        log_man << debug_dy.key << " : " << debug_dy.value << "\n";
        log_man << debug_dz.key << " : " << debug_dz.value << "\n";
        log_man << debug_level.key << " : " << (int)debug_level.value << "\n";

        log_man << cell_dimensions.key << " : " << cell_dimensions.value[0]
                << " " << cell_dimensions.value[1] << " "
                << cell_dimensions.value[2] << " " << cell_dimensions.value[3]
                << " " << cell_dimensions.value[4] << " "
                << cell_dimensions.value[5] << " " << cell_dimensions.value[6]
                << " " << cell_dimensions.value[7] << " "
                << cell_dimensions.value[8] << "\n";

        log_man << gp_sigma2.key << " : " << gp_sigma2.value << "\n";
        log_man << jitter_sigma2.key << " : " << jitter_sigma2.value << "\n";
        log_man << optimization_alg.key << " : " << optimization_alg.value
                << "\n";

        log_man << sigma2.key << " : " << sigma2.value << "\n";
        log_man << magnSigma2.key << " : " << magnSigma2.value << "\n";
        log_man << constSigma2.key << " : " << constSigma2.value << "\n";

        log_man << prior_mu.key << " : " << prior_mu.value << "\n";
        log_man << prior_nu.key << " : " << prior_nu.value << "\n";
        log_man << prior_s2.key << " : " << prior_s2.value << "\n";

        log_man << check_derivative.key << " : " << check_derivative.value
                << "\n";
        log_man << report_level.key << " : " << report_level.value << "\n";
        log_man << max_iter.key << " : " << max_iter.value << "\n";
        log_man << tolerance_func.key << " : " << tolerance_func.value << "\n";
        log_man << tolerance_sol.key << " : " << tolerance_sol.value << "\n";
        log_man << lambda_limit.key << " : " << lambda_limit.value << "\n";
        log_man << lambda.key << " : " << lambda.value << "\n";

        log_man << start_prune_at.key << " : " << start_prune_at.value << "\n";
        log_man << nprune_vals.key << " : " << nprune_vals.value << "\n";
        log_man << prune_threshold.key << " : " << prune_threshold.value
                << "\n";
    }
};

/**
 * @brief Structure of atoms positions and types.
 */
struct AtomsPositionAndType {
    Coord positions;
    Field<Index_t> type;

    void resize(const Index_t num_atoms)
    {
        positions.resize(1, num_atoms * 3);
        type.resize(1, num_atoms);
    }

    void clear()
    {
        positions.clear();
        type.clear();
    }
};

/**
 * @brief Stores configurations of all atoms in the system.
 */
struct AtomsConfiguration {
    Coord positions;           //!< Positions of all atoms in the system
    Field<uint8_t> is_frozen;  //!< Indicates wheather atom is frozen or not
    Field<Index_t> id;         //!< ID of each atom (we do not use it)
    Field<int> atomicNrs;      //!< Atomic Number of each atom
    Field<Index_t> type;       //!< Type of each atom

    AtomsPositionAndType atoms_froz_active;  //!< Info on active frozen atoms
    AtomsPositionAndType
        atoms_froz_inactive;         //!< Info on inactive frozen atoms
    AtomsPositionAndType atoms_mov;  //!< Info on moving atoms

    Field<int>
        pairtype;  //!< Pairtype indices for pairs of  atomtypes (n_at x n_at).
                   //!< Active pairtypes are indexed as 0,1,2,...,n_pt-1.
                   //!< Inactive pairtypes are given index -1.
    Index_t n_pt;  //!< Number of active pairtypes

    /**
     * @brief Extract data from a 5-column \e field.
     * The field should have 5 column: 3 coordinates, indicator if atom is
     * frozen, atom ID.
     */
    void assignFromField(const Field<double>& field)
    {
#ifndef NDEBUG
        assertMsg(field.getSize() != 0, "The reference field cannot be empty!");
        assertMsg(field.getNumCols() >= 5,
                  "The reference field should have 5 or more columns!");
#endif

        Index_t counter = 0;

        positions.resize(1, 3 * field.getNumRows());
        is_frozen.resize(1, field.getNumRows());
        id.resize(1, field.getNumRows());

        for (Index_t i = 0; i < field.getNumRows(); ++i) {
            for (Index_t j = 0; j < 3; ++j) {
                positions[counter++] = field(i, j);
            }
            if (field(i, 3) == 0)
                is_frozen[i] = MOVING_ATOM;
            else
                is_frozen[i] = FROZEN_ATOM;
            id[i] = field(i, 4);
        }
    }

    /**
     * @brief Clear the structure.
     */
    void clear()
    {
        positions.clear();
        is_frozen.clear();
        id.clear();
        atoms_froz_active.clear();
        atoms_froz_inactive.clear();
        atoms_mov.clear();
        pairtype.clear();
        n_pt = 0;
    }

    /**
     * @brief Clear distribution of active and inactive atoms.
     */
    void clearDistribution()
    {
        atoms_froz_active.clear();
        atoms_froz_inactive.clear();
        atoms_mov.clear();
        pairtype.clear();
        n_pt = 0;
    }

    /**
     * @brief Assign members of \e this structure to the \e other one.
     */
    void copy(const AtomsConfiguration& other)
    {
        positions = other.positions;
        is_frozen = other.is_frozen;
        id = other.id;
        atoms_froz_active = other.atoms_froz_active;
        atoms_froz_inactive = other.atoms_froz_inactive;
        atoms_mov = other.atoms_mov;
        pairtype = other.pairtype;
        n_pt = other.n_pt;
    }

    /**
     * @brief Return number of moving atoms.
     */
    Index_t countMovingAtoms()
    {
        return (Index_t)std::count(is_frozen.getInternalVector().begin(),
                                   is_frozen.getInternalVector().end(),
                                   MOVING_ATOM);
    }
};

/**
 * @brief Structure with information for the LBFGS method.
 */
struct LBFGSInfo {
    Field<double> deltaOrient;  //!< Change of orientation in m previous
                                //!< iterations
    Field<double> deltaF;       //!< Change of rotational force in m
                                //!< previous iterations excluding
                                //!< the last one

    Coord F_old;  //!< Force from the previous iteration

    Index_t num_cg_iter;     //!< Maximum number of conjugate gradient
                             //!< iterations before resetting the
                             //!< conjugate directions
    Index_t num_lbfgs_iter;  //!< Maximum number of previous LBFGS
                             //!< iterations kept in memory

    void clear()
    {
        deltaOrient.clear();
        deltaF.clear();
        F_old.clear();
        num_cg_iter = 0;
        num_lbfgs_iter = 0;
    }

    // For debugging
    void print()
    {
        io::LogManager log_man;

        log_man << "deltaOrient: \n";
        deltaOrient.print();

        log_man << "deltaF: \n";
        deltaF.print();

        log_man << "F_old: \n";
        F_old.print();

        log_man << "num_cg_iter: " << num_cg_iter << "\n";

        log_man << "num_lbfgs_iter: " << num_lbfgs_iter << "\n";
    }
};

/**
 * @brief Structure of translation parameters
 */
struct TransitionParameters {
    /**
     * Basic step length.
     */
    double step_length;

    /**
     * Maximum permissible step length.
     */
    double max_step_length;

    /**
     * Threshold at which to apply the rotation removal
     */
    double rotrem_thresh;
};

/**
 * @brief Stores stop criteria for dimer method (e.g. LBFGS).
 */
struct StopCriteriaForDimer {
    /**
     * @brief Final convergence threshold.
     */
    double force = 0.;

    /**
     * @brief Convergence threshold for rotation angle.
     */
    // used to be 2, also: T_anglerot_gp
    double angle_rotation = 0.;

    void clear()
    {
        force = angle_rotation = 0.;
    }
};

/**
 * @brief Observation data: energy, gradient and positions.
 */
struct Observation {
    Coord R;
    Field<double> E;
    Coord G;

    void clear()
    {
        R.clear();
        E.clear();
        G.clear();
    }

    /**
     * @brief Appends elements from \e other to correspoinding elemnets of \e
     * this.
     */
    void append(const Observation& other)
    {
        R.append(other.R);
        E.append(other.E);
        G.append(other.G);
    }

    // for debugging
    void printSizes()
    {
        io::LogManager log_man;

        log_man << "R size: " << R.getSize() << "\n";
        log_man << "R (_ni, _nj): " << R.getNumRows() << " " << R.getNumCols()
                << "\n\n";

        log_man << "E size: " << E.getSize() << "\n";
        log_man << "E (_ni, _nj): " << E.getNumRows() << " " << E.getNumCols()
                << "\n\n";

        log_man << "G size: " << G.getSize() << "\n";
        log_man << "G (_ni, _nj): " << G.getNumRows() << " " << G.getNumCols()
                << "\n\n";
    }
};

/**
 * @brief A structure with information used to write the energy surface to a
 * file.
 */
struct EnergySurfaceOutputInfo {
    std::string out_dir = "output";
    std::string file_name_R = "position";
    std::string file_name_E = "energy";
    std::string file_name_G = "gradient";
    std::string file_extension = "dat";
    double offset = 3.;
    double dy = 0.1;
    double dz = 0.1;
    uint8_t debug_level = DEBUG_L0;
};

} /* namespace gpr */

#endif /* STRUCTURES_STRUCTURES_H_ */
