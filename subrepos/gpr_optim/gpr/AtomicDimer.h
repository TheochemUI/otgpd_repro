/*
 * GPR.h
 *
 *  Created on: 23 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_GPR_H_
#define GPR_GPR_H_

#include <Eigen/Dense>
#include <chrono>
#include <cstdint>

#include "../data_types/Field.h"
#include "../managers/io/LogManager.h"
#include "../structures/Structures.h"
#include "auxiliary/ProblemSetUp.h"
#include "covariance_functions/ConstantCF.h"
#include "covariance_functions/SexpatCF.h"
#include "gpr/Enums.h"
#include "ml/GaussianProcessRegression.h"

namespace atmd {

/**
 * @brief Atomic dimer
 */
class AtomicDimer {
public:
    /**
     * @brief Default constructor.
     */
    AtomicDimer();
    virtual ~AtomicDimer();

    /* *********************************** */
    /* *** Initializers and executers **** */
    /* *********************************** */
    /**
     * @brief Initilalize with initial data.
     *
     * @param parameters Structure of parameters from the input file
     * @param initial_data Initial observation (energy, gradients, positions)
     * @param middle_point Initial observation at the middle point (energy,
     * gradients, positions)
     * @param initial_orientation Initial orientation of the dimer
     * @param init_atom_config Structure with initial configuration of atoms
     */
    void initialize(const gpr::InputParameters& parameters,
                    const gpr::Observation& initial_data,
                    const gpr::Observation& middle_point,
                    const gpr::Coord& initial_orientation,
                    const gpr::AtomsConfiguration& init_atom_config);

    /**
     * @brief Set configuration of the atoms.
     *
     * @param _atoms_config A reference configuration
     */
    inline void setAtomsConfiguration(
        const gpr::AtomsConfiguration& _atoms_config);

    /**
     * @brief This function uses the atomic GP-dimer method to converge to a
     * saddle point, starting from somewhere inside the convergence area.
     *
     * The relaxation of the dimer on the approximated energy surface
     * is done according to a dimer method, where a rotation step rotates
     * the dimer (a pair of images) towards its minimum energy orientation
     * to find the lowest curvature mode of the potential energy and
     * translation step moves the dimer towards the saddle point by
     * inverting the force component in the direction of the dimer.
     * After each relaxation phase, the energy and gradient are acquired at
     * the middle point of the dimer, and the GP hyperparameters are
     * reoptimized.
     *
     * @param general_potential Reference to a class of general potential or to
     * the object of \e Matter class (in a case of coupling with EON)
     */
    template <typename Pot>
    void execute(Pot& general_potential);

    /* *********************************** */
    /* ************* Getters ************* */
    /* *********************************** */
    /**
     * @brief Return final force at the middle point.
     */
    inline Eigen::VectorXd getFinalForceAtMidPoint();

    /**
     * @brief Return final coordinates of the middle point.
     */
    inline Eigen::VectorXd getFinalCoordOfMidPoint();

    /**
     * @brief Return the total number of force calls
     */
    inline gpr::Index_t getTotalForceCalls();

    /**
     * @brief Return the total number of force calls
     */
    inline gpr::Index_t getTotalGPRForceCalls();

    /**
     * @brief Return the total number of cycles
     *
     * This is essentially the number of relaxation phases
     */
    inline gpr::Index_t getIterations();

    /**
     * @brief Return final energy.
     */
    inline double getFinalEnergy();

    /**
     * @brief Return final curvature.
     */
    inline double getFinalCurvature();

    /**
     * @brief Return final orientation of the dimer.
     */
    inline gpr::Coord* getFinalOrientation();

    /**
     * @brief Return gpr::Coordinates of all observed points.
     */
    inline gpr::Coord& getMidPointCoordForAllObservationPointsRef();

    /**
     * @brief Return energies for all observed points.
     */
    inline gpr::Field<double>& getEnergyForAllObservationPointsRef();

    /**
     * @brief Return gradients for all observed points.
     */
    inline gpr::Coord& getGradientForAllObservationPointsRef();

    /**
     * @brief Get GPR model.
     */
    inline gpr::GaussianProcessRegression* getGPRModel();

    /**
     * @brief Call for a potential defined in a \e Matter object.
     *
     * @param cell_dimensions Dimensions of the cell in EON 9 element format
     * @param mid_point The point on which to obtain the potential
     * @param T The potential object pointer
     */
    template <typename T>
    void callGeneralPotentialFromEON(
        const gpr::vector3_reg (&cell_dimensions)[3],
        gpr::Observation& mid_point, T& potential);

    /**
     * @brief Final convergence test for force.
     *
     * @detail Function find maximum absolute component of the force acting on
     * the middle point of the dimer and returns \e true if the value is less
     * than \e T_dimer.
     *
     * @param middle_point The current mid-point of the configuration
     * @param stop_criteria The final convergence threshold
     * @returns \e true if final convergence is reached
     */
    bool isFinalConvergenceReached(const gpr::Observation& middle_point,
                                   const double stop_criteria);

    /**
     * @brief Rotational convergence test
     *
     * @param orient The orientation coordinates
     * @param stop_criteria The final convergence threshold
     * @returns \e true if rotational convergence is reached
     */
    bool isRotationalConvergenceReached(const gpr::Coord& orient,
                                        gpr::Observation& middle_point);

    /**
     * @brief Rotate the dimer.
     */
    bool rotateDimer(const double stop_criteria, gpr::Coord& orient,
                     gpr::LBFGSInfo& rot_info, gpr::Observation& middle_point,
                     gpr::Observation& approx_obs);

    /**
     * @brief Translate the dimer.
     */
    void translateDimer(const gpr::Observation& middle_point,
                        const gpr::Coord& orient,
                        const gpr::Observation& approx_obs,
                        const double curvature, gpr::LBFGSInfo& trans_info,
                        gpr::Coord& R_new);

private:
    /**
     * @brief Clears all internal data.
     */
    void clear();

    bool isAtomsConfigurationCorrect();

    /**
     * @brief Evaluate accurate energy and force at middle point of the dimer.
     */
    template <typename Pot>
    void evaluateAccurateEnergyAndForceAtMidPoint(
        gpr::Observation& middle_point, Pot& general_potential);

    /**
     * @brief Evaluate accurate energy and force at image 1 of the dimer.
     */
    template <typename Pot>
    void evaluateAccurateEnergyAndForce(gpr::Observation& obs,
                                        Pot& general_potential);

    void updateLocation(const gpr::Coord orient, const gpr::Coord& middle_point,
                        gpr::Coord& new_point);

    /**
     * @brief Update active frozen atoms.
     * @return True if new atoms are activated
     */
    bool updateActiveFrozenAtoms(const gpr::Coord& R_new,
                                 gpr::LBFGSInfo& rot_info,
                                 gpr::LBFGSInfo& trans_info);

    /**
     * @brief Limit translation of the atom.
     *
     * This function limits the move if any atom-wise step length is larger than
     * 99% of 0.5*(1-'ratio_at_limit') times the minimum inter-atomic distance.
     */
    bool limitTranslation(const gpr::Coord& R, gpr::LBFGSInfo& trans_info,
                          gpr::Coord& R_new);

    /**
     * @brief Return \e true if the inter atomic distance is too big.
     *
     * @param R_new New position of the middle point
     * @param R_all Old position of the middle point
     * @param num_es1 Number of estimations (not used)
     */
    bool isInterAtomicDistanceTooBig(const gpr::Coord& R_new,
                                     const gpr::Coord& R_all,
                                     gpr::Index_t& num_es1);

    /**
     * @brief Return total number of evaluations passed.
     */
    gpr::Index_t getNumEvaluations();

    /**
     * @brief Check if the iteration limit reached and print warning.
     * @param iter_counter Current iteration counter
     * @param iter_limit Iteration limit
     * @param message Message to be printed
     */
    void checkIfMaxIterReachedAndPrintWarning(const gpr::Index_t iter_counter,
                                              const gpr::Index_t iter_limit,
                                              const std::string& message);

    /**
     * @brief Check if the rotational convergence reached.
     *
     * @param orient Orientation of the dimer at the current step
     * @param orient_old Orientation of the dimer at the previous step
     * @param stop_criteria Stop criteria
     * @return True if the convergence reached, false otherwise
     */
    bool checkRotationalConvergence(const gpr::Coord& orient,
                                    const gpr::Coord& orient_old,
                                    const double stop_criteria);
    /**
     * @brief Define the start path for the relaxation phase
     */
    void defineStartPath(const gpr::Coord& R_latest_conv,
                         const gpr::Coord& orient_latest_conv,
                         const gpr::Coord& orient_init_gpr,
                         gpr::Coord& R_previous, gpr::Coord& orient_previous,
                         gpr::Observation& middle_point, gpr::Coord& orient);

    /**
     * @brief Perfom initial rotations of a dimer.
     *
     * This function corresponds to step 5 from section 3.2: "Minimum mode
     * saddle point searches using Gaussian process regression with
     * inverse-distance covariance function", Olli-Pekka Koistinen, Vilhjalmur
     * Asgeirsson, Aki Vehtari, Hannes J ́onsson
     */
    template <typename Pot>
    void performInitialRotations(gpr::Observation& middle_point,
                                 gpr::Coord& orient, Pot& potential);

    /**
     * @brief Perfom GPR iterations.
     *
     * This function corresponds to step 6 from section 3.2: "Minimum mode
     * saddle point searches using Gaussian process regression with
     * inverse-distance covariance function", Olli-Pekka Koistinen, Vilhjalmur
     * Asgeirsson, Aki Vehtari, Hannes J ́onsson
     */
    template <typename Pot>
    std::pair<double, gpr::Coord> performGPRIterations(
        gpr::Observation& middle_point, gpr::Coord& orient, Pot& potential);

    // Returns curvature
    double performRotationAndTranslationWithGPR(
        const double stop_dimer_gp, gpr::Observation& middle_point,
        gpr::Coord& orient, gpr::Coord& R_latest_conv, gpr::Coord& R_previous,
        gpr::Coord& orient_latest_conv, gpr::Coord& orient_previous,
        gpr::Coord& R_new, gpr::Index_t& num_esmax, gpr::Index_t& num_es1,
        const size_t& outer_iter_idx);

    // Returns curvature
    double performRotationAndTranslationPure(gpr::Observation& middle_point,
                                             gpr::Observation& approx_obs,
                                             gpr::Coord& orient,
                                             gpr::Coord& R_new,
                                             gpr::LBFGSInfo& rot_info,
                                             gpr::LBFGSInfo& trans_info);

    inline double calculateCurvature(const gpr::Coord& orient,
                                     const gpr::Observation& approx_obs);

    inline void assignFinalResults(
        const gpr::Observation& middle_point,
        const std::pair<double, gpr::Coord>& eigen_data);

    void writeLatticeData(gpr::Observation& middle_point);

#ifdef WITH_HDF5
    void writeHDF5(const gpr::Observation& middle_point,
                   const std::string& addl_group, const size_t& idx);
#endif

    void pruneObservedData();

private:
    /**
     * Configuration of atoms
     */
    gpr::AtomsConfiguration atoms_config;

    /**
     * Dimensions of the computational box.
     */
    gpr::vector3_reg cell_dimensions[3];

    /**
     * Unit vector along the initial direction
     */
    gpr::Coord orient_init;

    /**
     * Initila coordinates, energy and gradient at the middle point
     */
    gpr::Observation middle_point_init;

    /**
     * Coordinates, energy and gradient at the image 1
     */
    gpr::Observation image1;

    /**
     * Coordinates, energies and gradients of all observation points (N_obs x D)
     */
    gpr::Observation all_obs;

    /**
     * Final results
     */
    gpr::Observation middle_point_final;
    double curvature_final;
    gpr::Coord orient_final;

    /**
     * @brief Reference energy level
     */
    // FIXME: the reference energy level should be a scalar, not field!
    gpr::Field<double> E_all_init;

    /**
     * Zero level of biased potential
     */
    gpr::Field<double> E_zero_level;

    /**
     * Object of GPR model
     */
    gpr::GaussianProcessRegression* gpr_model;

    gpr::StopCriteriaForDimer stop_citeria_dimer;
    gpr::StopCriteriaForDimer stop_citeria_gpr;
    double divisor_stop_criteria_gpr;

    /**
     * Dimer separation (distance from the middle point of the dimer to the two
     * images)
     */
    double dimer_sep;

    /**
     * Maximum number of initial rotations (0 if initial rotations skipped)
     */
    gpr::IterationsLimits max_iter_init_rot;

    gpr::IterationsLimits max_iter_new_pairs;

    gpr::IterationsLimits max_iter_relax_rot;

    bool assume_many_iterations;

    /**
     * Limit for the ratio (< 1) of inter-atomic distances between image and its
     * "nearest" observed data point (the relaxation phase is stopped if
     * 'ratio_at_limit' is reached for any image)
     */
    double ratio_at_limit;

    /**
     * Activation distance for moving+frozen atom pairs (inf if all active)
     */
    double actdist_fro;

    /**
     * Number of outer iterations where the hyperparameter optimization is
     * started from values initialized based on the range of the current data
     * (after that, the optimization is started from the values of the previous
     * round)
     */
    gpr::Index_t num_bigiter_initparam;

    /**
     * Type of the rotation method
     */
    uint8_t method_rot;

    /**
     * Type of the translation method
     */
    uint8_t method_trans;

    /**
     * Early stopping criteria
     */
    struct early_stopping_t {
        /**
         * @brief The distance metric (e.g., EMD, RMSD) used to quantify the
         * change between two atomic configurations.
         */
        DistanceMetricType _dist_metric;

        /**
         * @brief A fixed distance threshold. This value is used directly when
         * the adaptive threshold mechanism is disabled.
         */
        double _threshold;

        /**
         * @brief A master switch to enable the size-dependent adaptive
         * threshold. If true, the threshold is calculated dynamically using
         * `_adaptive_A` and `_adaptive_floor`. If false, the fixed `_threshold`
         * is used.
         */
        bool _use_adaptive_threshold;

        /**
         * @brief The permissiveness constant for the adaptive threshold.
         *
         * This value controls the overall scale of the threshold, especially
         * for small-to-medium systems, following the relation `thresh ≈ A /
         * sqrt(N)`. A larger 'A' results in a looser (more permissive)
         * threshold. A good starting value for exploratory searches is ~1.3,
         * which yields a threshold of ~0.4 Å for a 10-atom system.
         */
        double _adaptive_A;

        /**
         * @brief The floor value, or minimum threshold, for the adaptive
         * mechanism.
         *
         * This is the most critical parameter for tailoring the check to the
         * specific scientific goal. It prevents the threshold from becoming
         * overly restrictive ("crippling") for very large systems.
         *
         */
        double _adaptive_floor;

    } early_params;

    /**
     * Parameters of the translation method.
     */
    gpr::TransitionParameters transition_param;

    /**
     * Parameters for Prune
     */
    bool use_prune;
    gpr::Index_t start_prune_at;
    gpr::Index_t nprune_vals;
    double prune_threshold;

    /**
     * Object for logging
     */
    gpr::io::LogManager log_man;

    /**
     * Just a counter to be returned to EON
     */
    gpr::Index_t totalIterations;

    gpr::Index_t num_of_gen_potential_calls;
    gpr::Index_t num_of_gpr_potential_calls;

    // Auxiliary counter
    gpr::Index_t file_counter = 0;
    gpr::EnergySurfaceOutputInfo debug_output_info;
    gpr::Index_t debug_level;
};

} /* namespace atmd */

#include "AtomicDimer.inl"

#endif /* GPR_GPR_H_ */
