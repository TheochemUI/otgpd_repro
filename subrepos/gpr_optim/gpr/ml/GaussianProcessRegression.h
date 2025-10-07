/*
 * GaussianProcessRegression.h
 *
 *  Created on: 3 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_ML_GAUSSIANPROCESSREGRESSION_H_
#define GPR_ML_GAUSSIANPROCESSREGRESSION_H_

#include <Eigen/Dense>
#include <cmath>
#include <deque>

#include "../../data_types/Field.h"
#include "../../managers/io/LogManager.h"
#include "../../structures/Structures.h"
#include "../covariance_functions/ConstantCF.h"
#include "../covariance_functions/SexpatCF.h"
#include "../observation_models/LikGaussian.h"
#include "gpr/Enums.h"
#include "gpr/ml/ParamOpt.h"

namespace gpr {

/**
 * @brief Gaussian Process Regression algorithm.
 */
class GaussianProcessRegression {
public:
    GaussianProcessRegression();
    virtual ~GaussianProcessRegression();

    /**
     * @brief Initialize GPR class.
     *
     * @param conf_info Information about the configurations necessary for the
     * GPR model
     * @param parameters Structure of parameters
     */
    void initialize(const InputParameters& parameters,
                    const AtomsConfiguration& conf_info);

    /**
     * @brief Set up hyperparameters for GPR.
     * @param all_obs Structure of all observations
     * @param conf_info Structure of atoms configuration
     * @param update_sexpat_cf_param Set to "true" to update parameters of the
     * SexpAt covariance function
     * @param update_const_cf_param Set to "true" to update parameters of the
     * Constant covariance function
     * @param update_sqrt_prior_param Set to "true" to update parameters of the
     * SQRT prior
     */
    void setHyperparameters(const Observation& all_obs,
                            const AtomsConfiguration& conf_info,
                            const bool update_sexpat_cf_param = true,
                            const bool update_const_cf_param = true,
                            const bool update_sqrt_prior_param = true);

    /**
     * @brief Set up default parameters.
     *
     * This function allocates memory for the private pointers and calls for \b
     * setUpDefault() method on each pointer.
     */
    void setUpDefault();

    /**
     * @brief Evaluate training covariance matrix (gp_cov + noise covariance).
     *
     * The function takes in matrix \e x that contains training input vectors to
     * GPR. Returns the noisy covariance matrix \e cov_matrix for observations
     * y, which is sum of noisless covariance matrix and diagonal term, for
     * example, from Gaussian noise.
     *
     * Every element ij of noisless covariance matrix contains covariance
     * between inputs i and j in \e x. The noisless covariance matrix is formed
     * with all covariance functions.
     *
     * This is essentially the application of the kernel function to the
     * training points
     *
     * @param x Training input vectors to GP
     * @param cov_matrix Noisy covariance matrix for observations y
     */
    // used to be trcov
    void evaluateTrainingCovarianceMatrix(const EigenMatrix& x,
                                          const Eigen::VectorXd& x_ind,
                                          EigenMatrix& cov_matrix);

    /**
     * @brief Evaluate covariance matrix between two input vectors.
     *
     * Function takes in Gaussian process GP and two matrixes TX and X that
     * contain input vectors to GP. Returns covariance matrix C. Every element
     * ij of C contains covariance between inputs i in TX and j in X. PREDCF is
     * an optional array specifying the indexes of covariance functions, which
     * are used for forming the matrix. If empty or not given, the matrix is
     * formed with all functions.
     *
     * @param x1 First input vector
     * @param x2 Second input vector
     * @param C Covariance matrix
     */
    void evaluateCovarianceMatrix(const EigenMatrix& x1, const EigenMatrix& x2,
                                  const Eigen::VectorXd& x1_ind,
                                  const Eigen::VectorXd& x2_ind,
                                  EigenMatrix& C);

    /**
     * @brief Extract unique indices from the input vector.
     */
    void extractUniqueIndices(const Eigen::VectorXd& input_ind,
                              Eigen::VectorXd& output_ind);

    /**
     * @brief Evaluate the energy function (un-normalized negative marginal log
     * posterior) and its gradient.
     *
     * This function takes a Gaussian process structure GP together with a
     * matrix X of input vectors and a matrix Y of targets, and evaluates the
     * energy function E and its gradient G. Each row of X corresponds to one
     * input vector and each row of Y corresponds to one target vector.
     *
     * The energy is minus log posterior cost function:
     *  E = EDATA + EPRIOR = - log p(Y|X, th) - log p(th),
     * where th represents the parameters (lengthScale, magnSigma2...), X is
     * inputs and Y is observations (regression) or latent values (non-Gaussian
     * likelihood).
     * @param w Vector of combined parameters from covariance functions and
     * likelyhood
     */
    // i.e. gp_eg
    void evaluateEnergyAndGradient(const Eigen::VectorXd& w,
                                   const EigenMatrix& x,
                                   const Eigen::VectorXd& x_ind,
                                   const Eigen::VectorXd& y,
                                   EnergyAndGradient& energy_and_gradient);

    /**
     * @brief Optimize paramaters of a Gaussian process
     *
     * This function optimizes the parameters of a GP structure given matrix X
     * of training inputs and vector Y of training targets.
     */
    void optimize(const Observation& observation);

    /**
     * @brief Calculate potential within GPR model.
     *
     * This auxiliary function calculates energy and gradient predictions for
     * given images \e R according to a GPR model.
     *
     * @param image1 Observation at image 1
     * @param num_of_calls Number of calls for GPR potential
     */
    void calculatePotential(Observation& image1);

    /**
     * @brief Set sigma2 jitter value.
     */
    inline void setJitterSigma2(const Index_t value);

    /**
     * @brief Set parameters of the GPR model.
     */
    inline void setParameters(const GPRSetup& parameters);

    /**
     * @brief Get pointer to Gaussian liklehood.
     */
    inline LikGaussian* getLikGaussian();

    /**
     * @brief Get pointer to SexpAt covariance function.
     */
    inline SexpatCF* getSexpAtCovarianceFunction();

    /**
     * @brief Get pointer to constant covariance function.
     */
    inline ConstantCF* getConstantCovarianceFunction();

    /**
     * @brief Returns number of GPR potential calls.
     */
    inline Index_t getNumberOfPotentialCalls();

protected:
    /**
     * @brief Calculate posterior mean prediction w.r.t energy.
     * This function mimics the following step from the MATLAB code:
     *   a = L'\(L\[E_all;G_all(:)]);
     */
    void calculatePosteriorMeanPrediction();

    /**
     * @brief Posterior mean prediction w.r.t energy.
     * This function mimics the following step from the MATLAB code:
     *   b = L\[E_all;G_all(:)];
     */
    void calculateMeanPrediction(const Eigen::VectorXd& y);

    /**
     * @brief Evaluate the gradient of energy (GP_E) for Gaussian Process.
     *
     * This function takes a full GP parameter vector W, GP structure GP, a
     * matrix X of input vectors and a matrix Y of target vectors, and evaluates
     * the gradient G of the energy function (gp_e). Each row of X corresponds
     * to one input vector and each row of Y corresponds to one target vector.
     *
     * @param x Matrix of coordinates
     * @param y Vector of energy and gradients
     * @param gradient Calculated new gradient
     */
    // i.e. gp_g
    void evaluateGradient(const EigenMatrix& x, const Eigen::VectorXd& x_ind,
                          const Eigen::VectorXd& y, Eigen::VectorXd& gradient);

    /**
     * @brief Evaluate the energy function (un-normalized negative log marginal
     * posterior)
     *
     * This function takes a Gaussian process structure GP together with a
     * matrix X of input vectors and a matrix Y of targets, and evaluates the
     * energy function E. Each row of X corresponds to one input vector and each
     * row of Y corresponds to one target vector.
     *
     * The energy is minus log posterior cost function:
     *  E = EDATA + EPRIOR = - log p(Y|X, th) - log p(th),
     * where th represents the parameters (lengthScale, magnSigma2...), X is
     * inputs and Y is observations (regression) or latent values (non-Gaussian
     * likelihood).
     *
     * Note that this implements Algo 2.1 of Rasmussen and Williams but the
     * variance is not computed, and the mean is not returned, so only the log
     * marginal liklihood is caclulated and used directly to predict.
     *
     * OPTIONS is optional parameter-value pair
     *  z - optional observed quantity in triplet (x_i,y_i,z_i)
     *  Some likelihoods may use this. For example, in case of Poisson
     * likelihood we have z_i=E_i, that is, expected value for ith case.
     *
     * @param x Matrix of coordinates
     * @param y Vector of energy and gradients
     * @return Energy Calculated new energy
     */
    // i.e. gp_e
    // FIXME: Consider breaking into a train and predict phase instead
    double evaluateEnergy(const EigenMatrix& x, const Eigen::VectorXd& x_ind,
                          const Eigen::VectorXd& y);

    /**
     * @brief Evaluate and decompose the covariance matrix using Cholesky
     * decomposition.
     */
    bool decomposeCovarianceMatrix(const EigenMatrix& x,
                                   const Eigen::VectorXd& x_ind);

private:
    /**
     * @brief Set GP parameters.
     * This function sets parameters in the covariance functions \e sexpat and
     * \e constant and in the Gaussian likelyhood \e lik_gaussian.
     * @param w Vector of combined parameters from covariance functions and
     * likelyhood
     */
    void setParameters(const Eigen::VectorXd& w);

    /**
     * @brief Calculate gradient with respect to covariance function.
     */
    template <typename CovFunc>
    void calculateGradientWithCovFunc(const EigenMatrix& x,
                                      const Eigen::VectorXd& ind_Ddim,
                                      const Eigen::VectorXd& uDdim,
                                      const Eigen::VectorXd& b,
                                      const EigenMatrix& invC,
                                      CovFunc& cov_func, Field<double>& gdata,
                                      Field<double>& gprior);

    /**
     *
     */
    template <typename CovFunc>
    void applyCovarianceFunction(const EigenMatrix& x,
                                 const Eigen::VectorXd& ind_Ddim,
                                 const Eigen::VectorXd& uDdim,
                                 CovFunc& sexpat_cov_func, EigenMatrix& Ktemp);

    // FIXME: rename!
    template <typename CovFunc>
    void evaluateCovarianceFunction(const EigenMatrix& x1,
                                    const EigenMatrix& x2,
                                    const Eigen::VectorXd& ind_Ddim,
                                    const Eigen::VectorXd& uDdim,
                                    const Eigen::VectorXd& ind_Ddim2,
                                    const Eigen::VectorXd& uDdim2,
                                    CovFunc& cov_func, EigenMatrix& Ktemp);

public:
    bool failedOptimizer;
    /**
     * @brief Extract coordinates from \e x using value \e ind from \e
     * ind_Ddim.
     */
    void extractCoordinatesByIndex(const EigenMatrix& x,
                                   const Eigen::VectorXd& ind_Ddim,
                                   const Index_t ind, Coord& x_loc);

    /**
     * @brief Assign \e field to \e matrix according to indices form \e
     * ind_Ddim1 and \e ind_Ddim2.
     *
     * This function performs the following operation (in MATLAB form):
     * \f[ matrix(ind_Ddim == row_val, ind_Ddim == col_val) = field; \f]
     *
     * \note \e matrix should be pre-allocated
     *
     * @param ind_Ddim1 Vector of dimensions used in rows
     * @param ind_Ddim2 Vector of dimensions used in columns
     * @param row_val Index of row
     * @param col_val Index of column
     * @param field Reference field
     * @param matrix Reference matrix
     * @param transpose_field True if the field should be transposed
     */
    void assignBlockToMatrix(const Eigen::VectorXd& ind_Ddim1,
                             const Eigen::VectorXd& ind_Ddim2,
                             const Index_t row_val, const Index_t col_val,
                             const Field<double>& field, EigenMatrix& matrix,
                             bool transpose_field = false);

    /**
     * @brief Remove column from the dense matrix.
     * @param column_ind Column index
     * @param matrix Dense matrix
     */
    void removeColumn(const Index_t column_ind, EigenMatrix& matrix);

    /**
     * @brief Count number of the first repetitive indices in \e ind_Ddim .
     */
    inline Index_t getNumberOfRepetitiveIndices(
        const Eigen::VectorXd& ind_Ddim);

    /**
     * TODO: document, move, get_variance.m
     * */
    void calculateVariance(Observation& image1);

    // NOTE(rg): should be private
    SexpatCF* sexpat_cov_func;
    // NOTE(rg): should be protected
    void updateModelWithFullData(const Observation& full_observation);
private:
    double sigma2;
    double jitter_sigma2;
    uint8_t optimization_alg;
    uint8_t report_level;

    LikGaussian* lik_gaussian;
    ConstantCF* const_cov_fun;

    io::LogManager log_man;

    funcmin::ParamOptimizer popt;
    OptimizationAlgorithmSettings opt_alg_settings;
    bool is_training_cov_matrix_evaluated;
    bool is_decomposed_succesfully;

    Index_t num_of_potential_calls;
    AtomsConfiguration atoms_config;
    struct farthest_point_sampling_t {
        int history;
        int latest_points;
        DistanceMetricType metric;
    } fps_options;

protected:
    EigenMatrix C;      // Auxiliary matrix
    EigenMatrix L;      // Low-triangle matrix of a decomposed covariance matrix
    Eigen::VectorXd b;  // Auxiliary vector
    Eigen::VectorXd a;  // Posterior mean prediction w.r.t energy
    EigenMatrix R_matrix;  // Matrix of coordinates assembled from gpr::Coord
    Eigen::VectorXd R_indices;  // Vector of indices that correspond to
                                // repetitive rows in R_matrix
    Eigen::VectorXd
        energy_and_gradient;  // Vector of combined energy and gradients
    // --- Members for Dynamic History Sizing ---
    std::deque<Eigen::VectorXd> hyperparameter_history_;

    // --- Tunable parameters for the heuristic ---
    const bool use_hod = true; // hyperparameter oscillation detection
    // How many past steps to look at.
    const int monitoring_window_ = 5;
    // If >this% of updates are oscillations, increase history.
    const double flip_threshold_ = 0.8;
    // How much to increase history by.
    const int history_increment_ = 2;
    // A sane upper limit to prevent performance issues.
    const int max_history_ = 30;
    Observation selectOptimizationSubset(
        const Observation& full_observation) const;
    Eigen::VectorXd performHyperparameterOptimization(
        const Observation& opt_observation);
    bool monitorAndAdaptHyperparameters(const Eigen::VectorXd& new_parameters);
};

} /* namespace gpr */

#include "GaussianProcessRegression.inl"

#endif /* GPR_ML_GAUSSIANPROCESSREGRESSION_H_ */
