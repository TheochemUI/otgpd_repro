/*
 * GaussianProcessRegression.cpp
 *
 *  Created on: 3 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "GaussianProcessRegression.h"

#include "gpr/Enums.h"
#include "gpr/auxiliary/AdditionalFunctionality.h"

#define EIGEN_NO_DEBUG
#include <Eigen/Core>
#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>

#include "../../backend/DistibutionFunctions.h"
#include "../auxiliary/Distance.h"

namespace gpr {

GaussianProcessRegression::GaussianProcessRegression()
{
    sigma2 = 1e-7;
    jitter_sigma2 = 1e-6;
    optimization_alg = SCG_opt;

    lik_gaussian = new LikGaussian;
    const_cov_fun = new ConstantCF;
    sexpat_cov_func = new SexpatCF;

    is_training_cov_matrix_evaluated = false;
    is_decomposed_succesfully = false;

    num_of_potential_calls = 0;
    failedOptimizer = false;
    fps_options.history = 5;
    fps_options.latest_points = 2;
    fps_options.metric = DistanceMetricType::EMD;
}

GaussianProcessRegression::~GaussianProcessRegression()
{
    if (lik_gaussian != nullptr) {
        delete lik_gaussian;
        lik_gaussian = nullptr;
    }
    if (const_cov_fun != nullptr) {
        delete const_cov_fun;
        const_cov_fun = nullptr;
    }
    if (sexpat_cov_func != nullptr) {
        delete sexpat_cov_func;
        sexpat_cov_func = nullptr;
    }
}

void GaussianProcessRegression::initialize(const InputParameters& parameters,
                                           const AtomsConfiguration& conf_info)
{
    this->atoms_config = conf_info;
    fps_options.history = parameters.fps_history.value;
    fps_options.metric = parameters.fps_metric.value;
    PriorBase prior_parameters;

    sigma2 = parameters.gp_sigma2.value;
    jitter_sigma2 = parameters.jitter_sigma2.value;
    report_level = parameters.report_level.value;

    if (parameters.optimization_alg.value == "SCG_opt") {
        popt.emplace<funcmin::SCG>();
    } else if (parameters.optimization_alg.value == "ADAM_opt") {
        popt.emplace<funcmin::ADAM>();
    } else if (parameters.optimization_alg.value == "Ceres_opt") {
        popt.emplace<funcmin::Ceres>();
    } else {
        popt.emplace<funcmin::SCG>();
    }

    // Initialize Gaussian likelihood
    prior_parameters.setMu(parameters.prior_mu.value);
    prior_parameters.setNu(parameters.prior_nu.value);
    prior_parameters.setS2(parameters.prior_s2.value);
    lik_gaussian->setSigma2(parameters.sigma2.value);

    // Initialize constant covariance function
    //    Field<double> lengthScale; <----- ? is it calculated or is in initial
    //    data?

    sexpat_cov_func->setMagnSigma2(parameters.magnSigma2.value);
    sexpat_cov_func->setConfInfo(conf_info);
    sexpat_cov_func->setPriorParametersSqrtt(prior_parameters);
    sexpat_cov_func->setPriorParametersGaussian(prior_parameters);

    // Initialize constant covariance function
    const_cov_fun->setConstSigma2(parameters.constSigma2.value);

    if (parameters.check_derivative.value == "true")
        opt_alg_settings.check_derivative = true;
    else
        opt_alg_settings.check_derivative = false;
    opt_alg_settings.report_level = parameters.report_level.value;
    opt_alg_settings.max_iter = parameters.max_iter.value;
    opt_alg_settings.tolerance_func = parameters.tolerance_func.value;
    opt_alg_settings.tolerance_sol = parameters.tolerance_sol.value;
    opt_alg_settings.lambda_limit = parameters.lambda_limit.value;
    opt_alg_settings.lambda = parameters.lambda.value;
}

// used to be update()
// calculation for R_all2 should be moved out
void GaussianProcessRegression::setHyperparameters(
    const Observation& all_obs, const AtomsConfiguration& conf_info,
    const bool update_sexpat_cf_param, const bool update_const_cf_param,
    const bool update_sqrt_prior_param)
{
    aux::Distance distance;
    Field<double> dist;
    Field<double> dummy(1, 1);
    double mean_y;
    double range_x;
    double range_y;
    double norm_inv;
    PriorBase prior_parameters;

    math::DistributionFunctions distr_func;
    aux::AuxiliaryFunctionality aux_func;
    Coord R_all2;

    mean_y = all_obs.E.getMean();
    range_y = all_obs.E.getMaxElt() - all_obs.E.getMinElt();

    // range_x = max(max(dist_at(R_all,R_all,conf_info,1)));
    dummy.set(1);
    distance.dist_at(all_obs.R, all_obs.R, conf_info, dummy, dist);
    range_x = dist.getMaxElt();

#ifndef NDEBUG
    assertMsg(const_cov_fun != nullptr,
              "Object of constant covariance function is not allocated!");
    assertMsg(sexpat_cov_func != nullptr,
              "Object of sexpAt covariance function is not allocated!");
#endif

    if (update_const_cf_param)
        const_cov_fun->setConstSigma2(std::max(1., mean_y * mean_y));

    if (update_sexpat_cf_param) {
        if (update_const_cf_param) {
            norm_inv = distr_func.normalCDFInverse(0.75, 0., range_y / 3.);
            sexpat_cov_func->setMagnSigma2(norm_inv * norm_inv);
        }

        norm_inv = distr_func.normalCDFInverse(0.75, 0., range_x / 3.);
        aux_func.repmatConst(1, conf_info.n_pt, norm_inv,
                             sexpat_cov_func->getLengthScaleRef());
    }

    if (update_sqrt_prior_param) {
        prior_parameters = sexpat_cov_func->getPriorParametersSqrtt();
        prior_parameters.setS2(std::max(1., (range_y / 3.) * (range_y / 3.)));
        sexpat_cov_func->setPriorParametersSqrtt(prior_parameters);
    }

    prior_parameters = sexpat_cov_func->getPriorParametersGaussian();
    prior_parameters.setS2(std::max(1., (range_x / 3.) * (range_x / 3.)));
    sexpat_cov_func->setPriorParametersGaussian(prior_parameters);
}

void GaussianProcessRegression::calculateVariance(Observation& image1)
{
    // NOTE: C is a truncated covariance!
    Index_t n = image1.R.getNumRows();
    EigenMatrix R_mod;
    Eigen::VectorXd R_mod_ind;
    EigenMatrix cov_matrix, KK;
    Eigen::VectorXd V;
    EigenMatrix v;
    Eigen::VectorXd VarEG_R;
    Index_t counter = 0;
    aux::AuxiliaryFunctionality aux_func;

    // R2 = [repmat(R,D+1,1),reshape(repmat(0:D,n,1),[],1)];
    aux_func.assembleMatrixOfRepetitiveCoordinates(image1.R, R_mod, R_mod_ind);

    // KK =
    // gp_cov(gp,R_all2,[repmat(R,D+1,1),reshape(repmat(0:D,n,1),[],1)]);
    evaluateCovarianceMatrix(R_matrix, R_mod, R_indices, R_mod_ind, KK);

    // Cov = gp_trcov(gp,R2);
    evaluateTrainingCovarianceMatrix(R_mod, R_mod_ind, cov_matrix);

    // V = diag(Cov);
    V = cov_matrix.diagonal();

    // v = L\KK;
    // v = L.triangularView<Eigen::Upper>().solve(KK);
    v = L.inverse() * KK;

    // VarEG_R = V - sum(v'.*v',2);
    // VarEG_R.resize(v.rows());
    // std::cout << "\nv is " << v.rows() << " " << v.cols() << std::endl;
    // std::cout << "\nL is " << L.rows() << " " << L.cols() << std::endl;
    // std::cout << "\nKK is " << KK.rows() << " " << KK.cols() << std::endl;
    // std::cout << "\ncov_matrix is " << cov_matrix.rows() << " "
    //           << cov_matrix.cols() << std::endl;
    // VarEG_R = V - sum(v'.*v',2);
    VarEG_R.resize(V.size());
    for (Index_t i = 0; i < V.size(); ++i) {
        double tmp = 0.;
        for (Index_t j = 0; j < v.rows(); ++j) {
            tmp += v(j, i) * v(j, i);
        }
        VarEG_R(i) = V(i) - tmp;
    }
    // std::cout << VarEG_R;
    // VarE_R = VarEG_R(1:n,:);
    image1.E.resize(1, n);
    for (Index_t idx = 0; idx < n; ++idx) {
        image1.E[idx] = VarEG_R(idx);
    }

    // G_R = reshape(EG_R((n+1):end,1), n, D);
    // G_R = reshape(EG_R((n+1):end,:), n, size(EG_R,1)/n-1);
    counter = n;
    image1.G.resize(n, image1.R.getNumCols());
    for (Index_t j = 0; j < image1.G.getNumCols(); ++j) {
        for (Index_t i = 0; i < image1.G.getNumRows(); ++i) {
            image1.G(i, j) = VarEG_R(counter++);
        }
    }
}

void GaussianProcessRegression::evaluateTrainingCovarianceMatrix(
    const EigenMatrix& x, const Eigen::VectorXd& x_ind, EigenMatrix& cov_matrix)
{
    Index_t n = (Index_t)x.rows();
    Eigen::VectorXd uDdim;  // unique derivatives' dimensions
    EigenMatrix Ktemp;

    cov_matrix.resize(n, n);
    cov_matrix.setZero();

    extractUniqueIndices(x_ind, uDdim);

    // Apply all covariance functions
    applyCovarianceFunction(x, x_ind, uDdim, *const_cov_fun, Ktemp);
    cov_matrix = cov_matrix + Ktemp;

    applyCovarianceFunction(x, x_ind, uDdim, *sexpat_cov_func, Ktemp);
    cov_matrix = cov_matrix + Ktemp;

    // TODO: check
    // The following lines from the MATLAB code make no sense!
    //    n = size(K,1);
    //    n1 = n+1;
    //    K(1:n1:end)=K(1:n1:end) + gp.jitterSigma2;

    //    C = C + gp.lik.fh.trcov(gp.lik, x1);
    Field<double> C;
    Coord x_coord;
    x_coord.resize((Index_t)x.rows(), (Index_t)x.cols());
    for (Index_t i = 0; i < x_coord.getNumRows(); ++i) {
        for (Index_t j = 0; j < x_coord.getNumCols(); ++j) {
            x_coord(i, j) = x(i, j);
        }
    }

#ifndef NDEBUG
    assertMsg(lik_gaussian != nullptr,
              "Object of Gausian likelihood is not allocated!");
#endif
    lik_gaussian->evaluateTrainingCovarianceMatrix(x_coord, C);

    for (Index_t i = 0; i < C.getNumRows(); ++i) {
        for (Index_t j = 0; j < C.getNumCols(); ++j) {
            cov_matrix(i, j) += C(i, j);
        }
    }

    // log_man << "Training covariance matrix is " << cov_matrix.size() << "
    // over ( " << cov_matrix.rows() << ", " << cov_matrix.cols() << " )" <<
    // "\n";
    is_training_cov_matrix_evaluated = true;
}

void GaussianProcessRegression::evaluateCovarianceMatrix(
    const EigenMatrix& x1, const EigenMatrix& x2, const Eigen::VectorXd& x1_ind,
    const Eigen::VectorXd& x2_ind, EigenMatrix& C)
{
    Index_t n = (Index_t)x1.rows();
    Index_t n1 = (Index_t)x2.rows();
    EigenMatrix Ktemp;
    Eigen::VectorXd uDdim, uDdim2;

    C.resize(n, n1);

    C.setZero();
    Ktemp.setZero();

    extractUniqueIndices(x1_ind, uDdim);
    extractUniqueIndices(x2_ind, uDdim2);

    /* Evaluate all covariance matrices */
    evaluateCovarianceFunction(x1, x2, x1_ind, uDdim, x2_ind, uDdim2,
                               *sexpat_cov_func, Ktemp);
    C += Ktemp;

    evaluateCovarianceFunction(x1, x2, x1_ind, uDdim, x2_ind, uDdim2,
                               *const_cov_fun, Ktemp);
    C += Ktemp;
}

// void GaussianProcessRegression::extractCoordinatesByIndex(
//    const EigenMatrix& x, const Eigen::VectorXd& ind_Ddim,
//    const Index_t ind, Coord& x_loc)
//{
//    // We know that ind_Ddim is a vector of repetitive indices, starting from
//    0,
//    // so we can use this information to determine the start/end rows of x_loc
//    Index_t Ddim_rows = (Index_t)ind_Ddim.rows();
//    Index_t num_rep_Ddim = Ddim_rows / (ind_Ddim(Ddim_rows - 1) + 1);
//
//    Index_t i_start = ind * num_rep_Ddim;
//    Index_t i_end = (ind + 1) * num_rep_Ddim;
//
//    Index_t counter_x_loc = 0;
////    Index_t counter_x = i_start * (Index_t)x.cols();
//
//    for (Index_t i = i_start; i < i_end; ++i) {
//        // Here we assume that the Eigen::Matrix is stored in a row-major
//        // order
//        std::copy(x.data() + i * x.cols(), x.data() + (i + 1) * x.cols(),
//                  x_loc.getInternalVector().data() + counter_x_loc *
//                  x_loc.getNj());
//        x_loc(counter_x_loc, x_loc.getNj() - 1) = ind;
//        ++counter_x_loc;
////        for (Index_t j = 0; j < x.cols(); ++j) {
////            // Here we assume that the Eigen::Matrix is stored in a
/// row-major /            // order /            x_loc[counter_x_loc++] =
/// x.data()[counter_x++]; // instead of x(i, j); /        } /
/// x_loc[counter_x_loc++] = ind;
//    }
//}
//
void GaussianProcessRegression::assignBlockToMatrix(
    const Eigen::VectorXd& ind_Ddim1, const Eigen::VectorXd& ind_Ddim2,
    const Index_t row_val, const Index_t col_val, const Field<double>& field,
    EigenMatrix& matrix, bool transpose_field)
{
    // We know that ind_Ddim1 and ind_Ddim2 are vectors of repetitive indices,
    // starting from 0, so we can use this information to determine the
    // start/end rows and columns of the matrix
    Index_t Ddim1_rows = (Index_t)ind_Ddim1.rows();
    Index_t Ddim2_rows = (Index_t)ind_Ddim2.rows();
    Index_t num_rep_Ddim1 = Ddim1_rows / (ind_Ddim1(Ddim1_rows - 1) + 1);
    Index_t num_rep_Ddim2 = Ddim2_rows / (ind_Ddim2(Ddim2_rows - 1) + 1);

    Index_t i_start = row_val * num_rep_Ddim1;
    Index_t i_end = (row_val + 1) * num_rep_Ddim1;
    Index_t j_start = col_val * num_rep_Ddim2;
    Index_t j_end = (col_val + 1) * num_rep_Ddim2;
    Index_t i_conter = 0;

    if (!transpose_field) {
        for (Index_t i = i_start; i < i_end; ++i) {
            Index_t j_conter = 0;
            for (Index_t j = j_start; j < j_end; ++j) {
                matrix(i, j) = field(i_conter, j_conter);
                ++j_conter;
            }
            ++i_conter;
        }
    } else {
        for (Index_t i = i_start; i < i_end; ++i) {
            Index_t j_conter = 0;
            for (Index_t j = j_start; j < j_end; ++j) {
                matrix(i, j) = field(j_conter, i_conter);
                ++j_conter;
            }
            ++i_conter;
        }
    }
}

void GaussianProcessRegression::extractCoordinatesByIndex(
    const EigenMatrix& x, const Eigen::VectorXd& ind_Ddim, const Index_t ind,
    Coord& x_loc)
{
    Index_t row_counter = 0;

    for (Index_t i = 0; i < x.rows(); ++i) {
        if (ind_Ddim[i] == ind) {
            for (Index_t j = 0; j < x.cols(); ++j) {
                x_loc(row_counter, j) = x(i, j);
            }
            ++row_counter;
        }
    }
}

// void GaussianProcessRegression::assignBlockToMatrix(
//    const Eigen::VectorXd& ind_Ddim1, const Eigen::VectorXd& ind_Ddim2,
//    const Index_t row_val, const Index_t col_val, const Field<double>& field,
//    EigenMatrix& matrix, bool transpose_field)
//{
//    // FIXME: this is a forehead approach to find and replace a block in
//    // a matrix. Should be changed!
//    Index_t i_conter = 0;
//    for (Index_t i = 0; i < matrix.rows(); ++i) {
//        if (ind_Ddim1(i) == row_val) {
//            Index_t j_conter = 0;
//            for (Index_t j = 0; j < matrix.cols(); ++j) {
//                if (ind_Ddim2(j) == col_val) {
//                    if (!transpose_field)
//                        matrix(i, j) = field(i_conter, j_conter);
//                    else
//                        matrix(i, j) = field(j_conter, i_conter);
//                    ++j_conter;
//                }
//            }
//            ++i_conter;
//        }
//    }
//}

void GaussianProcessRegression::removeColumn(const Index_t column_ind,
                                             EigenMatrix& matrix)
{
    Index_t num_rows = (Index_t)matrix.rows();
    Index_t num_cols = (Index_t)matrix.cols() - 1;

    if (column_ind < num_cols)
        matrix.block(0, column_ind, num_rows, num_cols - column_ind) =
            matrix.block(0, column_ind + 1, num_rows, num_cols - column_ind);

    matrix.conservativeResize(num_rows, num_cols);
}

void GaussianProcessRegression::extractUniqueIndices(
    const Eigen::VectorXd& input_ind, Eigen::VectorXd& output_ind)
{
    // TODO: check if this algorithm can be optimized
    std::vector<double> unique_ind_loc;

    unique_ind_loc.insert(unique_ind_loc.begin(), input_ind.data(),
                          input_ind.data() + input_ind.size());

    // erase all non positive values
    // FIXME: optimize
    auto it = std::remove_if(unique_ind_loc.begin(), unique_ind_loc.end(),
                             [](const double value) { return value <= 0; });
    unique_ind_loc.erase(it, unique_ind_loc.end());

    // erase all duplicates and sort
    sort(unique_ind_loc.begin(), unique_ind_loc.end());
    unique_ind_loc.erase(unique(unique_ind_loc.begin(), unique_ind_loc.end()),
                         unique_ind_loc.end());

    Eigen::Map<Eigen::VectorXd> map(unique_ind_loc.data(),
                                    unique_ind_loc.size());
    output_ind = map;
}

void GaussianProcessRegression::evaluateEnergyAndGradient(
    const Eigen::VectorXd& w, const EigenMatrix& x,
    const Eigen::VectorXd& x_ind, const Eigen::VectorXd& y,
    EnergyAndGradient& energy_and_gradient)
{
    setParameters(w);

    *energy_and_gradient.energy = evaluateEnergy(x, x_ind, y);

    if (fabs(*energy_and_gradient.energy) <= DBL_EPSILON) {
        io::ErrorManager err;
        err << "Energy field is empty!"
            << "\n";
        energy_and_gradient.gradient->setZero();
    } else {
        evaluateGradient(x, x_ind, y, *energy_and_gradient.gradient);
    }
}

double GaussianProcessRegression::evaluateEnergy(const EigenMatrix& x,
                                                 const Eigen::VectorXd& x_ind,
                                                 const Eigen::VectorXd& y)
{
    Index_t n = (Index_t)x.rows();
    double zc = 0., edata = 0., eprior = 0.;
    double energy = 0.;

    is_decomposed_succesfully = decomposeCovarianceMatrix(x, x_ind);

    // Test if the matrix is positive definite
    if (!is_decomposed_succesfully) {
        log_man << "Warning! Matrix L is not positive definite"
                << "\n";
        energy = std::numeric_limits<double>::quiet_NaN();
        ;
    } else {
        calculateMeanPrediction(y);

        zc = 0.;
        for (Index_t n = 0; n < L.rows(); ++n) {
            zc += log(L(n, n));
        }
        edata = 0.5 * n * log(2 * M_PI) + zc +
                0.5 * b.dot(b);  // Line 7, zc should be negative

        // Evaluate the prior contribution to the error from covariance
        // functions
#ifndef NDEBUG
        assertMsg(const_cov_fun != nullptr,
                  "Object of Gausian likelihood is not allocated!");
        assertMsg(sexpat_cov_func != nullptr,
                  "Object of Gausian likelihood is not allocated!");
#endif

        eprior -= const_cov_fun->calculateLogPrior();

        eprior -= sexpat_cov_func->calculateLogPrior();

        // Evaluate the prior contribution to the error from Gaussian likelihood
        eprior -= lik_gaussian->evaluateLogPrior();

        energy = edata + eprior;
    }
    //    std::cout << "\n in evaluateEnergy(): \n" << y << "\n\n";
    return energy;
}

void GaussianProcessRegression::evaluateGradient(const EigenMatrix& x,
                                                 const Eigen::VectorXd& x_ind,
                                                 const Eigen::VectorXd& y,
                                                 Eigen::VectorXd& gradient)
{
    EigenMatrix invC;
    Eigen::VectorXd b;
    Eigen::VectorXd uDdim;
    Field<double> gdata;
    Field<double> gprior;

    // Note: the calculateTrainingCovarianceMatrix() method is always called in
    // the evaluateEnergy() method, and we always call for the
    // evaluateGradient() method after evaluateEnergy() being called (see
    // evaluateEnergyAndGradient(). So, there is no point to re-evaluate
    // covariance metrix here. However, we keep it hear guarded by an
    // if-statement for sake of unit tests. Otherwise, the matrix C will not be
    // allocated.
    if (!is_training_cov_matrix_evaluated) {
        evaluateTrainingCovarianceMatrix(x, x_ind, C);
    }

    invC = C.inverse();

    b = C.lu().solve(y);

    // Gradient with respect to covariance function parameters
    extractUniqueIndices(x_ind, uDdim);

    calculateGradientWithCovFunc(x, x_ind, uDdim, b, invC, *const_cov_fun,
                                 gdata, gprior);

    calculateGradientWithCovFunc(x, x_ind, uDdim, b, invC, *sexpat_cov_func,
                                 gdata, gprior);

    // Gradient with respect to Gaussian likelihood function parameters
    // Evaluate the gradient from Gaussian likelihood
    Field<double> gprior_lik;

#ifndef NDEBUG
    assertMsg(lik_gaussian != nullptr,
              "Object of Gausian likelihood is not allocated!");
#endif
    lik_gaussian->evaluateLogPriorGradient(gprior_lik);

    if (!gprior_lik.isEmpty()) {
        gprior_lik *= -1.;
        gprior.append(gprior_lik);
    }

    if (!gprior.isEmpty()) {
        if (gradient.rows() != gprior.getSize())
            gradient.resize(gprior.getSize());
        for (Index_t n = 0; n < gradient.rows(); ++n)
            gradient[n] = gdata[n] + gprior[n];
    }
}

void GaussianProcessRegression::setParameters(const Eigen::VectorXd& w)
{
    if (w.rows() == 0) return;

#ifndef NDEBUG
    assertMsg(sexpat_cov_func != nullptr,
              "Object of Gausian likelihood is not allocated!");
    assertMsg(const_cov_fun != nullptr,
              "Object of Gausian likelihood is not allocated!");
    assertMsg(lik_gaussian != nullptr,
              "Object of Gausian likelihood is not allocated!");
#endif

    sexpat_cov_func->setParameters(w);
    const_cov_fun->setParameters(w);
    lik_gaussian->setParameters(w);
}

bool GaussianProcessRegression::decomposeCovarianceMatrix(
    const EigenMatrix& x, const Eigen::VectorXd& x_ind)
{
    Eigen::LLT<EigenMatrix> llt;

    const double minimum_jitter = 1e-8;
    C.diagonal().array() += minimum_jitter;
    evaluateTrainingCovarianceMatrix(x, x_ind, C);
    llt.compute(C);
    L = llt.matrixL();

    return llt.info() == Eigen::NumericalIssue ? false : true;
}

void GaussianProcessRegression::calculatePosteriorMeanPrediction()
{
    // This function mimics the following MATLAB code
    // a = L'\(L\[E_all;G_all(:)]);
    a = L.transpose().triangularView<Eigen::Upper>().solve(b);
}

void GaussianProcessRegression::calculateMeanPrediction(
    const Eigen::VectorXd& y)
{
    // This function mimics the following MATLAB code
    // b = L\[E_all;G_all(:)];
    b = L.triangularView<Eigen::Lower>().solve(y);
}

void GaussianProcessRegression::calculatePotential(Observation& image1)
{
    // NOTE: C is a truncated covariance!

    Index_t N_im = image1.R.getNumRows();
    EigenMatrix R_mod;
    Eigen::VectorXd R_mod_ind;
    EigenMatrix KK;
    Eigen::VectorXd EG_R;
    Index_t counter = 0;
    aux::AuxiliaryFunctionality aux_func;

    // [repmat(R,D+1,1),reshape(repmat(0:D,N_im,1),[],1)]
    aux_func.assembleMatrixOfRepetitiveCoordinates(image1.R, R_mod, R_mod_ind);

    // KK =
    // gp_cov(gp,R_all2,[repmat(R,D+1,1),reshape(repmat(0:D,N_im,1),[],1)]);
    evaluateCovarianceMatrix(R_matrix, R_mod, R_indices, R_mod_ind, KK);
    // log_man << "KK is " << KK.size() << " over ( " << KK.rows() << ", " <<
    // KK.cols() << " )" << "\n";

    // EG_R = KK'*a;
    EG_R = KK.transpose() * a;

    // E_R = EG_R(1:N_im,1);
    image1.E.resize(1, N_im);
    for (Index_t n = 0; n < N_im; ++n)
        image1.E[n] = EG_R(n);

    // G_R = reshape(EG_R((N_im+1):end,1),N_im,D);
    counter = N_im;
    image1.G.resize(N_im, image1.R.getNumCols());
    for (Index_t j = 0; j < image1.G.getNumCols(); ++j) {
        for (Index_t i = 0; i < image1.G.getNumRows(); ++i) {
            image1.G(i, j) = EG_R(counter++);
        }
    }

    ++num_of_potential_calls;
}

/**
 * @brief Main function to optimize the Gaussian Process hyperparameters.
 *
 * This function orchestrates the optimization process. It runs in a loop to
 * ensure that the resulting hyperparameters are stable.
 *
 * 1. Selects a representative subset of data.
 * 2. Performs numerical optimization of hyperparameters on this subset.
 * 3. Monitors the new parameters for instability.
 * 4. If instability is detected, the history size for the subset is increased,
 * and the entire process is repeated.
 * 5. If the parameters are stable, the loop terminates, and the final model is
 * updated with the full dataset.
 *
 * @param observation The full set of observed data (coordinates, energies,
 * gradients).
 */
void GaussianProcessRegression::optimize(const Observation& observation)
{
    const int max_retries = 3;  // To prevent potential infinite loops
    int retry_count = 0;

    while (retry_count <= max_retries) {
        // 1. Select a subset of data for efficient optimization based on the
        // current history size.
        Observation optimization_observation =
            selectOptimizationSubset(observation);

        // 2. Run the numerical optimization on the subset.
        Eigen::VectorXd optimized_parameters =
            performHyperparameterOptimization(optimization_observation);

        // If the numerical optimization itself fails, retrying with more data
        // is unlikely to help. We break the loop and proceed with the failed
        // state.
        if (failedOptimizer) {
            log_man << "Numerical optimization failed. Aborting optimization "
                       "attempts.\n";
            setParameters(
                optimized_parameters);  // Set the (likely poor) parameters
            break;
        }

        // 3. Set the new parameters so they can be monitored.
        setParameters(optimized_parameters);

        // 4. Monitor for instability. If it returns true, we need to retry.
        if (use_hod) {
            bool should_retry =
                monitorAndAdaptHyperparameters(optimized_parameters);

            if (should_retry) {
                retry_count++;
                log_man
                    << "Instability detected. Re-running optimization with an "
                       "increased dataset. Attempt "
                    << retry_count << "/" << max_retries << ".\n";
            } else {
                // Stability is acceptable, optimization is complete.
                log_man << "Hyperparameters appear stable. Finalizing "
                           "optimization.\n";
                break;  // Exit the while loop
            }
        } else {
            // no oscillation check, just exit
            break;  // Exit the while loop
        }
    }

    if (retry_count > max_retries) {
        log_man << "Warning: Reached max retries for optimization due to "
                   "persistent instability.\n";
    }

    // 5. Update the model's internal state using the full dataset and the last
    // successful set of hyperparameters.
    updateModelWithFullData(observation);
}

/**
 * @brief Selects a subset of observations for hyperparameter optimization using
 * Farthest Point Sampling (FPS), after filtering out high-force configurations.
 *
 * This function first filters the full dataset to exclude any points where the
 * maximum force component exceeds a defined threshold. It then uses a combined
 * strategy on the remaining valid points:
 * 1. Selects the most recent valid points.
 * 2. Selects the rest using Farthest Point Sampling (FPS) for diversity.
 *
 * @param full_observation The complete set of observations.
 * @return An Observation object containing the selected subset.
 */
Observation GaussianProcessRegression::selectOptimizationSubset(
    const Observation& full_observation) const
{
    Observation optimization_observation;
    const int num_optimization_points = fps_options.history;
    const int num_latest_points = fps_options.latest_points;
    const int num_total_points = full_observation.R.getNumRows();
    // TODO(rg): Make this a configurable parameter
    const double force_threshold = std::numeric_limits<double>::infinity();
    std::vector<int> valid_indices;
    valid_indices.reserve(num_total_points);

    for (int i = 0; i < num_total_points; ++i) {
        double max_force = 0.0;
        for (int j = 0; j < full_observation.G.getNumCols(); ++j) {
            max_force = std::max(max_force, std::abs(full_observation.G(i, j)));
        }
        if (max_force < force_threshold) {
            valid_indices.push_back(i);
        }
    }

    const int num_valid_points = valid_indices.size();

    if (num_valid_points <= num_optimization_points) {
        for (int idx: valid_indices) {
            optimization_observation.R.append(full_observation.R, idx);
            optimization_observation.E.append(full_observation.E, idx);
            optimization_observation.G.append(full_observation.G, idx);
        }
        return optimization_observation;
    }

    // --- Combined Latest + Farthest Point Sampling (FPS) on VALID points ---
    aux::Distance distance_calculator;
    std::vector<int> selected_indices;
    std::unordered_set<int> selected_indices_set;

    // 1. Start with the 'num_latest_points' most recent *valid* points.
    int actual_latest_points = std::min(num_latest_points, num_valid_points);
    int num_fps_points_to_select =
        std::max(0, num_optimization_points - actual_latest_points);

    for (int i = 0; i < actual_latest_points; ++i) {
        int current_idx = valid_indices[num_valid_points - 1 - i];
        selected_indices.push_back(current_idx);
        selected_indices_set.insert(current_idx);
    }

    if (num_fps_points_to_select > 0) {
        // --- FPS Implementation for the remaining points ---
        // Map keeps track of min distances for unselected valid points
        std::unordered_map<int, double> min_dists_map;

        // Initialize min_dists for all valid, unselected points
        for (int candidate_idx: valid_indices) {
            if (selected_indices_set.count(candidate_idx) == 0) {
                double current_min_dist = std::numeric_limits<double>::max();
                Coord candidate_coord;
                candidate_coord.append(full_observation.R, candidate_idx);

                for (int latest_idx: selected_indices) {
                    Coord latest_coord;
                    latest_coord.append(full_observation.R, latest_idx);
                    Field<double> dist_matrix;
                    switch (fps_options.metric) {
                        case DistanceMetricType::MAX_1D_LOG:
                            distance_calculator.dist_max1Dlog(
                                candidate_coord, latest_coord,
                                this->atoms_config, dist_matrix);
                            break;
                        case DistanceMetricType::EMD:
                            distance_calculator.dist_emd(
                                candidate_coord, latest_coord,
                                this->atoms_config, dist_matrix);
                            break;
                        case DistanceMetricType::RMSD:
                            distance_calculator.dist_rmsd(
                                candidate_coord, latest_coord,
                                this->atoms_config, dist_matrix);
                            break;
                        default:
                            throw("Invalid distance metric for FPS.");
                    }
                    current_min_dist =
                        std::min(current_min_dist, dist_matrix(0, 0));
                }
                min_dists_map[candidate_idx] = current_min_dist;
            }
        }

        // 2. Iteratively select the remaining FPS points.
        for (int i = 0; i < num_fps_points_to_select; ++i) {
            double max_min_dist = -1.0;
            int farthest_point_idx = -1;

            // Find the unselected point with the maximum minimum distance
            for (auto const& [idx, dist]: min_dists_map) {
                if (dist > max_min_dist) {
                    max_min_dist = dist;
                    farthest_point_idx = idx;
                }
            }

            if (farthest_point_idx != -1) {
                // Add the new point to our selection
                selected_indices.push_back(farthest_point_idx);
                selected_indices_set.insert(farthest_point_idx);
                min_dists_map.erase(
                    farthest_point_idx);  // Remove from candidates

                Coord farthest_coord;
                farthest_coord.append(full_observation.R, farthest_point_idx);

                // Update min_dists for all other unselected points
                for (auto& entry: min_dists_map) {
                    int candidate_idx = entry.first;
                    double& current_min_dist = entry.second;

                    Coord candidate_coord;
                    candidate_coord.append(full_observation.R, candidate_idx);

                    Field<double> dist_matrix;
                    switch (fps_options.metric) {
                        case DistanceMetricType::MAX_1D_LOG:
                            distance_calculator.dist_max1Dlog(
                                candidate_coord, farthest_coord,
                                this->atoms_config, dist_matrix);
                            break;
                        case DistanceMetricType::EMD:
                            distance_calculator.dist_emd(
                                candidate_coord, farthest_coord,
                                this->atoms_config, dist_matrix);
                            break;
                        case DistanceMetricType::RMSD:
                            distance_calculator.dist_rmsd(
                                candidate_coord, farthest_coord,
                                this->atoms_config, dist_matrix);
                            break;
                        default:
                            throw("Invalid distance metric for FPS.");
                    }
                    current_min_dist =
                        std::min(current_min_dist, dist_matrix(0, 0));
                }
            } else {
                break;  // No more points to select
            }
        }
    }

    // 3. Build the final optimization_observation from the selected indices
    std::sort(selected_indices.begin(), selected_indices.end());
    for (int idx: selected_indices) {
        optimization_observation.R.append(full_observation.R, idx);
        optimization_observation.E.append(full_observation.E, idx);
        optimization_observation.G.append(full_observation.G, idx);
    }

    return optimization_observation;
}

/**
 * @brief Performs the numerical optimization of hyperparameters.
 *
 * This function assembles the necessary matrices from the observation subset,
 * invokes the chosen optimization algorithm, and returns the optimized
 * parameters.
 *
 * @param opt_observation The subset of data to use for optimization.
 * @return An Eigen::VectorXd containing the optimized hyperparameters.
 */
Eigen::VectorXd GaussianProcessRegression::performHyperparameterOptimization(
    const Observation& opt_observation)
{
    aux::AuxiliaryFunctionality aux_func;
    Eigen::MatrixXd opt_R_matrix;
    Eigen::VectorXd opt_R_indices;
    Eigen::VectorXd opt_energy_and_gradient;

    aux_func.assembleMatrixOfRepetitiveCoordinates(opt_observation.R,
                                                   opt_R_matrix, opt_R_indices);
    aux_func.assembleVectorFromEnergyAndGradient(opt_observation,
                                                 opt_energy_and_gradient);

    // Get current parameters to serve as the initial guess
    Eigen::VectorXd parameters_cf = sexpat_cov_func->combineParameters();
    Eigen::VectorXd old_parameters_cf = parameters_cf;

    std::visit([&](auto& opt) { opt.setAlgorithmSettings(opt_alg_settings); },
               popt);

    if (this->report_level >= 2) {
        log_man << "Starting optimization with " << opt_R_matrix.rows()
                << " data points.\n";
    }

    double initial_strength = 1e-4;
    double growth_factor = 1e-3;
    double current_barrier_strength =
        initial_strength + growth_factor * opt_observation.R.getNumPoints();
    auto start = std::chrono::steady_clock::now();

    std::visit(
        [&](auto& _popt) {
            _popt.optimize(
                opt_R_matrix, opt_R_indices, opt_energy_and_gradient,
                parameters_cf,
                &gpr::GaussianProcessRegression::evaluateEnergyAndGradient,
                *this, current_barrier_strength);
            failedOptimizer =
                _popt.failedOptim;  // Set member flag based on result
        },
        popt);

    if (this->report_level >= 2) {
        auto elp_time = (std::chrono::steady_clock::now() - start);
        auto elp_seconds =
            std::chrono::duration_cast<std::chrono::duration<double>>(elp_time);

        log_man << "Optimization finished in: " << elp_seconds.count()
                << "s. Success: " << !failedOptimizer << "\n";
        log_man << "Parameter change norm: "
                << (old_parameters_cf - parameters_cf).norm() << "\n";
    }

    return parameters_cf;
}

/**
 * @brief Monitors hyperparameter history for oscillations and adapts the FPS
 * history size.
 *
 * If the ratio of sign flips in the recent history of hyperparameter updates
 * exceeds a threshold, it indicates instability. In response, the size of the
 * data subset for future optimizations is increased.
 *
 * @param new_parameters The latest vector of optimized hyperparameters.
 * @return bool Returns 'true' if the history size was increased (signaling a
 * need to retry), 'false' otherwise.
 */
bool GaussianProcessRegression::monitorAndAdaptHyperparameters(
    const Eigen::VectorXd& new_parameters)
{
    // 1. Store the latest hyperparameter vector and maintain history buffer
    // size.
    hyperparameter_history_.push_back(new_parameters);
    while (hyperparameter_history_.size() > monitoring_window_ + 1) {
        hyperparameter_history_.pop_front();
    }

    // 2. Check for instability only if we have enough data and history is not
    // maxed out.
    if (hyperparameter_history_.size() != monitoring_window_ + 1 ||
        fps_options.history >= max_history_) {
        return false;  // Not enough data or already maxed out, so no retry.
    }

    const int num_params = new_parameters.size();
    const int total_possible_flips = num_params * (monitoring_window_ - 1);
    if (total_possible_flips <= 0) {
        return false;  // Cannot compute flips, so no retry.
    }

    int sign_flips = 0;
    // Iterate through a 3-point sliding window to check for oscillations.
    for (size_t i = 0; i < monitoring_window_ - 1; ++i) {
        const Eigen::VectorXd& h_old = hyperparameter_history_[i];
        const Eigen::VectorXd& h_mid = hyperparameter_history_[i + 1];
        const Eigen::VectorXd& h_new = hyperparameter_history_[i + 2];

        Eigen::VectorXd delta1 = h_mid - h_old;
        Eigen::VectorXd delta2 = h_new - h_mid;

        // Check for sign flips for each hyperparameter component
        for (int j = 0; j < num_params; ++j) {
            if (delta1(j) * delta2(j) < -1e-9) {  // Using a small tolerance
                sign_flips++;
            }
        }
    }

    // 3. If instability is detected, increase the history size and signal for a
    // retry.
    double flip_ratio = static_cast<double>(sign_flips) / total_possible_flips;
    if (flip_ratio > flip_threshold_) {
        log_man << "Hyperparameter instability detected! Flip ratio: "
                << flip_ratio << " > " << flip_threshold_ << "\n";

        int old_history = fps_options.history;
        fps_options.history =
            std::min(old_history + history_increment_, max_history_);

        log_man << "Increasing FPS history size from " << old_history << " to "
                << fps_options.history << "\n";

        // Clear history after adjustment to start fresh monitoring on the next
        // attempt.
        hyperparameter_history_.clear();
        return true;  // Signal to retry
    }

    // No instability detected, no retry needed.
    return false;
}

/**
 * @brief Updates the GPR model's internal state using the full dataset.
 *
 * After hyperparameters are optimized, this function re-evaluates the
 * covariance matrix and its Cholesky decomposition using all available data
 * points.
 *
 * @param full_observation The complete set of observations.
 */
void GaussianProcessRegression::updateModelWithFullData(
    const Observation& full_observation)
{
    aux::AuxiliaryFunctionality aux_func;

    // Assemble matrices with the full dataset
    aux_func.assembleMatrixOfRepetitiveCoordinates(full_observation.R, R_matrix,
                                                   R_indices);
    aux_func.assembleVectorFromEnergyAndGradient(full_observation,
                                                 energy_and_gradient);

    // Re-decompose the covariance matrix with new hyperparameters and full data
    is_decomposed_succesfully = decomposeCovarianceMatrix(R_matrix, R_indices);

    if (!is_decomposed_succesfully) {
        log_man << "Warning! Matrix L is not positive definite after "
                   "re-evaluation with full dataset.\n";
    } else {
        calculateMeanPrediction(energy_and_gradient);
        calculatePosteriorMeanPrediction();
    }

    if (this->report_level >= 2) {
        log_man << "Final model updated with "
                << full_observation.R.getNumRows() << " data points.\n";
        log_man << "magnSigma2: " << sexpat_cov_func->getMagnSigma2() << "\n";
        log_man << "lengthScales:\n"
                << sexpat_cov_func->getLengthScaleRef().extractEigenVector()
                << "\n";
    }
}

} /* namespace gpr */
