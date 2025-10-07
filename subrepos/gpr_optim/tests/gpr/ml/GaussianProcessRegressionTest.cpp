/*
 * GaussianProcessRegressionTest.cpp
 *
 *  Created on: 16 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "GaussianProcessRegressionTest.h"

#include <cmath>

#include "../../../gpr/observation_models/LikGaussian.h"
#include "gpr/auxiliary/AdditionalFunctionality.h"

namespace gpr {
namespace tests {

GaussianProcessRegressionTest::GaussianProcessRegressionTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.

    gpr::Field<double> tmp_vec;
    gpr::EigenMatrix tmp_matrix;

    x1.resize(14, 6);
    x2.resize(7, 6);
    x1_ind.resize(x1.rows());
    x2_ind.resize(x2.rows());
    conf_info.atoms_froz_active.positions.resize(1, 26 * 3);
    conf_info.atoms_froz_active.type.resize(1, 26);
    conf_info.atoms_mov.type.resize(1, 2);
    conf_info.pairtype.resize(2, 2);

    tmp_vec.resize((gpr::Index_t)x1.rows(), (gpr::Index_t)x1.cols() + 1);
    io_manager.readFromPlainFile("tests/reference/gpr/ml/InputX1.dat",
                                 tmp_vec.getInternalVector());
    tmp_matrix = tmp_vec.extractEigenMatrix();
    x1 = tmp_matrix.block(0, 0, 14, 6);
    x1_ind = tmp_matrix.col(6);

    tmp_vec.resize((gpr::Index_t)x2.rows(), (gpr::Index_t)x2.cols() + 1);
    io_manager.readFromPlainFile("tests/reference/gpr/ml/InputX2.dat",
                                 tmp_vec.getInternalVector());
    tmp_matrix = tmp_vec.extractEigenMatrix();
    x2 = tmp_matrix.block(0, 0, 7, 6);
    x2_ind = tmp_matrix.col(6);

    io_manager.readFromPlainFile(
        "tests/reference/gpr/ml/ConfigurationFrozen.dat",
        conf_info.atoms_froz_active.positions.getInternalVector());

    conf_info.atoms_mov.type.set(0);
    conf_info.atoms_froz_active.type.set(1);

    conf_info.pairtype(0, 0) = 0;
    conf_info.pairtype(0, 1) = 1;
    conf_info.pairtype(1, 0) = 1;
    conf_info.pairtype(1, 1) = EMPTY;

    conf_info.n_pt = 2;
}

GaussianProcessRegressionTest::~GaussianProcessRegressionTest() { }

TEST_F(GaussianProcessRegressionTest, evaluateTrainingCovarianceMatrix)
{
    gpr::EigenMatrix K;
    gpr::Field<double> K_ref;
    gpr::GPRSetup gpr_parameters;

    this->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);

    K_ref.resize(14, 14);

    gpr_parameters.jitter_sigma2 = 0.;
    this->setParameters(gpr_parameters);

    this->getSexpAtCovarianceFunction()->setMagnSigma2(6.93874748072254e-009);
    this->getSexpAtCovarianceFunction()->setLengthScale(888.953211438594e-006);
    this->getSexpAtCovarianceFunction()->setConfInfo(conf_info);

    this->getConstantCovarianceFunction()->setConstSigma2(1.);

    io_manager.readFromPlainFile(
        "tests/reference/gpr/ml/ReferenceKevalTrCovMatrix.dat",
        K_ref.getInternalVector());

    this->evaluateTrainingCovarianceMatrix(x1, x1_ind, K);

    EXPECT_EQ(K.rows(), K_ref.getNumRows())
        << "Number of rows in the result matrix and the reference matrix "
           "are not equal.";

    EXPECT_EQ(K.cols(), K_ref.getNumCols())
        << "Number of columns in the result matrix and the reference matrix "
           "are not equal.";

    for (gpr::Index_t i = 0; i < K_ref.getNumRows(); ++i) {
        for (gpr::Index_t j = 0; j < K_ref.getNumCols(); ++j) {
            EXPECT_LE(std::fabs(K(i, j) - K_ref(i, j)), threshold)
                << "Element (" << i << ", " << j
                << ") of the result matrix "
                   "are not equal to the expected ones.";
        }
    }
}

TEST_F(GaussianProcessRegressionTest, evaluateCovarianceMatrix)
{
    gpr::EigenMatrix C;
    gpr::AtomsConfiguration conf_info;
    gpr::Field<double> K_ref;
    gpr::GPRSetup gpr_parameters;

    this->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
    conf_info.atoms_froz_active.positions.resize(1, 26 * 3);
    conf_info.atoms_froz_active.type.resize(1, 26);
    conf_info.atoms_mov.type.resize(1, 2);
    conf_info.pairtype.resize(2, 2);
    K_ref.resize(14, 7);

    gpr_parameters.jitter_sigma2 = 0.;
    this->setParameters(gpr_parameters);

    io_manager.readFromPlainFile(
        "tests/reference/gpr/ml/ConfigurationFrozen.dat",
        conf_info.atoms_froz_active.positions.getInternalVector());

    conf_info.atoms_mov.type.set(0);
    conf_info.atoms_froz_active.type.set(1);

    conf_info.pairtype(0, 0) = 0;
    conf_info.pairtype(0, 1) = 1;
    conf_info.pairtype(1, 0) = 1;
    conf_info.pairtype(1, 1) = EMPTY;

    conf_info.n_pt = 2;

    this->getSexpAtCovarianceFunction()->setMagnSigma2(8.17015647431529e-006);
    this->getSexpAtCovarianceFunction()->getLengthScaleRef()(0, 0) =
        1.01550784341955e+000;
    this->getSexpAtCovarianceFunction()->getLengthScaleRef()(0, 1) =
        31.8205633854612e-003;
    this->getSexpAtCovarianceFunction()->setConfInfo(conf_info);

    this->getConstantCovarianceFunction()->setConstSigma2(1.);

    // Reference C data stored in K_ref:
    io_manager.readFromPlainFile(
        "tests/reference/gpr/ml/ReferenceKevalCovMatrix.dat",
        K_ref.getInternalVector());

    this->evaluateCovarianceMatrix(x1, x2, x1_ind, x2_ind, C);

    EXPECT_EQ(C.rows(), K_ref.getNumRows())
        << "Number of rows in the result matrix and the reference matrix "
           "are not equal.";

    EXPECT_EQ(C.cols(), K_ref.getNumCols())
        << "Number of columns in the result matrix and the reference matrix "
           "are not equal.";

    for (gpr::Index_t i = 0; i < K_ref.getNumRows(); ++i) {
        for (gpr::Index_t j = 0; j < K_ref.getNumCols(); ++j) {
            EXPECT_LE(std::fabs(C(i, j) - K_ref(i, j)), threshold)
                << "Element (" << i << ", " << j
                << ") of the result matrix "
                   "are not equal to the expected ones.";
        }
    }
}

TEST_F(GaussianProcessRegressionTest, evaluateEnergyAndGradient)
{
    Eigen::VectorXd y;
    Eigen::VectorXd w;
    gpr::EigenMatrix K;
    gpr::EigenMatrix K_ref;
    gpr::PriorBase parameters;
    double energy;
    Eigen::VectorXd gradient;
    gpr::EnergyAndGradient energy_and_gradient;
    gpr::GPRSetup gpr_parameters;
    double energy_ref = -2.04481678030626e+000;
    gpr::vector3_reg gradient_ref;

    y.resize(14);
    w.resize(3);

    energy_and_gradient.energy = &energy;
    energy_and_gradient.gradient = &gradient;

    this->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);

    gradient_ref.set(3.02932909658129e+000, -2.72358921850395e+000,
                     -7.57715196734936e+000);

    w(0) = -18.7861445569979e+000;
    w(1) = -7.02546595438817e+000;
    w(2) = -7.02546595438817e+000;

    y(0) = 0.00000000000000e+000;
    y(1) = 370.498489473903e-006;
    y(2) = 22.2174753484639e-003;
    y(3) = 21.4452456338877e-003;
    y(4) = -22.2402530044595e-003;
    y(5) = -33.1231900862804e-003;
    y(6) = 18.5353722102027e-003;
    y(7) = 28.3157675762483e-003;
    y(8) = -17.3837548235882e-003;
    y(9) = -15.8030003032378e-003;
    y(10) = 21.9137885404279e-003;
    y(11) = 20.8628464314946e-003;
    y(12) = -25.8161071576406e-003;
    y(13) = -38.3251019392285e-003;

    this->getLikGaussian()->setSigma2(10.0000000000000e-009);

    parameters.setMu(0.);
    parameters.setS2(1.);
    parameters.setNu(20.);
    this->getSexpAtCovarianceFunction()->setPriorParametersGaussian(parameters);
    this->getSexpAtCovarianceFunction()->setPriorParametersSqrtt(parameters);
    this->getSexpAtCovarianceFunction()->setMagnSigma2(6.93874748072254e-009);
    this->getSexpAtCovarianceFunction()->setLengthScale(888.953211438594e-006);
    this->getSexpAtCovarianceFunction()->setConfInfo(conf_info);

    this->getConstantCovarianceFunction()->setConstSigma2(1.);

    gpr_parameters.jitter_sigma2 = 0.;
    this->setParameters(gpr_parameters);

    this->evaluateEnergyAndGradient(w, x1, x1_ind, y, energy_and_gradient);

    EXPECT_LE(std::fabs(*energy_and_gradient.energy - energy_ref), threshold)
        << "Calculated energy is not equal to the expected one.";

    EXPECT_EQ(energy_and_gradient.gradient->rows(), 3)
        << "Number of elements in gradient field is not equal to the "
           "expected one.";

    EXPECT_LE(
        std::fabs(energy_and_gradient.gradient->operator()(0) - gradient_ref.x),
        threshold * 1e3)
        << "X-component of the gradient is not equal "
           "to the expected one.";
    EXPECT_LE(
        std::fabs(energy_and_gradient.gradient->operator()(1) - gradient_ref.y),
        threshold * 1e3)
        << "Y-component of the gradient is not equal "
           "to the expected one.";
    EXPECT_LE(
        std::fabs(energy_and_gradient.gradient->operator()(2) - gradient_ref.z),
        threshold * 1e3)
        << "Z-component of the gradient is not equal "
           "to the expected one.";
}

TEST_F(GaussianProcessRegressionTest, calculatePotential)
{
    gpr::Observation observation_image1;
    gpr::Observation observation_all;
    gpr::Observation observation_ref;
    aux::AuxiliaryFunctionality aux_func;

    observation_all.R.resize(2, 6);
    observation_all.E.resize(1, 2);
    observation_all.G.resize(2, 6);

    observation_image1.R.resize(2, 6);
    observation_image1.E.clear();
    observation_image1.G.clear();

    observation_ref.E.resize(1, 2);
    observation_ref.G.resize(2, 6);

    observation_all.R(0, 0) = 8.98237316483057e+000;
    observation_all.R(0, 1) = 9.93723083577204e+000;
    observation_all.R(0, 2) = 7.89441632385049e+000;
    observation_all.R(0, 3) = 7.65248322727496e+000;
    observation_all.R(0, 4) = 9.95590549457398e+000;
    observation_all.R(0, 5) = 7.87787958998366e+000;

    observation_all.R(1, 0) = 8.97856277303058e+000;
    observation_all.R(1, 1) = 9.93211628067229e+000;
    observation_all.R(1, 2) = 7.89882761414426e+000;
    observation_all.R(1, 3) = 7.64888749663556e+000;
    observation_all.R(1, 4) = 9.95517051512886e+000;
    observation_all.R(1, 5) = 7.87274215046670e+000;

    observation_all.E(0, 0) = 0.00000000000000e+000;
    observation_all.E(0, 1) = 370.498489473903e-006;

    observation_all.G(0, 0) = 22.2174753484639e-003;
    observation_all.G(0, 1) = -22.2402530044595e-003;
    observation_all.G(0, 2) = 18.5353722102027e-003;
    observation_all.G(0, 3) = -17.3837548235882e-003;
    observation_all.G(0, 4) = 21.9137885404279e-003;
    observation_all.G(0, 5) = -25.8161071576406e-003;

    observation_all.G(1, 0) = 21.4452456338877e-003;
    observation_all.G(1, 1) = -33.1231900862804e-003;
    observation_all.G(1, 2) = 28.3157675762483e-003;
    observation_all.G(1, 3) = -15.8030003032378e-003;
    observation_all.G(1, 4) = 20.8628464314946e-003;
    observation_all.G(1, 5) = -38.3251019392285e-003;

    observation_image1.R(0, 0) = 8.98237316483057e+000;
    observation_image1.R(0, 1) = 9.93723083577204e+000;
    observation_image1.R(0, 2) = 7.89441632385049e+000;
    observation_image1.R(0, 3) = 7.65248322727496e+000;
    observation_image1.R(0, 4) = 9.95590549457398e+000;
    observation_image1.R(0, 5) = 7.87787958998366e+000;

    observation_image1.R(1, 0) = 8.97856277303058e+000;
    observation_image1.R(1, 1) = 9.93211628067229e+000;
    observation_image1.R(1, 2) = 7.89882761414426e+000;
    observation_image1.R(1, 3) = 7.64888749663556e+000;
    observation_image1.R(1, 4) = 9.95517051512886e+000;
    observation_image1.R(1, 5) = 7.87274215046670e+000;

    this->getLikGaussian()->setSigma2(10.0000000000000e-009);

    this->getConstantCovarianceFunction()->setConstSigma2(1.);

    this->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
    this->getSexpAtCovarianceFunction()->getLengthScaleRef()(0, 0) =
        1.01550784341955e+000;
    this->getSexpAtCovarianceFunction()->getLengthScaleRef()(0, 1) =
        31.8205633854612e-003;
    this->getSexpAtCovarianceFunction()->setMagnSigma2(8.17015647431529e-006);
    this->getSexpAtCovarianceFunction()->setConfInfo(conf_info);

    this->setJitterSigma2(0.);

    observation_ref.E(0, 0) = -704.595901167195e-009;
    observation_ref.E(0, 1) = 371.203048873380e-006;

    observation_ref.G(0, 0) = 22.2187412441741e-003;
    observation_ref.G(0, 1) = -22.2419538432585e-003;
    observation_ref.G(0, 2) = 18.5375109120299e-003;
    observation_ref.G(0, 3) = -17.3807794897030e-003;
    observation_ref.G(0, 4) = 21.9136490152291e-003;
    observation_ref.G(0, 5) = -25.8188654019664e-003;

    observation_ref.G(1, 0) = 21.4437403317027e-003;
    observation_ref.G(1, 1) = -33.1212926208074e-003;
    observation_ref.G(1, 2) = 28.3134580120646e-003;
    observation_ref.G(1, 3) = -15.8057733045760e-003;
    observation_ref.G(1, 4) = 20.8628503344299e-003;
    observation_ref.G(1, 5) = -38.3221119846374e-003;

    aux_func.assembleMatrixOfRepetitiveCoordinates(
        observation_all.R, this->R_matrix, this->R_indices);
    aux_func.assembleVectorFromEnergyAndGradient(observation_all,
                                                 this->energy_and_gradient);
    this->decomposeCovarianceMatrix(this->R_matrix, this->R_indices);
    this->calculateMeanPrediction(this->energy_and_gradient);
    this->calculatePosteriorMeanPrediction();
    this->calculatePotential(observation_image1);

    EXPECT_EQ(observation_image1.E.getNumRows(), observation_ref.E.getNumRows())
        << "Number of rows in the field of energy is not equal to the "
           "expected one.";

    EXPECT_EQ(observation_image1.E.getNumCols(), observation_ref.E.getNumCols())
        << "Number of columns in the field of energy is not equal to the "
           "expected one.";

    EXPECT_EQ(observation_image1.G.getNumRows(), observation_ref.G.getNumRows())
        << "Number of rows in the field of gradients is not equal to the "
           "expected one.";

    EXPECT_EQ(observation_image1.G.getNumCols(), observation_ref.G.getNumCols())
        << "Number of columns in the field of gradients is not equal to the "
           "expected one.";

    for (gpr::Index_t i = 0; i < observation_image1.E.getNumRows(); ++i) {
        for (gpr::Index_t j = 0; j < observation_image1.E.getNumCols(); ++j) {
            EXPECT_LE(
                std::fabs(observation_image1.E(i, j) - observation_ref.E(i, j)),
                threshold)
                << "Element (" << i << ", " << j
                << ") in the field of energy "
                   "is not equal to the expected one.";
        }
    }

    for (gpr::Index_t i = 0; i < observation_image1.G.getNumRows(); ++i) {
        for (gpr::Index_t j = 0; j < observation_image1.G.getNumCols(); ++j) {
            EXPECT_LE(
                std::fabs(observation_image1.G(i, j) - observation_ref.G(i, j)),
                threshold)
                << "Element (" << i << ", " << j
                << ") in the field of gradients "
                   "is not equal to the expected one.";
        }
    }
}

TEST_F(GaussianProcessRegressionTest, calculatePosteriorMeanPrediction)
{
    gpr::Field<double> tmp;
    gpr::Observation obs_all;
    aux::AuxiliaryFunctionality aux_func;

    gpr::Coord conf_fro;
    Eigen::VectorXd a_ref;

    R_matrix.resize(14, 7);
    tmp.resize((gpr::Index_t)R_matrix.rows(), (gpr::Index_t)R_matrix.cols());
    obs_all.E.resize(1, 2);
    obs_all.G.resize(2, 6);
    obs_all.R.resize(2, 6);
    a_ref.resize(14);
    conf_info.atoms_froz_active.positions.resize(1, 26 * 3);
    conf_info.atoms_froz_active.type.resize(1, 26);
    conf_info.atoms_mov.type.resize(1, 2);
    conf_info.pairtype.resize(2, 2);

    io_manager.readFromPlainFile(
        "tests/reference/gpr/ml/InputMatrixRprepAuxVector.dat",
        tmp.getInternalVector());
    R_matrix = tmp.extractEigenMatrix();

    obs_all.E(0, 0) = 0.00000000000000e+000;
    obs_all.E(0, 1) = 370.498489473903e-006;

    obs_all.G(0, 0) = 22.2174753484639e-003;
    obs_all.G(0, 1) = -22.2402530044595e-003;
    obs_all.G(0, 2) = 18.5353722102027e-003;
    obs_all.G(0, 3) = -17.3837548235882e-003;
    obs_all.G(0, 4) = 21.9137885404279e-003;
    obs_all.G(0, 5) = -25.8161071576406e-003;

    obs_all.G(1, 0) = 21.4452456338877e-003;
    obs_all.G(1, 1) = -33.1231900862804e-003;
    obs_all.G(1, 2) = 28.3157675762483e-003;
    obs_all.G(1, 3) = -15.8030003032378e-003;
    obs_all.G(1, 4) = 20.8628464314946e-003;
    obs_all.G(1, 5) = -38.3251019392285e-003;

    obs_all.R(0, 0) = 8.98237316483057e+000;
    obs_all.R(0, 1) = 9.93723083577204e+000;
    obs_all.R(0, 2) = 7.89441632385049e+000;
    obs_all.R(0, 3) = 7.65248322727496e+000;
    obs_all.R(0, 4) = 9.95590549457398e+000;
    obs_all.R(0, 5) = 7.87787958998366e+000;

    obs_all.R(1, 0) = 8.97856277303058e+000;
    obs_all.R(1, 1) = 9.93211628067229e+000;
    obs_all.R(1, 2) = 7.89882761414426e+000;
    obs_all.R(1, 3) = 7.64888749663556e+000;
    obs_all.R(1, 4) = 9.95517051512886e+000;
    obs_all.R(1, 5) = 7.87274215046670e+000;

    a_ref(0) = 70.4595903686042e+000;
    a_ref(1) = -70.4559397222489e+000;
    a_ref(2) = -126.589571024787e+000;
    a_ref(3) = 150.530218499774e+000;
    a_ref(4) = 170.083879912637e+000;
    a_ref(5) = -189.746547307559e+000;
    a_ref(6) = -213.870182712076e+000;
    a_ref(7) = 230.956418369173e+000;
    a_ref(8) = -297.533388526113e+000;
    a_ref(9) = 277.300133822936e+000;
    a_ref(10) = 13.9525198853451e+000;
    a_ref(11) = -390.293526144286e-003;
    a_ref(12) = 275.824432578006e+000;
    a_ref(13) = -298.995459107190e+000;

    io_manager.readFromPlainFile(
        "tests/reference/gpr/ml/ConfigurationFrozen.dat",
        conf_info.atoms_froz_active.positions.getInternalVector());

    conf_info.atoms_mov.type.set(0);
    conf_info.atoms_froz_active.type.set(1);

    conf_info.pairtype(0, 0) = 0;
    conf_info.pairtype(0, 1) = 1;
    conf_info.pairtype(1, 0) = 1;
    conf_info.pairtype(1, 1) = EMPTY;

    conf_info.n_pt = 2;

    this->getLikGaussian()->setSigma2(10.0000000000000e-009);

    this->getConstantCovarianceFunction()->setConstSigma2(1.);

    this->getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
    this->getSexpAtCovarianceFunction()->getLengthScaleRef()(0, 0) =
        1.01550784341955e+000;
    this->getSexpAtCovarianceFunction()->getLengthScaleRef()(0, 1) =
        31.8205633854612e-003;
    this->getSexpAtCovarianceFunction()->setMagnSigma2(8.17015647431529e-006);
    this->getSexpAtCovarianceFunction()->setConfInfo(conf_info);

    //    this->setDeriv(7);
    this->setJitterSigma2(0.);

    aux_func.assembleMatrixOfRepetitiveCoordinates(obs_all.R, this->R_matrix,
                                                   this->R_indices);
    aux_func.assembleVectorFromEnergyAndGradient(obs_all,
                                                 this->energy_and_gradient);

    this->decomposeCovarianceMatrix(this->R_matrix, this->R_indices);
    this->calculateMeanPrediction(this->energy_and_gradient);
    this->calculatePosteriorMeanPrediction();

    EXPECT_EQ(this->a.rows(), a_ref.rows())
        << "Number of rows in the 'a' vector is not equal to the "
           "expected one.";

    for (gpr::Index_t n = 0; n < a.rows(); ++n) {
        EXPECT_LE(std::fabs(this->a(n) - a_ref(n)), threshold * 1e6)
            << "Element (" << n
            << ") in the field of energy "
               "is not equal to the expected one.";
    }
}

} /* namespace tests */
} /* namespace gpr */
