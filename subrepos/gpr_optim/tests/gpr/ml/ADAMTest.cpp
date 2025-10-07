#include "ADAMTest.h"

#include "../../gpr/ml/GaussianProcessRegressionTest.h"

namespace gpr {
namespace tests {

ADAMTest::ADAMTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

ADAMTest::~ADAMTest() { }

TEST_F(ADAMTest, optimization)
{
    gpr::GaussianProcessRegression gpr_model;
    Eigen::VectorXd w(3);
    Eigen::VectorXd y(14);
    gpr::EigenMatrix x(14, 6);
    gpr::EigenMatrix tmp_matrix(14, 7);
    Eigen::VectorXd x_ind(14, 1);
    Eigen::VectorXd w_ref(3);
    gpr::Field<double> tmp;
    gpr::AtomsConfiguration conf_info;
    gpr::PriorBase parameters;
    gpr::GPRSetup gpr_parameters;
    gpr::OptimizationAlgorithmSettings settings;
    io::FileManager io_manager;

    gpr_model.getSexpAtCovarianceFunction()->getLengthScaleRef().resize(1, 2);
    conf_info.atoms_froz_active.positions.resize(1, 26 * 3);
    conf_info.atoms_froz_active.type.resize(1, 26);
    conf_info.atoms_mov.type.resize(1, 2);
    conf_info.pairtype.resize(2, 2);

    // TODO(rg): These are just different convergence regimes
    w_ref(0) = -11.971961999501948;
    w_ref(1) = 0.0038825279841619817;
    w_ref(2) = -3.5177815029865771;

    w(0) = -14.7861445569979e+000;
    w(1) = 10.3888257106399e-003;
    w(2) = -8.02546595438817e+000;

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

    tmp.resize((gpr::Index_t)tmp_matrix.rows(),
               (gpr::Index_t)tmp_matrix.cols());
    io_manager.readFromPlainFile("tests/reference/gpr/ml/InputX1.dat",
                                 tmp.getInternalVector());
    tmp_matrix = tmp.extractEigenMatrix();
    x = tmp_matrix.block(0, 0, 14, 6);
    x_ind = tmp_matrix.col(6);

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

    gpr_model.getLikGaussian()->setSigma2(10.0000000000000e-009);

    parameters.setMu(0.);
    parameters.setS2(1.);
    parameters.setNu(20.);

    gpr_model.getSexpAtCovarianceFunction()->setPriorParametersGaussian(
        parameters);
    gpr_model.getSexpAtCovarianceFunction()->setPriorParametersSqrtt(
        parameters);
    gpr_model.getSexpAtCovarianceFunction()->setMagnSigma2(
        6.93874748072254e-009);
    gpr_model.getSexpAtCovarianceFunction()->setLengthScale(
        888.953211438594e-006);
    gpr_model.getSexpAtCovarianceFunction()->setConfInfo(conf_info);

    gpr_model.getConstantCovarianceFunction()->setConstSigma2(1.);

    gpr_parameters.jitter_sigma2 = 0.;
    gpr_model.setParameters(gpr_parameters);

    // Set up solver
    settings.max_iter = 8000;
    settings.learning_rate = 0.08;
    settings.learning_rate_decay = 0.999;
    settings.tolerance_func = 1e-8;
    settings.tolerance_sol = 1e-4;
    settings.report_level = 2;
    settings.beta1 = 0.9;
    settings.beta2 = 0.94;
    settings.epsilon = 1.0e-8;
    settings.amsgrad = true;
    settings.weight_decay = 0.0;

    adam.setAlgorithmSettings(settings);

    adam.optimize(x, x_ind, y, w,
                  &gpr::GaussianProcessRegression::evaluateEnergyAndGradient,
                  gpr_model);

    const double absolute_tolerance = 1e-3;

    for (gpr::Index_t n = 0; n < w_ref.rows(); ++n) {
        EXPECT_NEAR(w(n), w_ref(n), absolute_tolerance)
            << "Component " << n << " of the result vector is not "
            << "close enough to the expected value.";
    }
}

} /* namespace tests */
} /* namespace gpr */
