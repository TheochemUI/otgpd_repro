//
//  SCGTest.cpp
//  gpr_dimer
//
//  Created by Maxim Masterov on 25/11/2020.
//

#include "SCGTest.h"

#include "../../gpr/ml/GaussianProcessRegressionTest.h"

namespace gpr {
namespace tests {

SCGTest::SCGTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

SCGTest::~SCGTest() { }

TEST_F(SCGTest, optimization)
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

    w_ref(0) = -11.7150224969728e+000;
    w_ref(1) = 15.3888257106399e-003;
    w_ref(2) = -3.44764255084076e+000;

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
    settings.check_derivative = false;
    settings.lambda_limit = 1e20;
    settings.max_iter = 400;
    settings.tolerance_func = 1e-4;
    settings.tolerance_sol = 1e-4;
    settings.report_level = 0;

    scg.setAlgorithmSettings(settings);

    scg.optimize(x, x_ind, y, w,
                 &gpr::GaussianProcessRegression::evaluateEnergyAndGradient,
                 gpr_model);

    for (gpr::Index_t n = 0; n < w_ref.rows(); ++n)
        EXPECT_LE(std::fabs(w(0) - w_ref(0)), settings.tolerance_func)
            << "Components of the result vector are not "
               "equal to the expected ones.";
}

} /* namespace tests */
} /* namespace gpr */
