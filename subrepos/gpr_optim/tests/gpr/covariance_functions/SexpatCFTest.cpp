/*
 * SexpAtTest.cpp
 *
 *  Created on: 14 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "SexpatCFTest.h"

#include <cmath>

namespace gpr {
namespace tests {

SexpAtTest::SexpAtTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.

    x1.resize(2, 2 * 3);
    x2.resize(2, 2 * 3);
    conf_info.atoms_froz_active.positions.resize(1, 26 * 3);
    conf_info.atoms_froz_active.type.resize(1, 26);
    conf_info.atoms_mov.type.resize(1, 2);
    conf_info.pairtype.resize(2, 2);

    x1(0, 0) = 8.98237316483057e+000;
    x1(0, 1) = 9.93723083577204e+000;
    x1(0, 2) = 7.89441632385049e+000;
    x1(0, 3) = 7.65248322727496e+000;
    x1(0, 4) = 9.95590549457398e+000;
    x1(0, 5) = 7.87787958998366e+000;

    x1(1, 0) = 8.97856277303058e+000;
    x1(1, 1) = 9.93211628067229e+000;
    x1(1, 2) = 7.89882761414426e+000;
    x1(1, 3) = 7.64888749663556e+000;
    x1(1, 4) = 9.95517051512886e+000;
    x1(1, 5) = 7.87274215046670e+000;

    x2(0, 0) = 8.98237316483057e+000;
    x2(0, 1) = 9.93723083577204e+000;
    x2(0, 2) = 7.89441632385049e+000;
    x2(0, 3) = 7.65248322727496e+000;
    x2(0, 4) = 9.95590549457398e+000;
    x2(0, 5) = 7.87787958998366e+000;

    x2(1, 0) = 8.97856277303058e+000;
    x2(1, 1) = 9.93211628067229e+000;
    x2(1, 2) = 7.89882761414426e+000;
    x2(1, 3) = 7.64888749663556e+000;
    x2(1, 4) = 9.95517051512886e+000;
    x2(1, 5) = 7.87274215046670e+000;

    io_manager.readFromPlainFile(
        "tests/reference/gpr/covariance_functions/ConfigurationFrozen.dat",
        conf_info.atoms_froz_active.positions.getInternalVector());

    conf_info.atoms_mov.type.set(0);
    conf_info.atoms_froz_active.type.set(1);

    conf_info.pairtype(0, 0) = 0;
    conf_info.pairtype(0, 1) = 1;
    conf_info.pairtype(1, 0) = 1;
    conf_info.pairtype(1, 1) = EMPTY;

    conf_info.n_pt = 2;
}

SexpAtTest::~SexpAtTest() { }

TEST_F(SexpAtTest, sexpat_lp)
{
    gpr::PriorBase parameters;
    double lp;
    double lp_ref = -26.9064625685046e+000;

    parameters.setMu(0.);
    parameters.setS2(1.);
    parameters.setNu(20.);
    sexpat.setMagnSigma2(6.93874748072254e-009);

    sexpat.getLengthScaleRef().resize(1, 2);
    sexpat.setLengthScale(888.953211438594e-006);

    sexpat.setPriorParametersGaussian(parameters);
    sexpat.setPriorParametersSqrtt(parameters);

    lp = sexpat.calculateLogPrior();

    EXPECT_LE(fabs(lp - lp_ref), threshold)
        << "Value of `lp` is not equal to the expected one.";
}

TEST_F(SexpAtTest, sexpat_lpg)
{
    gpr::PriorBase parameters;
    gpr::Field<double> lpg;
    double lpg_ref[3] = {499.999996357158e-003, 999.999209762188e-003,
                         999.999209762188e-003};

    parameters.setMu(0.);
    parameters.setS2(1.);
    parameters.setNu(20.);
    sexpat.setMagnSigma2(6.93874748072254e-009);

    sexpat.getLengthScaleRef().resize(1, 2);
    sexpat.setLengthScale(888.953211438594e-006);
    sexpat.setPriorParametersGaussian(parameters);
    sexpat.setPriorParametersSqrtt(parameters);

    lpg = sexpat.calculateLogPriorGradient();

    EXPECT_LE(fabs(lpg[0] - lpg_ref[0]), threshold)
        << "Value of `lpg[0]` is not equal to the expected one.";

    EXPECT_LE(fabs(lpg[1] - lpg_ref[1]), threshold)
        << "Value of `lpg[1]` is not equal to the expected one.";

    EXPECT_LE(fabs(lpg[2] - lpg_ref[2]), threshold)
        << "Value of `lpg[2]` is not equal to the expected one.";
}

TEST_F(SexpAtTest, sexpat_cfg)
{
    gpr::PriorBase parameters;
    std::vector<gpr::Field<double> > DKff;
    std::vector<gpr::Field<double> > DKff_ref;

    DKff_ref.resize(3);
    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n) {
        DKff_ref[n].resize(1, 4);
    }

    sexpat.getLengthScaleRef().resize(1, 2);
    sexpat.setLengthScale(888.953211438594e-006);

    sexpat.setMagnSigma2(6.93874748072254e-009);

    sexpat.setConfInfo(conf_info);

    parameters.setMu(0.);
    parameters.setS2(1.);
    parameters.setNu(20.);

    DKff_ref[0](0, 0) = 6.93874748072254e-009;
    DKff_ref[0](0, 1) = 351.124132742883e-015;
    DKff_ref[0](0, 2) = 351.124132742883e-015;
    DKff_ref[0](0, 3) = 6.93874748072254e-009;

    DKff_ref[1](0, 0) = 0.00000000000000e+000;
    DKff_ref[1](0, 1) = 14.2348374481591e-018;
    DKff_ref[1](0, 2) = 14.2348374481591e-018;
    DKff_ref[1](0, 3) = 0.00000000000000e+000;

    DKff_ref[2](0, 0) = 0.00000000000000e+000;
    DKff_ref[2](0, 1) = 6.94626888098018e-012;
    DKff_ref[2](0, 2) = 6.94626888098018e-012;
    DKff_ref[2](0, 3) = 0.00000000000000e+000;

    sexpat.setPriorParametersGaussian(parameters);
    sexpat.setPriorParametersSqrtt(parameters);

    sexpat.calculateGradOfCovMatrix(x1, x2, DKff);

    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n) {
        for (gpr::Index_t m = 0; m < DKff_ref[n].getSize(); ++m) {
            EXPECT_LE(fabs(DKff[n][m] - DKff_ref[n][m]), threshold)
                << "Value of `DKff['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }
}

TEST_F(SexpAtTest, sexpat_cfdg)
{
    gpr::Field<gpr::Index_t> dims;
    std::vector<gpr::Field<double> > DKff;
    std::vector<gpr::Field<double> > DKff_ref;

    dims.resize(1, 1);

    DKff_ref.resize(3);
    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n)
        DKff_ref[n].resize(2, 2);

    dims(0, 0) = 1;

    sexpat.getLengthScaleRef().resize(1, 2);
    sexpat.setLengthScale(888.953211438594e-006);

    sexpat.setMagnSigma2(6.93874748072254e-009);

    sexpat.setConfInfo(conf_info);

    DKff_ref[0](0, 0) = 0.00000000000000e+000;
    DKff_ref[0](0, 1) = -169.946842311334e-012;
    DKff_ref[0](1, 0) = 170.285367533359e-012;
    DKff_ref[0](1, 1) = 0.00000000000000e+000;

    DKff_ref[1](0, 0) = 0.00000000000000e+000;
    DKff_ref[1](0, 1) = -4.02676104048135e-012;
    DKff_ref[1](1, 0) = 4.02606171683247e-012;
    DKff_ref[1](1, 1) = 0.00000000000000e+000;

    DKff_ref[2](0, 0) = 0.00000000000000e+000;
    DKff_ref[2](0, 1) = -3.01813522425624e-009;
    DKff_ref[2](1, 0) = 3.02415591220194e-009;
    DKff_ref[2](1, 1) = 0.00000000000000e+000;

    sexpat.calculateGradOfCovMatrixWithDerivatives(x1, x2, dims, DKff);

    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n) {
        for (gpr::Index_t m = 0; m < DKff_ref[n].getSize(); ++m) {
            EXPECT_LE(fabs(DKff[n][m] - DKff_ref[n][m]), threshold)
                << "Value of `DKff['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }
}

TEST_F(SexpAtTest, sexpat_ginput4)
{
    gpr::Field<gpr::Index_t> dims;
    std::vector<gpr::Field<double> > DKff;
    std::vector<gpr::Field<double> > DKff_ref;

    dims.resize(1, 1);

    DKff_ref.resize(1);
    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n)
        DKff_ref[n].resize(2, 2);

    dims(0, 0) = 1;

    sexpat.getLengthScaleRef().resize(1, 2);
    sexpat.setLengthScale(888.953211438594e-006);

    sexpat.setMagnSigma2(6.93874748072254e-009);

    sexpat.setConfInfo(conf_info);

    DKff_ref[0](0, 0) = 0.00000000000000e+000;
    DKff_ref[0](0, 1) = -169.946842311334e-012;
    DKff_ref[0](1, 0) = 170.285367533359e-012;
    DKff_ref[0](1, 1) = 0.00000000000000e+000;

    sexpat.ginput4(x1, x2, dims, DKff);

    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n) {
        for (gpr::Index_t m = 0; m < DKff_ref[n].getSize(); ++m) {
            EXPECT_LE(fabs(DKff[n][m] - DKff_ref[n][m]), threshold)
                << "Value of `DKff['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }
}

TEST_F(SexpAtTest, sexpat_ginput2)
{
    gpr::Field<gpr::Index_t> dims;
    std::vector<gpr::Field<double> > DKff;
    std::vector<gpr::Field<double> > DKff_ref;

    dims.resize(1, 1);

    DKff_ref.resize(1);
    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n)
        DKff_ref[n].resize(2, 2);

    dims(0, 0) = 1;

    sexpat.getLengthScaleRef().resize(1, 2);
    sexpat.setLengthScale(888.953211438594e-006);

    sexpat.setMagnSigma2(6.93874748072254e-009);

    sexpat.setConfInfo(conf_info);

    DKff_ref[0](0, 0) = 6.55347337726538e-003;
    DKff_ref[0](0, 1) = 249.037058460746e-009;
    DKff_ref[0](1, 0) = 249.037058460746e-009;
    DKff_ref[0](1, 1) = 6.54671861564645e-003;

    sexpat.ginput2(x1, x2, dims, DKff);

    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n) {
        for (gpr::Index_t m = 0; m < DKff_ref[n].getSize(); ++m) {
            EXPECT_LE(fabs(DKff[n][m] - DKff_ref[n][m]), threshold)
                << "Value of `DKff['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }
}

TEST_F(SexpAtTest, sexpat_ginput3)
{
    gpr::Field<gpr::Index_t> dims1, dims2;
    std::vector<gpr::Field<double> > DKff;
    std::vector<gpr::Field<double> > DKff_ref;

    dims1.resize(1, 1);
    dims2.resize(1, 1);

    DKff_ref.resize(1);
    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n)
        DKff_ref[n].resize(2, 2);

    dims1(0, 0) = 1;
    dims2(0, 0) = 2;

    sexpat.getLengthScaleRef().resize(1, 2);
    sexpat.setLengthScale(888.953211438594e-006);

    sexpat.setMagnSigma2(6.93874748072254e-009);

    sexpat.setConfInfo(conf_info);

    DKff_ref[0](0, 0) = -72.4168570312556e-006;
    DKff_ref[0](0, 1) = -215.345156286196e-009;
    DKff_ref[0](1, 0) = -214.991981108543e-009;
    DKff_ref[0](1, 1) = -87.6234764908369e-006;

    sexpat.ginput3(x1, x2, dims1, dims2, DKff);

    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n) {
        for (gpr::Index_t m = 0; m < DKff_ref[n].getSize(); ++m) {
            EXPECT_LE(fabs(DKff[n][m] - DKff_ref[n][m]), threshold)
                << "Value of `DKff['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }
}

TEST_F(SexpAtTest, sexpat_cfdg2)
{
    gpr::Field<gpr::Index_t> dims1, dims2;
    std::vector<gpr::Field<double> > DKff;
    std::vector<gpr::Field<double> > DKff_ref;

    dims1.resize(1, 1);
    dims2.resize(1, 1);

    DKff_ref.resize(3);
    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n)
        DKff_ref[n].resize(2, 2);

    dims1(0, 0) = 1;
    dims2(0, 0) = 1;

    sexpat.getLengthScaleRef().resize(1, 2);
    sexpat.setLengthScale(888.953211438594e-006);

    sexpat.setMagnSigma2(6.93874748072254e-009);

    sexpat.setConfInfo(conf_info);

    DKff_ref[0](0, 0) = 6.55347337726538e-003;
    DKff_ref[0](0, 1) = 249.037058460746e-009;
    DKff_ref[0](1, 0) = 249.037058460746e-009;
    DKff_ref[0](1, 1) = 6.54671861564645e-003;

    DKff_ref[1](0, 0) = -11.2166246269546e-003;
    DKff_ref[1](0, 1) = -571.383262689728e-009;
    DKff_ref[1](1, 0) = -571.383262689729e-009;
    DKff_ref[1](1, 1) = -11.2126457536494e-003;

    DKff_ref[2](0, 0) = -1.89032212757614e-003;
    DKff_ref[2](0, 1) = 5.16484421469534e-006;
    DKff_ref[2](1, 0) = 5.16484421469534e-006;
    DKff_ref[2](1, 1) = -1.88079147764354e-003;

    sexpat.calculateGradOfCovMatrixWithDerivatives2(x1, x2, dims1, dims2, DKff);

    for (gpr::Index_t n = 0; n < DKff_ref.size(); ++n) {
        for (gpr::Index_t m = 0; m < DKff_ref[n].getSize(); ++m) {
            EXPECT_LE(fabs(DKff[n][m] - DKff_ref[n][m]), threshold)
                << "Value of `DKff['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }
}

TEST_F(SexpAtTest, sexpat_cov)
{
    gpr::PriorBase parameters;
    gpr::Field<double> C;
    gpr::Field<double> C_ref;

    C_ref.resize(2, 2);

    sexpat.getLengthScaleRef().resize(1, 2);
    sexpat.getLengthScaleRef()(0, 0) = 1.00138353271710e+000;
    sexpat.getLengthScaleRef()(0, 1) = 20.5486885227074e-003;

    sexpat.setMagnSigma2(2.09859544785255e-006);

    sexpat.setConfInfo(conf_info);

    C_ref[0] = 2.09859544785255e-006;
    C_ref[1] = 2.06010387765538e-006;
    C_ref[2] = 2.06010387765538e-006;
    C_ref[3] = 2.09859544785255e-006;

    sexpat.calculateCovarianceMatrix(x1, x2, C);

    for (gpr::Index_t n = 0; n < C_ref.getSize(); ++n) {
        EXPECT_LE(fabs(C[n] - C_ref[n]), threshold)
            << "Value of `C['" << n << "]` is not equal to the expected one.";
    }
}

} /* namespace tests */
} /* namespace gpr */
