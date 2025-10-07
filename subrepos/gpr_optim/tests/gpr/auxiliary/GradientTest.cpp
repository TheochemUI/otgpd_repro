/*
 * GradientTest.cpp
 *
 *  Created on: 23 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "GradientTest.h"

#include "../../../gpr/Enums.h"

namespace gpr {
namespace tests {

GradientTest::GradientTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.

    conf_info.atoms_froz_active.positions.resize(1, 26 * 3);
    conf_info.atoms_froz_active.type.resize(1, 26);
    conf_info.atoms_mov.type.resize(1, 2);
    conf_info.pairtype.resize(2, 2);

    io_manager.readFromPlainFile(
        "tests/reference/gpr/auxiliary/ConfigurationFrozen.dat",
        conf_info.atoms_froz_active.positions.getInternalVector());

    conf_info.atoms_mov.type.set(0);
    conf_info.atoms_froz_active.type.set(1);

    conf_info.pairtype(0, 0) = 0;
    conf_info.pairtype(0, 1) = 1;
    conf_info.pairtype(1, 0) = 1;
    conf_info.pairtype(1, 1) = EMPTY;

    conf_info.n_pt = 2;
}

GradientTest::~GradientTest() { }

TEST_F(GradientTest, calculateDerivativesX1X2DiffDimV_vec)
{
    gpr::Coord x1, x2;
    gpr::Field<gpr::Index_t> dims1, dims2;
    gpr::Derivatives<std::vector<gpr::Field<double> > > derivatives;
    gpr::Derivatives<std::vector<gpr::Field<double> > > derivatives_pt;
    gpr::Derivatives<std::vector<gpr::Field<double> > > derivatives_ref;
    gpr::Derivatives<std::vector<gpr::Field<double> > > derivatives_pt_ref;
    gpr::Field<double> lengthScale;
    uint8_t calc_options =
        OptionsForGradCalculation::D1 | OptionsForGradCalculation::D2 |
        OptionsForGradCalculation::D12 | OptionsForGradCalculation::D1_pt |
        OptionsForGradCalculation::D2_pt | OptionsForGradCalculation::D12_pt;

    x1.resize(2, 2 * 3);
    x2.resize(2, 2 * 3);

    dims1.resize(1, 1);
    dims2.resize(1, 1);

    lengthScale.resize(1, 2);

    derivatives_ref.D1.resize(1);
    derivatives_ref.D2.resize(1);
    derivatives_ref.D12.resize(1);
    derivatives_pt_ref.D1.resize(2);
    derivatives_pt_ref.D2.resize(2);
    derivatives_pt_ref.D12.resize(2);

    for (gpr::Index_t n = 0; n < derivatives_ref.D1.size(); ++n) {
        derivatives_ref.D1[n].resize(2, 2);
        derivatives_ref.D2[n].resize(2, 2);
        derivatives_ref.D12[n].resize(2, 2);
    }

    for (gpr::Index_t n = 0; n < derivatives_pt_ref.D1.size(); ++n) {
        derivatives_pt_ref.D1[n].resize(2, 2);
        derivatives_pt_ref.D2[n].resize(2, 2);
        derivatives_pt_ref.D12[n].resize(2, 2);
    }

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

    dims1(0, 0) = 1;
    dims2(0, 0) = 4;

    lengthScale.set(888.953211438594e-006);

    derivatives_ref.D1[0](0, 0) = 0.00000000000000e+000;
    derivatives_ref.D1[0](0, 1) = 968.015732691209e+000;
    derivatives_ref.D1[0](1, 0) = -969.943969405563e+000;
    derivatives_ref.D1[0](1, 1) = 0.00000000000000e+000;

    derivatives_ref.D2[0](0, 0) = 0.00000000000000e+000;
    derivatives_ref.D2[0](0, 1) = -937.586882744879e+000;
    derivatives_ref.D2[0](1, 0) = 936.714863624322e+000;
    derivatives_ref.D2[0](1, 1) = 0.00000000000000e+000;

    derivatives_ref.D12[0](0, 0) = 1.61652007918101e+006;
    derivatives_ref.D12[0](0, 1) = 1.61623333966243e+006;
    derivatives_ref.D12[0](1, 0) = 1.61623333966243e+006;
    derivatives_ref.D12[0](1, 1) = 1.61594665100592e+006;

    derivatives_pt_ref.D1[0](0, 0) = 0.00000000000000e+000;
    derivatives_pt_ref.D1[0](0, 1) = 11.4485758520161e+000;
    derivatives_pt_ref.D1[0](1, 0) = -11.4465450952254e+000;
    derivatives_pt_ref.D1[0](1, 1) = 0.00000000000000e+000;

    derivatives_pt_ref.D1[1](0, 0) = 0.00000000000000e+000;
    derivatives_pt_ref.D1[1](0, 1) = -979.464308543225e+000;
    derivatives_pt_ref.D1[1](1, 0) = 981.390514500789e+000;
    derivatives_pt_ref.D1[1](1, 1) = 0.00000000000000e+000;

    derivatives_pt_ref.D2[0](0, 0) = 0.00000000000000e+000;
    derivatives_pt_ref.D2[0](0, 1) = 11.4465450952254e+000;
    derivatives_pt_ref.D2[0](1, 0) = -11.4485758520161e+000;
    derivatives_pt_ref.D2[0](1, 1) = 0.00000000000000e+000;

    derivatives_pt_ref.D2[1](0, 0) = 0.00000000000000e+000;
    derivatives_pt_ref.D2[1](0, 1) = 926.140337649654e+000;
    derivatives_pt_ref.D2[1](1, 0) = -925.266287772305e+000;
    derivatives_pt_ref.D2[1](1, 1) = 0.00000000000000e+000;

    derivatives_pt_ref.D12[0](0, 0) = -1.61652007918101e+006;
    derivatives_pt_ref.D12[0](0, 1) = -1.61623333966243e+006;
    derivatives_pt_ref.D12[0](1, 0) = -1.61623333966243e+006;
    derivatives_pt_ref.D12[0](1, 1) = -1.61594665100592e+006;

    derivatives_pt_ref.D12[1](0, 0) = 0.00000000000000e+000;
    derivatives_pt_ref.D12[1](0, 1) = 0.00000000000000e+000;
    derivatives_pt_ref.D12[1](1, 0) = 0.00000000000000e+000;
    derivatives_pt_ref.D12[1](1, 1) = 0.00000000000000e+000;

    gradient.calculateDerivativesDiffDim(x1, x2, conf_info, lengthScale, dims1,
                                         dims2, calc_options, Directions::x,
                                         derivatives, derivatives_pt);

    for (gpr::Index_t n = 0; n < derivatives_ref.D1.size(); ++n) {
        for (gpr::Index_t m = 0; m < derivatives_ref.D1[n].getSize(); ++m) {
            EXPECT_LE(fabs(derivatives.D1[n][m] - derivatives_ref.D1[n][m]),
                      threshold * 1e3)
                << "Value of `D1['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }

    for (gpr::Index_t n = 0; n < derivatives_ref.D2.size(); ++n) {
        for (gpr::Index_t m = 0; m < derivatives_ref.D2[n].getSize(); ++m) {
            EXPECT_LE(fabs(derivatives.D2[n][m] - derivatives_ref.D2[n][m]),
                      threshold * 1e3)
                << "Value of `D2['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }

    for (gpr::Index_t n = 0; n < derivatives_ref.D12.size(); ++n) {
        for (gpr::Index_t m = 0; m < derivatives_ref.D12[n].getSize(); ++m) {
            EXPECT_LE(fabs(derivatives.D12[n][m] - derivatives_ref.D12[n][m]),
                      threshold * 1e5)
                << "Value of `D12['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }

    for (gpr::Index_t n = 0; n < derivatives_pt_ref.D1.size(); ++n) {
        for (gpr::Index_t m = 0; m < derivatives_pt_ref.D1[n].getSize(); ++m) {
            EXPECT_LE(
                fabs(derivatives_pt.D1[n][m] - derivatives_pt_ref.D1[n][m]),
                threshold * 1e3)
                << "Value of `D1_pt['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }

    for (gpr::Index_t n = 0; n < derivatives_pt_ref.D2.size(); ++n) {
        for (gpr::Index_t m = 0; m < derivatives_pt_ref.D2[n].getSize(); ++m) {
            EXPECT_LE(
                fabs(derivatives_pt.D2[n][m] - derivatives_pt_ref.D2[n][m]),
                threshold * 1e3)
                << "Value of `D2_pt['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }

    for (gpr::Index_t n = 0; n < derivatives_pt_ref.D12.size(); ++n) {
        for (gpr::Index_t m = 0; m < derivatives_pt_ref.D12[n].getSize(); ++m) {
            EXPECT_LE(
                fabs(derivatives_pt.D12[n][m] - derivatives_pt_ref.D12[n][m]),
                threshold * 1e5)
                << "Value of `D12_pt['" << n << "][" << m
                << "]` is not equal to the expected one.";
        }
    }
}

} /* namespace tests */
} /* namespace gpr */
