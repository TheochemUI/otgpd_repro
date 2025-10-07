/*
 * PriorLogUnifTest.cpp
 *
 *  Created on: 13 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "PriorLogUnifTest.h"

#include <cmath>

namespace gpr {
namespace tests {

PriorLogUnifTest::PriorLogUnifTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

PriorLogUnifTest::~PriorLogUnifTest() { }

TEST_F(PriorLogUnifTest, lp)
{
    gpr::Field<double> x;
    double lp = 0.;
    double lp_ref = 14.0509319087763e+000;

    x.resize(1, 2);

    prior.setMu(0.);
    prior.setNu(20.);
    prior.setS2(1.);

    x(0, 0) = 888.953211438594e-006;
    x(0, 1) = 888.953211438594e-006;

    lp = prior.calculateLogPrior(x);

    EXPECT_LE(fabs(lp - lp_ref), threshold)
        << "Value of `lp` is not equal to the expected one.";
}

TEST_F(PriorLogUnifTest, lpg)
{
    gpr::Field<double> x;
    gpr::Field<double> lpg;
    double lpg_ref[2] = {-1.12491859766354e+003, -1.12491859766354e+003};

    x.resize(1, 2);

    prior.setMu(0.);
    prior.setNu(20.);
    prior.setS2(1.);

    x(0, 0) = 888.953211438594e-006;
    x(0, 1) = 888.953211438594e-006;

    lpg = prior.calculateLogPriorGradient(x);

    EXPECT_LE(fabs(lpg[0] - lpg_ref[0]), threshold * 1e3)
        << "Value of `lpg[0]` is not equal to the expected one.";
    EXPECT_LE(fabs(lpg[1] - lpg_ref[1]), threshold * 1e3)
        << "Value of `lpg[1]` is not equal to the expected one.";
}

} /* namespace tests */
} /* namespace gpr */
