/*
 * PriorSqrttTest.cpp
 *
 *  Created on: 13 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "PriorSqrttTest.h"

#include <cmath>

namespace gpr {
namespace tests {

PriorSqrttTest::PriorSqrttTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

PriorSqrttTest::~PriorSqrttTest() { }

TEST_F(PriorSqrttTest, lp)
{
    gpr::Field<double> x;
    double lp = 0.;
    double lp_ref = 7.76849175391677e+000;

    x.resize(1, 1);

    prior.setMu(0.);
    prior.setNu(20.);
    prior.setS2(1.);

    x(0, 0) = 6.93874748072254e-009;

    lp = prior.calculateLogPrior(x);

    EXPECT_LE(fabs(lp - lp_ref), threshold)
        << "Value of `lp` is not equal to the expected one.";
}

TEST_F(PriorSqrttTest, lpg)
{
    gpr::Field<double> x;
    gpr::Field<double> lpg;
    double lpg_ref = -72.0591151402986e+006;

    x.resize(1, 1);

    prior.setMu(0.);
    prior.setNu(20.);
    prior.setS2(1.);

    x(0, 0) = 6.93874748072254e-009;

    lpg = prior.calculateLogPriorGradient(x);

    EXPECT_LE(fabs(lpg[0] - lpg_ref), threshold * 1e5)
        << "Value of `lpg[0]` is not equal to the expected one.";
}

} /* namespace tests */
} /* namespace gpr */
