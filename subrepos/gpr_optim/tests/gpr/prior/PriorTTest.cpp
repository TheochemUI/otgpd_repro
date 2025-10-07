/*
 * PriorTTest.cpp
 *
 *  Created on: 13 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "PriorTTest.h"

#include <cmath>

namespace gpr {
namespace tests {

PriorTTest::PriorTTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

PriorTTest::~PriorTTest() { }

TEST_F(PriorTTest, lp)
{
    gpr::Field<double> x;
    double lp = 0.;
    double lp_ref = -931.433340379401e-003;

    x.resize(1, 1);

    prior.setMu(0.);
    prior.setNu(20.);
    prior.setS2(1.);

    x(0, 0) = 6.93874748072254e-009;

    lp = prior.calculateLogPrior(x);

    EXPECT_LE(fabs(lp - lp_ref), threshold)
        << "Value of `lp` is not equal to the expected one.";
}

TEST_F(PriorTTest, lpg)
{
    gpr::Field<double> x;
    gpr::Field<double> lpg;
    double lpg_ref = -7.28568485475866e-009;

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
