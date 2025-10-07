/*
 * LikGaussianTest.cpp
 *
 *  Created on: 23 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "LikGaussianTest.h"

namespace gpr {
namespace tests {

LikGaussianTest::LikGaussianTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

LikGaussianTest::~LikGaussianTest() { }

TEST_F(LikGaussianTest, trcov)
{
    gpr::LikGaussian lik_gaussian;
    gpr::Coord x;
    gpr::Field<double> C;
    double sigma2 = 10.0000000000000e-009;

    x.resize(14, 2 * 3);

    lik_gaussian.setSigma2(sigma2);

    lik_gaussian.evaluateTrainingCovarianceMatrix(x, C);

    for (gpr::Index_t i = 0; i < C.getNumRows(); ++i) {
        for (gpr::Index_t j = 0; j < C.getNumCols(); ++j) {
            if (i == j) {
                EXPECT_LE(fabs(C(i, j) - sigma2), threshold)
                    << "Value of `C` is not equal to the expected one.";
            } else {
                EXPECT_LE(fabs(C(i, j)), threshold)
                    << "Value of `C` is not equal to the expected one.";
            }
        }
    }
}

} /* namespace tests */
} /* namespace gpr */
