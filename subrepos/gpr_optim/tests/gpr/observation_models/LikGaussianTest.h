/*
 * LikGaussianTest.h
 *
 *  Created on: 23 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_LIKGAUSSIANTEST_H_
#define TESTS_LIKGAUSSIANTEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/observation_models/LikGaussian.h"

namespace gpr {
namespace tests {

class LikGaussianTest : public ::testing::Test {
public:
    LikGaussianTest();
    virtual ~LikGaussianTest();

    gpr::LikGaussian lik_gaussian;

    double threshold;  // Epsilon for comparison operators.
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_LIKGAUSSIANTEST_H_ */
