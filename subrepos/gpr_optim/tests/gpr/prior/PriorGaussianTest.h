/*
 * PriorGaussianTest.h
 *
 *  Created on: 13 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_PRIORGAUSSIANTEST_H_
#define TESTS_PRIORGAUSSIANTEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/prior/PriorGaussian.h"

namespace gpr {
namespace tests {

class PriorGaussianTest : public ::testing::Test {
protected:
    PriorGaussianTest();
    ~PriorGaussianTest();

    gpr::PriorGaussian prior;

    double threshold;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_PRIORGAUSSIANTEST_H_ */
