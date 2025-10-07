/*
 * PriorLogUnifTest.h
 *
 *  Created on: 13 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_PRIORLOGUNIFTEST_H_
#define TESTS_PRIORLOGUNIFTEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/prior/PriorLogUnif.h"

namespace gpr {
namespace tests {

class PriorLogUnifTest : public ::testing::Test {
protected:
    PriorLogUnifTest();
    ~PriorLogUnifTest();

    gpr::PriorLogUnif prior;

    double threshold;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_PRIORLOGUNIFTEST_H_ */
