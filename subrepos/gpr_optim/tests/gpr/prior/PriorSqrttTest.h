/*
 * PriorSqrttTest.h
 *
 *  Created on: 13 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_PRIORSQRTTTEST_H_
#define TESTS_PRIORSQRTTTEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/prior/PriorSqrtt.h"

namespace gpr {
namespace tests {

class PriorSqrttTest : public ::testing::Test {
protected:
    PriorSqrttTest();
    ~PriorSqrttTest();

    gpr::PriorSqrtt prior;

    double threshold;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_PRIORSQRTTTEST_H_ */
