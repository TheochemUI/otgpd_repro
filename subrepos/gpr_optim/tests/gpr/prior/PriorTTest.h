/*
 * PriorTTest.h
 *
 *  Created on: 13 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_PRIORTTEST_H_
#define TESTS_PRIORTTEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/prior/PriorT.h"

namespace gpr {
namespace tests {

class PriorTTest : public ::testing::Test {
protected:
    PriorTTest();
    ~PriorTTest();

    gpr::PriorT prior;

    double threshold;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_PRIORTTEST_H_ */
