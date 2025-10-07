/*
 * LBFGSTest.h
 *
 *  Created on: 30 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_LBFGSTEST_H_
#define TESTS_LBFGSTEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/dimer/Dimer.h"

namespace gpr {
namespace tests {

class DimerTest : public ::testing::Test {
public:
    DimerTest();
    virtual ~DimerTest();

    dimer::Dimer dimer;

    double threshold;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_LBFGSTEST_H_ */
