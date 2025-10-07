/*
 * AtomicDimerTest.h
 *
 *  Created on: 18 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_GPR_ATOMICDIMERTEST_H_
#define TESTS_GPR_ATOMICDIMERTEST_H_

#include <gtest/gtest.h>

#include "../../gpr/AtomicDimer.h"

namespace gpr {
namespace tests {

class AtomicDimerTest : public ::testing::Test {
public:
    AtomicDimerTest();
    virtual ~AtomicDimerTest();

    atmd::AtomicDimer dimer;

    double threshold;  // Epsilon for comparison operators.
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_GPR_ATOMICDIMERTEST_H_ */
