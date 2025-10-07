/*
 * GPDimerTest.h
 *
 *  Created on: 24 Sep 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_GPDIMERTEST_H_
#define TESTS_GPDIMERTEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/auxiliary/ProblemSetUp.h"
#include "../../../managers/io/FileManager.h"

namespace gpr {
namespace tests {

class ProblemSetUpTest : public ::testing::Test {
public:
    ProblemSetUpTest();
    virtual ~ProblemSetUpTest();

    aux::ProblemSetUp setup;
    io::FileManager io_manager;

    double threshold;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_GPDIMERTEST_H_ */
