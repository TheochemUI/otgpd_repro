/*
 * GradientTest.h
 *
 *  Created on: 23 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_GRADIENTTEST_H_
#define TESTS_GRADIENTTEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/auxiliary/Gradient.h"
#include "../../../managers/io/FileManager.h"

namespace gpr {
namespace tests {

class GradientTest : public ::testing::Test {
public:
    GradientTest();
    virtual ~GradientTest();

    aux::Gradient gradient;
    gpr::AtomsConfiguration conf_info;
    io::FileManager io_manager;

    double threshold;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_GRADIENTTEST_H_ */
