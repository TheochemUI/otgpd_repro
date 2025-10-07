/*
 * GPCFTest.h
 *
 *  Created on: 14 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_SEXPATTEST_H_
#define TESTS_SEXPATTEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/covariance_functions/SexpatCF.h"
#include "../../../managers/io/FileManager.h"

namespace gpr {
namespace tests {

class SexpAtTest : public ::testing::Test {
protected:
    SexpAtTest();
    ~SexpAtTest();

    gpr::SexpatCF sexpat;
    io::FileManager io_manager;
    gpr::Coord x1, x2;
    gpr::AtomsConfiguration conf_info;
    double threshold;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_SEXPATTEST_H_ */
