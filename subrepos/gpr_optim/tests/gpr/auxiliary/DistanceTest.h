/*
 * DistanceTest.h
 *
 *  Created on: 30 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_DISTANCETEST_H_
#define TESTS_DISTANCETEST_H_

#include <gtest/gtest.h>

#include "../../../gpr/auxiliary/Distance.h"
#include "../../../managers/io/FileManager.h"

namespace gpr {
namespace tests {

class DistanceTest : public ::testing::Test {
public:
    DistanceTest();
    virtual ~DistanceTest();

    aux::Distance distance;
    gpr::AtomsConfiguration conf_info;
    io::FileManager io_manager;

    double threshold;  // Epsilon for comparison operators.
protected:
    void setup_local_conf(AtomsConfiguration& conf, Index_t num_atoms,
                          const std::vector<Index_t>& types);
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_DISTANCETEST_H_ */
