/*
 * GaussianProcessRegressionTest.h
 *
 *  Created on: 16 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_GAUSSIANPROCESSREGRESSIONTEST_H_
#define TESTS_GAUSSIANPROCESSREGRESSIONTEST_H_

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include "../../../gpr/ml/GaussianProcessRegression.h"
#include "../../../managers/io/FileManager.h"

namespace gpr {
namespace tests {

class GaussianProcessRegressionTest : public gpr::GaussianProcessRegression,
                                      public ::testing::Test {
public:
    GaussianProcessRegressionTest();
    virtual ~GaussianProcessRegressionTest();

    gpr::EigenMatrix x1;
    gpr::EigenMatrix x2;
    Eigen::VectorXd x1_ind;
    Eigen::VectorXd x2_ind;
    gpr::AtomsConfiguration conf_info;
    io::FileManager io_manager;

    double threshold;  // Epsilon for comparison operators.
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_GAUSSIANPROCESSREGRESSIONTEST_H_ */
