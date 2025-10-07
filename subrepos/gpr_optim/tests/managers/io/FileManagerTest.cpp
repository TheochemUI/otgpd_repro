/*
 * FileManagerTest.cpp
 *
 *  Created on: 17 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "FileManagerTest.h"

#include <cmath>

namespace gpr {
namespace tests {

FileManagerTest::FileManagerTest()
{
    // TODO Auto-generated constructor stub
}

FileManagerTest::~FileManagerTest()
{
    // TODO Auto-generated destructor stub
}

TEST_F(FileManagerTest, ReadInputFile)
{
    gpr::InputParameters parameters;

    io_manager.readInputFile("input/input.dat", parameters);

    EXPECT_EQ(parameters.i_dist.value, 0)
        << "Parameter 'i_dist' from the input file was not read correctly.";

    EXPECT_LE(std::fabs(parameters.dist_sp.value[0] - 0.02), DBL_EPSILON)
        << "Parameter 'dist_sp[0]' from the input file was not read correctly.";
    EXPECT_LE(std::fabs(parameters.dist_sp.value[1] - 0.3), DBL_EPSILON)
        << "Parameter 'dist_sp[1]' from the input file was not read correctly.";
    EXPECT_LE(std::fabs(parameters.dist_sp.value[2] - 1.), DBL_EPSILON)
        << "Parameter 'dist_sp[2]' from the input file was not read correctly.";
    EXPECT_LE(std::fabs(parameters.dist_sp.value[3] - 3.), DBL_EPSILON)
        << "Parameter 'dist_sp[3]' from the input file was not read correctly.";

    EXPECT_EQ(parameters.method_rot.value, "LBFGS_alg")
        << "Parameter 'method_rot' from the input file was not read correctly.";

    EXPECT_LE(std::fabs(parameters.ratio_at_limit.value - 0.666666666666667),
              DBL_EPSILON)
        << "Parameter 'ratio_at_limit' from the input file was not read "
           "correctly.";

    EXPECT_EQ(parameters.num_bigiter.value, 300)
        << "Parameter 'num_bigiter' from the input file was not read "
           "correctly.";

    EXPECT_EQ(parameters.start_prune_at.value, 8)
        << "Parameter 'start_prune_at' from the input file was not read "
           "correctly.";

    EXPECT_EQ(parameters.nprune_vals.value, 3)
        << "Parameter 'nprune_vals' from the input file was not read "
           "correctly.";

    EXPECT_EQ(parameters.prune_threshold.value, 0.3)
        << "Parameter 'prune_threshold' from the input file was not read "
           "correctly.";
}

TEST_F(FileManagerTest, ReadDataFile)
{
    gpr::Field<double> orient_init, orient_start;
    std::map<std::string, gpr::Field<double>*> dict;
    std::vector<double> orient_init_ref(60);
    std::vector<double> orient_start_ref(60);

    dict["orient_init"] = &orient_init;
    dict["orient_start"] = &orient_start;

    io_manager.readDataFile("input/orients_CuH2.dat", dict);

    // Reference values
    if (!io_manager.readFromPlainFile(
            "tests/reference/managers/io/ReferenceOrientInit.dat",
            orient_init_ref)) {
        return;
    }

    if (!io_manager.readFromPlainFile(
            "tests/reference/managers/io/ReferenceStartInit.dat",
            orient_start_ref)) {
        return;
    }

    for (std::pair<std::string, gpr::Field<double>*> elt: dict) {
        gpr::Index_t counter = 0;
        for (gpr::Index_t n = 0; n < elt.second->getNumRows(); ++n) {
            for (gpr::Index_t m = 0; m < elt.second->getNumCols(); ++m) {
                double value = elt.second->operator()(n, m);
                double ref_value = 0.;

                if (elt.first == "orient_init") {
                    ref_value = orient_init_ref[counter];
                }
                if (elt.first == "orient_start") {
                    ref_value = orient_start_ref[counter];
                }
                ++counter;

                EXPECT_LE(std::fabs(value - ref_value), DBL_EPSILON)
                    << "Field " << elt.first
                    << " was not read correctly from "
                       "the data file ";
            }
        }
    }
}

} /* namespace tests */
} /* namespace gpr */
