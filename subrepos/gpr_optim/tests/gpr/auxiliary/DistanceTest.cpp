/*
 * DistanceTest.cpp
 *
 *  Created on: 30 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "DistanceTest.h"

#include <vector>

#include "data_types/Vector3_reg.h"

#ifdef USE_HIGHS
#include "Highs.h"
#endif

namespace gpr {
namespace tests {
using aux::calculate_aligned_rmsd;

DistanceTest::DistanceTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.

    conf_info.atoms_froz_active.positions.resize(1, 26 * 3);
    conf_info.atoms_froz_active.type.resize(1, 26);
    conf_info.atoms_mov.type.resize(1, 2);
    conf_info.pairtype.resize(2, 2);

    io_manager.readFromPlainFile(
        "tests/reference/gpr/auxiliary/ConfigurationFrozen.dat",
        conf_info.atoms_froz_active.positions.getInternalVector());

    conf_info.atoms_mov.type.set(0);
    conf_info.atoms_froz_active.type.set(1);

    conf_info.pairtype(0, 0) = 0;
    conf_info.pairtype(0, 1) = 1;
    conf_info.pairtype(1, 0) = 1;
    conf_info.pairtype(1, 1) = EMPTY;

    conf_info.n_pt = 2;
}

void DistanceTest::setup_local_conf(AtomsConfiguration& conf, Index_t num_atoms,
                                    const std::vector<Index_t>& types)
{
    conf.atoms_mov.type.resize(1, num_atoms);
    for (Index_t i = 0; i < num_atoms; ++i) {
        conf.atoms_mov.type(0, i) = types[i];
    }
}

DistanceTest::~DistanceTest() { }

TEST_F(DistanceTest, DistMax1Dlog)
{
    gpr::Coord x1, x2;
    gpr::Field<double> dist;
    std::vector<double> dist_ref;

    x1.resize(1, 2 * 3);
    x2.resize(6, 2 * 3);
    dist_ref.resize(6);

    x1(0, 0) = 8.984799577839603;  // .x
    x1(0, 1) = 9.946773537381196;  // .y
    x1(0, 2) = 7.883158760864689;  // .z
    x1(0, 3) = 7.648053020353811;  // .x
    x1(0, 4) = 9.947018098632240;  // .y
    x1(0, 5) = 7.884341447664946;  // .z

    x2(0, 0) = 8.982373164830571;  // .x
    x2(0, 1) = 9.937230835772036;  // .y
    x2(0, 2) = 7.894416323850492;  // .z
    x2(0, 3) = 7.652483227274959;  // .x
    x2(0, 4) = 9.955905494573976;  // .y
    x2(0, 5) = 7.877879589983660;  // .z
                                   //
    x2(1, 0) = 8.978562773030582;  // .x
    x2(1, 1) = 9.932116280672286;  // .y
    x2(1, 2) = 7.898827614144257;  // .z
    x2(1, 3) = 7.648887496635559;  // .x
    x2(1, 4) = 9.955170515128863;  // .y
    x2(1, 5) = 7.872742150466696;  // .z
                                   //
    x2(2, 0) = 8.976993347159784;  // .x
    x2(2, 1) = 9.935874361177897;  // .y
    x2(2, 2) = 7.894533919329630;  // .z
    x2(2, 3) = 7.644176926015588;  // .x
    x2(2, 4) = 9.955573130196294;  // .y
    x2(2, 5) = 7.878193601315998;  // .z
                                   //
    x2(3, 0) = 8.982700234844764;  // .x
    x2(3, 1) = 9.938213927286730;  // .y
    x2(3, 2) = 7.892109355815186;  // .z
    x2(3, 3) = 7.643018372071518;  // .x
    x2(3, 4) = 9.956281713287362;  // .y
    x2(3, 5) = 7.875909963059283;  // .z
                                   //
    x2(4, 0) = 8.989491283507103;  // .x
    x2(4, 1) = 9.936555712126957;  // .y
    x2(4, 2) = 7.892532513641982;  // .z
    x2(4, 3) = 7.645972793452610;  // .x
    x2(4, 4) = 9.955825824010615;  // .y
    x2(4, 5) = 7.876166184685688;  // .z
                                   //
    x2(5, 0) = 8.983330780846810;  // .x
    x2(5, 1) = 9.944154132051644;  // .y
    x2(5, 2) = 7.875729042237787;  // .z
    x2(5, 3) = 7.650585914617185;  // .x
    x2(5, 4) = 9.929441987442432;  // .y
    x2(5, 5) = 7.882144032513800;  // .z

    dist_ref = {0.005675600658215, 0.008024540882219, 0.005889796208082,
                0.004873773574595, 0.006608214559330, 0.007841185655827};

    distance.dist_max1Dlog(x1, x2, conf_info, dist);

    for (gpr::Index_t n = 0; n < dist.getSize(); ++n) {
        EXPECT_LE(std::fabs(dist(0, n) - dist_ref[n]), threshold)
            << "Elements of the matrix are not equal to the expected ones.";
    }
}

TEST_F(DistanceTest, DistAt)
{
    gpr::Coord x1;
    gpr::Coord x2;
    gpr::Field<double> lengthscale;
    gpr::Field<double> dist;
    gpr::Field<double> dist_ref(6, 6);

    x1.resize(6, 2 * 3);
    x2.resize(6, 2 * 3);
    dist.resize(5, 5);
    lengthscale.resize(1, 1);

    x1.set(0, 0, {8.982373164830571, 9.937230835772036, 7.894416323850492});
    x1.set(1, 0, {8.978562773030582, 9.932116280672286, 7.898827614144257});
    x1.set(2, 0, {8.976993347159784, 9.935874361177897, 7.894533919329630});
    x1.set(3, 0, {8.982700234844764, 9.938213927286730, 7.892109355815186});
    x1.set(4, 0, {8.989491283507103, 9.936555712126957, 7.892532513641982});
    x1.set(5, 0, {8.983330780846810, 9.944154132051644, 7.875729042237787});

    x1.set(0, 1, {7.652483227274959, 9.955905494573976, 7.877879589983660});
    x1.set(1, 1, {7.648887496635559, 9.955170515128863, 7.872742150466696});
    x1.set(2, 1, {7.644176926015588, 9.955573130196294, 7.878193601315998});
    x1.set(3, 1, {7.643018372071518, 9.956281713287362, 7.875909963059283});
    x1.set(4, 1, {7.645972793452610, 9.955825824010615, 7.876166184685688});
    x1.set(5, 1, {7.650585914617185, 9.929441987442432, 7.882144032513800});

    x2.set(0, 0, {8.982373164830571, 9.937230835772036, 7.894416323850492});
    x2.set(1, 0, {8.978562773030582, 9.932116280672286, 7.898827614144257});
    x2.set(2, 0, {8.976993347159784, 9.935874361177897, 7.894533919329630});
    x2.set(3, 0, {8.982700234844764, 9.938213927286730, 7.892109355815186});
    x2.set(4, 0, {8.989491283507103, 9.936555712126957, 7.892532513641982});
    x2.set(5, 0, {8.983330780846810, 9.944154132051644, 7.875729042237787});

    x2.set(0, 1, {7.652483227274959, 9.955905494573976, 7.877879589983660});
    x2.set(1, 1, {7.648887496635559, 9.955170515128863, 7.872742150466696});
    x2.set(2, 1, {7.644176926015588, 9.955573130196294, 7.878193601315998});
    x2.set(3, 1, {7.643018372071518, 9.956281713287362, 7.875909963059283});
    x2.set(4, 1, {7.645972793452610, 9.955825824010615, 7.876166184685688});
    x2.set(5, 1, {7.650585914617185, 9.929441987442432, 7.882144032513800});

    lengthscale(0, 0) = 1.;

    dist_ref(0, 0) = 0.;
    dist_ref(0, 1) = 0.003953892010280;
    dist_ref(0, 2) = 0.004070202425876;
    dist_ref(0, 3) = 0.008460855559739;
    dist_ref(0, 4) = 0.011285766759530;
    dist_ref(0, 5) = 0.014672259642108;
    dist_ref(1, 0) = 0.003953892010280;
    dist_ref(1, 1) = 0.;
    dist_ref(1, 2) = 0.004423202945557;
    dist_ref(1, 3) = 0.009080320379056;
    dist_ref(1, 4) = 0.011918448534057;
    dist_ref(1, 5) = 0.016476334213804;
    dist_ref(2, 0) = 0.004070202425876;
    dist_ref(2, 1) = 0.004423202945557;
    dist_ref(2, 2) = 0.;
    dist_ref(2, 3) = 0.005989481708574;
    dist_ref(2, 4) = 0.009472766608845;
    dist_ref(2, 5) = 0.014771352424969;
    dist_ref(3, 0) = 0.008460855559739;
    dist_ref(3, 1) = 0.009080320379056;
    dist_ref(3, 2) = 0.005989481708574;
    dist_ref(3, 3) = 0.;
    dist_ref(3, 4) = 0.003964897696479;
    dist_ref(3, 5) = 0.015433995102774;
    dist_ref(4, 0) = 0.011285766759530;
    dist_ref(4, 1) = 0.011918448534057;
    dist_ref(4, 2) = 0.009472766608845;
    dist_ref(4, 3) = 0.003964897696479;
    dist_ref(4, 4) = 0.;
    dist_ref(4, 5) = 0.016861382029173;
    dist_ref(5, 0) = 0.014672259642108;
    dist_ref(5, 1) = 0.016476334213804;
    dist_ref(5, 2) = 0.014771352424969;
    dist_ref(5, 3) = 0.015433995102774;
    dist_ref(5, 4) = 0.016861382029173;
    dist_ref(5, 5) = 0.;

    distance.dist_at(x1, x2, conf_info, lengthscale, dist);

    for (gpr::Index_t n = 0; n < dist.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < dist.getNumCols(); ++m) {
            EXPECT_LE(std::fabs(dist(n, m) - dist_ref(n, m)), threshold)
                << "Elements of the matrix are not equal to the expected ones.";
        }
    }
}

TEST_F(DistanceTest, MinDistInteratomic)
{
    gpr::Coord x;
    gpr::Field<double> dist;
    gpr::Field<double> dist_ref(1, 2);

    x.resize(1, 3 * 2);

    x.set(0, 0, {8.982538896311098, 9.937443825779654, 7.894092264564180});
    x.set(0, 1, {7.652281581062163, 9.955691736800297, 7.877995340784432});

    dist_ref(0, 0) = 1.330479846515945;
    dist_ref(0, 1) = 1.330479846515945;

    distance.mindist_interatomic(x, conf_info, dist);

    for (gpr::Index_t n = 0; n < dist.getNumRows(); ++n) {
        for (gpr::Index_t m = 0; m < dist.getNumCols(); ++m) {
            EXPECT_LE(std::fabs(dist(n, m) - dist_ref(n, m)), threshold)
                << "Elements of the matrix are not equal to the expected ones.";
        }
    }
}

TEST_F(DistanceTest, RotationalInvarianceNoPermute)
{
    gpr::Coord x1, x2;
    x1.resize(1, 3 * 3);
    x2.resize(1, 3 * 3);
    // Define the original atomic coordinates
    x1.set(0, 0, {0.0, 0.0, 0.0});
    x1.set(0, 1, {1.0, 0.0, 0.0});
    x1.set(0, 2, {0.0, 1.0, 0.0});
    // Define the rotated atomic coordinates
    x2.set(0, 0, {0.0, 0.0, 0.0});   // Rotated Atom 1
    x2.set(0, 1, {0.0, 1.0, 0.0});   // Rotated Atom 2
    x2.set(0, 2, {-1.0, 0.0, 0.0});  // Rotated Atom 3
    // Define the atom types and configuration
    AtomsConfiguration conf_info;
    conf_info.atoms_mov.type.resize(1, 3);
    conf_info.atoms_mov.type(0, 0) = 0;  // Type of Atom 1
    conf_info.atoms_mov.type(0, 1) = 1;  // Type of Atom 2
    conf_info.atoms_mov.type(0, 2) = 0;  // Type of Atom 3
    conf_info.n_pt = 3;
    conf_info.pairtype.resize(2, 2);
    conf_info.pairtype.set(EMPTY);
    conf_info.pairtype(0, 0) = 0;
    conf_info.pairtype(0, 1) = 1;
    conf_info.pairtype(1, 0) = 1;
    conf_info.pairtype(1, 1) = 2;
    gpr::Field<double> lengthscale(1, 1);
    lengthscale.set(1.0);
    Field<double> dist;
    distance.dist_at(x1, x2, conf_info, lengthscale, dist);
    ASSERT_NEAR(dist(0, 0), 0.0, 1e-9);
}

TEST_F(DistanceTest, RotationalInvarianceWithPermutations)
/* nice idea, but breaks the GP pretty badly, need to do more */
{
    gpr::Coord x1_orig, x2_orig;
    x1_orig.resize(1, 3 * 3);
    x2_orig.resize(1, 3 * 3);

    // Define the original atomic coordinates.
    x1_orig.set(0, 0, {0.0, 0.0, 0.0});  // Atom 0 (type 0)
    x1_orig.set(0, 1, {1.0, 0.0, 0.0});  // Atom 1 (type 1)
    x1_orig.set(0, 2, {0.0, 1.0, 0.0});  // Atom 2 (type 0)

    // Define the rotated AND permuted atomic coordinates.
    x2_orig.set(0, 0, {-1.0, 0.0, 0.0});  // Rotated position of original Atom 2
    x2_orig.set(0, 1, {0.0, 1.0, 0.0});   // Rotated position of original Atom 1
    x2_orig.set(0, 2, {0.0, 0.0, 0.0});   // Rotated position of original Atom 0

    // Define the atom types and configuration.
    AtomsConfiguration conf_info;
    conf_info.atoms_mov.type.resize(1, 3);
    conf_info.atoms_mov.type(0, 0) = 0;  // Type of Atom 1
    conf_info.atoms_mov.type(0, 1) = 1;  // Type of Atom 2
    conf_info.atoms_mov.type(0, 2) = 0;  // Type of Atom 3

    auto canonical1 =
        aux::get_canonical_configuration(x1_orig, conf_info.atoms_mov.type);
    auto canonical2 =
        aux::get_canonical_configuration(x2_orig, conf_info.atoms_mov.type);

    gpr::Coord& x1_canonical = canonical1.first;
    gpr::Coord& x2_canonical = canonical2.first;

    // The canonicalization should reorder the types as well.
    // Update conf_info to reflect the new canonical order for atom types.
    conf_info.atoms_mov.type = canonical1.second;

    conf_info.n_pt = 3;
    conf_info.pairtype.resize(2, 2);
    conf_info.pairtype.set(EMPTY);
    conf_info.pairtype(0, 0) = 0;
    conf_info.pairtype(0, 1) = 1;
    conf_info.pairtype(1, 0) = 1;
    conf_info.pairtype(1, 1) = 2;

    gpr::Field<double> lengthscale(1, 1);
    lengthscale.set(1.0);

    Field<double> dist;
    distance.dist_at(x1_canonical, x2_canonical, conf_info, lengthscale, dist);

    ASSERT_NEAR(dist(0, 0), 0.0, 1e-9);
}

TEST_F(DistanceTest, RmsdIdentity)
{
    gpr::Coord x;
    x.resize(1, 3 * 3);
    x.set(0, 0, {0.0, 0.0, 0.0});
    x.set(0, 1, {1.5, 0.0, 0.0});
    x.set(0, 2, {0.0, 1.2, 0.5});

    double rmsd = calculate_aligned_rmsd(x, x);
    EXPECT_NEAR(rmsd, 0.0, 1e-9);
}

TEST_F(DistanceTest, RmsdTranslationalInvariance)
{
    gpr::Coord x1, x2;
    x1.resize(1, 3 * 3);
    x2.resize(1, 3 * 3);

    // Original configuration
    x1.set(0, 0, {0.0, 0.0, 0.0});
    x1.set(0, 1, {1.5, 0.0, 0.0});
    x1.set(0, 2, {0.0, 1.2, 0.5});

    // Translated configuration
    const gpr::vector3_reg translation_vec = {10.0, -5.0, 2.5};
    x2.set(0, 0, x1.at(0, 0) + translation_vec);
    x2.set(0, 1, x1.at(0, 1) + translation_vec);
    x2.set(0, 2, x1.at(0, 2) + translation_vec);

    double rmsd = calculate_aligned_rmsd(x1, x2);
    EXPECT_NEAR(rmsd, 0.0, 1e-9);
}

TEST_F(DistanceTest, RmsdTranslationalInvarianceWithTypes)
{
    gpr::Coord x1, x2;
    x1.resize(1, 3 * 3);
    x2.resize(1, 3 * 3);

    // Original configuration
    x1.set(0, 0, {0.0, 0.0, 0.0});
    x1.set(0, 1, {1.5, 0.0, 0.0});
    x1.set(0, 2, {0.0, 1.2, 0.5});

    // Translated configuration
    const gpr::vector3_reg translation_vec = {10.0, -5.0, 2.5};
    x2.set(0, 0, x1.at(0, 0) + translation_vec);
    x2.set(0, 1, x1.at(0, 1) + translation_vec);
    x2.set(0, 2, x1.at(0, 2) + translation_vec);

    AtomsConfiguration local_conf;
    setup_local_conf(local_conf, 3, {0, 1, 2});

    gpr::Field<double> dist;
    distance.dist_rmsd(x1, x2, local_conf, dist);

    ASSERT_EQ(dist.getNumRows(), 1);
    ASSERT_EQ(dist.getNumCols(), 1);
    EXPECT_NEAR(dist(0, 0), 0.0, 1e-9);
}

TEST_F(DistanceTest, RmsdRotationalInvariance)
{
    gpr::Coord x1, x2;
    x1.resize(1, 3 * 3);
    x2.resize(1, 3 * 3);

    // Original configuration
    x1.set(0, 0, {1.0, 2.0, 3.0});
    x1.set(0, 1, {4.0, 5.0, 6.0});
    x1.set(0, 2, {7.0, 8.0, 9.0});

    // Rotated configuration (90 degrees around z-axis)
    // x' = -y, y' = x, z' = z
    x2.set(0, 0, {-2.0, 1.0, 3.0});
    x2.set(0, 1, {-5.0, 4.0, 6.0});
    x2.set(0, 2, {-8.0, 7.0, 9.0});

    double rmsd = calculate_aligned_rmsd(x1, x2);
    EXPECT_NEAR(rmsd, 0.0, 1e-9);
}

TEST_F(DistanceTest, RmsdKnownValue)
{
    gpr::Coord x1, x2;
    x1.resize(1, 3 * 2);
    x2.resize(1, 3 * 2);

    // Simple 2-atom system
    x1.set(0, 0, {0.0, 0.0, 0.0});
    x1.set(0, 1, {1.0, 0.0, 0.0});

    // Same system, but with the second atom moved
    x2.set(0, 0, {0.0, 0.0, 0.0});
    x2.set(0, 1, {1.0, 1.0, 0.0});

    // Expected RMSD calculation:
    // Centered x1: {-0.5, 0, 0}, {0.5, 0, 0}
    // Centered x2: {-0.5, -0.5, 0}, {0.5, 0.5, 0}
    // After optimal rotation (Kabsch), the deviation is known.
    // For this specific case, the aligned coordinates of x1 will be
    // rotated to match x2 as closely as possible. The resulting
    // RMSD should be sqrt(2 - sqrt(2)) / sqrt(2)...
    double expected_rmsd = 0.2071067811865475;

    double rmsd = calculate_aligned_rmsd(x1, x2);
    EXPECT_NEAR(rmsd, expected_rmsd, 1e-9);
}

TEST_F(DistanceTest, DistRmsdMultiConfig)
{
    // This test uses the main dist_rmsd function
    gpr::Coord x1, x2;
    gpr::Field<double> dist;

    x1.resize(1, 3 * 3);  // One test configuration
    x2.resize(3, 3 * 3);  // Three reference configurations

    // Define the base configuration
    x1.set(0, 0, {1.0, 2.0, 3.0});
    x1.set(0, 1, {4.0, 5.0, 6.0});
    x1.set(0, 2, {7.0, 8.0, 9.0});

    // Config 1: A rotated version of x1 (should have RMSD ~ 0)
    x2.set(0, 0, {-2.0, 1.0, 3.0});
    x2.set(0, 1, {-5.0, 4.0, 6.0});
    x2.set(0, 2, {-8.0, 7.0, 9.0});

    // Config 2: A translated version of x1 (should have RMSD ~ 0)
    const gpr::vector3_reg T = {1, 1, 1};
    x2.set(1, 0, x1.at(0, 0) + T);
    x2.set(1, 1, x1.at(0, 1) + T);
    x2.set(1, 2, x1.at(0, 2) + T);

    // Config 3: A distorted version (stretching one bond)
    x2.set(2, 0, {1.0, 2.0, 3.0});
    x2.set(2, 1, {4.0, 5.0, 7.0});  // z coord changed from 6.0 to 7.0
    x2.set(2, 2, {7.0, 8.0, 9.0});
    // Expected RMSD for this case is sqrt( (1.0^2)/3 ) = 1/sqrt(3)
    double expected_rmsd_distorted = 0.4714045207910317;

    AtomsConfiguration local_conf;
    // NOTE(rg): This shouldn't be breaking with the conf_info.. but it does..
    // TBD
    setup_local_conf(local_conf, 3, {0, 1, 2});
    distance.dist_rmsd(x1, x2, local_conf, dist);

    ASSERT_EQ(dist.getNumRows(), 1);
    ASSERT_EQ(dist.getNumCols(), 3);
    EXPECT_NEAR(dist(0, 0), 0.0, 1e-9) << "Rotation failed";
    EXPECT_NEAR(dist(0, 1), 0.0, 1e-9) << "Translation failed";
    EXPECT_NEAR(dist(0, 2), expected_rmsd_distorted, 1e-9)
        << "Distortion failed";
}

TEST_F(DistanceTest, RmsdPermutationalInvariance)
{
    gpr::Coord x1, x2;
    x1.resize(1, 3 * 3);
    x2.resize(1, 3 * 3);
    gpr::Field<double> dist;

    // A-B-A type molecule
    // Atom 0 (Type 0) at {0,0,0}
    // Atom 1 (Type 1) at {1,0,0}
    // Atom 2 (Type 0) at {1,1,0}
    x1.set(0, 0, {0.0, 0.0, 0.0});
    x1.set(0, 1, {1.0, 0.0, 0.0});
    x1.set(0, 2, {1.0, 1.0, 0.0});

    // Same molecule, but rotated 90 deg around Z and with atoms 0 and 2
    // permuted. Rotated positions: Atom 0 -> {0,0,0} Atom 1 -> {0,1,0} Atom 2
    // -> {-1,1,0} Permuted x2:
    x2.set(0, 0, {-1.0, 1.0, 0.0});  // Original Atom 2's new position
    x2.set(0, 1, {0.0, 1.0, 0.0});   // Original Atom 1's new position
    x2.set(0, 2, {0.0, 0.0, 0.0});   // Original Atom 0's new position

    // Call the actual class method
    distance.dist_rmsd(x1, x2, conf_info, dist);

    ASSERT_EQ(dist.getNumRows(), 1);
    ASSERT_EQ(dist.getNumCols(), 1);
    EXPECT_NEAR(dist(0, 0), 0.0, 1e-9) << "RMSD failed to handle permutation";
}

TEST_F(DistanceTest, EmdIdentity)
{
    gpr::Coord x;
    x.resize(1, 3 * 3);
    x.set(0, 0, {0.0, 0.0, 0.0});
    x.set(0, 1, {1.5, 0.0, 0.0});
    x.set(0, 2, {0.0, 1.2, 0.5});

    AtomsConfiguration local_conf;
    setup_local_conf(local_conf, 3, {0, 1, 2});

    gpr::Field<double> dist;
    distance.dist_emd(x, x, local_conf, dist);

    ASSERT_EQ(dist.getNumRows(), 1);
    ASSERT_EQ(dist.getNumCols(), 1);
    EXPECT_NEAR(dist(0, 0), 0.0, 1e-9);
}

TEST_F(DistanceTest, EmdSimpleTranslation)
{
    gpr::Coord x1, x2;
    x1.resize(1, 3 * 2);
    x2.resize(1, 3 * 2);

    x1.set(0, 0, {0.0, 0.0, 0.0});
    x1.set(0, 1, {1.0, 0.0, 0.0});

    const gpr::vector3_reg T = {3.0, 4.0, 0.0};  // length = 5.0
    x2.set(0, 0, x1.at(0, 0) + T);
    x2.set(0, 1, x1.at(0, 1) + T);

    AtomsConfiguration local_conf;
    setup_local_conf(local_conf, 2, {0, 1});

    gpr::Field<double> dist;
    distance.dist_emd(x1, x2, local_conf, dist);

    // EMD is the sum of distances each atom travels.
    // Each of the 2 atoms travels a distance of 5.0.
    double expected_emd = 5.0 + 5.0;

    ASSERT_EQ(dist.getNumRows(), 1);
    ASSERT_EQ(dist.getNumCols(), 1);
    EXPECT_NEAR(dist(0, 0), expected_emd, 1e-9);
}

TEST_F(DistanceTest, EmdPermutationInvariance)
{
    gpr::Coord x1, x2;
    x1.resize(1, 3 * 3);
    x2.resize(1, 3 * 3);
    gpr::Field<double> dist;

    // A-B-A type molecule, e.g. water
    // Atom 0 (Type 0) at {0,0,0}
    // Atom 1 (Type 1) at {1,0,0}
    // Atom 2 (Type 0) at {1,1,0}
    x1.set(0, 0, {0.0, 0.0, 0.0});
    x1.set(0, 1, {1.0, 0.0, 0.0});
    x1.set(0, 2, {1.0, 1.0, 0.0});

    // Same molecule, but with atoms 0 and 2 permuted in the list.
    // The EMD should be 0 because the Hungarian algorithm finds the
    // optimal zero-cost assignment for the identical atoms.
    x2.set(0, 0, {1.0, 1.0, 0.0});  // Original Atom 2's position
    x2.set(0, 1, {1.0, 0.0, 0.0});  // Original Atom 1's position
    x2.set(0, 2, {0.0, 0.0, 0.0});  // Original Atom 0's position

    AtomsConfiguration local_conf;
    setup_local_conf(local_conf, 3, {0, 1, 0});

    distance.dist_emd(x1, x2, local_conf, dist);

    ASSERT_EQ(dist.getNumRows(), 1);
    ASSERT_EQ(dist.getNumCols(), 1);
    EXPECT_NEAR(dist(0, 0), 0.0, 1e-9) << "EMD failed to handle permutation";
}

TEST_F(DistanceTest, EmdMultiTypeAndDistortion)
{
    gpr::Coord x1, x2;
    x1.resize(1, 4 * 3);
    x2.resize(1, 4 * 3);
    gpr::Field<double> dist;

    // Two atoms of type 0, two of type 1
    x1.set(0, 0, {0.0, 0.0, 0.0});  // Type 0, idx 0
    x1.set(0, 1, {1.0, 0.0, 0.0});  // Type 0, idx 1
    x1.set(0, 2, {5.0, 0.0, 0.0});  // Type 1, idx 2
    x1.set(0, 3, {6.0, 0.0, 0.0});  // Type 1, idx 3

    // Distort one atom (idx 1 -> moved by 2.0)
    // and permute the type 1 atoms (idx 2 and 3)
    x2.set(0, 0, {0.0, 0.0, 0.0});  // Original pos of idx 0
    x2.set(0, 1, {1.0, 2.0, 0.0});  // Original pos of idx 1 is moved
    x2.set(0, 2, {6.0, 0.0, 0.0});  // Original pos of idx 3
    x2.set(0, 3, {5.0, 0.0, 0.0});  // Original pos of idx 2

    AtomsConfiguration local_conf;
    setup_local_conf(local_conf, 4, {0, 0, 1, 1});

    distance.dist_emd(x1, x2, local_conf, dist);

    // The optimal assignment for type 0 will be (0->0, 1->1), cost = 2.0
    // The optimal assignment for type 1 will be (2->3, 3->2), cost = 0.0
    // Total EMD should be the sum of costs = 2.0
    double expected_emd = 2.0;

    ASSERT_EQ(dist.getNumRows(), 1);
    ASSERT_EQ(dist.getNumCols(), 1);
    EXPECT_NEAR(dist(0, 0), expected_emd, 1e-9);
}

TEST_F(DistanceTest, EmdFailsRotationInvariance)
{
    gpr::Coord x1, x2;
    x1.resize(1, 3 * 2);
    x2.resize(1, 3 * 2);

    x1.set(0, 0, {0.0, 0.0, 0.0});
    x1.set(0, 1, {1.0, 0.0, 0.0});

    // Create a 90-degree rotation matrix around the Z-axis
    Eigen::Matrix3d R;
    R << cos(M_PI / 2), -sin(M_PI / 2), 0, sin(M_PI / 2), cos(M_PI / 2), 0, 0,
        0, 1;

    // Convert gpr::vector3_reg to Eigen::Vector3d to rotate
    Eigen::Vector3d p0(x1.at(0, 0).x, x1.at(0, 0).y, x1.at(0, 0).z);
    Eigen::Vector3d p1(x1.at(0, 1).x, x1.at(0, 1).y, x1.at(0, 1).z);

    Eigen::Vector3d p0_rot = R * p0;  // Stays at origin
    Eigen::Vector3d p1_rot = R * p1;  // Moves to {0, 1, 0}

    x2.set(0, 0, {p0_rot.x(), p0_rot.y(), p0_rot.z()});
    x2.set(0, 1, {p1_rot.x(), p1_rot.y(), p1_rot.z()});

    AtomsConfiguration local_conf;
    setup_local_conf(local_conf, 2, {0, 0});

    gpr::Field<double> dist;
    distance.dist_emd(x1, x2, local_conf, dist);

    // An ideal rotation-invariant distance should be 0.
    // The current EMD will be non-zero.
    // Expected cost: dist(p0 to p0_rot) + dist(p1 to p1_rot) = 0 + sqrt((1-0)^2
    // + (0-1)^2) = sqrt(2) NOTE(rg): Given that the inverse distance kernel is
    // not invariant to rotations, this shouldn't be either
    EXPECT_GT(dist(0, 0), 1.0);
    EXPECT_FALSE(std::abs(dist(0, 0)) < 1e-9)
        << "Metric should NOT be rotation invariant";
}

#ifdef USE_HIGHS
TEST_F(DistanceTest, AssignmentSolversGiveSameResult)
{
    Eigen::MatrixXd cost_matrix(4, 4);
    cost_matrix << 10, 19, 8, 15, 10, 18, 7, 17, 13, 16, 9, 14, 12, 19, 8, 18;
    double custom_solver_result = aux::solve_assignment_problem(cost_matrix);
    double highs_solver_result =
        aux::solve_assignment_problem_highs(cost_matrix);

    std::cout << "Custom solver result: " << custom_solver_result << std::endl;
    std::cout << "HiGHS solver result:  " << highs_solver_result << std::endl;

    EXPECT_NEAR(custom_solver_result, highs_solver_result, 1e-9)
        << "The custom solver and the HiGHS-based solver produced different "
           "minimum costs.";

    // See scripts/gen_constraint_test.py to verify
    EXPECT_NEAR(custom_solver_result, 49.0, 1e-9);
}
#endif  // USE_HIGHS

} /* namespace tests */
} /* namespace gpr */
