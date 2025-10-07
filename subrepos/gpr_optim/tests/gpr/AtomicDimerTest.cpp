/*
 * AtomicDimerTest.cpp
 *
 *  Created on: 18 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "AtomicDimerTest.h"

#include "../../managers/io/FileManager.h"
#include "EONPotential/Morse.h"

namespace gpr {
namespace tests {

AtomicDimerTest::AtomicDimerTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

AtomicDimerTest::~AtomicDimerTest() { }

TEST_F(AtomicDimerTest, InterfaceToEONPotential)
{
    gpr::Observation middle_point;
    gpr::AtomsConfiguration atom_config;
    gpr::vector3_reg box[3];
    io::FileManager io_manager;
    gpr::Index_t num_atoms = 2;
    Morse morse;
    double E_ref;
    gpr::Coord G_ref;

    middle_point.E.resize(1);
    middle_point.R.resize(1, num_atoms * 3);
    G_ref.resize(1, num_atoms * 3);

    atom_config.positions.resize(1, num_atoms * 3);
    atom_config.is_frozen.resize(1, num_atoms);
    atom_config.id.resize(1, num_atoms);

    box[0].set(10., 0., 0.);
    box[1].set(0., 10., 0.);
    box[2].set(0., 0., 10.);

    middle_point.R(0, 0) = 0.;
    middle_point.R(0, 1) = 0.;
    middle_point.R(0, 2) = 0.40150089545845091;
    middle_point.R(0, 3) = 0.;
    middle_point.R(0, 4) = 0.;
    middle_point.R(0, 5) = 3.29849910454155237;

    atom_config.is_frozen.set(MOVING_ATOM);

    for (gpr::Index_t n = 0; n < num_atoms; ++n)
        atom_config.id[n] = n;

    io_manager.readFromPlainFile(
        "tests/reference/gpr/auxiliary/InputAtomsPositions.dat",
        atom_config.positions.getInternalVector());

    dimer.setAtomsConfiguration(atom_config);
    dimer.callGeneralPotentialFromEON(box, middle_point, morse);

    E_ref = -0.71016446199686;
    G_ref(0, 0) = 0.;
    G_ref(0, 1) = 0.;
    G_ref(0, 2) = -6.5505184668007e-06;
    G_ref(0, 3) = 0.;
    G_ref(0, 4) = 0.;
    G_ref(0, 5) = 6.5505184668007e-06;

    EXPECT_LE(fabs(middle_point.E[0] - E_ref), threshold)
        << "Energy is not equal to the expected one.";
    for (gpr::Index_t n = 0; n < G_ref.getSize(); ++n)
        EXPECT_LE(fabs(middle_point.G[n] - G_ref[n]), threshold)
            << "A component of the gradient is not equal to the expected one.";
}

} /* namespace tests */
} /* namespace gpr */
