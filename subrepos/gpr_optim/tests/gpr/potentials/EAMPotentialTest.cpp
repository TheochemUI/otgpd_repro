//
//  EAMPotentialTest.cpp
//  gpr_dimer
//
//  Created by Maxim Masterov on 02/02/2021.
//

#include "EAMPotentialTest.h"

namespace gpr {
namespace tests {

EAMPotentialTest::EAMPotentialTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

EAMPotentialTest::~EAMPotentialTest() { }

TEST_F(EAMPotentialTest, CuH2)
{
    gpr::Coord R;
    gpr::Coord F;
    gpr::Coord F_ref;
    double E;
    double E_ref;
    double celldim[9];
    pot::EAMPotential eam_potential;
    gpr::Index_t num_atoms = 218;

    R.resize(1, 3 * num_atoms);
    F_ref.resize(1, 3 * num_atoms);
    F.resize(1, 3 * num_atoms);

    celldim[0] = 15.3455999999999992;
    celldim[1] = 0.;
    celldim[2] = 0.;

    celldim[3] = 0.;
    celldim[4] = 21.7020000000000017;
    celldim[5] = 0.;

    celldim[6] = 0.;
    celldim[7] = 0.;
    celldim[8] = 100.0000000000000000;

    io_manager.readFromPlainFile(
        "tests/reference/gpr/potentials/InputRCuH2.dat", R.getInternalVector());
    // Reference
    E_ref = -696.58719149384660;

    io_manager.readFromPlainFile(
        "tests/reference/gpr/potentials/ReferenceFCuH2.dat",
        F_ref.getInternalVector());

    eam_potential.force(num_atoms, R.getInternalVector().data(), nullptr,
                        F.getInternalVector().data(), &E, nullptr, celldim);

    EXPECT_LE(std::fabs(E - E_ref), threshold)
        << "Energy value is greater than the expected one.";

    for (gpr::Index_t n = 0; n < R.getNumPoints(); ++n) {
        EXPECT_LE(std::fabs(F.at(n).x - F_ref.at(n).x), threshold)
            << "X-component of the force field is greater than the expected "
               "one."
            << n;
        EXPECT_LE(std::fabs(F.at(n).y - F_ref.at(n).y), threshold)
            << "Y-component of the force field is greater than the expected "
               "one."
            << n;
        EXPECT_LE(std::fabs(F.at(n).z - F_ref.at(n).z), threshold)
            << "Z-component of the force field is greater than the expected "
               "one."
            << n;
    }
}

} /* namespace tests */
} /* namespace gpr */
