/*
 * ProblemSetUpTest.cpp
 *
 *  Created on: 24 Sep 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "ProblemSetUpTest.h"

#include <cmath>

#include "../../../gpr/potentials/EAMPotential.h"

namespace gpr {
namespace tests {

ProblemSetUpTest::ProblemSetUpTest()
{
    threshold = 1e3 * DBL_EPSILON;  // The reference values in F_rot_ref
                                    // are taken from the Matlab code.
                                    // Due to difference in the precision
                                    // between C++ and Matlab, we are
                                    // scaling the threshold by 1e3.
}

ProblemSetUpTest::~ProblemSetUpTest() { }

TEST_F(ProblemSetUpTest, activateFrozenAtoms)
{
    gpr::Coord R_new;
    gpr::Index_t actdist_fro;
    bool atoms_were_changed;
    gpr::AtomsConfiguration atoms_conf;
    gpr::Coord conf_fro_ref;

    R_new.resize(1, 3 * 2);
    atoms_conf.atoms_froz_inactive.resize(216);
    atoms_conf.atoms_froz_active.resize(216);
    atoms_conf.pairtype.resize(2, 2);
    atoms_conf.atoms_mov.resize(2);
    conf_fro_ref.resize(1, 3 * 190);

    R_new.set(0, 0, {8.982373164830571, 9.937230835772036, 7.894416323850492});
    R_new.set(0, 1, {7.652483227274959, 9.955905494573976, 7.877879589983660});

    actdist_fro = 5;

    atoms_conf.atoms_mov.positions.set(0, 0, R_new.at(0));
    atoms_conf.atoms_mov.positions.set(0, 1, R_new.at(1));
    atoms_conf.atoms_mov.type.set(0);
    atoms_conf.atoms_froz_inactive.type.set(1);
    atoms_conf.pairtype.set(EMPTY);
    atoms_conf.pairtype(0, 0) = 0;
    atoms_conf.n_pt = 1;

    io_manager.readFromPlainFile(
        "tests/reference/gpr/auxiliary/ConfigurationInactiveFrozen.dat",
        atoms_conf.atoms_froz_inactive.positions.getInternalVector());

    io_manager.readFromPlainFile(
        "tests/reference/gpr/auxiliary/ReferenceConfigurationFrozen.dat",
        conf_fro_ref.getInternalVector());

    atoms_were_changed =
        setup.activateFrozenAtoms(R_new, actdist_fro, atoms_conf);

    EXPECT_EQ(atoms_were_changed, true)
        << "Variable `atoms_were_changed` is not equal to 'true'.";

    EXPECT_EQ(atoms_conf.pairtype(0, 0), 0) << "Wrong value of pairtype(0, 0).";
    EXPECT_EQ(atoms_conf.pairtype(0, 1), 1) << "Wrong value of pairtype(0, 1).";
    EXPECT_EQ(atoms_conf.pairtype(1, 0), 1) << "Wrong value of pairtype(1, 0).";
    EXPECT_EQ(atoms_conf.pairtype(1, 1), EMPTY)
        << "Wrong value of pairtype(1, 1).";

    EXPECT_EQ(atoms_conf.atoms_froz_inactive.positions.getNumRows(), 1)
        << "Wrong number of elements in i-th direction in field "
           "atoms_conf.conf_fro.";
    EXPECT_EQ(atoms_conf.atoms_froz_inactive.positions.getNumCols(), 3 * 190)
        << "Wrong number of elements in j-th direction in field "
           "atoms_conf.conf_fro.";

    EXPECT_EQ(atoms_conf.atoms_froz_inactive.type.getNumRows(), 1)
        << "Wrong number of elements in i-th direction in field "
           "atoms_conf.atomtype_fro.";
    EXPECT_EQ(atoms_conf.atoms_froz_inactive.type.getNumCols(), 190)
        << "Wrong number of elements in j-th direction in field "
           "atoms_conf.atomtype_fro.";

    for (gpr::Index_t i = 0;
         i < atoms_conf.atoms_froz_inactive.positions.getNumRows(); ++i) {
        for (gpr::Index_t j = 0;
             j < atoms_conf.atoms_froz_inactive.positions.getNumCols(); ++j) {
            EXPECT_LE(fabs(atoms_conf.atoms_froz_inactive.positions(i, j) -
                           conf_fro_ref(i, j)),
                      DBL_EPSILON)
                << "Elements of the `atoms_conf.conf_fro` "
                   "are not equal to the expected ones.";
        }
    }

    for (gpr::Index_t n = 0; n < atoms_conf.atoms_froz_inactive.type.getSize();
         ++n) {
        EXPECT_EQ(atoms_conf.atoms_froz_inactive.type[n], 1)
            << "Elements of the `atoms_conf.atomtype_fro` are not "
               "equal to the expected ones.";
    }
}

TEST_F(ProblemSetUpTest, cutOffEnergy)
{
    gpr::Observation middle_point;
    gpr::Coord G_R_ref;
    gpr::Field<double> E_level;
    gpr::Field<double> E_R_ref;
    gpr::vector3_reg box[3];
    gpr::AtomsConfiguration atom_config;
    pot::EAMPotential eam_potential;
    gpr::Index_t num_of_potential_calls;

    middle_point.R.resize(1, 2 * 3);
    atom_config.positions.resize(1, 218 * 3);
    atom_config.is_frozen.resize(1, 218);
    E_level.resize(1, 1);
    E_R_ref.resize(1, 1);
    G_R_ref.resize(1, 2 * 3);

    box[0].set(15.3455999999999992, 0., 0.);
    box[1].set(0., 21.7020000000000017, 0.);
    box[2].set(0., 0., 100.0000000000000000);

    E_level(0, 0) = -696.587191493847e+000;

    middle_point.R(0, 0) = 8.97856277303058e+000;
    middle_point.R(0, 1) = 9.93211628067229e+000;
    middle_point.R(0, 2) = 7.89882761414426e+000;
    middle_point.R(0, 3) = 7.64888749663556e+000;
    middle_point.R(0, 4) = 9.95517051512886e+000;
    middle_point.R(0, 5) = 7.87274215046670e+000;

    atom_config.is_frozen.set(FROZEN_ATOM);
    atom_config.is_frozen(0, 216) = MOVING_ATOM;
    atom_config.is_frozen(0, 217) = MOVING_ATOM;

    io_manager.readFromPlainFile(
        "tests/reference/gpr/auxiliary/InputAtomsPositions.dat",
        atom_config.positions.getInternalVector());

    E_R_ref(0, 0) = 370.498489473903e-006;

    G_R_ref(0, 0) = 21.4452456338877e-003;
    G_R_ref(0, 1) = -33.1231900862804e-003;
    G_R_ref(0, 2) = 28.3157675762483e-003;
    G_R_ref(0, 3) = -15.8030003032378e-003;
    G_R_ref(0, 4) = 20.8628464314946e-003;
    G_R_ref(0, 5) = -38.3251019392285e-003;

    setup.initializeEnergyAndGradient(E_level, atom_config, box, middle_point,
                                      eam_potential, num_of_potential_calls);
    setup.cutOffEnergy(E_level, middle_point.E);

    EXPECT_LE(std::fabs(middle_point.E(0, 0) - E_R_ref(0, 0)), 10 * threshold)
        << "Energy value is greater than the expected one.";

    for (gpr::Index_t n = 0; n < middle_point.G.getNumPoints(); ++n) {
        EXPECT_LE(std::fabs(middle_point.G.at(n).x - G_R_ref.at(n).x),
                  threshold)
            << "X-component of the gradient is greater than the expected one."
            << n;
        EXPECT_LE(std::fabs(middle_point.G.at(n).y - G_R_ref.at(n).y),
                  threshold)
            << "Y-component of the gradient is greater than the expected one."
            << n;
        EXPECT_LE(std::fabs(middle_point.G.at(n).z - G_R_ref.at(n).z),
                  threshold)
            << "Z-component of the gradient is greater than the expected one."
            << n;
    }
}

TEST_F(ProblemSetUpTest, UpdateConf)
{
    gpr::Coord unfreezed_atoms;
    gpr::AtomsConfiguration atoms_config;
    gpr::Coord coord_atom_new;
    gpr::Coord coord_atom_ref;
    gpr::Index_t size = 10;

    unfreezed_atoms.resize(1, 3 * 2);
    atoms_config.positions.resize(1, 3 * size);
    atoms_config.is_frozen.resize(1, size);
    atoms_config.id.resize(1, size);
    coord_atom_ref.resize(1, 3 * size);

    unfreezed_atoms.set(
        0, 0, {8.978562773030582, 9.932116280672286, 7.898827614144257});
    unfreezed_atoms.set(
        0, 1, {7.648887496635559, 9.955170515128863, 7.872742150466696});

    atoms_config.positions.set(
        0, 0, {5.7546000000000, 8.1385070000000, 6.9753600000000});
    atoms_config.positions.set(
        0, 1, {8.3122000000000, 8.1385070000000, 6.9753600000000});
    atoms_config.positions.set(
        0, 2, {10.869800000000, 8.1385070000000, 6.9753600000000});
    atoms_config.positions.set(
        0, 3, {5.7546000000000, 11.755500000000, 6.9753590000000});
    atoms_config.positions.set(
        0, 4, {8.3122000000000, 11.755500000000, 6.9753590000000});
    atoms_config.positions.set(
        0, 5, {10.869800000000, 11.755500000000, 6.9753590000000});
    atoms_config.positions.set(
        0, 6, {7.0334000000000, 9.9470030000000, 5.7487170000000});
    atoms_config.positions.set(
        0, 7, {9.5910000000000, 9.9470030000000, 5.7487170000000});
    atoms_config.positions.set(
        0, 8, {8.6822270000000, 9.9470040000000, 11.732978000000});
    atoms_config.positions.set(
        0, 9, {7.9421730000000, 9.9470040000000, 11.732978000000});

    atoms_config.is_frozen(0, 0) = FROZEN_ATOM;
    atoms_config.is_frozen(0, 1) = FROZEN_ATOM;
    atoms_config.is_frozen(0, 2) = FROZEN_ATOM;
    atoms_config.is_frozen(0, 3) = FROZEN_ATOM;
    atoms_config.is_frozen(0, 4) = FROZEN_ATOM;
    atoms_config.is_frozen(0, 5) = FROZEN_ATOM;
    atoms_config.is_frozen(0, 6) = FROZEN_ATOM;
    atoms_config.is_frozen(0, 7) = FROZEN_ATOM;
    atoms_config.is_frozen(0, 8) = MOVING_ATOM;
    atoms_config.is_frozen(0, 9) = MOVING_ATOM;

    atoms_config.id(0, 0) = 208;
    atoms_config.id(0, 1) = 209;
    atoms_config.id(0, 2) = 210;
    atoms_config.id(0, 3) = 211;
    atoms_config.id(0, 4) = 212;
    atoms_config.id(0, 5) = 213;
    atoms_config.id(0, 6) = 214;
    atoms_config.id(0, 7) = 215;
    atoms_config.id(0, 8) = 216;
    atoms_config.id(0, 9) = 217;

    coord_atom_ref.set(0, 0,
                       {5.7546000000000, 8.1385070000000, 6.9753600000000});
    coord_atom_ref.set(0, 1,
                       {8.3122000000000, 8.1385070000000, 6.9753600000000});
    coord_atom_ref.set(0, 2,
                       {10.869800000000, 8.1385070000000, 6.9753600000000});
    coord_atom_ref.set(0, 3,
                       {5.7546000000000, 11.755500000000, 6.9753590000000});
    coord_atom_ref.set(0, 4,
                       {8.3122000000000, 11.755500000000, 6.9753590000000});
    coord_atom_ref.set(0, 5,
                       {10.869800000000, 11.755500000000, 6.9753590000000});
    coord_atom_ref.set(0, 6,
                       {7.0334000000000, 9.9470030000000, 5.7487170000000});
    coord_atom_ref.set(0, 7,
                       {9.5910000000000, 9.9470030000000, 5.7487170000000});
    coord_atom_ref.set(
        0, 8, {8.978562773030582, 9.932116280672286, 7.898827614144257});
    coord_atom_ref.set(
        0, 9, {7.648887496635559, 9.955170515128863, 7.872742150466696});

    setup.updateConf(unfreezed_atoms, atoms_config, coord_atom_new);

    for (gpr::Index_t n = 0; n < coord_atom_new.getNumCols() / 3; ++n) {
        EXPECT_LE((coord_atom_new.at(0, n) - coord_atom_ref.at(0, n)).length(),
                  DBL_EPSILON)
            << "Elements of the field are not equal to the expected ones.";
    }
}

TEST_F(ProblemSetUpTest, SetPairTypeMov)
{
    gpr::Field<gpr::Index_t> atomtype_mov;
    gpr::Index_t n_pt;
    gpr::Field<int> pairtype;

    atomtype_mov.resize(1, 2);
    pairtype.resize(2, 2);

    atomtype_mov.set(0);
    pairtype.set(EMPTY);
    n_pt = 0;

    setup.setPairtypeForMovingAtoms(atomtype_mov, n_pt, pairtype);

    EXPECT_EQ(n_pt, 1) << "Value of n_pt is not equal to one.";

    EXPECT_EQ(pairtype(0, 0), 0)
        << "Value of pairtype(0, 0) is not equal to one.";
    EXPECT_EQ(pairtype(0, 1), EMPTY)
        << "Value of pairtype(0, 1) is not equal to zero.";
    EXPECT_EQ(pairtype(1, 0), EMPTY)
        << "Value of pairtype(1, 0) is not equal to zero.";
    EXPECT_EQ(pairtype(1, 1), EMPTY)
        << "Value of pairtype(1, 1) is not equal to zero.";
}

TEST_F(ProblemSetUpTest, generalPotential)
{
    gpr::Observation middle_point;
    gpr::vector3_reg box[3];
    double E_ref;
    gpr::vector3_reg G_ref[2];
    gpr::AtomsConfiguration atom_config;
    pot::EAMPotential eam_potential;
    gpr::Index_t num_atoms = 218;
    gpr::Index_t num_of_gen_potential_calls = 0;

    middle_point.E.resize(1);
    middle_point.R.resize(1, 2 * 3);

    atom_config.positions.resize(1, num_atoms * 3);
    atom_config.is_frozen.resize(1, num_atoms);
    atom_config.id.resize(1, num_atoms);

    box[0].set(15.3455999999999992, 0., 0.);
    box[1].set(0., 21.7020000000000017, 0.);
    box[2].set(0., 0., 100.0000000000000000);

    middle_point.R(0, 0) = 8.98237316483057e+000;
    middle_point.R(0, 1) = 9.93723083577204e+000;
    middle_point.R(0, 2) = 7.89441632385049e+000;
    middle_point.R(0, 3) = 7.65248322727496e+000;
    middle_point.R(0, 4) = 9.95590549457398e+000;
    middle_point.R(0, 5) = 7.87787958998366e+000;

    atom_config.is_frozen.set(FROZEN_ATOM);
    atom_config.is_frozen[num_atoms - 2] = MOVING_ATOM;
    atom_config.is_frozen[num_atoms - 1] = MOVING_ATOM;

    for (gpr::Index_t n = 0; n < num_atoms; ++n)
        atom_config.id[n] = n;

    io_manager.readFromPlainFile(
        "tests/reference/gpr/auxiliary/InputAtomsPositions.dat",
        atom_config.positions.getInternalVector());

    E_ref = -696.58719149384660;
    G_ref[0].set(22.2174753484639e-003, -22.2402530044595e-003,
                 18.5353722102027e-003);
    G_ref[1].set(-17.3837548235882e-003, 21.9137885404279e-003,
                 -25.8161071576406e-003);

    setup.calculateGeneralPotential(atom_config, box, middle_point,
                                    eam_potential, num_of_gen_potential_calls);

    EXPECT_LE(std::fabs(middle_point.E[0] - E_ref), threshold)
        << "Energy value is greater than the expected one.";

    for (gpr::Index_t n = 0; n < middle_point.G.getNumPoints(); ++n) {
        EXPECT_LE(std::fabs(middle_point.G.at(n).x - G_ref[n].x), threshold)
            << "X-component of the gradient is greater than the expected one."
            << n;
        EXPECT_LE(std::fabs(middle_point.G.at(n).y - G_ref[n].y), threshold)
            << "Y-component of the gradient is greater than the expected one."
            << n;
        EXPECT_LE(std::fabs(middle_point.G.at(n).z - G_ref[n].z), threshold)
            << "Z-component of the gradient is greater than the expected one."
            << n;
    }
}

} /* namespace tests */
} /* namespace gpr */
