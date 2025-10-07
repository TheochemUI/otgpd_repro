/*
 * main.cpp
 *
 *  Created on: 28 May 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include <unordered_map>

#include "gpr/AtomicDimer.h"
#include "gpr/Enums.h"
#include "gpr/auxiliary/ProblemSetUp.h"
#include "gpr/potentials/EAMPotential.h"
#include "managers/io/FileManager.h"

#ifdef USE_CAPNP
#include <capnp/message.h>
#include <capnp/serialize.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include "managers/io/CapnProtoManager.h"
#endif


/**
 * @brief Fill in the structure of \e atoms_config.
 *
 * This function will populate a structure of \e atoms_config and determine
 * frozen active, frozen inactive and moving atoms. It will also set up the
 * vector of pairtypes.
 *
 * @param R_init Initial coordinates of moving atoms
 * @param parameters Structure of input parameters
 * @param atomtypes Map of atom types, where key represents the atom names
 * @param atoms_config Structure of all atoms configuration. It should contain
 * information on positions of all atoms, their IDs and an idicator if atom is
 * frozen or moving
 */
void initializeAtomsConfiguration(
    const gpr::Coord& R_init, const gpr::InputParameters& parameters,
    const std::unordered_map<std::string, gpr::Index_t> atomtypes,
    gpr::AtomsConfiguration& atoms_config)
{
    aux::ProblemSetUp problem_setup;
    gpr::Index_t number_of_mov_atoms = atoms_config.countMovingAtoms();
    gpr::Index_t number_of_fro_atoms =
        atoms_config.is_frozen.getSize() - number_of_mov_atoms;

    // Resize structures for moving and frozen atoms
    atoms_config.atoms_mov.resize(number_of_mov_atoms);
    atoms_config.atoms_froz_inactive.resize(number_of_fro_atoms);

    atoms_config.atoms_mov.type.set(atomtypes.at("H2"));
    atoms_config.atoms_froz_inactive.type.set(atomtypes.at("Cu"));

    // Assign moving and frozen atoms and list all frozen atoms as inactive
    gpr::Index_t counter_f = 0, counter_m = 0;
    for (gpr::Index_t n = 0; n < atoms_config.is_frozen.getSize(); ++n) {
        if (atoms_config.is_frozen[n] == MOVING_ATOM)
            atoms_config.atoms_mov.positions.set(0, counter_m++,
                                                 atoms_config.positions.at(n));
        else
            atoms_config.atoms_froz_inactive.positions.set(
                0, counter_f++, atoms_config.positions.at(n));
    }

    // Pairtype indices for pairs of atomtypes (n_at x n_at)
    // Active pairtypes are indexed as 0,1,...,n_pt-1. Inactive pairtypes are
    // given index EMPTY.
    atoms_config.pairtype.resize((gpr::Index_t)atomtypes.size(),
                                 (gpr::Index_t)atomtypes.size());
    atoms_config.pairtype.set(EMPTY);

    // Set pairtype indices for moving+moving atom pairs (and update number of
    // active pairtypes)
    problem_setup.setPairtypeForMovingAtoms(
        atoms_config.atoms_mov.type, atoms_config.n_pt, atoms_config.pairtype);

    // Activate frozen atoms within activation distance
    problem_setup.activateFrozenAtoms(R_init, parameters.actdist_fro.value,
                                      atoms_config);
}

void setUpProblem(atmd::AtomicDimer& atomic_dimer)
{
    typedef gpr::Field<double> FieldDbl;

    std::map<std::string, gpr::Coord*> dict_coord;
    std::map<std::string, FieldDbl*> dict_field;

    gpr::io::FileManager fm;

    aux::ProblemSetUp problem_setup;
    gpr::Coord R_all_init;  // no initial data points
    FieldDbl E_all_init;    // no initial data
    gpr::Coord G_all_init;  // no initial data
    gpr::Coord R_init;
    FieldDbl E_init;
    gpr::Coord G_init;
    gpr::Coord orient_init, orient_init_tmp, orient_start;
    gpr::Coord R_sp, R;
    FieldDbl conf;
    gpr::InputParameters parameters;
    gpr::Observation init_observations;
    gpr::Observation init_middle_point;
    gpr::AtomsConfiguration atoms_config;

    atoms_config.clear();

    // Read input parameters
#ifdef USE_CAPNP
    try {
        gpr::io::loadParametersFromCapnp("input/capnp_params.bin", parameters);
    } catch (const std::exception& e) {
        std::cerr << "Fatal Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
#else
    fm.readInputFile("input/input.dat", parameters);
#endif

    // Read the initial orientation from the input file
    dict_coord["orient_init"] = &orient_init;
    dict_coord["orient_start"] = &orient_start;
    fm.readDataFile("input/orients_CuH2.dat", dict_coord);
    dict_coord.clear();

    // Read the initial middle point from the input file
    dict_coord["R_sp"] = &R_sp;
    fm.readDataFile("input/displacedH2.dat", dict_coord);
    dict_coord.clear();

    // Read the initial configuration data from the input file
    dict_field["conf"] = &conf;
    fm.readDataFile("input/CuH2-init_PES_point.dat", dict_field);
    dict_field.clear();

    // Set up a structure with configuration of all atoms
    atoms_config.assignFromField(conf);

    // Define the initial middle point of the dimer
    // Initial middle point not observed
    R_init.resize(1, R_sp.getNumCols());
    double dist_sp = parameters.dist_sp.value[parameters.i_dist.value];
    for (gpr::Index_t n = 0; n < R_sp.getNumCols(); ++n) {
        R_init[n] = R_sp[n] + dist_sp * orient_start(parameters.i_run.value, n);
    }
    E_init.clear();
    G_init.clear();

    /* ********************************************************* */
    /* ********************************************************* */
    /* ********************************************************* */
    std::unordered_map<std::string, gpr::Index_t> atomtypes = {{"H2", 0},
                                                               {"Cu", 1}};
    // We know that only the last two atoms are H2
    atoms_config.type.resize(1, atoms_config.id.getSize());
    atoms_config.type.set(atomtypes["Cu"]);
    atoms_config.type[atoms_config.type.getSize() - 2] = atomtypes["H2"];
    atoms_config.type[atoms_config.type.getSize() - 1] = atomtypes["H2"];
    initializeAtomsConfiguration(R_init, parameters, atomtypes, atoms_config);
    /* ********************************************************* */
    /* ********************************************************* */
    /* ********************************************************* */

    init_observations.clear();
    init_middle_point.clear();
    init_middle_point.R = R_init;

    atomic_dimer.initialize(parameters, init_observations, init_middle_point,
                            orient_init, atoms_config);
}

int main(int argc, char** argv)
{
    int error = 0;

    //    ::testing::InitGoogleTest(&argc, argv);
    //    error = RUN_ALL_TESTS();

    // Declare main objects
    atmd::AtomicDimer atomic_dimer;
    pot::EAMPotential* eam_potential = new pot::EAMPotential;

    // Set up the problem
    // Define (default) properties of the internal algorithms
    // TODO: just make a separate class to initialize everything
    //    atomic_dimer.setUpDefault();
    setUpProblem(atomic_dimer);

    // Execute the dimer method with GPR
    atomic_dimer.execute(*eam_potential);

    delete eam_potential;

    return error;
}
