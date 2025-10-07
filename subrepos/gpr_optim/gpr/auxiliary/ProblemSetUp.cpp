/*
 * ProblemSetUpimer.cpp
 *
 *  Created on: 24 Sep 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "ProblemSetUp.h"

#include <algorithm>
#include <cfloat>
#include <iostream>
#include <set>

namespace aux {

ProblemSetUp::ProblemSetUp() { }

ProblemSetUp::~ProblemSetUp() { }

// FIXME: change arguments' names (see header)
void ProblemSetUp::cutOffEnergy(const gpr::Field<double>& E_reference,
                                gpr::Field<double>& E_R)
{
    if (E_R.getSize() != 0) E_R = E_R - E_reference;
}

bool ProblemSetUp::activateFrozenAtoms(const gpr::Coord& R,
                                       const double activation_distance,
                                       gpr::AtomsConfiguration& atoms_config)
{
    gpr::AtomsPositionAndType* inactive_atoms =
        &atoms_config.atoms_froz_inactive;
    gpr::AtomsPositionAndType* active_atoms = &atoms_config.atoms_froz_active;
    gpr::AtomsPositionAndType* moving_atoms = &atoms_config.atoms_mov;
    std::vector<int> activated_atoms(inactive_atoms->positions.getNumCols() / 3,
                                     EMPTY);
    std::vector<gpr::Index_t> activated_atomtypes;  // auxiliary vector
    bool new_active_atoms =
        false;  // Indicator of occurrence of new active frozen atoms

    //!> We let it go on if all atoms are moving; but if all the atoms are not
    //! moving and there are no inactive atoms, we throw an error
    if (moving_atoms->positions.getSize() != atoms_config.positions.getSize()) {
        gpr::assertMsg(
            !inactive_atoms->type.isEmpty(),
            "Error! The vector 'atomtype_froz_inactive' is empty when it "
            "shouldn't "
            "be!\n");
    }

    // First, check that there are no clones of moving atoms among inactive
    // frozen atoms. If such a clone exists - remove it from the set of inactive
    // frozen atoms. We are checking atoms by calculating the distance between
    // them.
    removeDuplicatedAtoms(*moving_atoms, *active_atoms, *inactive_atoms);

    // Do it only if inactive atoms exist
    if (!inactive_atoms->positions.isEmpty()) {
        if (activation_distance > 0.) {
            // Activate inactive frozen atoms (for every omage) within the
            // radius of 'activation_distance' from some moving atom
            gpr::Index_t num_inactive_atoms =
                inactive_atoms->positions.getNumPoints();
            std::vector<int> activated_atoms(num_inactive_atoms, EMPTY);
            for (gpr::Index_t i = 0; i < R.getNumRows(); ++i) {
                for (gpr::Index_t n = 0; n < num_inactive_atoms; ++n) {
                    gpr::vector3_reg atom_fro_inactive =
                        inactive_atoms->positions.at(n);
                    for (gpr::Index_t m = 0, end = R.getNumCols() / 3; m < end;
                         ++m) {
                        if ((atom_fro_inactive - R.at(i, m)).length() <
                            activation_distance)
                            activated_atoms[n] = n;
                    }
                }
            }

            // Check how many atoms were activated
            long num_activated_atoms = activated_atoms.size() -
                                       std::count(activated_atoms.begin(),
                                                  activated_atoms.end(), EMPTY);
            // std::cout <<num_activated_atoms << " were activated\n" ;

            if (num_activated_atoms != 0) {
                gpr::Index_t tot_activated_atoms =
                    active_atoms->type.getSize() +
                    (gpr::Index_t)num_activated_atoms;
                std::cout << "\n Now we have " << tot_activated_atoms
                          << " total active atoms\n";
                gpr::Index_t counter = 0;
                active_atoms->resize(tot_activated_atoms);
                activated_atomtypes.resize(num_activated_atoms);

                // First reassign inactive atoms to active ones
                for (gpr::Index_t n = 0; n < activated_atoms.size(); ++n) {
                    if (activated_atoms[n] != EMPTY) {
                        active_atoms->positions.set(
                            0, counter, inactive_atoms->positions.at(n));
                        active_atoms->type[counter] = inactive_atoms->type[n];
                        activated_atomtypes[counter] =
                            active_atoms->type[counter];
                        ++counter;
                    }
                }

                // Now remove reassigned inactive atoms
                for (int n = (int)activated_atoms.size() - 1; n >= 0; --n) {
                    if (activated_atoms[n] != EMPTY) {
                        inactive_atoms->positions.deleteColumnAt(
                            activated_atoms[n]);
                        inactive_atoms->type.deleteColumn(activated_atoms[n]);
                    }
                }

                new_active_atoms = true;
            }
        } else {
            // Activate all inactive frozen atoms
            new_active_atoms = true;
            atoms_config.atoms_froz_active = atoms_config.atoms_froz_inactive;
            inactive_atoms->clear();
        }

        // Activate new pairtypes if necessary
        if (new_active_atoms) {
            // Get unique indices from `atoms_config.atomtype_mov`.
            std::set<gpr::Index_t> unique_atomtype_mov(
                moving_atoms->type.getInternalVector().begin(),
                moving_atoms->type.getInternalVector().end());
            std::set<gpr::Index_t> unique_atomtype_activated(
                activated_atomtypes.begin(), activated_atomtypes.end());

            // TODO: check, this looks like setPairtypeForMovingAtoms()
            for (auto& at_i: unique_atomtype_mov) {
                for (auto& at_j: unique_atomtype_activated) {
                    if (atoms_config.pairtype(at_i, at_j) == EMPTY) {
                        ++atoms_config.n_pt;
                        atoms_config.pairtype(at_i, at_j) =
                            atoms_config.n_pt + EMPTY;
                        atoms_config.pairtype(at_j, at_i) =
                            atoms_config.n_pt + EMPTY;
                    }
                }
            }
        }
    }
    return new_active_atoms;
}

void ProblemSetUp::updateConf(const gpr::Coord& unfreezed_atoms,
                              const gpr::AtomsConfiguration& atoms_config,
                              gpr::Coord& coord_atom_new)
{
    gpr::Index_t counter = 0;
    gpr::Index_t num_atoms = atoms_config.is_frozen.getSize();

    // FIXME: no need to make deep copy, keeping the allocation pattern should
    // be enough.
    coord_atom_new = atoms_config.positions;

    // Note, in the MATLAB version `conf_atom` and `coord_atom_old` are stored
    // in a single matrix with 4 columns
    for (gpr::Index_t n = 0; n < num_atoms; ++n) {
        if (atoms_config.is_frozen[n] == MOVING_ATOM) {
            coord_atom_new[3 * n] = unfreezed_atoms[3 * counter];
            coord_atom_new[3 * n + 1] = unfreezed_atoms[3 * counter + 1];
            coord_atom_new[3 * n + 2] = unfreezed_atoms[3 * counter + 2];
            ++counter;
        }
    }
}

void ProblemSetUp::setPairtypeForMovingAtoms(
    const gpr::Field<gpr::Index_t>& atomtype_mov, gpr::Index_t& n_pt,
    gpr::Field<int>& pairtype)
{
    gpr::Index_t N_mov = atomtype_mov.getSize();
    for (gpr::Index_t i = 0; i < N_mov - 1; ++i) {
        gpr::Index_t at_i = atomtype_mov(0, i);
        for (gpr::Index_t j = i + 1; j < N_mov; ++j) {
            gpr::Index_t at_j = atomtype_mov(0, j);
            if (pairtype(at_i, at_j) == EMPTY) {
                ++n_pt;
                pairtype(at_i, at_j) = n_pt + EMPTY;
                pairtype(at_j, at_i) = n_pt + EMPTY;
            }
        }
    }
}

void ProblemSetUp::removeDuplicatedAtoms(
    const gpr::AtomsPositionAndType& moving_atoms,
    gpr::AtomsPositionAndType& active_atoms,
    gpr::AtomsPositionAndType& inactive_atoms)
{
    std::vector<gpr::Index_t> clones_act;
    std::vector<gpr::Index_t> clones_inact;
    gpr::Index_t num_active_atoms = active_atoms.type.getSize();
    gpr::Index_t num_inactive_atoms = inactive_atoms.type.getSize();
    gpr::Index_t num_moving_atoms = moving_atoms.type.getSize();

    if (num_inactive_atoms == 0 || num_moving_atoms == 0 ||
        num_active_atoms == 0)
        return;

    for (gpr::Index_t m = 0; m < num_moving_atoms; ++m) {
        gpr::vector3_reg mov_at_pos = moving_atoms.positions.at(m);
        for (gpr::Index_t n = 0; n < num_active_atoms; ++n) {
            gpr::vector3_reg act_at_pos = active_atoms.positions.at(n);
            if ((mov_at_pos - act_at_pos).length() <= DBL_EPSILON) {
                clones_act.push_back(n);
            }
        }
        for (gpr::Index_t n = 0; n < num_inactive_atoms; ++n) {
            gpr::vector3_reg inact_at_pos = inactive_atoms.positions.at(n);
            if ((mov_at_pos - inact_at_pos).length() <= DBL_EPSILON) {
                clones_inact.push_back(n);
            }
        }
    }

    eraseClones(clones_act, active_atoms);
    eraseClones(clones_inact, inactive_atoms);
}

void ProblemSetUp::eraseClones(const std::vector<gpr::Index_t>& clones,
                               gpr::AtomsPositionAndType& atoms)
{
    if (!clones.empty()) {
        gpr::io::LogManager log_man;
        log_man << "Warning! Some (in)active frozen atoms were discovered as "
                   "duplicates of the moving atoms and will be removed from "
                   "the data sets..\n";
        for (int n = (int)clones.size() - 1; n >= 0; --n) {
            atoms.positions.deleteColumnAt(clones[n]);
            atoms.type.deleteColumn(clones[n]);
        }
    }
}

} /* namespace aux */
