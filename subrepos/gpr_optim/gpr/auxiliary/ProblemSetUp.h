/*
 * GPDimer.h
 *
 *  Created on: 24 Sep 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_GPDIMER_H_
#define GPR_GPDIMER_H_

#include <vector>

#include "../../data_types/Coord.h"
#include "../../data_types/Field.h"
#include "../../structures/Structures.h"

namespace aux {

// TODO: rename or move to AtomicDimer class
// used to be GP_Dimer

/**
 * @brief Set up the GP-dimer problem.
 */
class ProblemSetUp {
public:
    ProblemSetUp();
    virtual ~ProblemSetUp();

    /**
     * @brief Activate inactive frozen atoms within the given radius.
     *
     * Activates inactive frozen atoms within the radius \e actdist_fro from
     * some moving atom in configurations \e R_new. These frozen atoms are then
     * taken into account in the covariance function.
     *
     * When a frozen atoms is activated, its coordinates and atomtype index are
     * added to \e conf_info.conf_fro and \e conf_info.atomtype_fro,
     * respectively, and removed from  \e conf_info_inactive.conf_fro and
     * \e conf_info_inactive.atomtype_fro. If the new frozen atom activates new
     * pairtypes, also \e conf_info.pairtype and \e conf_info.n_pt are updated.
     *
     * @param R Coordinates of moving atoms
     * @param activation_distance Activation distance for moving+frozen atom
     * pairs (pass a negative value if you want to activate all atoms)
     * @param atoms_config Structure including information about the
     *                     configurations necessary for the GP model
     * @return True if new active frozen atoms included, False if not
     */
    // used to be update_active_fro
    bool activateFrozenAtoms(const gpr::Coord& R,
                             const double activation_distance,
                             gpr::AtomsConfiguration& atoms_config);

    /**
     * @brief Subtract reference energy value from \e E_R.
     *
     * @param E_reference A field with the reference enerfy
     * @param E_R In/out energy field
     */
    void cutOffEnergy(const gpr::Field<double>& E_reference,
                      gpr::Field<double>& E_R);

    /**
     * @brief Initialize energy and gradient.
     */
    template <typename Pot>
    void initializeEnergyAndGradient(
        const gpr::Field<double>& E_reference,
        const gpr::AtomsConfiguration& atoms_config,
        const gpr::vector3_reg (&system_size)[3],
        gpr::Observation& middle_point, Pot& general_potential,
        gpr::Index_t& num_of_calls);

    /**
     * @brief Invoke potential calculation.
     *
     * Each row of 'R' represents one configuration including the coordinates of
     * the moving atoms: [x_1,y_1,z_1,x_2,y_2,z_2,...].
     *
     * @param atoms_config Full configuration of the system (i.e. both
     * constrained and unconstrained atoms)
     * @param cell_dimensions Cell dimension
     * @param mid_point Coordinates, energy and gradients of the middle point of
     * the dimer (1 x D)
     * @param general_potential Reference to a class with general potential
     * @param num_of_gen_potential_calls Number of calls for general potential
     */
    template <typename Pot>
    void calculateGeneralPotential(const gpr::AtomsConfiguration& atoms_config,
                                   const gpr::vector3_reg (&cell_dimensions)[3],
                                   gpr::Observation& mid_point,
                                   Pot& general_potential,
                                   gpr::Index_t& num_of_calls);

    /**
     * @brief Update an atomic configuration with new coordinates for the
     * unfreezed atoms.
     *
     * @param unfreezed_atoms The coordinates of the unfreezed atoms
     * @param atoms_config The initial configuration of atoms.
     * @param coord_atom_new  The new coordinates of all atoms.
     */
    void updateConf(const gpr::Coord& unfreezed_atoms,
                    const gpr::AtomsConfiguration& atoms_config,
                    gpr::Coord& coord_atom_new);

    /**
     * @brief Set pairtype indices for moving+moving atom pairs and
     * gives the updated number of active pairtypes.
     *
     * The provided \e n_pt should be equal to the number of active pairtypes
     * before the update.
     *
     * @param atomtype_mov Atomtype indices for moving atoms (1 x N_mov)
     * @param pairtype     Pairtype indices for pairs of atomtypes before the
     * update (n_at x n_at)
     * @param n_pt         Number of active pairtypes before and after the
     * update
     */
    void setPairtypeForMovingAtoms(const gpr::Field<gpr::Index_t>& atomtype_mov,
                                   gpr::Index_t& n_pt,
                                   gpr::Field<int>& pairtype);

    /**
     * @brief Remove clones of the moving atoms from the set of inactive frozen
     * atoms.
     */
    void removeDuplicatedAtoms(const gpr::AtomsPositionAndType& moving_atoms,
                               gpr::AtomsPositionAndType& active_atoms,
                               gpr::AtomsPositionAndType& inactive_atoms);

private:
    void eraseClones(const std::vector<gpr::Index_t>& clones,
                     gpr::AtomsPositionAndType& atoms);
};

} /* namespace aux */

#include "ProblemSetUp.inl"

#endif /* GPR_GPDIMER_H_ */
