//
//  ProblemSetUp.inl
//  gpr_dimer
//
//  Created by Maxim Masterov on 02/02/2021.
//

#ifndef ProblemSetUp_inl
#define ProblemSetUp_inl

namespace aux {

template <typename Pot>
void ProblemSetUp::calculateGeneralPotential(
    const gpr::AtomsConfiguration& atoms_config,
    const gpr::vector3_reg (&cell_dimensions)[3],
    gpr::Observation& middle_point, Pot& general_potential,
    gpr::Index_t& num_of_calls)
{
    gpr::Index_t num_im = middle_point.R.getNumRows();
    gpr::Index_t dim = middle_point.R.getNumCols();
    gpr::Coord coord_atom_new;
    gpr::Coord G_local;
    double E_local = 0.;
    gpr::Index_t num_atoms = atoms_config.is_frozen.getSize();
    double box[9];

    box[0] = cell_dimensions[0].x;
    box[1] = cell_dimensions[0].y;
    box[2] = cell_dimensions[0].z;

    box[3] = cell_dimensions[1].x;
    box[4] = cell_dimensions[1].y;
    box[5] = cell_dimensions[1].z;

    box[6] = cell_dimensions[2].x;
    box[7] = cell_dimensions[2].y;
    box[8] = cell_dimensions[2].z;

    middle_point.E.resize(num_im, 1);
    middle_point.G.resize(num_im, dim);

    updateConf(middle_point.R, atoms_config, coord_atom_new);

    G_local.resize(coord_atom_new.getNumRows() * coord_atom_new.getNumCols());

    for (gpr::Index_t n = 0; n < num_im; ++n) {
        general_potential.force(num_atoms,
                                coord_atom_new.getInternalVector().data(),
                                atoms_config.atomicNrs.getData(),
                                G_local.getInternalVector().data(), &E_local,
                                nullptr /*variance*/, box);

        middle_point.E(n, 0) = E_local;

        gpr::Index_t counter = 0;
        for (gpr::Index_t m = 0; m < num_atoms; ++m) {
            if (atoms_config.is_frozen[m] == MOVING_ATOM) {
                middle_point.G.set(n, counter, G_local.at(n, m));
                counter++;
                if (!((counter) % (dim / 3))) counter = 0;
            }
        }
        middle_point.G *= -1.;

        ++num_of_calls;
        // FIXME: check if transposition and reshaping are actually needed
        // Gi = reshape(Gi',1,size(Gi,1)*3);
        // G(i,:) = Gi;
    }
}

template <typename Pot>
void ProblemSetUp::initializeEnergyAndGradient(
    const gpr::Field<double>& E_reference,
    const gpr::AtomsConfiguration& atoms_config,
    const gpr::vector3_reg (&system_size)[3], gpr::Observation& middle_point,
    Pot& general_potential, gpr::Index_t& num_of_calls)
{
    middle_point.E.resize(E_reference.getNumRows(), E_reference.getNumCols());
    calculateGeneralPotential(atoms_config, system_size, middle_point,
                              general_potential, num_of_calls);
}

} /* namespace aux */

#endif /* ProblemSetUp_h */
