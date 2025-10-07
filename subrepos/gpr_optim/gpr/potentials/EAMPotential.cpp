//
//  EAMPotential.cpp
//  gpr_dimer
//
//  Created by Maxim Masterov on 30/01/2021.
//

#include "EAMPotential.h"

#include <fstream>

namespace pot {

void EAMPotential::force(long N, const double *R, const int * /*atomicNrs*/,
                         double *F, double *U, double *variance,
                         const double *box)
{
    if (variance != nullptr) {
        *variance = 0.0;
    }
    gpr::vector3_reg cell_dimensions[3];
    std::string con_file_name = "tmp.con";
    std::string forces_file_name = "forces.out";
    std::string energy_file_name = "energy.out";
    std::string energy_excutable_name = "./mEAMCUH2";
    std::ofstream os;
    std::ifstream is;
    gpr::Index_t num_atoms = (gpr::Index_t)N;

    cell_dimensions[0].x = box[0];
    cell_dimensions[0].y = box[1];
    cell_dimensions[0].z = box[2];

    cell_dimensions[1].x = box[3];
    cell_dimensions[1].y = box[4];
    cell_dimensions[1].z = box[5];

    cell_dimensions[2].x = box[6];
    cell_dimensions[2].y = box[7];
    cell_dimensions[2].z = box[8];

    gpr::vector3_reg celldim = {cell_dimensions[0].x, cell_dimensions[1].y,
                                cell_dimensions[2].z};

    /*
     * THIS BLOCK WRITES PROPER INPUT (.con) FOR THE POTENTIAL
     * note: rueit only works if two H atoms are included and placed at the end
     * of the Cu gpr::Coordinates
     */
    os.open(con_file_name.c_str(), std::ios::out);
    assert(((void)("Output file cannot be opened"), os.is_open()));

    os << "Random Number Seed"
       << "\n";
    os << "Time"
       << "\n";

    os << celldim << "\n";

    // FIXME: should not be hard-coded, or there should be a comment
    os << "90.0000000000000000 90.0000000000000000 90.0000000000000000"
       << "\n";
    os << "0 0"
       << "\n";

    os << num_atoms << " " << 0 << " " << 1 << "\n";
    os << 2 << "\n";
    assert(((void)("The number of points in vector R is too small"),
            num_atoms * 3 >= 2));
    os << num_atoms - 2 << " " << 2 << "\n";
    os << "63.5459999999999994 1.00792999999999999"
       << "\n";

    os.setf(std::ios::fixed);
    os.precision(14);
    os << "Cu"
       << "\n";
    os << "Components of Type 1"
       << "\n";
    for (gpr::Index_t n = 0; n < num_atoms - 2; ++n) {
        os << R[3 * n] << " " << R[3 * n + 1] << " " << R[3 * n + 2] << " " << 0
           << " " << n << "\n";
    }

    os << "H"
       << "\n";
    os << "Components of Type 2"
       << "\n";
    for (gpr::Index_t n = num_atoms - 2; n < num_atoms; ++n) {
        os << R[3 * n] << " " << R[3 * n + 1] << " " << R[3 * n + 2] << " " << 0
           << " " << n << "\n";
    }

    os.close();

    // Call the executable (compiled from a fortran code)
    // Note, the executable for the energy estimation should be placed in the
    // same folder as the executable for the GPR code.
    int error = std::system(energy_excutable_name.c_str());
    gpr::assertMsg(error == 0, "Error! Binary " + energy_excutable_name +
                                   " exited with the error code " +
                                   std::to_string(error));

    // Get forces
    is.open(forces_file_name.c_str(), std::ios::in);
    gpr::Index_t dummy_id;
    gpr::vector3_reg dummy_vec;
    for (gpr::Index_t n = 0; n < num_atoms; ++n) {
        is >> dummy_id >> dummy_vec.x >> dummy_vec.y >> dummy_vec.z;
        F[3 * n] = dummy_vec.x;
        F[3 * n + 1] = dummy_vec.y;
        F[3 * n + 2] = dummy_vec.z;
    }
    is.close();

    // Get energy
    is.open(energy_file_name.c_str(), std::ios::in);
    is >> *U;
    is.close();
}

} /* namespace pot */
