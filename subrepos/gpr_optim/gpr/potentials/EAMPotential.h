//
//  EAMPotential.hpp
//  gpr_dimer
//
//  Created by Maxim Masterov on 30/01/2021.
//

#ifndef EAMPotential_hpp
#define EAMPotential_hpp

#include <Eigen/Dense>

#include "../../data_types/Coord.h"
#include "../../structures/Structures.h"

namespace pot {
/**
 * @brief EAM potential. Compatible with EON.
 */
class EAMPotential {
public:
    /**
     * @brief Invoke EAM potential calculation.
     *
     * Each row of 'R' represents one configuration including the coordinates of
     * the moving atoms: [x_1,y_1,z_1,x_2,y_2,z_2,...]. Note, \e F should be
     * pre-allocated. The size of \e F should be equal to 3*N, same as the size
     * of \e R.
     *
     * @param N Number of atoms
     * @param R Pointer to array of positions
     * @param atomicNrs ???
     * @param F Pointer to array of forces
     * @param U Pointer to internal energy
     * @param box Adress to supercell size
     * @param nImages Number of images (always 1 for EAM potential)
     */
    void force(long N, const double* R, const int* atomicNrs, double* F,
               double* U, double *variance, const double* box);
};

} /* End namespace pot */

#endif /* EAMPotential_hpp */
