//
//  AdditionalFunctionality.h
//  gpr_dimer
//
//  Created by Maxim Masterov on 04/12/2020.
//

#ifndef AdditionalFunctionality_h
#define AdditionalFunctionality_h

#include "../../data_types/Coord.h"
#include "../../structures/Structures.h"

namespace aux {

/**
 * @brief Some auxiliary functions.
 */
class AuxiliaryFunctionality {
public:
    AuxiliaryFunctionality() { }
    virtual ~AuxiliaryFunctionality() { }

    /**
     * @brief Assemble Eigen matrix with repetitive rows based on provided set
     * of coordinates \e R.
     *
     * This function mimics the following MATAB code:
     * \code{.m}
     * [N_im,D] = size(R);
     * matrix = [repmat(R,D+1,1),reshape(repmat(0:D,N_im,1),[],1)];
     * \endcode
     * The last column in the resulting \e matrix consists of repetitive integer
     * indices.
     *
     * \note The resulting \e matrix will have size [(size(R, 2) + 1) * size(R,
     * 1); size(R, 2) + 1].
     *
     * @param coord Original set of coordinates
     * @param matrix Matrix assembled using \e coord
     * @param ind Vector of indices
     */
    inline void assembleMatrixOfRepetitiveCoordinates(const gpr::Coord& coord,
                                                      gpr::EigenMatrix& matrix,
                                                      Eigen::VectorXd& ind);

    /**
     * @brief Combine Energy and gradient into a single Eigen vector.
     */
    inline void assembleVectorFromEnergyAndGradient(
        const gpr::Observation& observation, Eigen::VectorXd& vector);

    /**
     * @brief Simple replica of the repmat from MATLAB.
     * \b repmat(ref_value, Ni, Nj)
     */
    template <typename D, template <typename> class T>
    inline void repmatConst(const gpr::Index_t Ni, const gpr::Index_t Nj,
                            const D value, T<D>& field);
};

} /* namespace aux */

#include "AdditionalFunctionality.inl"

#endif /* AdditionalFunctionality_h */
