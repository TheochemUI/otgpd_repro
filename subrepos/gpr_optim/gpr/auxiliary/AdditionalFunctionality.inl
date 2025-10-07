//
//  AdditionalFunctionality.inl
//  gpr_dimer
//
//  Created by Maxim Masterov on 21/01/2021.
//

#ifndef AdditionalFunctionality_inl
#define AdditionalFunctionality_inl

namespace aux {

inline void AuxiliaryFunctionality::assembleMatrixOfRepetitiveCoordinates(
    const gpr::Coord& coord, gpr::EigenMatrix& matrix, Eigen::VectorXd& ind)
{
    // [repmat(R,D+1,1),reshape(repmat(0:D,N_im,1),[],1)]
    gpr::Index_t N_im = coord.getNumRows();

    matrix.resize(coord.getNumRows() * (coord.getNumCols() + 1),
                  coord.getNumCols() + 1);
    ind.resize(coord.getNumRows() * (coord.getNumCols() + 1));

    for (gpr::Index_t i = 0; i < matrix.rows(); ++i) {
        for (gpr::Index_t j = 0; j < matrix.cols() - 1; ++j) {
            matrix(i, j) = coord(i % N_im, j);
        }
    }

    gpr::Index_t counter = 0;
    for (gpr::Index_t i = 0; i < matrix.rows(); ++i) {
        matrix(i, matrix.cols() - 1) = counter;
        ind(i) = counter;
        if (((i + 1) % N_im) == 0 || N_im == 1) ++counter;
    }
}

inline void AuxiliaryFunctionality::assembleVectorFromEnergyAndGradient(
    const gpr::Observation& observation, Eigen::VectorXd& vector)
{
    gpr::Index_t counter = 0;

    vector.resize(observation.E.getSize() + observation.G.getSize());

    for (gpr::Index_t n = 0; n < observation.E.getSize(); ++n)
        vector(counter++) = observation.E[n];

    for (gpr::Index_t j = 0; j < observation.G.getNumCols(); ++j) {
        for (gpr::Index_t i = 0; i < observation.G.getNumRows(); ++i) {
            vector(counter++) = observation.G(i, j);
        }
    }
}

template <typename D, template <typename> class T>
inline void AuxiliaryFunctionality::repmatConst(const gpr::Index_t Ni,
                                                const gpr::Index_t Nj,
                                                const D value, T<D>& field)
{
    field.resize(Ni, Nj);
    field.set(value);
}

} /* namespace aux */

#endif /* AdditionalFunctionality_inl */
