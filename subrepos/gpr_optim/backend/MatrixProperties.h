//
//  MatrixProperties.hpp
//  gpr_dimer
//
//  Created by Maxim Masterov on 18/01/2021.
//

#ifndef MatrixProperties_hpp
#define MatrixProperties_hpp

#include "../backend/Macros.h"

namespace math {
/**
 * @brief Contains methods for checking basic matrix properties.
 */
class MatrixProperties {
public:
    /**
     * @brief Return true if the given matrix is Diagonally Dominant (DD).
     */
    bool isDD(const gpr::EigenMatrix& Matrix);

    /**
     * @brief Return true if the given matrix is Positive Definite (PD).
     */
    bool isPD(const gpr::EigenMatrix& Matrix);

    /**
     * @brief Return true if the given matrix is Positive Semi-Definite (PSD).
     */
    bool isPSD(const gpr::EigenMatrix& Matrix);

private:
    /**
     * @brief Calculate x^T A x, where 'x' is an arbitrary vector.
     */
    double xTAx(const gpr::EigenMatrix& Matrix);
};
}  // namespace math

#endif /* MatrixProperties_hpp */
