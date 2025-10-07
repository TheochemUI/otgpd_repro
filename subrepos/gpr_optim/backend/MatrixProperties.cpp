//
//  MatrixProperties.cpp
//  gpr_dimer
//
//  Created by Maxim Masterov on 18/01/2021.
//

#include "MatrixProperties.h"

namespace math {

bool MatrixProperties::isDD(const gpr::EigenMatrix& Matrix)
{
    bool res = true;

    for (long i = 0; i < Matrix.rows(); ++i) {
        double off_diag = 0.;
        double diag = 0.;
        for (long j = 0; j < Matrix.cols(); ++j) {
            double coeff = Matrix(i, j);
            if (j != i)
                off_diag += fabs(coeff);
            else
                diag = fabs(coeff);
        }

        /*
         * If absolute value of diagonal element is smaller than a sum of
         * absolute values of off-diagonal elements matrix is not DD
         */
        if (fabs(diag) < fabs(off_diag)) {
            res = false;
            break;
        }
    }

    return res;
}

bool MatrixProperties::isPD(const gpr::EigenMatrix& Matrix)
{
    return xTAx(Matrix) > 0 ? true : false;
}

bool MatrixProperties::isPSD(const gpr::EigenMatrix& Matrix)
{
    return xTAx(Matrix) >= 0 ? true : false;
}

double MatrixProperties::xTAx(const gpr::EigenMatrix& Matrix)
{
    Eigen::VectorXd z(Matrix.rows());

    z.setOnes();

    // Just add some zeros at some positions
    for (int i = 0; i < z.rows(); i += 4) {
        z(i) = 0.;
    }

    return z.dot(Matrix * z);
}
}  // namespace math
