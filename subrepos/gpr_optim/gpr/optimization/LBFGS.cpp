/*
 * LBFGS.cpp
 *
 *  Created on: 31 Mar 2021
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "LBFGS.h"

namespace gpr {
namespace optim {

LBFGS::LBFGS()
{
    default_scaling = 0.01;
}

LBFGS::~LBFGS() { }

void LBFGS::execute(const Coord& F, const LBFGSInfo& info,
                    Eigen::VectorXd& rot_orient)
{
    Index_t dof = F.getNumCols();
    Index_t m = info.deltaOrient.getNumRows();  // History size
    Eigen::VectorXd q(dof);
    Eigen::VectorXd s(dof);
    Eigen::VectorXd y(dof);
    Eigen::VectorXd* z = &rot_orient;
    Eigen::VectorXd alpha(m);
    double beta;
    double rho = 0.;
    double gamma = 0.;

    // should take: - positions
    //              - rotational forces
    z->resize(dof);
    alpha.setZero();

    // The negative sign is due to conversion of forces to gradients
    q = -F.extractEigenVector();

    for (int i = m - 1; i >= 0; --i) {
        s = info.deltaOrient.extractRowAsEigenVector(i);
        y = -info.deltaF.extractRowAsEigenVector(i);
        rho = 1. / y.dot(s);
        alpha(i) = rho * s.dot(q);
        q = q - alpha(i) * y;
    }

    if (m > 0) {
        s = info.deltaOrient.extractRowAsEigenVector(m - 1);
        y = -info.deltaF.extractRowAsEigenVector(m - 1);
        gamma = s.dot(y) / y.dot(y);
    } else {
        gamma = default_scaling;
    }

    *z = q * gamma;

    for (gpr::Index_t i = 0; i < m; ++i) {
        s = info.deltaOrient.extractRowAsEigenVector(i);
        y = -info.deltaF.extractRowAsEigenVector(i);
        rho = 1. / y.dot(s);
        beta = rho * y.dot(*z);
        *z = *z + s * (alpha(i) - beta);
    }
}

} /* namespace optim */
} /* namespace gpr */
