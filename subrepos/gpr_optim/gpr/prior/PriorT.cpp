/*
 * PriorT.cpp
 *
 *  Created on: 10 Jul 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "PriorT.h"

#include <cmath>

namespace gpr {

double PriorT::calculateLogPrior(const Field<double>& x)
{
    double lp = 0.;
    double nu_rec = 1. / nu;
    double s2_rec = 1. / s2;

    // TODO: optimize it
    for (Index_t n = 0; n < x.getSize(); ++n) {
        lp += lgamma(0.5 * (nu + 1)) - lgamma(0.5 * nu) -
              0.5 * log(nu * M_PI * s2) -
              0.5 * (nu + 1) *
                  log(1. + (x[n] - mu) * (x[n] - mu) * nu_rec * s2_rec);
    }

    return lp;
}

Field<double> PriorT::calculateLogPriorGradient(const Field<double>& x)
{
    Field<double> lpg;

    lpg.resize(x.getNumRows(), x.getNumCols());

    for (Index_t n = 0; n < lpg.getSize(); ++n) {
        double tmp = x[n] - mu;
        lpg[n] = -(nu + 1) * tmp / (nu * s2 + tmp * tmp);
    }

    return lpg;
}

} /* namespace gpr */
