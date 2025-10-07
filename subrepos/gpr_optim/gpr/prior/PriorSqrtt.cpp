/*
 * PriorSqrtt.cpp
 *
 *  Created on: 9 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "PriorSqrtt.h"

#include <cmath>

namespace gpr {

double PriorSqrtt::calculateLogPrior(const Field<double>& x)
{
    double lp = 0.;
    double nu_rec = 1. / nu;
    double s2_rec = 1. / s2;

    // TODO: optimize it
    for (Index_t n = 0; n < x.getSize(); ++n) {
        double sqrt_x = sqrt(x[n]);
        lp += lgamma(0.5 * (nu + 1)) - lgamma(0.5 * nu) -
              0.5 * log(nu * M_PI * s2) -
              0.5 * (nu + 1) *
                  log(1. + (sqrt_x - mu) * (sqrt_x - mu) * nu_rec * s2_rec) -
              log(2. * sqrt_x);
    }

    return lp;
}

Field<double> PriorSqrtt::calculateLogPriorGradient(const Field<double>& x)
{
    Field<double> lpg;

    lpg.resize(x.getNumRows(), x.getNumCols());

    // TODO: optimize it
    for (Index_t n = 0; n < x.getSize(); ++n) {
        double sqrt_x = sqrt(x[n]);
        double sqrt_x_2_rec = 0.5 / sqrt_x;
        lpg[n] = sqrt_x_2_rec * (-(nu + 1) * (sqrt_x - mu) /
                                 (nu * s2 + (sqrt_x - mu) * (sqrt_x - mu))) -
                 0.5 / x[n];
    }

    return lpg;
}

} /* namespace gpr */
