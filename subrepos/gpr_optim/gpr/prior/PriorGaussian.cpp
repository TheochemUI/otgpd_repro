/*
 * PriorGaussian.cpp
 *
 *  Created on: 9 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "PriorGaussian.h"

#include <cmath>

namespace gpr {

double PriorGaussian::calculateLogPrior(const Field<double>& x)
{
    double lp = 0.;
    double log_2_pi = log(2 * M_PI);
    double s2_rec = 1. / s2;
    Field<double> sum_x_min_mu;
    double sum_all = 0.;

    sum_x_min_mu.resize(1, x.getNumCols());
    sum_x_min_mu.setZero();

    for (Index_t i = 0; i < x.getNumRows(); ++i) {
        for (Index_t j = 0; j < x.getNumCols(); ++j) {
            sum_x_min_mu(0, j) += (x(i, j) - mu) * (x(i, j) - mu);
        }
    }

    for (Index_t n = 0; n < sum_x_min_mu.getSize(); ++n) {
        sum_all += -log_2_pi - log(s2) - s2_rec * sum_x_min_mu[n];
    }

    lp = 0.5 * sum_all;

    return lp;
}

Field<double> PriorGaussian::calculateLogPriorGradient(const Field<double>& x)
{
    Field<double> lpg;

    lpg.resize(x.getNumRows(), x.getNumCols());

    double s2_rec = 1. / s2;
    for (Index_t n = 0; n < lpg.getSize(); ++n) {
        lpg[n] = s2_rec * (mu - x[n]);
    }

    return lpg;
}

} /* namespace gpr */
