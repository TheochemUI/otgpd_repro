/*
 * PriorLogUnif.cpp
 *
 *  Created on: 10 Jul 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "PriorLogUnif.h"

#include <cmath>

namespace gpr {

double PriorLogUnif::calculateLogPrior(const Field<double>& x)
{
    double res = 0.;

    for (Index_t n = 0; n < x.getSize(); ++n) {
        res -= log(x[n]);
    }

    return res;
}

Field<double> PriorLogUnif::calculateLogPriorGradient(const Field<double>& x)
{
    Field<double> lpg;

    lpg.resize(x.getNumRows(), x.getNumCols());

    for (Index_t n = 0; n < x.getSize(); ++n) {
        lpg[n] = -1. / x[n];
    }

    return lpg;
}

} /* namespace gpr */
