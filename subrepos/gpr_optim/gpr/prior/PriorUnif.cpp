/*
 * PriorUnif.cpp
 *
 *  Created on: 10 Jul 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "PriorUnif.h"

namespace gpr {

double PriorUnif::calculateLogPrior(const Field<double>& x)
{
    return 0.;
}

Field<double> PriorUnif::calculateLogPriorGradient(const Field<double>& x)
{
    Field<double> lpg;

    lpg.resize(x.getNumRows(), x.getNumCols());
    lpg.setZero();

    return lpg;
}

} /* namespace gpr */
