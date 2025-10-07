/*
 * LikGaussian.cpp
 *
 *  Created on: 22 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "LikGaussian.h"

namespace gpr {

double LikGaussian::evaluateLogPrior()
{
    return 0.;
}

void LikGaussian::evaluateLogPriorGradient(Field<double>& lpg)
{
    lpg.clear();
}

void LikGaussian::evaluateTrainingCovarianceMatrix(const Coord& x,
                                                   Field<double>& C)
{
    Index_t n = x.getNumRows();

    C.resize(n, n);
    for (Index_t i = 0; i < C.getNumRows(); ++i) {
        C(i, i) = sigma2;
    }
}

} /* namespace gpr */
