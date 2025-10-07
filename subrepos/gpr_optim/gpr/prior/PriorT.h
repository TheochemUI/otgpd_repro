/*
 * PriorT.h
 *
 *  Created on: 10 Jul 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_PRIOR_PRIORT_H_
#define GPR_PRIOR_PRIORT_H_

#include <string>

#include "../../data_types/Field.h"
#include "PriorBase.h"

namespace gpr {

/**
 * @brief Student-t prior.
 *
 * Creates Student's t-distribution prior in which the named parameters
 * have the specified values. Any unspecified parameters are set to
 * default values.
 *
 * The parameterization is as in Gelman, Carlin, Stern, Dunson, Vehtari,
 * and Rubin (2013). Bayesian Data Analysis, third edition.

 */
class PriorT : public PriorBase {
public:
    PriorT() { }
    virtual ~PriorT() { }

    double calculateLogPrior(const Field<double>& x) override;

    Field<double> calculateLogPriorGradient(const Field<double>& x) override;
};

} /* namespace gpr */

#endif /* GPR_PRIOR_PRIORT_H_ */
