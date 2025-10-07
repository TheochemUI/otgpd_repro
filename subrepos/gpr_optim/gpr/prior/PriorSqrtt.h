/*
 * PriorSqrtt.h
 *
 *  Created on: 9 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_PRIOR_PRIORSQRTT_H_
#define GPR_PRIOR_PRIORSQRTT_H_

#include <string>

#include "../../data_types/Field.h"
#include "PriorBase.h"

namespace gpr {

/**
 * @brief Student-t prior for the square root of the parameter.
 *
 * Creates for the square root of the parameter Student's
 * t-distribution prior structure in which the named parameters
 * have the specified values. Any unspecified parameters are set
 * to default values.
 */
class PriorSqrtt : public PriorBase {
public:
    PriorSqrtt() { }
    virtual ~PriorSqrtt() { }

    double calculateLogPrior(const Field<double>& x) override;

    Field<double> calculateLogPriorGradient(const Field<double>& x) override;
};

} /* namespace gpr */

#endif /* GPR_PRIOR_PRIORSQRTT_H_ */
