/*
 * PriorLogUnif.h
 *
 *  Created on: 10 Jul 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_PRIOR_PRIORLOGUNIF_H_
#define GPR_PRIOR_PRIORLOGUNIF_H_

#include <string>

#include "../../data_types/Field.h"
#include "PriorBase.h"

namespace gpr {

/**
 * @brief Uniform prior for the logarithm of the parameter.
 */
class PriorLogUnif : public PriorBase {
public:
    PriorLogUnif() { }
    virtual ~PriorLogUnif() { }

    double calculateLogPrior(const Field<double>& x) override;

    Field<double> calculateLogPriorGradient(const Field<double>& x) override;
};

} /* namespace gpr */

#endif /* GPR_PRIOR_PRIORLOGUNIF_H_ */
