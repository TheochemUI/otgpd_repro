/*
 * PriorUnif.h
 *
 *  Created on: 10 Jul 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_PRIOR_PRIORUNIF_H_
#define GPR_PRIOR_PRIORUNIF_H_

#include <string>

#include "../../data_types/Field.h"
#include "PriorBase.h"

namespace gpr {

/**
 * @brief Uniform prior.
 */
class PriorUnif : public PriorBase {
public:
    PriorUnif() { }
    virtual ~PriorUnif() { }

    double calculateLogPrior(const Field<double>& x) override;

    Field<double> calculateLogPriorGradient(const Field<double>& x) override;
};

} /* namespace gpr */

#endif /* GPR_PRIOR_PRIORUNIF_H_ */
