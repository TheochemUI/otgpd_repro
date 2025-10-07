/*
 * PriorGaussian.h
 *
 *  Created on: 9 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_PRIOR_PRIORGAUSSIAN_H_
#define GPR_PRIOR_PRIORGAUSSIAN_H_

#include <string>

#include "../../data_types/Field.h"
#include "PriorBase.h"

namespace gpr {

/**
 * @brief Gaussian prior.
 *
 * Creates Gaussian prior in which parameters have the specified values. Any
 * unspecified parameters are set to default values.
 */
class PriorGaussian : public PriorBase {
public:
    PriorGaussian() { }
    virtual ~PriorGaussian() { }

    double calculateLogPrior(const Field<double>& x) override;

    Field<double> calculateLogPriorGradient(const Field<double>& x) override;
};

} /* namespace gpr */

#endif /* GPR_PRIOR_PRIORGAUSSIAN_H_ */
