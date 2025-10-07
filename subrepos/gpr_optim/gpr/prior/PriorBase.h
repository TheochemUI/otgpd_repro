/*
 * PriorBase.h
 *
 *  Created on: 10 Jul 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_PRIOR_PRIORPARAMETERS_H_
#define GPR_PRIOR_PRIORPARAMETERS_H_

#include "../../data_types/Field.h"

namespace gpr {

/**
 * @brief Class of prior parameters.
 */
class PriorBase {
public:
    PriorBase()
    {
        clear();
    }

    virtual double calculateLogPrior(const Field<double>& x)
    {
        return 0.;
    }

    virtual Field<double> calculateLogPriorGradient(const Field<double>& x)
    {
        Field<double> lpg;
        lpg.resize(x.getNumRows(), x.getNumCols(), true);
        return lpg;
    }

    inline void clear()
    {
        mu = s2 = nu = 0.;
    }

    inline double getMu() const
    {
        return mu;
    }

    inline double getS2() const
    {
        return s2;
    }

    inline double getNu() const
    {
        return nu;
    }

    inline void setMu(const double value)
    {
        mu = value;
    }

    inline void setS2(const double value)
    {
        s2 = value;
    }

    inline void setNu(const double value)
    {
        nu = value;
    }

    inline void setParameters(const PriorBase& other)
    {
        mu = other.getMu();
        s2 = other.getS2();
        nu = other.getNu();
    }

protected:
    double mu;  // location
    double s2;  // scale squared (variance)
    double nu;  // degrees of freedom
};

} /* namespace gpr */

#endif /* GPR_PRIOR_PRIORPARAMETERS_H_ */
