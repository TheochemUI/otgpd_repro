/*
 * LikGaussian.h
 *
 *  Created on: 22 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_OBSERVATION_MODELS_LIKGAUSSIAN_H_
#define GPR_OBSERVATION_MODELS_LIKGAUSSIAN_H_

#include <Eigen/Dense>
#include <cmath>

#include "../../data_types/Coord.h"
#include "../../data_types/Field.h"
#include "../prior/PriorBase.h"

namespace gpr {

/**
 * @brief Gaussian noise with variance.
 */
class LikGaussian {
public:
    LikGaussian() : sigma2(0.) { }
    virtual ~LikGaussian() { }

    /**
     * \breif Evaluate the log prior of covariance function parameters.
     *
     * This function returns \f[ log(p(th)) \f], where th collects the
     * parameters. This function is needed when there are likelihood parameters.
     */
    double evaluateLogPrior();

    /**
     * @brief Evaluate gradient of the log prior with respect to the parameters.
     *
     * This function takes a Gaussian likelihood and returns
     * \f[ LPG = d log (p(th))/dth \f], where th is the vector of parameters.
     * This function is needed when there are likelihood parameters.
     *
     * @param lpg Result
     */
    void evaluateLogPriorGradient(Field<double>& lpg);

    /**
     * @brief Evaluate training covariance matrix of inputs.
     *
     * @param x Set of coordinates
     * @param C Covariance matrix
     */
    void evaluateTrainingCovarianceMatrix(const Coord& x, Field<double>& C);

    /**
     * @brief Set sigma2 value.
     */
    inline void setSigma2(const double value);

    inline void setPriorParametersGaussian(const PriorBase& prior);

    inline void setPriorParametersSqrtt(const PriorBase& prior);

    inline void setParameters(const Eigen::VectorXd& w);

    /**
     * @brief Combine parameters of the gaussian likelihood into one Eigen
     * vector.
     */
    inline Eigen::VectorXd combineParameters();

private:
    double sigma2;
    PriorBase prior_parameters_gaussian;
    PriorBase prior_parameters_sqrtt;
};

} /* namespace gpr */

#include "LikGaussian.inl"

#endif /* GPR_OBSERVATION_MODELS_LIKGAUSSIAN_H_ */
