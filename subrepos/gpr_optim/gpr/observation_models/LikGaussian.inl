//
//  LikGaussian.inl
//  gpr_dimer
//
//  Created by Maxim Masterov on 21/01/2021.
//

#ifndef LikGaussian_h
#define LikGaussian_h

namespace gpr {

inline void LikGaussian::setSigma2(const double value)
{
    sigma2 = value;
}

inline void LikGaussian::setPriorParametersGaussian(const PriorBase& prior)
{
    prior_parameters_gaussian.setParameters(prior);
}

inline void LikGaussian::setPriorParametersSqrtt(const PriorBase& prior)
{
    prior_parameters_sqrtt.setParameters(prior);
}

inline void LikGaussian::setParameters(const Eigen::VectorXd& w)
{
    // FIXME: this function does nothing in the MATLAB code
    // see `lik_gaussian_unpak`
}

/**
 * @brief Combine parameters of the gaussian likelihood into one Eigen vector.
 */
inline Eigen::VectorXd LikGaussian::combineParameters()
{
    Eigen::VectorXd res;
    return res;
}

} /* namespace gpr */

#endif /* LikGaussian_h */
