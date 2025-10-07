//
//  SexpatCF.inl
//  gpr_dimer
//
//  Created by Maxim Masterov on 21/01/2021.
//

#ifndef SexpatCF_inl
#define SexpatCF_inl

namespace gpr {

inline void SexpatCF::setMagnSigma2(double value)
{
    magnSigma2 = value;
}

inline void SexpatCF::setLengthScale(double value)
{
    lengthScale.set(value);
}

inline void SexpatCF::setConfInfo(const AtomsConfiguration& value)
{
    conf_info = value;
}

inline void SexpatCF::setParameters(const Eigen::VectorXd& w)
{
    magnSigma2 = std::exp(w(0));
    for (Index_t n = 0; n < lengthScale.getSize(); ++n)
        lengthScale[n] = std::exp(w(n + 1));
}

inline double SexpatCF::getMagnSigma2()
{
    return magnSigma2;
}

inline Field<double>& SexpatCF::getLengthScaleRef()
{
    return lengthScale;
}

inline void SexpatCF::clear()
{
    initialized = false;
    magnSigma2 = 0.;
    lengthScale.set(0.);
    conf_info.clearDistribution();
}

inline void SexpatCF::setPriorParametersGaussian(const PriorBase& prior)
{
    prior_parameters_gaussian.setParameters(prior);
}

inline void SexpatCF::setPriorParametersSqrtt(const PriorBase& prior)
{
    prior_parameters_sqrtt.setParameters(prior);
}

inline PriorBase SexpatCF::getPriorParametersGaussian()
{
    return prior_parameters_gaussian;
}

inline PriorBase SexpatCF::getPriorParametersSqrtt()
{
    return prior_parameters_sqrtt;
}

inline Eigen::VectorXd SexpatCF::combineParameters()
{
    Eigen::VectorXd res(1 + lengthScale.getSize());

    res(0) = log(magnSigma2);

    for (Index_t n = 0; n < lengthScale.getSize(); ++n)
        res(n + 1) = log(lengthScale[n]);

    // TODO: check if parameters or prior functions should also be appended
    return res;
}

} /* namespace gpr */

#endif /* SexpatCF_inl */
