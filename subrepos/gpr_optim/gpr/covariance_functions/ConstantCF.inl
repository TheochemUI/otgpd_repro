//
//  ConstantCF.inl
//  gpr_dimer
//
//  Created by Maxim Masterov on 21/01/2021.
//

#ifndef ConstantCF_inl
#define ConstantCF_inl

namespace gpr {

inline void ConstantCF::setConstSigma2(const double _constSigma2)
{
    constSigma2 = _constSigma2;
}

inline void ConstantCF::setParameters(const Eigen::VectorXd &w)
{
    // FIXME: this function does nothing in MATLAB code (see
    // gpcf_constant_unpak)
    //        constSigma2 = std::exp(w(0));
}

inline double ConstantCF::getConstSigma2()
{
    return constSigma2;
}

inline Eigen::VectorXd ConstantCF::combineParameters()
{
    Eigen::VectorXd res;
    return res;
}

} /* namespace gpr */

#endif /* ConstantCF_inl */
