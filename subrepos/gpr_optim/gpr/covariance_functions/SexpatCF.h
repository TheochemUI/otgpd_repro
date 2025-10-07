/*
 * GPCF.h
 *
 *  Created on: 8 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_SEXPAT_H_
#define GPR_SEXPAT_H_

#include <cmath>
#include <vector>

#include "../../data_types/Coord.h"
#include "../../data_types/Field.h"
#include "../../structures/Structures.h"
#include "../prior/PriorBase.h"

namespace gpr {

/**
 * @brief Squared exponential covariance function.
 */
class SexpatCF {
public:
    SexpatCF();
    virtual ~SexpatCF();

    /**
     * @brief Evaluate the log prior density of the parameters.
     *
     * This function takes a covariance function structure GPCF and returns
     * log(p(th)), where th collects the parameters. This is a mandatory
     * subfunction used for example in energy computations.
     *
     * Evaluate the prior contribution to the error. The parameters that are
     * sampled are transformed, e.g., W = log(w) where w is all the "real"
     * samples. On the other hand errors are evaluated in the W-space so we
     * need take into account also the Jacobian of transformation, e.g.,
     * W -> w = exp(W). See Gelman et al. (2013), Bayesian Data Analysis,
     * third edition, p. 21.
     */
    double calculateLogPrior();

    /**
     * @brief Evaluate the gradient of the log prior density of the parameters
     *        with respect to parameters.
     *
     * This function takes a covariance function structure GPCF and returns
     * LPG = d log (p(th))/dth, where th is the vector of parameters. This is a
     * mandatory subfunction used in gradient computations.
     *
     * @param lpg
     */
    Field<double> calculateLogPriorGradient();

    /**
     * @brief Evaluate covariance matrix between two input vectors.
     *
     * The function takes in covariance function of a Gaussian process GP and
     * two matrixes TX and X that contain input vectors to GP. Returns
     * covariance matrix C. Every element ij of C contains covariance between
     * inputs i in TX and j in X. This is a mandatory subfunction used for
     * example in prediction and energy computations.
     * @param x1 ???
     * @param x2 ???
     * @param C Covariance matrix
     */
    void calculateCovarianceMatrix(const Coord& x1, Coord& x2,
                                   Field<double>& C);

    /**
     * @brief Evaluate training variance vector.
     *
     * Function takes in covariance function of a Gaussian process GPCF and
     * matrix TX that contains training inputs. Returns variance vector C. Every
     * element i of C contains variance of input i in TX. This is a mandatory
     * subfunction used for example in prediction and energy computations.
     * @param x ???
     * @param C ???
     */
    void calculateTrainingCovarianceMatrix(const Coord& x, Field<double>& C);

    /**
     * @brief Evaluate gradient of covariance function with respect to the
     *        parameters.
     *
     * The function takes a covariance function structure GPCF, a matrix X of
     * input vectors and returns DKff, the gradients of covariance matrix Kff =
     * k(X,X) with respect to th (cell array with matrix elements). This is a
     * mandatory subfunction used for example in gradient computations.
     *
     * @param x2
     * @param x2
     * @param DKff
     */
    void calculateGradOfCovMatrix(const Coord& x1, Coord& x2,
                                  std::vector<Field<double> >& DKff);

    /**
     * @brief Evaluate gradient of covariance function, of which has been taken
     *        partial derivative with respect to x1, with respect to parameters.
     *
     * The function takes a covariance function structure GPCF, a matrix X of
     * input vectors and returns DKff, the gradients of derivatived covariance
     * matrix dK(df,f)/dhyp = d(d k(X,X)/dx)/dhyp, with respect to the
     * parameters.
     *
     * Evaluate: DKff{1:m} = d Kff / d magnSigma2
     *           DKff{m+1:2m} = d Kff / d lengthScale_m
     * m is the dimension of inputs. If ARD is used, then multiple
     * lengthScales. This subfunction is needed when using derivative
     * observations.
     *
     *      dims - is a vector of input dimensions with respect to which the
     *             derivatives of the covariance function have been calculated
     *             [by default dims=1:size(x,2)]
     *
     * Note! When coding the derivatives of the covariance function, remember
     * to double check them. See gp_cov for lines of code to check the
     * matrices
     *
     * @param x1
     * @param x2
     * @param dims
     * @param DKff
     */
    void calculateGradOfCovMatrixWithDerivatives(
        const Coord& x1, Coord& x2, Field<Index_t>& dims,
        std::vector<Field<double> >& DKff);

    /**
     * @brief Evaluate gradient of covariance function, of which has been taken
     * partial derivatives with respect to both input variables x1 and x2 with
     * respect to parameters.
     *
     * The function takes a covariance function structure GPCF, a matrix X of
     * input vectors and returns DKff, the gradients of derivative covariance
     * matrix dK(df,df)/dhyp = d(d^2 k(X1,X2)/dX1dX2)/dhyp with respect to the
     * parameters.
     *
     * Evaluate: DKff{1-m} = d Kff / d magnSigma2
     *           DKff{m+1-2m} = d Kff / d lengthScale_m
     * m is the dimension of inputs. If ARD is used, then multiple
     * lengthScales. This subfunction is needed when using derivative
     * observations.
     *
     * @param x1
     * @param x2
     * @param dims1
     * @param dims2
     * @param DKff
     */
    void calculateGradOfCovMatrixWithDerivatives2(
        const Coord& x1, Coord& x2, Field<Index_t>& dims1,
        Field<Index_t>& dims2, std::vector<Field<double> >& DKff);

    // TODO: ginput2, ginput3, ginput4 can be merged
    /**
     * @brief Calculates covariances between the derivatives
     *
     * @param x1 Matrix of coordinates
     * @param x2 Matrix of coordinates
     * @param dims
     * @param DKff Covariances between the derivatives
     */
    void ginput2(const Coord& x1, Coord& x2, Field<Index_t>& dims,
                 std::vector<Field<double> >& DKff);

    /**
     * @brief Evaluate gradient of covariance function with respect to both
     * input variables x1 and x2 (in different dimensions).
     *
     * This function takes a covariance function structure GPCF, a matrix X of
     * input vectors and returns DKff, the gradients of twice derivatived
     * covariance matrix K(df,df) = dk(X1,X2)/dX1dX2 (cell array with matrix
     * elements). The derivative is calculated in multidimensional problem
     * between input's observation dimensions which are not same. This
     * subfunction is needed when using derivative observations.
     *
     * @param x1
     * @param x2
     * @param dims1
     * @param dims2
     * @param DKff
     */
    void ginput3(const Coord& x1, Coord& x2, Field<Index_t>& dims1,
                 Field<Index_t>& dims2, std::vector<Field<double> >& DKff);

    /**
     * @brief Evaluate gradient of covariance function with respect to x.
     *
     * Simplified and faster version of sexp_ginput, returns full matrices.
     *
     * Returns DKff, the gradients of covariance matrix Kff = k(X,X2) with
     * respect to dimensions DIMS of X.
     *
     * @param x1
     * @param x2
     * @param dims
     * @param DKff
     */
    // TODO: check if x1 is always equal to x2!
    void ginput4(const Coord& x1, Coord& x2, Field<Index_t>& dims,
                 std::vector<Field<double> >& DKff);

    inline void setMagnSigma2(double value);

    inline void setLengthScale(double value);

    inline void setConfInfo(const AtomsConfiguration& value);

    /**
     * @brief Set  parameters of \e this covariance function  using the provided
     * vector.
     * @param w Vector of combined parameters from covariance functions and
     * likelyhood
     */
    inline void setParameters(const Eigen::VectorXd& w);

    inline double getMagnSigma2();

    inline Field<double>& getLengthScaleRef();

    inline void clear();

    inline void setPriorParametersGaussian(const PriorBase& prior);

    inline void setPriorParametersSqrtt(const PriorBase& prior);

    inline PriorBase getPriorParametersGaussian();

    inline PriorBase getPriorParametersSqrtt();

    /**
     * @brief Combine parameters of the covariance function into one Eigen
     * vector.
     */
    inline Eigen::VectorXd combineParameters();

private:
    void truncateValue(double& value);

private:
    bool initialized;
    double magnSigma2;
    Field<double> lengthScale;
    AtomsConfiguration conf_info;

    PriorBase prior_parameters_gaussian;
    PriorBase prior_parameters_sqrtt;
};

} /* namespace gpr */

#include "SexpatCF.inl"

#endif /* GPR_SEXPAT_H_ */
