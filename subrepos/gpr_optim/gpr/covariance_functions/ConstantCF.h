#ifndef GPR_GPCF_CONSTANT_H
#define GPR_GPCF_CONSTANT_H

#include <stdio.h>

#include <vector>

#include "../../data_types/Coord.h"
#include "../../data_types/Field.h"
#include "../../structures/Structures.h"

namespace gpr {

/**
 * @brief Constant  covariance function.
 */
class ConstantCF {
public:
    ConstantCF();
    ~ConstantCF() { }

    /**
     * @brief Evaluate the log prior of covariance function parameters.
     *
     * @return Result
     */
    double calculateLogPrior();

    /**
     * @brief Evaluate gradient of the log prior with respect to the parameters.
     *
     * @return Result
     */
    Field<double> calculateLogPriorGradient();

    /**
     * @brief The function takes a matrix x of input vectors and returns DKff,
     * the gradients of covariance matrix Kff = k(X,X2) with respect to th (cell
     * array with matrix elements).
     *
     * This is a mandatory subfunction used in gradient computations.
     *
     * @param x
     * @param x2
     * @param DKff
     */
    void calculateGradOfCovMatrix(const Coord &x, Coord &x2,
                                  std::vector<Field<double>> &DKff);

    /**
     * @brief Evaluate gradient of covariance function, of which has been taken
     * partial derivative with respect to x, with respect to parameters.
     *
     * @param x
     * @param x2
     * @param dims
     * @param DKff
     */
    void calculateGradOfCovMatrixWithDerivatives(
        const Coord &x, Coord &x2, Field<Index_t> &dims,
        std::vector<Field<double>> &DKff);

    /**
     * @brief Evaluate gradient of covariance function, of which has been taken
     * partial derivatives with respect to both input variables x, with respect
     * to parameters.
     *
     * Returns DKff, the gradients of derivative covariance matrix
     * dK(df,df)/dhyp = d(d^2 k(X1,X2)/dX1dX2)/dhyp with respect to the
     * parameters.
     *
     * @param x
     * @param x2
     * @param dims1
     * @param dims2
     * @param DKff
     */
    void calculateGradOfCovMatrixWithDerivatives2(
        const Coord &x, Coord &x2, Field<Index_t> &dims1, Field<Index_t> &dims2,
        std::vector<Field<double>> &DKff);

    /**
     * @brief Evaluate gradient of covariance function with respect to both
     * input variables x and x2 (in same dimension).
     *
     * Returns DKff, the gradients of twice derivatived covariance matrix
     * K(df,df) = dk(X1,X2)/dX1dX2 (cell array with matrix elements). Input
     * variable's dimensions are expected to be same.
     *
     * @param x
     * @param x2
     * @param dims
     * @param DKff
     */
    void ginput2(const Coord &x, Coord &x2, Field<Index_t> &dims,
                 std::vector<Field<double>> &DKff);

    /**
     * @brief  Evaluate gradient of covariance function with
     * respect to both input variables x and x2 (in
     * different dimensions).
     *
     * DKff = GPCF_CONSTANT_GINPUT3(GPCF, X, X2) takes a covariance
     * function structure GPCF, a matrix X of input vectors and
     * returns DKff, the gradients of twice derivatived covariance
     * matrix K(df,df) = dk(X1,X2)/dX1dX2 (cell array with matrix
     * elements). The derivative is calculated in multidimensional
     * problem between input's observation dimensions which are not
     * same. This subfunction is needed when using derivative
     * observations.
     *
     *
     * @param x
     * @param x2
     * @param dims1
     * @param dims2
     * @param DKff
     */
    void ginput3(const Coord &x, Coord &x2, Field<Index_t> &dims1,
                 Field<Index_t> &dims2, std::vector<Field<double>> &DKff);

    /**
     * @brief Evaluate gradient of covariance function with respect to x.
     * Simplified and faster version of constant_ginput, returns full matrices.
     *
     * Returns DKff, the gradients of covariance matrix Kff = k(X,X2) with
     * respect to X (whole matrix). This function is needed when using
     * derivative observations.
     *
     * @param x
     * @param x2
     * @param dims
     * @param DKff
     */
    void ginput4(const Coord &x, Coord &x2, Field<Index_t> &dims,
                 std::vector<Field<double>> &DKff);

    /**
     * @brief Evaluate covariance matrix between two input vectors
     *
     *    C = GP_CONSTANT_COV(GP, TX, X) takes in covariance function
     *    of a Gaussian process GP and two matrixes TX and X that
     *    contain input vectors to GP. Returns covariance matrix C.
     *    Every element ij of C contains covariance between inputs i
     *    in TX and j in X. This is a mandatory subfunction used for
     *    example in prediction and energy computations.
     *
     * @param x
     * @param x2
     * @param C
     */
    void calculateCovarianceMatrix(const Coord &x, Coord &x2, Field<double> &C);

    /**
     * @brief Set the constSigma2 value.
     */
    inline void setConstSigma2(const double _constSigma2);

    /**
     * @brief Set  parameters of \e this covariance function using the provided
     * vector.
     * @param w Vector of hyperparameters
     */
    inline void setParameters(const Eigen::VectorXd &w);

    inline double getConstSigma2();

    /**
     * @brief Combine parameters of the covariance function into one Eigen
     * vector.
     */
    inline Eigen::VectorXd combineParameters();

private:
    double constSigma2;
};

} /* namespace gpr */

#include "ConstantCF.inl"

#endif /* GPR_GPCF_CONSTANT_H */
