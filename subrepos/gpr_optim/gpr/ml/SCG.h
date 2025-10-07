/*
 * SCG.h
 *
 *  Created on: 20 Nov 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_ML_SCG_H_
#define GPR_ML_SCG_H_

#include <Eigen/Dense>

#include "../../backend/Macros.h"
#include "../../data_types/Coord.h"
#include "../../structures/Structures.h"

/**
 * \namespace fmin
 * @brief Function minimization algorithms
 */
namespace funcmin {

/**
 * @brief Scaled Conjugate Gradient method.
 */
class SCG {
public:
    SCG();
    bool failedOptim;
    virtual ~SCG();

    /**
     * @brief Optimization function.
     * Incorporates Squared Conjugate Gradient method. See "A Scaled Conjugate
     * Gradient Algorithm for Fast Supervised Learning", MARTIN FODSLETTE
     * MEILLER, 1993, Neural Networks, vol. 6, pp 525-533
     * @param x Training inputs
     * @param y Training targets
     * @param w Combination of all parameters of covariance functions and
     * likelihood
     * @param func_to_min Function pointer. Should consist of a full name of the
     * function that will be used during the minimization process (e.g.
     * ClassName::FunctionName)
     * @param holder Object of the class that owns the function \e func_to_min
     */
    template <typename ClassName, typename FuncName>
    void optimize(const gpr::EigenMatrix& x, const Eigen::VectorXd& x_ind,
                  const Eigen::VectorXd& y, Eigen::VectorXd& w,
                  FuncName func_to_min, ClassName& holder, double);

    /**
     * @brief Set parameters of the optimization algorithm.
     */
    inline void setAlgorithmSettings(
        const gpr::OptimizationAlgorithmSettings& _settings);

private:
    /**
     * @brief Look for a maximum absolute value in the vector.
     */
    inline double findMaxAbsValue(const Eigen::VectorXd& vec);

    /**
     * @brief Chack is the vector has any NaN coefficients.
     */
    inline bool isNanCoeff(const Eigen::VectorXd& vec);

    /**
     * @brief Check if the vector ha any Inf coefficients.
     */
    inline bool isInfCoeff(const Eigen::VectorXd& vec);

private:
    gpr::OptimizationAlgorithmSettings settings;
};

} /* namespace funcmin */

#include "SCG.inl"

#endif /* GPR_ML_SCG_H_ */
