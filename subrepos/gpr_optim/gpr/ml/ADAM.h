#ifndef GPR_ML_ADAM_H_
#define GPR_ML_ADAM_H_

#include <Eigen/Dense>

#include "../../backend/Macros.h"
#include "../../data_types/Coord.h"
#include "../../structures/Structures.h"

namespace funcmin {

/**
 * @brief Adaptive Moment Estimation (ADAM) optimizer.
 */
class ADAM {
public:
    ADAM();
    virtual ~ADAM();

    // This flag will be set to true if optimization fails to converge.
    bool failedOptim;

    /**
     * @brief Optimization function using the ADAM algorithm.
     * @param x Training inputs (passed to func_to_min)
     * @param x_ind Training input indices (passed to func_to_min)
     * @param y Training targets (passed to func_to_min)
     * @param w On input, the initial parameters; on output, the optimized
     * parameters.
     * @param func_to_min The function to minimize. It must calculate energy and
     * gradient.
     * @param holder Object of the class that owns the function func_to_min.
     */
    template <typename ClassName, typename FuncName>
    void optimize(const gpr::EigenMatrix& x, const Eigen::VectorXd& x_ind,
                  const Eigen::VectorXd& y, Eigen::VectorXd& w,
                  FuncName func_to_min, ClassName& holder, double);

    /**
     * @brief Set parameters of the optimization algorithm.
     * For ADAM, this includes learning_rate, beta1, beta2, epsilon, etc.
     */
    inline void setAlgorithmSettings(
        const gpr::OptimizationAlgorithmSettings& _settings);

private:
    gpr::OptimizationAlgorithmSettings settings;
};

} /* namespace funcmin */

#include "ADAM.inl"

#endif /* GPR_ML_ADAM_H_ */
