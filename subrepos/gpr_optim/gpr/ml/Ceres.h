#pragma once

#include <Eigen/Dense>

#include "../../backend/Macros.h"
#include "../../data_types/Coord.h"
#include "../../structures/Structures.h"

namespace funcmin {

/**
 * @brief Wrapper for the Ceres Solver for robust nonlinear optimization.
 *
 * This class uses Ceres's trust-region minimizer, which natively
 * handles the bounds constraints required for stable hyperparameter
 * optimization.
 */
class Ceres {
public:
    Ceres();
    virtual ~Ceres();
    bool failedOptim;

    template <typename ClassName, typename FuncName>
    void optimize(const gpr::EigenMatrix& x_data,
                         const Eigen::VectorXd& x_ind, const Eigen::VectorXd& y,
                         Eigen::VectorXd& w, FuncName func_to_min,
                         ClassName& holder, double);

    inline void setAlgorithmSettings(
        const gpr::OptimizationAlgorithmSettings& _settings);

private:
    gpr::OptimizationAlgorithmSettings settings;
};

}  // namespace funcmin

#include "Ceres.inl"
