#include <ceres/ceres.h>

#include <cfloat>
#include <vector>

namespace funcmin {

template <typename ClassName, typename FuncName>
class GprCostFunction : public ceres::CostFunction {
public:
    GprCostFunction(const gpr::EigenMatrix& x_data,
                    const Eigen::VectorXd& x_ind_data,
                    const Eigen::VectorXd& y_data, FuncName func,
                    ClassName& gpr_instance)
        : x_(x_data),
          x_ind_(x_ind_data),
          y_(y_data),
          func_to_min_(func),
          holder_(gpr_instance),
          normalization_factor_(-1.0)
    {
        const int n_vars = holder_.sexpat_cov_func->combineParameters().size();
        mutable_parameter_block_sizes()->push_back(n_vars);
        set_num_residuals(1);
    }

    virtual ~GprCostFunction() { }

    virtual bool Evaluate(double const* const* parameters, double* residuals,
                          double** jacobians) const override
    {
        const int n_vars = parameter_block_sizes()[0];
        Eigen::Map<const Eigen::VectorXd> w_current(parameters[0], n_vars);

        gpr::EnergyAndGradient eg;
        double raw_energy;
        eg.energy = &raw_energy;

        if (jacobians != nullptr && jacobians[0] != nullptr) {
            Eigen::VectorXd raw_gradient(n_vars);
            eg.gradient = &raw_gradient;

            // 1. Call your function to get the raw, unscaled energy and
            // gradient.
            (holder_.*func_to_min_)(w_current, x_, x_ind_, y_, eg);

            // 2. On the first evaluation, compute the normalization factor.
            if (normalization_factor_ < 0.0) {
                normalization_factor_ = raw_gradient.norm();
                if (normalization_factor_ < DBL_EPSILON) {
                    normalization_factor_ = 1.0;
                }
            }

            // 3. Scale BOTH the energy (residual) and the gradient (Jacobian).
            residuals[0] = raw_energy / normalization_factor_;
            Eigen::Map<Eigen::VectorXd> jacobian_map(jacobians[0], n_vars);
            jacobian_map = raw_gradient / normalization_factor_;

            // --- THE HACK HAS BEEN REMOVED ---
            // const double sign_correction = (residuals[0] < 0.0) ? -1.0 : 1.0;
            // jacobian_map *= sign_correction;

        } else {
            // Gradient is not requested.
            Eigen::VectorXd grad_dummy;
            eg.gradient = &grad_dummy;
            (holder_.*func_to_min_)(w_current, x_, x_ind_, y_, eg);

            // Still scale the residual even if jacobian is not requested.
            if (normalization_factor_ > 0.0) {
                residuals[0] = raw_energy / normalization_factor_;
            } else {
                residuals[0] = raw_energy;  // Not initialized yet
            }
        }

        return std::isfinite(residuals[0]);
    }

private:
    const gpr::EigenMatrix& x_;
    const Eigen::VectorXd& x_ind_;
    const Eigen::VectorXd& y_;
    FuncName func_to_min_;
    ClassName& holder_;
    mutable double normalization_factor_;
};

template <typename ClassName, typename FuncName>
void Ceres::optimize(const gpr::EigenMatrix& x_data,
                     const Eigen::VectorXd& x_ind, const Eigen::VectorXd& y,
                     Eigen::VectorXd& w, FuncName func_to_min,
                     ClassName& holder, double)
{
    const int n_vars = w.size();
    this->failedOptim = true;
    if (n_vars == 0) {
        this->failedOptim = false;
        return;
    }

    Eigen::VectorXd lower_bounds(n_vars);
    Eigen::VectorXd upper_bounds(n_vars);
    for (int i = 0; i < n_vars; ++i) {
        lower_bounds(i) = std::log(0.2);
        upper_bounds(i) = std::log(3.0);
    }
    upper_bounds(0) = std::log(10);

    ceres::Problem problem;
    ceres::CostFunction* cost_function =
        new GprCostFunction<ClassName, FuncName>(x_data, x_ind, y, func_to_min,
                                                 holder);
    problem.AddResidualBlock(cost_function, new ceres::TrivialLoss(), w.data());

    for (int i = 0; i < n_vars; ++i) {
        problem.SetParameterLowerBound(w.data(), i, lower_bounds(i));
        problem.SetParameterUpperBound(w.data(), i, upper_bounds(i));
    }

    ceres::Solver::Options options;
    options.max_num_iterations = settings.max_iter;
    options.function_tolerance = settings.tolerance_func;
    options.gradient_tolerance = settings.tolerance_sol;
    options.linear_solver_type = ceres::DENSE_QR;
    options.line_search_direction_type = ceres::LBFGS;
    options.minimizer_progress_to_stdout = (settings.report_level >= 2);

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    if (summary.IsSolutionUsable()) {
        this->failedOptim = false;
    }
}

inline void Ceres::setAlgorithmSettings(
    const gpr::OptimizationAlgorithmSettings& _settings)
{
    settings = _settings;
}

}  // namespace funcmin
