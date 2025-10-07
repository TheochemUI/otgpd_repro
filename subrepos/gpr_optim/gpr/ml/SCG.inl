//
//  SCG.inl
//  gpr_dimer
//
//  Created by Maxim Masterov on 24/11/2020.
//

#ifndef SCG_h
#define SCG_h

#include <cfloat>
#include <cmath>

namespace funcmin {

template <typename ClassName, typename FuncName>
void SCG::optimize(const gpr::EigenMatrix& x, const Eigen::VectorXd& x_ind,
                   const Eigen::VectorXd& y, Eigen::VectorXd& w,
                   FuncName func_to_min, ClassName& holder,
                   double current_barrier_strength)
{
    double sigma0 = 1e-4;
    gpr::Index_t func_count = 1;
    gpr::Index_t grad_count = 1;
    gpr::Index_t fail_count = 0;
    double lambdamin = 1.0e-15;
    double lambdamax = 1.0e100;
    gpr::Index_t iter_counter = 0;
    bool success = true;
    this->failedOptim = false;
    gpr::Index_t nsuccess = 0;
    gpr::Index_t nparams = (gpr::Index_t)w.rows();
    gpr::io::LogManager log_man;

    Eigen::VectorXd r, r_new, r_old;  //, r_plus;
    Eigen::VectorXd p;
    Eigen::VectorXd x_new, x_plus;

    double delta = 0., gamma = 0., kappa = 0., Delta = 0.;
    double alpha = 0., mu = 0., beta = 0., sigma = 0.;

    double f_new = 0., f_old = 0.;
    double tmp = 0.;

    uint8_t exit_code = 0;

    // --- Lambda wrapper for adaptive barriers ---
    const double max_log_magnSigma2 = std::log(2.0);
    const double barrier_strength = std::min(current_barrier_strength, 0.5);
    if (settings.report_level >= 1)
        log_man << "Effective barrier strength " << barrier_strength << "\n";
    if (settings.report_level >= 1)
        log_man << "Barrier at " << max_log_magnSigma2 << "\n";

    // This lambda calls the original objective function and applies the barrier
    // penalty.
    auto evaluate_penalized = [&](const Eigen::VectorXd& w_current,
                                  double& energy, Eigen::VectorXd& gradient) {
        // 1. Call the original function to get the base energy and gradient
        gpr::EnergyAndGradient eg_wrapper;
        eg_wrapper.energy = &energy;
        eg_wrapper.gradient = &gradient;
        (holder.*func_to_min)(w_current, x, x_ind, y, eg_wrapper);
        func_count++;
        grad_count++;

        // 2. Apply the barrier penalty
        const double boundary_gap = max_log_magnSigma2 - w_current(0);

        if (boundary_gap <= 1e-9) {  // Use a small tolerance
            // Hit or crossed the boundary. This is an invalid state for the
            // optimizer.
            energy = std::numeric_limits<double>::infinity();
            // if (settings.report_level >= 2)
            //     log_man << "Invalid, hit barrier\n";
        } else {
            // Add the barrier's contribution to the energy and gradient
            energy -= barrier_strength * std::log(boundary_gap);
            gradient(0) += barrier_strength / boundary_gap;
            // if (settings.report_level >= 2)
            //     log_man << "Barrier contribution (-e,g) " << energy << ", "
            //             << gradient(0) << "\n";
        }
    };

    // Initial function value and gradient
    evaluate_penalized(w, f_old, r_old);

    r = r_old;
    p = -r;

    if (settings.report_level >= 1) log_man << "\n";

    // Main optimization loop
    while (iter_counter < settings.max_iter) {
        if (settings.report_level >= 1)
            log_man << " SCG iteration: " << iter_counter << "\n";

        // Step 2
        if (success) {
            mu = p.dot(r);
            if (mu >= 0.) {
                p = -r;
                mu = p.dot(r);
            }

            kappa = p.dot(p);

            if (kappa < DBL_EPSILON) {
                if (settings.report_level >= 2) {
                    log_man << " Gradient smaller than machine precission"
                            << "\n";
                }
                this->failedOptim = true;
                exit_code = 1;
                break;
            }

            sigma = sigma0 / sqrt(kappa);
            x_plus = w + sigma * p;
            evaluate_penalized(x_plus, tmp, r_new);
            while ((isInfCoeff(r_new) || isNanCoeff(r_new)) &&
                   !std::isnan(f_old)) {
                sigma = 2. * sigma;
                kappa = 0.25 * kappa;
                x_plus = w + sigma * p;
                // [tmp,gplus] = fun(xplus);
                evaluate_penalized(x_plus, tmp, r_new);
            }
            gamma = (p.dot(r_new - r)) / sigma;
        }

        // Increase effective curvature and evaluate step size alpha
        delta = gamma + settings.lambda * kappa;

        if (delta <= 0.) {
            delta = settings.lambda * kappa;
            settings.lambda -= gamma / kappa;
        }

        // Step 5
        alpha = -mu / delta;

        // Calculate the comparison ratio
        x_new = w + alpha * p;
        evaluate_penalized(x_new, f_new, r_new);

        while (std::isinf(f_new) || std::isnan(f_new)) {
            log_man << " Warning! Function value at xnew not finite "
                       "or a number"
                    << "\n";

            settings.lambda = std::min(4. * settings.lambda, lambdamax);
            delta = gamma + settings.lambda * kappa;
            if (delta <= 0) {
                delta = settings.lambda * kappa;
                settings.lambda = settings.lambda - gamma / kappa;
            }
            alpha = -mu / delta;
            x_new = w + alpha * p;
            evaluate_penalized(x_new, f_new, r_new);
            ++func_count;
            ++grad_count;
            fail_count++;
            if (fail_count > 100) {
                log_man << " Critical! Too many failures"
                        << "\n";
                this->failedOptim = true;
                break;
            }
        }

        // Step 6
        // check data types of f_new and f_old
        Delta = 2 * (f_new - f_old) / (alpha * mu);

        // Step 7
        if (Delta >= 0) {
            success = true;
            ++nsuccess;
            w = x_new;
        } else {
            success = false;
        }

        if (success) {
            if (findMaxAbsValue(p * alpha) < settings.tolerance_sol) {
                exit_code = 2;
                this->failedOptim = false;
                break;
            } else if (std::fabs(f_new - f_old) < settings.tolerance_func) {
                exit_code = 3;
                this->failedOptim = false;
                break;
            } else {
                // Update variables for new position
                f_old = f_new;

                r_old.swap(r);
                r.swap(r_new);

                // If the gradient is zero then we are done.
                if (r.dot(r) < DBL_EPSILON) {  //  && all(isreal(grad_new)
                    this->failedOptim = false;
                    exit_code = 1;
                    break;
                }
            }
        }

        // Step 8
        // Adjust lambda according to comparison ratio.
        if (Delta < 0.25) {
            settings.lambda = std::min(4.0 * settings.lambda, lambdamax);
        }
        if (Delta > 0.75) {
            settings.lambda = std::max(0.5 * settings.lambda, lambdamin);
        }

        // If scale parameter is at its limit, stop optimization
        if (settings.lambda >= settings.lambda_limit) {
            this->failedOptim = true;
            exit_code = 0;
            if (settings.report_level >= 1) {
                log_man << " Warning: Optimization stopped because lambda "
                           "parameter reached limit. Check that the analytic "
                           "gradients are correct!"
                        << "\n";
            }
            break;
        }

        // Update search direction using Polak-Ribiere formula, or re-start
        // in direction of negative gradient after nparams steps.
        if (nsuccess == nparams) {
            p = -r;
            nsuccess = 0;
        } else {
            if (success) {
                beta = (r_old - r).dot(r) / mu;
                p = beta * p - r;
            }
        }

        ++iter_counter;
    }

    // If we get here, then we haven't terminated in the given number of
    // iterations.
    exit_code = 0;
    if (settings.report_level >= 1 && iter_counter == settings.max_iter) {
        this->failedOptim = true;
        log_man << " Maximum number of iterations has been exceeded."
                << "\n";
    }

    if (settings.report_level >= 2) {
        log_man << " Func-count " << func_count << ". Final f(x)=" << f_new
                << "."
                << "\n";
    }

    if (settings.report_level >= 1) log_man << "\n";
}

inline double SCG::findMaxAbsValue(const Eigen::VectorXd& vec)
{
    return vec.lpNorm<Eigen::Infinity>();
}

inline bool SCG::isNanCoeff(const Eigen::VectorXd& vec)
{
    return vec.array().isNaN().any();
}

inline bool SCG::isInfCoeff(const Eigen::VectorXd& vec)
{
    return vec.array().isInf().any();
}

inline void SCG::setAlgorithmSettings(
    const gpr::OptimizationAlgorithmSettings& _settings)
{
    settings = _settings;
}

} /* namespace funcmin */

#endif /* SCG_h */
