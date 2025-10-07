#ifndef ADAM_inl
#define ADAM_inl

#include <cmath>
#include <cfloat>

namespace funcmin {

inline void ADAM::setAlgorithmSettings(
    const gpr::OptimizationAlgorithmSettings& _settings)
{
    settings = _settings;
}

template <typename ClassName, typename FuncName>
void ADAM::optimize(const gpr::EigenMatrix& x, const Eigen::VectorXd& x_ind,
                    const Eigen::VectorXd& y, Eigen::VectorXd& w,
                    FuncName func_to_min, ClassName& holder, double)
{
    gpr::Index_t nparams = (gpr::Index_t)w.rows();
    gpr::io::LogManager log_man;

    Eigen::VectorXd m = Eigen::VectorXd::Zero(nparams);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(nparams);
    Eigen::VectorXd g(nparams);

    Eigen::VectorXd v_max = Eigen::VectorXd::Zero(nparams);

    double current_energy = 0.0;
    double last_energy = std::numeric_limits<double>::max();
    double current_lr = settings.learning_rate;

    double grad_norm = 1.0;

    this->failedOptim = true;
    if (settings.report_level >= 1) log_man << "\n";

    const double one_minus_beta1 = 1.0 - settings.beta1;
    const double one_minus_beta2 = 1.0 - settings.beta2;
    double beta1_t = 1.0;
    double beta2_t = 1.0;

    for (gpr::Index_t t = 1; t <= settings.max_iter; ++t) {
        gpr::EnergyAndGradient energy_and_gradient;
        energy_and_gradient.energy = &current_energy;
        energy_and_gradient.gradient = &g;
        (holder.*func_to_min)(w, x, x_ind, y, energy_and_gradient);

        if (g.array().isNaN().any() || std::isnan(current_energy)) {
            if (settings.report_level >= 1)
                log_man << "ADAM Error: NaN encountered. Stopping.\n";
            return;
        }

        // --- NORMALIZATION START ---
        // On the first step, determine the gradient's characteristic scale.
        if (t == 1) {
            grad_norm = g.norm();
            // If the initial gradient is effectively zero, avoid normalization.
            if (grad_norm < DBL_EPSILON) {
                grad_norm = 1.0;
            }
        }
        // Normalize the current gradient by the initial gradient's norm.
        g /= grad_norm;
        // --- NORMALIZATION END ---

        // Apply weight decay (L2 regularization) to the normalized gradient
        if (settings.weight_decay > 0.0) {
            g += settings.weight_decay * w;
        }

        // Standard ADAM moment updates with the normalized gradient
        m = settings.beta1 * m + one_minus_beta1 * g;
        v = settings.beta2 * v + one_minus_beta2 * g.array().square().matrix();

        beta1_t *= settings.beta1;
        beta2_t *= settings.beta2;

        Eigen::VectorXd m_hat = m / (1.0 - beta1_t);
        Eigen::VectorXd v_hat;

        if (settings.amsgrad) {
            v_max = v_max.array().max(v.array()).matrix();
            v_hat = v_max / (1.0 - beta2_t);
        } else {
            v_hat = v / (1.0 - beta2_t);
        }

        Eigen::VectorXd w_update = (current_lr * m_hat.array() /
                                    (v_hat.array().sqrt() + settings.epsilon))
                                       .matrix();
        w -= w_update;

        current_lr *= settings.learning_rate_decay;

        if (settings.report_level >= 2 && (t % 10 == 0 || t == 1)) {
            log_man << " ADAM iteration: " << t
                    << ", Energy: " << current_energy
                    << ", Update norm: " << w_update.norm()
                    << ", LR: " << current_lr << "\n";
        }

        if (w_update.norm() < settings.tolerance_sol) {
            if (settings.report_level >= 1)
                log_man
                    << "ADAM Converged: Parameter change below tolerance.\n";
            this->failedOptim = false;
            break;
        }

        if (t > 1 &&
            std::abs(last_energy - current_energy) < settings.tolerance_func) {
            if (settings.report_level >= 1)
                log_man << "ADAM Converged: Function value change below "
                           "tolerance.\n";
            this->failedOptim = false;
            break;
        }

        last_energy = current_energy;
    }

    if (this->failedOptim && settings.report_level >= 1) {
        log_man << "ADAM Warning: Maximum number of iterations has been "
                   "exceeded.\n";
    }

    if (settings.report_level >= 1) log_man << "\n";
}

} /* namespace funcmin */

#endif /* ADAM_inl */
