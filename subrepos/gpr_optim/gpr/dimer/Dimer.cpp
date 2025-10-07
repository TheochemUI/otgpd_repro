/*
 * LBFGS.cpp
 *
 *  Created on: 30 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "Dimer.h"

#include "../../managers/io/ErrorManager.h"
#include "../optimization/LBFGS.h"

namespace dimer {
using gpr::Index_t;

double project_out_rot_trans_with_feedback(const gpr::Coord& R,
                                           Eigen::VectorXd& step)
{
    const gpr::Index_t num_atoms = R.getNumCols() / 3;
    const gpr::Index_t dof = num_atoms * 3;

    if (num_atoms == 0) {
        throw;
    }

    // 1. Calculate the Center of Mass (COM)
    gpr::vector3_reg com(0.0, 0.0, 0.0);
    for (Index_t i = 0; i < num_atoms; ++i) {
        com += R.at(0, i);
    }
    com /= static_cast<double>(num_atoms);

    // 2. Construct the 6 basis vectors for translation and rotation
    std::vector<Eigen::VectorXd> basis;
    basis.reserve(6);

    // Translational basis vectors
    for (int i = 0; i < 3; ++i) {
        Eigen::VectorXd t = Eigen::VectorXd::Zero(dof);
        for (Index_t j = 0; j < num_atoms; ++j) {
            t(3 * j + i) = 1.0;
        }
        basis.push_back(t);
    }

    // Rotational basis vectors (using coordinates relative to COM)
    Eigen::VectorXd r_x = Eigen::VectorXd::Zero(dof);
    Eigen::VectorXd r_y = Eigen::VectorXd::Zero(dof);
    Eigen::VectorXd r_z = Eigen::VectorXd::Zero(dof);

    for (Index_t i = 0; i < num_atoms; ++i) {
        gpr::vector3_reg pos = R.at(0, i) - com;
        // Rotation around x-axis: dR = a * (i x R) = a * (0, -z, y)
        r_x(3 * i + 1) = -pos.z;
        r_x(3 * i + 2) = pos.y;
        // Rotation around y-axis: dR = a * (j x R) = a * (z, 0, -x)
        r_y(3 * i + 0) = pos.z;
        r_y(3 * i + 2) = -pos.x;
        // Rotation around z-axis: dR = a * (k x R) = a * (-y, x, 0)
        r_z(3 * i + 0) = -pos.y;
        r_z(3 * i + 1) = pos.x;
    }
    basis.push_back(r_x);
    basis.push_back(r_y);
    basis.push_back(r_z);

    // 3. Orthonormalize the basis using Gram-Schmidt
    std::vector<Eigen::VectorXd> ortho_basis;
    ortho_basis.reserve(6);

    for (auto& v: basis) {
        Eigen::VectorXd u = v;
        for (const auto& ou: ortho_basis) {
            u -= v.dot(ou) * ou;
        }

        // This check correctly handles linear molecules, where one rotational
        // basis vector will be linearly dependent on the others.
        if (u.norm() > 1e-9) {
            u.normalize();
            ortho_basis.push_back(u);
        }
    }

    // Calculate the projection and its magnitude before modifying the step
    Eigen::VectorXd removed_component = Eigen::VectorXd::Zero(dof);
    for (const auto& u: ortho_basis) {
        removed_component += step.dot(u) * u;
    }

    // Now, modify the step
    step -= removed_component;

    // Return the norm of what was removed
    return removed_component.norm();
}

void Dimer::translate(const gpr::Coord& R, const gpr::Coord& orient,
                      const gpr::Coord& F, const double curv,
                      const gpr::TransitionParameters& param_trans,
                      gpr::LBFGSInfo& transinfo, gpr::Coord& R_new)
{
    R_new.clear();
    R_new = R;

    if (curv < 0.) {
        translateDimerForNegCurv(orient, F, param_trans, transinfo, R_new);
    } else {
        translateDimerForPosCurv(orient, F, param_trans, transinfo, R_new);
    }
}

void Dimer::translateDimerForNegCurv(
    const gpr::Coord& orient, const gpr::Coord& F,
    const gpr::TransitionParameters& param_trans, gpr::LBFGSInfo& transinfo,
    gpr::Coord& R)
{
    gpr::Index_t m = transinfo.deltaOrient.getNumCols();
    gpr::Coord F_trans;
    gpr::optim::LBFGS lbfgs;
    Eigen::VectorXd r_raw;
    translateForce(F, orient, F_trans);

    if (m > 0) {
        transinfo.deltaF.append(F_trans - transinfo.F_old);
    }

    lbfgs.execute(F_trans, transinfo, r_raw);
    Eigen::VectorXd r_projected = r_raw;
    Eigen::VectorXd r_final;

    if (this->is_subsystem) {
        r_final = r_raw;
    } else {
        Eigen::VectorXd r_projected = r_raw;
        double removed_magnitude =
            project_out_rot_trans_with_feedback(R, r_projected);

        if (removed_magnitude > param_trans.rotrem_thresh) {
            // The rotational component is large. So the GPR model is likely
            // producing a large, unphysical torque. Use the original,
            // unprojected step. This large rotation will be caught by the
            // existing step length checks, triggering the original resampling
            // mechanism.
            // NOTE(rg): In practice this should never happen, just remove all
            // rotations with a gaint threshold
            r_final = r_raw;
        } else {
            // The rotational component is small. The GPR is behaving.
            // We use the cleaned-up, projected step for a more precise
            // movement.
            r_final = r_projected;
        }
    }

    updateLBFGSInfoStructTrans(param_trans, F_trans, r_final, transinfo);
    R -= r_final;
}

void Dimer::translateDimerForPosCurv(
    const gpr::Coord& orient, const gpr::Coord& F,
    const gpr::TransitionParameters& param_trans, gpr::LBFGSInfo& transinfo,
    gpr::Coord& R)
{
    gpr::Coord F_trans;
    gpr::Coord orient_search;

    F_trans.resize(1, R.getNumCols());
    F_trans.setZero();

    F_trans = orient * F.dot(orient);
    F_trans *= -1.;
    orient_search = F_trans / F_trans.norm();
    R = R + orient_search * param_trans.step_length;

    transinfo.F_old.setZero();
    transinfo.deltaOrient.clear();
    transinfo.deltaF.clear();
}

void Dimer::rotate(const gpr::Coord& R, const gpr::Coord& orient,
                   const gpr::Coord& G01, const uint8_t potential,
                   const double T_anglerot, const bool estim_Curv,
                   const gpr::Observation& all_obs,
                   gpr::GaussianProcessRegression* gpr_model,
                   gpr::LBFGSInfo& rotinfo, gpr::Coord& orient_new,
                   double& Curv, gpr::Observation& image1)
{
    gpr::Index_t m = rotinfo.deltaOrient.getNumRows();
    gpr::Coord F_rot;
    gpr::optim::LBFGS lbfgs;
    Eigen::VectorXd r;

    rotateForce(G01, orient, F_rot);

    if (m > 0) {
        rotinfo.deltaF.append(F_rot - rotinfo.F_old);
    }

    // TODO(rg): As per the brms_dimer_rot paper this ought to be the CG
    // esp. without rotation removal..
    lbfgs.execute(F_rot, rotinfo, r);

    optimizeRotationAngle(R, orient, G01, potential, T_anglerot, estim_Curv,
                          false, r, F_rot, gpr_model, orient_new, Curv, image1);

    updateLBFGSInfoStructRot(image1, orient, orient_new, F_rot, rotinfo);
}

void Dimer::updateLBFGSInfoStructRot(const gpr::Observation& observation,
                                     const gpr::Coord& orient,
                                     const gpr::Coord& orient_new,
                                     gpr::Coord& F_rot, gpr::LBFGSInfo& rotinfo)
{
    gpr::Index_t m = rotinfo.deltaOrient.getNumRows();
    gpr::Index_t num_lbfgsiter_rot = rotinfo.num_lbfgs_iter;

    if (observation.R.isEmpty()) {
        rotinfo.deltaOrient.clear();
        rotinfo.deltaF.clear();
    } else {
        rotinfo.deltaOrient.append(orient_new - orient);
        if (m >= num_lbfgsiter_rot) {
            rotinfo.deltaOrient.deleteRow(0);
            rotinfo.deltaF.deleteRow(0);
        }
    }

    rotinfo.F_old = F_rot;
}

void Dimer::updateLBFGSInfoStructTrans(
    const gpr::TransitionParameters& param_trans, const gpr::Coord& F_trans,
    Eigen::VectorXd& r, gpr::LBFGSInfo& transinfo)
{
    gpr::Index_t m = transinfo.deltaOrient.getNumRows();
    double steplength = r.norm();

    if (steplength > param_trans.max_step_length) {
        r = r * param_trans.max_step_length / steplength;
        transinfo.F_old.setZero();
        transinfo.deltaOrient.clear();
        transinfo.deltaF.clear();
    } else {
        Eigen::VectorXd tmp(F_trans.getNumCols());
        tmp = -r;
        transinfo.deltaOrient.appendVector(tmp);
        if (m >= transinfo.num_lbfgs_iter) {
            transinfo.deltaOrient.deleteRow(0);
            transinfo.deltaF.deleteRow(0);
        }
        transinfo.F_old = F_trans;
    }
}

void Dimer::optimizeRotationAngle(
    const gpr::Coord& R, const gpr::Coord& orient, const gpr::Coord& G01,
    const uint8_t potential, const double T_anglerot, const bool estim_Curv,
    const bool estim_G1, const Eigen::VectorXd& r, const gpr::Coord& F_rot,
    gpr::GaussianProcessRegression* gpr_model, gpr::Coord& orient_new,
    double& Curv, gpr::Observation& image1)
{
    gpr::Coord F_rot_oriented;
    gpr::Coord orient_rot_new;
    gpr::Coord G1;
    double F_0 = 0.;
    double C_0 = 0.;
    double omega_est = 0.;  // Estimated rotational angle

    calculateOrientedRotationalForce(orient, r, F_rot, G01, F_rot_oriented);

    F_0 = F_rot_oriented.norm();
    C_0 = calculateCurvature(G01, orient);

    omega_est = estimateRotationalAngle(orient, G01, F_rot_oriented);

    if (omega_est < T_anglerot) {
        orient_new = orient;
        orient_rot_new.gpr::Field<double>::set(0.);
        Curv = 0.;
        G1.clear();
        image1.R.clear();
        image1.E.clear();
        image1.G.clear();
        std::cout << "Emptying R_obs\n";
    } else {
        gpr::Coord orient_rot;
        gpr::Coord orient_est;
        gpr::Coord orient_rot_est;
        gpr::Coord F_rot_est;
        gpr::Coord tmp_G;
        double a1;
        double b1;
        double omega;  // Rotational angle

        orient_rot = F_rot_oriented / F_0;

        // Estimate orientation and rotation direction
        calculateOrientationVector(orient, orient_rot, omega_est, orient_est);
        calculateRotationDirection(orient, orient_rot, omega_est,
                                   orient_rot_est);

        // New location of image 1
        image1.R = R + orient_est * dimer_sep;

        calculatePotential(potential, gpr_model, image1);

        tmp_G = G01;
        for (gpr::Index_t j = 0; j < tmp_G.getNumCols(); ++j)
            tmp_G(1, j) = image1.G(0, j);

        rotateForce(tmp_G, orient_est, F_rot_est);

        omega = calculateRotationAngle(F_rot_est, orient_rot_est,
                                       F_rot_oriented, omega_est, a1, b1);

        // Calculate orientation and rotation direction
        calculateOrientationVector(orient, orient_rot, omega, orient_new);
        calculateRotationDirection(orient, orient_rot, omega, orient_rot_new);

        // Calculate curvature energy along the new orientation
        if (estim_Curv)
            Curv = C_0 + a1 * (cos(2. * omega) - 1.) + b1 * sin(2. * omega);
        else
            Curv = 0.;

        if (estim_G1) {
            // Eq. 19 from "Minimum mode saddle point searches using Gaussian
            // process regression with inverse-distance covariance function",
            // Olli-Pekka Koistinen, Vilhjalmur Asgeirsson, Aki Vehtari,
            // Hannes Jonsson
            double tmp_1 = sin(omega_est - omega) / sin(omega_est);
            double tmp_2 = sin(omega) / sin(omega_est);
            double tmp_3 =
                (1. - cos(omega) - sin(omega) * tan(omega_est * 0.5));

            G1.resize(1, G01.getNumCols());

            for (gpr::Index_t j = 0; j < G1.getNumCols(); ++j) {
                G1(0, j) = tmp_1 * G01(1, j) + tmp_2 * image1.G(0, j) +
                           tmp_3 * G01(0, j);
            }
        } else {
            G1.clear();
        }
    }
}

void Dimer::rotateForce(const gpr::Coord& G01, const gpr::Coord& orient,
                        gpr::Coord& F_rot)
{
    gpr::Index_t D = G01.getNumCols();
    gpr::Field<double> F1(1, D);
    gpr::Field<double> F2(1, D);
    gpr::Field<double> F1_rot(1, D);
    gpr::Field<double> F2_rot(1, D);

    F_rot.resize(1, G01.getNumCols());

    F1 = G01.extractRowAsVector(1);
    F1 *= -1.;

    F2 = G01.extractRowAsVector(0);
    F2 *= -2.;
    F2 += G01.extractRowAsVector(1);

    F1_rot = F1 - orient * F1.dot(orient);
    F2_rot = F2 - orient * F2.dot(orient);

    F_rot = (F1_rot - F2_rot) / dimer_sep;
}

void Dimer::translateForce(const gpr::Coord& F, const gpr::Coord& orient,
                           gpr::Coord& F_trans)
{
    F_trans = F - orient * (F.dot(orient)) * 2.;
}

double Dimer::calculateCurvature(const gpr::Coord& G01,
                                 const gpr::Coord& orient)
{
    double C = 0.;

    for (gpr::Index_t n = 0; n < G01.getNumCols(); ++n) {
        C += (G01(1, n) - G01(0, n)) * orient(0, n);
    }
    C /= dimer_sep;

    return C;
}

void Dimer::calculatePotential(const uint8_t potential,
                               gpr::GaussianProcessRegression* gpr_model,
                               gpr::Observation& image1)
{
    gpr::io::ErrorManager err;

    switch (potential) {
        case POTENTIAL_GP:
            gpr::assertMsg(gpr_model != nullptr,
                           "Error! NULL pointer is passed!");
            gpr_model->calculatePotential(image1);
            break;

        // This case is used in unit tests
        case POTENTIAL_TEST1:
            image1.E.resize(1, 1);
            image1.G.resize(1, 2 * 3);

            image1.E(0, 0) = 140.981934401792e-006;

            image1.G(0, 0) = 22.0169292639464e-003;
            image1.G(0, 1) = -25.3799869163130e-003;
            image1.G(0, 2) = 20.4383773983493e-003;
            image1.G(0, 3) = -17.5501833843735e-003;
            image1.G(0, 4) = 21.9491835909580e-003;
            image1.G(0, 5) = -28.6553826196841e-003;

            break;

        // This case is used in unit tests
        case POTENTIAL_TEST2:
            image1.E.resize(1, 1);
            image1.G.resize(1, 2 * 3);

            image1.E(0, 0) = 268.422097239501e-006;

            image1.G(0, 0) = 21.4263025463240e-003;
            image1.G(0, 1) = -30.4813824740716e-003;
            image1.G(0, 2) = 26.0833530773849e-003;
            image1.G(0, 3) = -16.2806875633080e-003;
            image1.G(0, 4) = 21.0852193621926e-003;
            image1.G(0, 5) = -35.3650094447835e-003;

            break;

        default:
            err << "Error! The type of the potential is unknown!";

            break;
    }
}

void Dimer::calculateOrientedRotationalForce(const gpr::Coord& orient,
                                             const Eigen::VectorXd& r,
                                             const gpr::Coord& F_rot,
                                             const gpr::Coord& G01,
                                             gpr::Coord& F_rot_oriented)
{
    Eigen::VectorXd orient_rot(r.rows());
    double tmp = 0.;

    for (gpr::Index_t n = 0; n < r.rows(); ++n) {
        tmp += orient[n] * r(n);
    }
    for (gpr::Index_t n = 0; n < r.rows(); ++n) {
        orient_rot(n) = r(n) - tmp * orient[n];
    }

    orient_rot /= orient_rot.norm();
    F_rot_oriented.resize(1, (gpr::Index_t)orient_rot.rows());

    tmp = 0.;
    for (gpr::Index_t n = 0; n < orient_rot.rows(); ++n) {
        tmp += F_rot[n] * orient_rot(n);
    }
    for (gpr::Index_t n = 0; n < F_rot_oriented.getSize(); ++n) {
        F_rot_oriented[n] = tmp * orient_rot[n];
    }
}

double Dimer::estimateRotationalAngle(const gpr::Coord& orient,
                                      const gpr::Coord& G01,
                                      const gpr::Coord& F_rot_oriented)
{
    double F_0 = F_rot_oriented.norm();
    double C_0 = calculateCurvature(G01, orient);

    return 0.5 * atan(0.5 * F_0 / fabs(C_0));
}

double Dimer::calculateRotationAngle(const gpr::Coord& F_rot_omega_est,
                                     const gpr::Coord& orient_rot_omega_est,
                                     const gpr::Coord& F_rot_oriented,
                                     const double omega_est, double& a1,
                                     double& b1)
{
    double F_0 = F_rot_oriented.norm();
    double F_omega_est;
    double omega = 0.;

    F_omega_est = F_rot_omega_est.dot(orient_rot_omega_est);

    a1 = (F_omega_est - F_0 * cos(2. * omega_est)) / (2. * sin(2. * omega_est));
    b1 = -0.5 * F_0;
    omega = 0.5 * atan(b1 / a1);

    if (omega < 0.) omega = M_PI_2 + omega;

    return omega;
}

void Dimer::calculateOrientationVector(const gpr::Coord& orient,
                                       const gpr::Coord orient_rot,
                                       const double omega,
                                       gpr::Coord& orient_new)
{
    orient_new = orient * cos(omega) + orient_rot * sin(omega);
    orient_new /= orient_new.norm();
}

void Dimer::calculateRotationDirection(const gpr::Coord& orient,
                                       const gpr::Coord orient_rot,
                                       const double omega,
                                       gpr::Coord& orient_rot_new)
{
    orient_rot_new = orient_rot * cos(omega) - orient * sin(omega);
    orient_rot_new /= orient_rot_new.norm();
}

} /* namespace dimer */
