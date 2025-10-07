/*
 * LBFGS.h
 *
 *  Created on: 30 Oct 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_DIMER_H_
#define GPR_DIMER_H_

#include <Eigen/Dense>

#include "../../data_types/Coord.h"
#include "../../structures/Structures.h"
#include "../ml/GaussianProcessRegression.h"

namespace dimer {

/**
 * @brief Projects rigid-body translations and rotations out of a Cartesian
 * step vector using QR decomposition.
 *
 * This function constructs a basis for the 3 translational and 3 rotational
 * degrees of freedom. It then uses Eigen's numerically stable Householder
 * QR decomposition to find an orthonormal basis for this subspace. Finally,
 * it subtracts the component of the step vector that lies within this
 * subspace.
 *
 * @param R The current atomic coordinates, used to define the rotational
 * basis.
 * @param step The 3N-dimensional optimization step vector to be modified
 * in-place.
 * @return The L2 norm of the rotational/translational component that was
 * removed.
 */
double project_out_rot_trans_with_feedback(const gpr::Coord& R,
                                           Eigen::VectorXd& step);

/**
 * @brief Dimer method.
 */
class Dimer {
public:
    /**
     * @brief Default constructor.
     *
     * @param separation Dimer separation
     * @param has_frozen_atoms Flag to indicate if the system is a subset of a
     * larger system
     */
    Dimer(double separation, bool is_subsystem = false)
        : dimer_sep(separation), is_subsystem(is_subsystem) { };

    virtual ~Dimer() { }

    /**
     * @brief Translate dimer.
     *
     * The dimer is translated one step towards saddle point according to the
     * L-BFGS method.
     *
     * @param R Coordinates of the middle point  of the dimer.
     * @param orient Unit vector along the direction of the dimer.
     * @param F Force at the middle point  of the dimer.
     * @param curv Curvature of energy along the direction of the dimer.
     * @param param_trans Predefined step length for convex regions, maximum
     * step length.
     * @param transinfo Structure with information on previous external
     * iterations.
     * @param R_new Coordinates of the new middle point  of the dimer.
     */
    void translate(const gpr::Coord& R, const gpr::Coord& orient,
                   const gpr::Coord& F, const double curv,
                   const gpr::TransitionParameters& param_trans,
                   gpr::LBFGSInfo& transinfo, gpr::Coord& R_new);

    /**
     * @brief Rotate dimer.
     *
     * The dimer is rotated one step towards its minimum energy orientation
     * according to the modified Newton method on a rotation plane chosen based
     * on the L-BFGS method. The rotation angle is then optimized according to
     * the modified Newton method based on a finite difference estimation of the
     * scalar rotational force along the direction of the rotation.
     *
     * @param R Coordinates of the middle point  of the dimer.
     * @param orient Unit vector along the direction of the dimer.
     * @param G01 Gradient at the middle point and image 1 of the dimer.
     * @param potential Potential and gradient function.
     * @param T_anglerot Convergence threshold for the rotation angle.
     * @param estim_Curv If \e true, an estimate for the curvature along the
     * direction of the dimer after the rotation is calculated.
     * @param all_obs Object of all observed points.
     * @param gpr_model Pointer to the GPR model.
     * @param rotinfo Structure with information on previous external
     * iterations.
     * @param orient_new Unit vector along the direction of the dimer after
     * optimal rotation.
     * @param Curv Estimate for the curvature along the direction of the dimer
     * after the rotation (0 if estim_Curv = 0 or no rotation).
     * @param image1 Observation of the new observed location.
     */
    void rotate(const gpr::Coord& R, const gpr::Coord& orient,
                const gpr::Coord& G01, const uint8_t potential,
                const double T_anglerot, const bool estim_Curv,
                const gpr::Observation& all_obs,
                gpr::GaussianProcessRegression* gpr_model,
                gpr::LBFGSInfo& rotinfo, gpr::Coord& orient_new, double& Curv,
                gpr::Observation& image1);

    /**
     * @brief Calculate rotational force acting on image 1 of the dimer.
     *
     * @param G01 Gradient vectors at the middle point and image 1 of the dimer.
     * @param orient Unit vector along the direction of the dimer.
     * @param F_rot Rotational force acting on image 1 of the dimer.
     */
    void rotateForce(const gpr::Coord& G01, const gpr::Coord& orient,
                     gpr::Coord& F_rot);

    /**
     * @brief Calculate translational force acting on image 1 of the dimer.
     *
     * @param F Force at the middle point of the dimer.
     * @param orient Unit vector along the direction of the dimer.
     * @param F_trans Translational force acting on image 1 of the dimer.
     */
    void translateForce(const gpr::Coord& F, const gpr::Coord& orient,
                        gpr::Coord& F_trans);

    /**
     * @brief Calculate a rough estimate of the optimal Rotation angle.
     *
     * @param orient Unit vector along the direction of the dimer (before
     * rotation).
     * @param G01 Gradient at the middle point and image 1 of the dimer (before
     * rotation). Each point occupies one row in the field.
     * @param F_rot_oriented Oriented rotational force acting on image 1 of the
     * dimer.
     * @return Rough estimate of the optimal Rotation angle.
     */
    double estimateRotationalAngle(const gpr::Coord& orient,
                                   const gpr::Coord& G01,
                                   const gpr::Coord& F_rot_oriented);

private:
    /**
     * @brief Translate dimer for negative curvature.
     *
     * @param orient Unit vector along the direction of the dimer.
     * @param F Force at the middle point of the dimer.
     * @param param_trans Predefined step length for convex regions, maximum
     * step length.
     * @param transinfo Structure with information on previous external
     * iterations.
     * @param R Coordinates of the new middle point of the dimer.
     */
    void translateDimerForNegCurv(const gpr::Coord& orient, const gpr::Coord& F,
                                  const gpr::TransitionParameters& param_trans,
                                  gpr::LBFGSInfo& transinfo, gpr::Coord& R);

    /**
     * @brief Translate dimer for positive curvature.
     *
     * @param orient Unit vector along the direction of the dimer.
     * @param F Force at the middle point of the dimer.
     * @param param_trans Predefined step length for convex regions, maximum
     * step length.
     * @param transinfo Structure with information on previous external
     * iterations.
     * @param R Coordinates of the new middle point of the dimer.
     */
    void translateDimerForPosCurv(const gpr::Coord& orient, const gpr::Coord& F,
                                  const gpr::TransitionParameters& param_trans,
                                  gpr::LBFGSInfo& transinfo, gpr::Coord& R);

    /**
     * @brief Update \e LBFGSInfo structure for rotational force.
     *
     * @param observation Observation of the new observed location.
     * @param orient Previous unit vector along the direction of the dimer.
     * @param orient_new New unit vector along the direction of the dimer.
     * @param F_rot Rotational force acting on image 1 of the dimer.
     * @param rotinfo Updated structure with information on previous external
     * iterations.
     */
    void updateLBFGSInfoStructRot(const gpr::Observation& observation,
                                  const gpr::Coord& orient,
                                  const gpr::Coord& orient_new,
                                  gpr::Coord& F_rot, gpr::LBFGSInfo& rotinfo);

    /**
     * @brief Update \e LBFGSInfo structure for translational force.
     *
     * @param param_trans Predefined step length for convex regions, maximum
     * step length.
     * @param F_trans Translational force acting on image 1 of the dimer.
     * @param r New unit vector along the direction of the dimer.
     * @param rotinfo Updated structure with information on previous external
     * iterations.
     */
    void updateLBFGSInfoStructTrans(
        const gpr::TransitionParameters& param_trans, const gpr::Coord& F_trans,
        Eigen::VectorXd& r, gpr::LBFGSInfo& transinfo);

    /**
     * @brief Calculate curvature of the energy along the dimer.
     *
     * @param G01 Gradient at the middle point and image 1 of the dimer.
     * @param orient Unit vector along the direction of the dimer.
     * @return Curvature.
     */
    double calculateCurvature(const gpr::Coord& G01, const gpr::Coord& orient);

    /**
     * @brief Optimize rotation angle with modified Newton method.
     *
     * The rotation angle is optimized according to the modified Newton method
     * based on a finite difference estimation of the scalar rotational force
     * along the direction of the rotation.
     *
     * @param R Coordinates of the middle point of the dimer.
     * @param orient Unit vector along the direction of the dimer (before
     * rotation).
     * @param G01 Gradient at the middle point and image 1 of the dimer (before
     * rotation). Each point occupies one row in the field.
     * @param potential Potential and gradient function.
     * @param T_anglerot Convergence threshold for the rotation angle.
     * @param estim_Curv If \e true, an estimate for the curvature along the
     * direction of the dimer after the rotation is calculated
     * @param estim_G1 If \e true, an estimate for the gradient at image 1 after
     * the rotation is calculated
     * @param r New unit vector along the direction of the dimer.
     * @param F_rot Rotational force acting on image 1 of the dimer.
     * @param gpr_model Pointer to the GPR model.
     * @param orient_new Unit vector along the direction of the dimer after
     * optimal rotation.
     * @param Curv Estimate for the curvature along the direction of the dimer
     * after the rotation (0 if estim_Curv is false)
     * @param image1 Observation of image 1 of the dimer after a small test
     * rotation
     */
    void optimizeRotationAngle(const gpr::Coord& R, const gpr::Coord& orient,
                               const gpr::Coord& G01, const uint8_t potential,
                               const double T_anglerot, const bool estim_Curv,
                               const bool estim_G1, const Eigen::VectorXd& r,
                               const gpr::Coord& F_rot,
                               gpr::GaussianProcessRegression* gpr_model,
                               gpr::Coord& orient_new, double& Curv,
                               gpr::Observation& image1);

    /**
     * @brief Calculate potential.
     *
     * @param potential Potential type.
     * @param gpr_model Pointer to the GPR model (or nullptr).
     * @param dtheta_obs Observed energy and gradients.
     */
    void calculatePotential(const uint8_t potential,
                            gpr::GaussianProcessRegression* gpr_model,
                            gpr::Observation& image1);

    void calculateOrientedRotationalForce(const gpr::Coord& orient,
                                          const Eigen::VectorXd& r,
                                          const gpr::Coord& F_rot,
                                          const gpr::Coord& G01,
                                          gpr::Coord& F_rot_oriented);

    /**
     * @brief Calculate optimal Rotation angle.
     *
     * @param F_rot Rotational force acting on image 1 of the dimer.
     * @param orient_rot Orientation of the rotation on image 1 of the dimer.
     * @param F_rot_oriented Oriented rotational force acting on image 1 of the
     * dimer.
     * @param omega_est Estimated optimal Rotation angle.
     * @return Optimal Rotation angle.
     */
    double calculateRotationAngle(const gpr::Coord& F_rot_omega_est,
                                  const gpr::Coord& orient_rot_omega_est,
                                  const gpr::Coord& F_rot_oriented,
                                  const double omega_est, double& a1,
                                  double& b1);

    /**
     * @brief Calculate orientation vector of the dimer after the rotation.
     *
     * @param orient Old orientation vector.
     * @param orient_rot Old rotation direction.
     * @param omega Rotation angle.
     * @param orient_new New orientation.
     */
    void calculateOrientationVector(const gpr::Coord& orient,
                                    const gpr::Coord orient_rot,
                                    const double omega, gpr::Coord& orient_new);

    /**
     * @brief Calculate the rotation direction in the end of the rotation.
     *
     * @param orient Old orientation vector.
     * @param orient_rot Old rotation direction.
     * @param omega Rotation angle.
     * @param orient_new New rotation direction.
     */
    void calculateRotationDirection(const gpr::Coord& orient,
                                    const gpr::Coord orient_rot,
                                    const double omega,
                                    gpr::Coord& orient_rot_new);

private:
    /**
     * Dimer separation (distance from the middle point of the dimer to the two
     * images)
     */
    double dimer_sep;
    /**
     * Flag for deactivating the rotation removal, since it is not robust when
     * applied to parts of a larger system
     * Practically this is just a check for if frozen atoms exist
     */
    bool is_subsystem;
};

} /* namespace dimer */

#endif /* GPR_DIMER_H_ */
