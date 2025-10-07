#include "../../gpr/dimer/Dimer.h"
#include "gtest/gtest.h"

namespace gpr {
namespace tests {

class DimerProjectionTest : public ::testing::Test { };

// Test that a pure translational vector is completely removed.
TEST_F(DimerProjectionTest, PureTranslation)
{
    gpr::Coord R;
    R.resize(1, 3 * 3);  // A 3-atom system (e.g., water)
    R.set(0, 0, {0.0, 0.0, 0.0});
    R.set(0, 1, {0.95, 0.0, 0.0});
    R.set(0, 2, {0.0, 0.95, 0.0});

    // A step where every atom moves 0.5 Angstrom in the x-direction.
    Eigen::VectorXd step(9);
    step << 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0;

    double original_norm = step.norm();
    double removed_magnitude =
        dimer::project_out_rot_trans_with_feedback(R, step);

    EXPECT_NEAR(removed_magnitude, original_norm, 1e-9);
    EXPECT_NEAR(step.norm(), 0.0, 1e-9);
}

// Test that a pure rotational vector is completely removed.
TEST_F(DimerProjectionTest, PureRotation)
{
    gpr::Coord R;
    R.resize(1, 3 * 3);  // A 3-atom system
    R.set(0, 0, {0.0, 0.0, 0.0});
    R.set(0, 1, {1.0, 0.0, 0.0});
    R.set(0, 2, {0.0, 1.0, 0.0});

    // A step representing a small rotation around the z-axis
    Eigen::VectorXd step(9);
    step << 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, -0.1, 0.0, 0.0;

    double original_norm = step.norm();
    double removed_magnitude =
        dimer::project_out_rot_trans_with_feedback(R, step);

    EXPECT_NEAR(removed_magnitude, original_norm, 1e-9);
    EXPECT_NEAR(step.norm(), 0.0, 1e-9);
}

// Test that a pure internal (vibrational) motion is completely unaffected.
TEST_F(DimerProjectionTest, PureInternalMotion)
{
    gpr::Coord R;
    R.resize(1, 3 * 2);             // A 2-atom system (diatomic)
    R.set(0, 0, {0.0, 0.0, -0.5});  // Atom 1
    R.set(0, 1, {0.0, 0.0, 0.5});   // Atom 2

    // A step representing a pure bond stretch.
    // This motion conserves center of mass and has zero angular momentum.
    Eigen::VectorXd step(6);
    step << 0.0, 0.0, -0.1,  // Atom 1 moves away
        0.0, 0.0, 0.1;       // Atom 2 moves away

    Eigen::VectorXd original_step = step;
    double original_norm = step.norm();

    double removed_magnitude =
        dimer::project_out_rot_trans_with_feedback(R, step);

    // The magnitude of the removed component should be zero.
    EXPECT_NEAR(removed_magnitude, 0.0, 1e-9);
    // The step vector should be unchanged.
    EXPECT_NEAR(step.norm(), original_norm, 1e-9);
    ASSERT_EQ(step.size(), original_step.size());
    for (int i = 0; i < step.size(); ++i) {
        EXPECT_NEAR(step(i), original_step(i), 1e-9);
    }
}

// Test that for a mixed motion, only the rot/trans components are removed.
TEST_F(DimerProjectionTest, MixedMotion)
{
    gpr::Coord R;
    R.resize(1, 3 * 2);  // A 2-atom system
    R.set(0, 0, {0.0, 0.0, -0.5});
    R.set(0, 1, {0.0, 0.0, 0.5});

    // A mixed step: translation in x + bond stretch
    Eigen::VectorXd translational_part(6);
    translational_part << 0.5, 0.0, 0.0, 0.5, 0.0, 0.0;

    Eigen::VectorXd internal_part(6);
    internal_part << 0.0, 0.0, -0.1, 0.0, 0.0, 0.1;

    Eigen::VectorXd step = translational_part + internal_part;

    double removed_magnitude =
        dimer::project_out_rot_trans_with_feedback(R, step);

    // The magnitude removed should match the norm of the translational part.
    EXPECT_NEAR(removed_magnitude, translational_part.norm(), 1e-9);
    // The final step vector should be equal to the internal part.
    ASSERT_EQ(step.size(), internal_part.size());
    for (int i = 0; i < step.size(); ++i) {
        EXPECT_NEAR(step(i), internal_part(i), 1e-9);
    }
}

// Test that the projection works correctly for a linear molecule
// (rank-deficient case).
TEST_F(DimerProjectionTest, LinearMolecule)
{
    gpr::Coord R;
    R.resize(1, 3 * 3);
    R.set(0, 0, {0.0, 0.0, 0.0});
    R.set(0, 1, {0.0, 0.0, 1.2});
    R.set(0, 2, {0.0, 0.0, -1.2});

    // A mixed step: translation in x + asymmetric stretch in z
    Eigen::VectorXd translational_part(9);
    translational_part << 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0;

    Eigen::VectorXd internal_part(9);
    internal_part << 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, -0.1;

    Eigen::VectorXd step = translational_part + internal_part;

    double removed_magnitude =
        dimer::project_out_rot_trans_with_feedback(R, step);

    // The magnitude removed should match the norm of the translational part.
    EXPECT_NEAR(removed_magnitude, translational_part.norm(), 1e-9);
    // The final step vector should be equal to the internal part.
    ASSERT_EQ(step.size(), internal_part.size());
    for (int i = 0; i < step.size(); ++i) {
        EXPECT_NEAR(step(i), internal_part(i), 1e-9);
    }
}

}  // namespace tests
}  // namespace gpr
