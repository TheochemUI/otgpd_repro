#include "CapnpParamsTest.h"

#include <capnp/message.h>
#include <capnp/serialize.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../../structures/Structures.h"
#include "managers/io/FileManager.h"
#include "managers/io/CapnProtoManager.h"

namespace gpr {
namespace tests {

CapnpParamsTest::CapnpParamsTest()
{
    // TODO Auto-generated constructor stub
}

CapnpParamsTest::~CapnpParamsTest()
{
    // TODO Auto-generated destructor stub
}

TEST(ParameterLoading, LegacyVsCapnp)
{
    // --- Phase 1: Load parameters using the legacy FileManager ---
    gpr::InputParameters params_legacy;
    gpr::io::FileManager fm;
    ASSERT_NO_THROW({ fm.readInputFile("input/input.dat", params_legacy); })
        << "Test prerequisite failed: could not parse legacy input.dat";

    // --- Phase 2: Load parameters using the new Cap'n Proto loader ---
    gpr::InputParameters params_capnp;
    ASSERT_NO_THROW({
        gpr::io::loadParametersFromCapnp("input/capnp_params.bin", params_capnp);
    }) << "The new Cap'n Proto loader threw an exception.";

    // --- Phase 3: Assert that the two structures are identical ---

    // Integers & Bools
    EXPECT_EQ(params_legacy.i_dist.value, params_capnp.i_dist.value);
    EXPECT_EQ(params_legacy.i_run.value, params_capnp.i_run.value);
    EXPECT_EQ(params_legacy.actdist_fro.value, params_capnp.actdist_fro.value);
    EXPECT_EQ(params_legacy.num_iter_initrot.value,
              params_capnp.num_iter_initrot.value);
    EXPECT_EQ(params_legacy.num_iter_rot_gp.value,
              params_capnp.num_iter_rot_gp.value);
    EXPECT_EQ(params_legacy.num_bigiter.value, params_capnp.num_bigiter.value);
    EXPECT_EQ(params_legacy.num_iter.value, params_capnp.num_iter.value);
    EXPECT_EQ(params_legacy.islarge_num_iter.value,
              params_capnp.islarge_num_iter.value);
    EXPECT_EQ(params_legacy.max_iter.value, params_capnp.max_iter.value);
    EXPECT_EQ(params_legacy.start_prune_at.value,
              params_capnp.start_prune_at.value);
    EXPECT_EQ(params_legacy.nprune_vals.value, params_capnp.nprune_vals.value);
    EXPECT_EQ(params_legacy.divisor_T_dimer_gp.value,
              params_capnp.divisor_T_dimer_gp.value);
    EXPECT_EQ(params_legacy.report_level.value,
              params_capnp.report_level.value);
    EXPECT_EQ(params_legacy.debug_level.value, params_capnp.debug_level.value);
    EXPECT_EQ(static_cast<bool>(params_legacy.initrot_nogp.value),
              static_cast<bool>(params_capnp.initrot_nogp.value));
    EXPECT_EQ(params_legacy.use_prune.value, params_capnp.use_prune.value);

    // Doubles
    EXPECT_DOUBLE_EQ(params_legacy.dimer_sep.value,
                     params_capnp.dimer_sep.value);
    EXPECT_DOUBLE_EQ(params_legacy.T_dimer.value, params_capnp.T_dimer.value);
    EXPECT_DOUBLE_EQ(params_legacy.T_anglerot_init.value,
                     params_capnp.T_anglerot_init.value);
    EXPECT_DOUBLE_EQ(params_legacy.T_anglerot_gp.value,
                     params_capnp.T_anglerot_gp.value);
    EXPECT_DOUBLE_EQ(params_legacy.disp_max.value, params_capnp.disp_max.value);
    EXPECT_DOUBLE_EQ(params_legacy.ratio_at_limit.value,
                     params_capnp.ratio_at_limit.value);
    EXPECT_DOUBLE_EQ(params_legacy.gp_sigma2.value,
                     params_capnp.gp_sigma2.value);
    EXPECT_DOUBLE_EQ(params_legacy.jitter_sigma2.value,
                     params_capnp.jitter_sigma2.value);
    EXPECT_DOUBLE_EQ(params_legacy.sigma2.value, params_capnp.sigma2.value);
    EXPECT_DOUBLE_EQ(params_legacy.prior_mu.value, params_capnp.prior_mu.value);
    EXPECT_DOUBLE_EQ(params_legacy.prior_nu.value, params_capnp.prior_nu.value);
    EXPECT_DOUBLE_EQ(params_legacy.prior_s2.value, params_capnp.prior_s2.value);
    EXPECT_DOUBLE_EQ(params_legacy.lambda.value, params_capnp.lambda.value);
    EXPECT_DOUBLE_EQ(params_legacy.lambda_limit.value,
                     params_capnp.lambda_limit.value);
    EXPECT_DOUBLE_EQ(params_legacy.tolerance_func.value,
                     params_capnp.tolerance_func.value);
    EXPECT_DOUBLE_EQ(params_legacy.tolerance_sol.value,
                     params_capnp.tolerance_sol.value);
    EXPECT_DOUBLE_EQ(params_legacy.prune_threshold.value,
                     params_capnp.prune_threshold.value);
    EXPECT_DOUBLE_EQ(params_legacy.rotation_removal_projection_threshold.value,
                     params_capnp.rotation_removal_projection_threshold.value);

    // Strings
    EXPECT_STREQ(params_legacy.method_rot.value.c_str(),
                 params_capnp.method_rot.value.c_str());
    EXPECT_STREQ(params_legacy.method_trans.value.c_str(),
                 params_capnp.method_trans.value.c_str());
    EXPECT_STREQ(params_legacy.optimization_alg.value.c_str(),
                 params_capnp.optimization_alg.value.c_str());
    EXPECT_STREQ(params_legacy.check_derivative.value.c_str(),
                 params_capnp.check_derivative.value.c_str());
    EXPECT_STREQ(params_legacy.debug_output_dir.value.c_str(),
                 params_capnp.debug_output_dir.value.c_str());

    // Arrays
    for (int i = 0; i < 4; ++i) {
        EXPECT_DOUBLE_EQ(params_legacy.dist_sp.value[i],
                         params_capnp.dist_sp.value[i]);
    }
    for (int i = 0; i < 2; ++i) {
        EXPECT_DOUBLE_EQ(params_legacy.param_trans.value[i],
                         params_capnp.param_trans.value[i]);
    }
    for (int i = 0; i < 9; ++i) {
        EXPECT_DOUBLE_EQ(params_legacy.cell_dimensions.value[i],
                         params_capnp.cell_dimensions.value[i]);
    }
}

} /* namespace tests */
} /* namespace gpr */
