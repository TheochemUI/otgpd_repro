#include "CapnProtoManager.h"

#include "gpr/Enums.h"
#include "schema/gprd_params.capnp.h"

namespace gpr::io {

void loadParametersFromCapnp(const std::string& binary_filename,
                             gpr::InputParameters& parameters)
{
    int fd = open(binary_filename.c_str(), O_RDONLY);
    if (fd < 0) {
        throw std::runtime_error("Error: Could not open file: " +
                                 binary_filename);
    }

    struct stat file_stats;
    if (fstat(fd, &file_stats) < 0) {
        close(fd);
        throw std::runtime_error("Error: Could not get file stats for: " +
                                 binary_filename);
    }
    off_t size = file_stats.st_size;

    if (size == 0) {
        close(fd);
        throw std::runtime_error("Error: Input file is empty: " +
                                 binary_filename);
    }
    const void* mapped_ptr = mmap(nullptr, size, PROT_READ, MAP_SHARED, fd, 0);
    if (mapped_ptr == MAP_FAILED) {
        close(fd);
        throw std::runtime_error("Error: Could not mmap file: " +
                                 binary_filename);
    }

    kj::ArrayPtr<const capnp::word> wordPtr(
        reinterpret_cast<const capnp::word*>(mapped_ptr),
        size / sizeof(capnp::word));
    capnp::FlatArrayMessageReader reader(wordPtr);
    auto paramsReader = reader.getRoot<::InputParameters>();

    // --- Adapt from Cap'n Proto Reader to legacy gpr::InputParameters struct
    // ---

    // --- Run Selector ---
    auto runSelectorReader = paramsReader.getRunSelector();
    parameters.i_dist.value = runSelectorReader.getIDist();
    parameters.i_run.value = runSelectorReader.getIRun();
    auto sepDistsReader = runSelectorReader.getSepDists();
    for (size_t i = 0; i < sepDistsReader.size() && i < 4; ++i) {
        parameters.dist_sp.value[i] = sepDistsReader[i];
    }

    // --- Problem Group ---
    auto problemReader = paramsReader.getProblem();
    parameters.actdist_fro.value = problemReader.getActdistFro();
    parameters.dimer_sep.value = problemReader.getDimerSep();
    if (problemReader.getMethodRot() == DimerOptimizer::LBFGS)
        parameters.method_rot.value = "LBFGS_alg";
    if (problemReader.getMethodTrans() == DimerOptimizer::LBFGS)
        parameters.method_trans.value = "LBFGS_alg";

    auto paramTransReader = problemReader.getParamTrans();
    parameters.param_trans.value[0] = paramTransReader.getStepLength();
    parameters.param_trans.value[1] = paramTransReader.getMaxStepLength();
    parameters.rotation_removal_projection_threshold.value =
        paramTransReader.getRotremThresh();

    // Final convergence on the true PES
    parameters.T_dimer.value = problemReader.getFinalConvergenceForce();

    // Initial rotations phase
    auto initialRotationsReader = problemReader.getInitialRotations();
    parameters.num_iter_initrot.value =
        initialRotationsReader.getMaxIterations();
    // Note the inversion: initrot_nogp=1 means NO GP, so useGpr=false
    parameters.initrot_nogp.value = initialRotationsReader.getUseGpr() ? 0 : 1;
    parameters.T_anglerot_init.value =
        initialRotationsReader.getConvergenceAngle();

    // Main relaxation phase
    auto relaxationReader = problemReader.getRelaxation();
    parameters.num_bigiter.value = relaxationReader.getMaxOuterIterations();
    parameters.num_iter.value = relaxationReader.getMaxInnerIterations();
    parameters.num_iter_rot_gp.value =
        relaxationReader.getMaxRotationIterations();
    parameters.islarge_num_iter.value =
        relaxationReader.getAssumeManyIterations() ? 1 : 0;
    // Note: 'inittrans_nogp' is not present in the legacy struct, but if it
    // were, the logic would be:
    // !relaxationReader.getUseGprForInitialTranslation()

    auto relaxationConvergenceReader = relaxationReader.getConvergence();
    parameters.T_anglerot_gp.value = relaxationConvergenceReader.getAngle();
    // The legacy code re-used T_dimer for the relaxation force criterion.
    // If you had a separate legacy variable, you'd map `getForce()` to it here.

    // Other problem parameters
    parameters.divisor_T_dimer_gp.value = problemReader.getDivisorTdimerGP();
    parameters.disp_max.value = problemReader.getMidpointMaxDisp();
    parameters.ratio_at_limit.value = problemReader.getRatioAtLimit();

    // --- GPR Model ---
    auto gprReader = paramsReader.getGpr();
    parameters.gp_sigma2.value = gprReader.getGpSigma2();
    parameters.jitter_sigma2.value = gprReader.getJitterSigma2();
    parameters.sigma2.value = gprReader.getNoiseSigma2();
    auto priorReader = gprReader.getPrior();
    parameters.prior_mu.value = priorReader.getMu();
    parameters.prior_nu.value = priorReader.getNu();
    parameters.prior_s2.value = priorReader.getS2();

    // --- Early stopping
    auto earlyStoppingReader = paramsReader.getEarlyStopping();
    switch (earlyStoppingReader.getDistMetric()) {
        case (::DistanceMetric::MAX1_D_LOG): {
            parameters.es_dist_metric.value = DistanceMetricType::MAX_1D_LOG;
            break;
        }
        case (::DistanceMetric::RMSD): {
            parameters.es_dist_metric.value = DistanceMetricType::RMSD;
            break;
        }
        case (::DistanceMetric::EMD): {
            parameters.es_dist_metric.value = DistanceMetricType::EMD;
            break;
        }
        default:
            throw("Shouldn't be here");
    }
    parameters.es_threshold.value = earlyStoppingReader.getThreshold();

    // --- Hyperparameter Optimizer ---
    auto optimizerReader = paramsReader.getOptimizer();
    parameters.check_derivative.value =
        optimizerReader.getCheckDerivative() ? "true" : "false";
    parameters.max_iter.value = optimizerReader.getMaxIter();
    parameters.tolerance_func.value = optimizerReader.getTolFunc();
    parameters.tolerance_sol.value = optimizerReader.getTolSol();

    switch (optimizerReader.getAlgorithm()) {
        case HyperparameterOptimizer::SCG: {
            parameters.optimization_alg.value = "SCG_opt";
            auto scgReader = optimizerReader.getScg();
            parameters.lambda_limit.value = scgReader.getLambdaLimit();
            parameters.lambda.value = scgReader.getLambda();
            break;
        }
        case HyperparameterOptimizer::ADAM: {
            parameters.optimization_alg.value =
                "ADAM_opt";  // Assuming this legacy string
            auto adamReader = optimizerReader.getAdam();
            parameters.learning_rate.value = adamReader.getLr();
            parameters.learning_rate_decay.value = adamReader.getLrd();
            parameters.beta1.value = adamReader.getB1();
            parameters.beta2.value = adamReader.getB2();
            parameters.epsilon.value = adamReader.getEps();
            parameters.weight_decay.value = adamReader.getWeightDecay();
            parameters.amsgrad.value = adamReader.getAmsgrad();
            break;
        }
    }

    // --- Pruning ---
    auto pruningReader = paramsReader.getPruning();
    parameters.use_prune.value = pruningReader.getUsePrune();
    parameters.start_prune_at.value = pruningReader.getPruneBegin();
    parameters.nprune_vals.value = pruningReader.getPruneNVals();
    parameters.prune_threshold.value = pruningReader.getPruneThreshold();

    // --- Debug ---
    auto debugReader = paramsReader.getDebug();
    parameters.report_level.value = debugReader.getReportLevel();
    parameters.debug_level.value = debugReader.getDebugLevel();
    parameters.debug_output_dir.value = debugReader.getDebugOutDir().cStr();
    parameters.debug_output_file_R.value = debugReader.getDebugPosFile().cStr();
    parameters.debug_output_file_E.value =
        debugReader.getDebugEnergyFile().cStr();
    parameters.debug_output_file_G.value =
        debugReader.getDebugGradFile().cStr();
    parameters.debug_output_file_extension.value =
        debugReader.getDebugOutExt().cStr();
    parameters.debug_offset_from_mid_point.value =
        debugReader.getDebugOffsetMidPoint();
    parameters.debug_dy.value = debugReader.getDebugDy();
    parameters.debug_dz.value = debugReader.getDebugDz();

    // --- Cell Dimensions ---
    auto cellDimsReader = paramsReader.getCellDimensions();
    for (size_t i = 0; i < cellDimsReader.size() && i < 9; ++i) {
        parameters.cell_dimensions.value[i] = cellDimsReader[i];
    }

    // Cleanup
    munmap(const_cast<void*>(mapped_ptr), size);
    close(fd);
}
}  // namespace gpr::io
