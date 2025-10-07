//
//  AtomicDimerInit.cpp
//  gpr_dimer
//
//  Created by Maxim Masterov on 16/12/2020.
//

#include "AtomicDimer.h"
#include "Enums.h"

namespace atmd {

void AtomicDimer::initialize(const gpr::InputParameters& parameters,
                             const gpr::Observation& init_observations,
                             const gpr::Observation& init_middle_point,
                             const gpr::Coord& init_orientation,
                             const gpr::AtomsConfiguration& init_atoms_config)
{
    if (parameters.method_rot.value == "LBFGS_alg") method_rot = LBFGS_alg;
    if (parameters.method_trans.value == "LBFGS_alg") method_trans = LBFGS_alg;

    stop_citeria_dimer.force = parameters.T_dimer.value;
    stop_citeria_dimer.angle_rotation = parameters.T_anglerot_init.value;

    stop_citeria_gpr.angle_rotation = parameters.T_anglerot_gp.value;
    stop_citeria_gpr.force = parameters.T_dimer.value;
    divisor_stop_criteria_gpr = parameters.divisor_T_dimer_gp.value;

    max_iter_init_rot.outer = parameters.num_iter_initrot.value;
    max_iter_init_rot.inner = parameters.num_iter.value;

    max_iter_new_pairs.outer = parameters.num_bigiter.value;
    max_iter_new_pairs.inner = parameters.num_iter.value;

    max_iter_relax_rot.outer = 1;
    max_iter_relax_rot.inner = parameters.num_iter_rot_gp.value;

    if (parameters.islarge_num_iter.value == 1)
        assume_many_iterations = true;
    else
        assume_many_iterations = false;

    dimer_sep = parameters.dimer_sep.value;
    ratio_at_limit = parameters.ratio_at_limit.value;
    actdist_fro = parameters.actdist_fro.value;

    // ====================================================================== //
    // ====================================================================== //
    gpr::Coord atoms_coords;
    gpr::Field<gpr::Index_t> freezed_atoms;
    cell_dimensions[0].set(parameters.cell_dimensions.value[0],
                           parameters.cell_dimensions.value[1],
                           parameters.cell_dimensions.value[2]);
    cell_dimensions[1].set(parameters.cell_dimensions.value[3],
                           parameters.cell_dimensions.value[4],
                           parameters.cell_dimensions.value[5]);
    cell_dimensions[2].set(parameters.cell_dimensions.value[6],
                           parameters.cell_dimensions.value[7],
                           parameters.cell_dimensions.value[8]);

    orient_init = init_orientation;

    middle_point_init = init_middle_point;
    all_obs = init_observations;

    atoms_config = init_atoms_config;

    gpr_model->initialize(parameters, atoms_config);

    transition_param.step_length = parameters.param_trans.value[0];
    transition_param.max_step_length = parameters.param_trans.value[1];
    transition_param.rotrem_thresh =
        parameters.rotation_removal_projection_threshold.value;

    use_prune = parameters.use_prune.value;
    start_prune_at = parameters.start_prune_at.value;
    nprune_vals = parameters.nprune_vals.value;
    prune_threshold = parameters.prune_threshold.value;

    early_params._dist_metric = parameters.es_dist_metric.value;
    early_params._threshold = parameters.es_threshold.value;

    debug_output_info.out_dir = parameters.debug_output_dir.value;
    debug_output_info.file_name_E = parameters.debug_output_file_E.value;
    debug_output_info.file_name_R = parameters.debug_output_file_R.value;
    debug_output_info.file_name_G = parameters.debug_output_file_G.value;
    debug_output_info.file_extension =
        parameters.debug_output_file_extension.value;
    debug_output_info.offset = parameters.debug_offset_from_mid_point.value;
    debug_output_info.dy = parameters.debug_dy.value;
    debug_output_info.dz = parameters.debug_dz.value;
    debug_level = parameters.debug_level.value;
}

} /* namespace atmd */
