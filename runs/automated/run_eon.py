#!/usr/bin/env python3

from io import StringIO
import uuid
import shutil
from pathlib import Path
from zipfile import ZipFile
import logging
import subprocess
import os

import ase.io
from ase.io import read as aseread
from ase import Atoms
from ase.constraints import FixAtoms
import ase.data

from rgpycrumbs._aux import getstrform, get_gitroot, switchdir

import eon.akmc
from eon.explorer import MinModeExplorer
from eon.config import ConfigClass as eonConf
import eon.fileio as eio

import click

# XXX: eon.config.config is a leaky global state.. why.


def setup_logger():
    """Sets up the logger."""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


def generate_hyperparameter_labels(atoms: Atoms) -> list[str]:
    """Generates hyperparameter labels from atom pairs (e.g., H-H, H-C)."""
    fixed_indices = []
    if atoms.constraints:
        for constr in atoms.constraints:
            if isinstance(constr, FixAtoms):
                fixed_indices.extend(constr.get_indices())

    all_numbers = atoms.get_atomic_numbers()
    moving_numbers = [n for i, n in enumerate(all_numbers) if i not in fixed_indices]

    if not moving_numbers:
        return []

    unique_moving_numbers = sorted(list(set(moving_numbers)))
    symbols = [ase.data.chemical_symbols[n] for n in unique_moving_numbers]

    labels = []
    for i in range(len(symbols)):
        for j in range(i, len(symbols)):
            labels.append(f"{symbols[i]}-{symbols[j]}")
    return labels


def prepare_con_file(run_dir, spin, mol_idx_str, gitroot) -> Atoms:
    """Prepares the pos.con file and returns the ase.Atoms object."""
    with ZipFile(gitroot / "data" / "sella_si_data.zip", "r") as zdat:
        with zdat.open(f"sella_si/{spin}/{mol_idx_str}.xyz", "r") as atmdat:
            atoms = aseread(StringIO(atmdat.read().decode()), format="xyz")
            # Centering to handle gh-188 in eOn
            atoms.set_cell([25, 25, 25])
            atoms.center()
            ase.io.write(getstrform(run_dir / "pos.con"), atoms)
    return atoms


def prepare_nwchem_settings(run_dir: Path, spin: str) -> Path:
    """Writes the appropriate nwchem settings file based on spin and returns its path."""
    if spin == "singlet":
        filename = "nwchem_singlet.nwi"
        content = """#-*- mode: conf -*-

basis noprint
  * library 3-21G
end

scf
  nopen 0
  thresh 1e-8
  maxiter 200
end
"""
    elif spin == "doublet":
        filename = "nwchem_doublet.nwi"
        content = """#-*- mode: conf -*-

basis noprint
  * library 3-21G
end

scf
  uhf
  nopen 1
  thresh 1e-8
  maxiter 200
end
"""
    else:
        raise ValueError(f"Invalid spin '{spin}' provided for NWChem settings.")

    settings_path = run_dir / filename
    with open(settings_path, "w") as f:
        f.write(content)
    return settings_path


def prepare_config_ini(
    run_dir: Path, mol_idx: int, atoms: Atoms, mult: int, nwchem_settings_path: Path
):
    """
    Prepares the config.ini file using the SocketNWChem potential.
    """
    unique_tag = f"{uuid.uuid4().hex[:6]}"
    socket_name = f"eon_{unique_tag}"

    # Calculate prior_degrees_of_freedom based on the number of unique pair types + 1
    # XXX(rg): Needs to be higher for convergence !???? o.O xref S034
    # NOTE(rg): The actual prior is drawn from (Ediff/3)**2 as per the OP thesis
    # Taking the MATLAB value as the minimum dof here
    pair_labels = generate_hyperparameter_labels(atoms)
    prior_dof = max(28, 2 * (len(pair_labels) + 1))

    config_content = f"""
[Main]
job = saddle_search
temperature = 300
random_seed = 1995
finite_difference = 0.01

[Potential]
potential = SocketNWChem

[SocketNWChemPot]
unix_socket_path = {socket_name}
unix_socket_mode = true
nwchem_settings = {nwchem_settings_path.absolute()}
make_template_input = false

[Optimizer]
opt_method = lbfgs
converged_force = 0.01
max_iterations = 1000
max_move = 0.05
convergence_metric = norm

[LBFGS]
lbfgs_memory = 25
lbfgs_inverse_curvature = 0.01
lbfgs_auto_scale = true

[Saddle Search]
displace_softest_mode_weight = 1.0
displace_magnitude = 0.01
min_mode_method = gprdimer
max_energy = 10.0

[GPR Dimer]
finite_angle = 0.05
relaxation_converged_angle = 0.1
max_outer_iterations = 300
max_midpoint_displacement = 0.5
rotation_opt_method = lbfgs
translation_opt_method = lbfgs
active_radius = 30.3
dimer_separation = 0.01
convex_region_step_size = 0.1
max_step_size = 0.05
has_many_iterations = true
hyperparameter_opt_method = SCG_opt
gpr_noise_variance = 1e-5
prior_mean = 0.0
prior_variance = 1.0
opt_max_iterations = 400
# Set these to 1e-3 for more robust convergence
# xref S089, no longer required, with barrier
opt_tol_sol = 1e-2
opt_tol_func = 1e-2
scg_lambda_limit = 1e16
scg_lambda_init = 100
gpr_jitter_variance = 0
report_level = 2
debug_level = 0
converged_angle = 0.0873
ratio_at_limit = 0.666666666667
divisor_t_dimer = 10
# This defines how close to a Gaussian the distribution's prior on magsigma2 is
# MATLAB uses 20, JCTC uses 28
prior_degrees_of_freedom = {prior_dof}
nogp_initial_rotations = true
gpr_variance = 1e-7
max_inner_iterations = 200

max_initial_rotation_iterations = 6
max_relaxation_rotation_iterations = 10

rotation_removal_projection_threshold = 50000
# this is actually overridden
es_threshold = 0.4
es_dist_metric = emd

fps_metric = emd
fps_history = 10

[Debug]
write_movies = True
"""
    with open(run_dir / "config.ini", "w") as f:
        f.write(config_content)


def setup_initial_displacement(run_dir, logger, input_dir):
    """Sets up the initial displacement by running eOn's MinModeExplorer or
    by using files from an input directory."""

    # Check if a valid input directory is provided
    if input_dir and Path(input_dir).exists():
        pos_path = Path(input_dir) / "pos.con"
        dir_path = Path(input_dir) / "direction.dat"
        disp_path = Path(input_dir) / "displacement.con"

        if pos_path.exists() and dir_path.exists() and disp_path.exists():
            shutil.copy(pos_path, run_dir / "pos.con")
            shutil.copy(dir_path, run_dir / "direction.dat")
            shutil.copy(disp_path, run_dir / "displacement.con")
            logger.info(f"Using displacement files from {input_dir}")
            return

    # Fallback to MinModeExplorer if no valid input_dir is found
    logger.warning("No valid input directory with all required files found. Falling back to MinModeExplorer.")
    with switchdir(run_dir):
        econf = eonConf()
        econf.init("config.ini")

        kT = econf.main_temperature / 11604.5  # in eV
        states = eon.akmc.get_statelist(kT, econf)
        start_state_num = 0
        current_state = states.get_state(start_state_num)
        previous_state = current_state
        explore_state = current_state
        state_explorer = MinModeExplorer(
            states, previous_state, explore_state, superbasin=None, config=econf
        )
        displacement, mode, disp_type = state_explorer.generate_displacement()
        eio.savecon(Path.cwd() / "displacement.con", displacement)
        eio.save_mode(Path.cwd() / "direction.dat", mode)

        # Delete eOn's generated files that we don't need here
        if (Path.cwd() / "jobs").exists():
            shutil.rmtree(Path.cwd() / "jobs")
        if (Path.cwd() / "states").exists():
            shutil.rmtree(Path.cwd() / "states")

        logger.info(f"Initial displacement setup complete in {run_dir}")


@click.command()
@click.argument("mol_idx", type=int)
@click.option(
    "--spin",
    required=True,
    type=click.Choice(["singlet", "doublet"]),
    help="Spin multiplicity.",
)
@click.option("--rdir", required=True, type=str, help="Name of the run directory.")
@click.option(
    "--input-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help="Optional directory containing pre-generated displacement files (pos.con, direction.dat, displacement.con).",
)
def main(mol_idx, spin, rdir, input_dir):
    """
    Prepares all necessary input files for an eOn GPRD benchmark run.
    """
    logger = setup_logger()
    gitroot = get_gitroot()
    mult = 1 if spin == "singlet" else 2
    mol_idx_str = f"{mol_idx:03d}"

    run_dir = Path.cwd() / "snake_runs" / rdir / spin / mol_idx_str
    run_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Preparing files for {spin} {mol_idx_str} in {run_dir}")

    atoms = prepare_con_file(run_dir, f"{spin}s", mol_idx_str, gitroot)
    settings_path = prepare_nwchem_settings(run_dir, spin)
    prepare_config_ini(run_dir, mol_idx, atoms, mult, settings_path)

    # Check for the environment variable and use it if set
    env_input_dir = os.getenv("EON_DISPLACEMENT_DIR")
    if env_input_dir:
        # Construct the full path to the specific run directory
        # NOTE(rg): JCTC used singlets and doublets, so {spin}s
        input_dir = Path(env_input_dir) / f"{spin}s" / mol_idx_str
        logger.info(f"Using input directory from EON_DISPLACEMENT_DIR environment variable: {input_dir}")

    setup_initial_displacement(run_dir, logger, input_dir)

    logger.info("Running rgpycrumbs CLI to generate final NWChem socket input...")
    command = ["python", "-m", "rgpycrumbs.cli", "eon", "generate_nwchem_input"]
    try:
        result = subprocess.run(
            command, cwd=run_dir, check=True, capture_output=True, text=True
        )
        logger.info("rgpycrumbs CLI finished successfully.")
        logger.debug(result.stdout)
    except subprocess.CalledProcessError as e:
        logger.error("The rgpycrumbs CLI tool failed to run.")
        logger.error(f"Return code: {e.returncode}")
        logger.error(f"Stdout:\n{e.stdout}")
        logger.error(f"Stderr:\n{e.stderr}")
        exit(1)

    logger.info("Preparation complete.")


if __name__ == "__main__":
    main()
