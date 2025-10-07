#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "tomli==2.0.1",
#   "pycapnp==2.1.0",
#   "rich==13.7.1",
# ]
# ///

"""
Generates a Cap'n Proto binary buffer from a TOML configuration file.

This script parses a TOML input file containing simulation parameters,
populates a predefined Cap'n Proto schema with these parameters, and
serializes the message to a binary output file.
"""

# Standard library imports
import argparse
import logging
import sys
from typing import Any, Dict

import tomli
from rich.logging import RichHandler
from rich.console import Console


try:
    import capnp
    import gprd_params_capnp as schema
except ImportError:
    print(
        "Error: Could not import the compiled schema 'gprd_params_capnp.py'.",
        file=sys.stderr,
    )
    print(
        "Please generate it first by running: capnp compile -o python gpr_dimer.capnp",
        file=sys.stderr,
    )
    sys.exit(1)


# --- Global Configuration ---
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[
        RichHandler(
            console=Console(stderr=True),
            rich_tracebacks=True,
            show_path=False,
            markup=True,
        )
    ],
)
log = logging.getLogger("rich")


# --- Helper Functions for Populating Schema Sections ---


def _populate_problem_group(cfg: Dict[str, Any], problem):
    """Populates the ProblemGroup section of the schema."""
    problem.actdistFro = cfg.get("actdistFro", 0.0)
    problem.dimerSep = cfg.get("dimerSep", 0.0)
    problem.methodRot = cfg.get("methodRot", "lbfgs")
    problem.methodTrans = cfg.get("methodTrans", "lbfgs")
    problem.divisorTdimerGP = cfg.get("divisorTdimerGP", 0)
    problem.midpointMaxDisp = cfg.get("midpointMaxDisp", 0.0)
    problem.ratioAtLimit = cfg.get("ratioAtLimit", 0.0)
    problem.finalConvergenceForce = cfg.get("finalConvergenceForce", 0.01)

    if "paramTrans" in cfg:
        problem.paramTrans.stepLength = cfg["paramTrans"].get("stepLength", 0.0)
        problem.paramTrans.maxStepLength = cfg["paramTrans"].get("maxStepLength", 0.0)
        problem.paramTrans.rotremThresh = cfg["paramTrans"].get("rotremThresh", 1.0e-6)

    if "initialRotations" in cfg:
        rot_cfg = cfg["initialRotations"]
        problem.initialRotations.maxIterations = rot_cfg.get("maxIterations", 0)
        problem.initialRotations.useGpr = rot_cfg.get("useGpr", True)
        problem.initialRotations.convergenceAngle = rot_cfg.get("convergenceAngle", 0.0)

    if "relaxation" in cfg:
        relax_cfg = cfg["relaxation"]
        relaxation = problem.relaxation
        relaxation.maxOuterIterations = relax_cfg.get("maxOuterIterations", 0)
        relaxation.maxInnerIterations = relax_cfg.get("maxInnerIterations", 0)
        relaxation.maxRotationIterations = relax_cfg.get("maxRotationIterations", 0)
        relaxation.assumeManyIterations = relax_cfg.get("assumeManyIterations", False)
        relaxation.useGprForInitialTranslation = relax_cfg.get(
            "useGprForInitialTranslation", True
        )
        if "convergence" in relax_cfg:
            relaxation.convergence.force = relax_cfg["convergence"].get("force", 0.0)
            relaxation.convergence.angle = relax_cfg["convergence"].get("angle", 0.0)


def _populate_optimizer(cfg: Dict[str, Any], optimizer: schema.HyperparameterOptimizer):
    """Populates the HyperparameterOptimizer section of the schema."""
    optimizer.algorithm = cfg.get("algorithm", "scg")
    optimizer.checkDerivative = cfg.get("checkDerivative", False)
    optimizer.maxIter = cfg.get("maxIter", 0)
    optimizer.tolFunc = cfg.get("tolFunc", 0.0)
    optimizer.tolSol = cfg.get("tolSol", 0.0)

    # Union logic
    if optimizer.algorithm == "scg" and "scg" in cfg:
        scg_opts = optimizer.init("scg")
        scg_opts.lambdaLimit = cfg["scg"].get("lambdaLimit", 1e20)
        setattr(scg_opts, "lambda", cfg["scg"].get("lambda", 0.0))
    elif optimizer.algorithm == "adam" and "adam" in cfg:
        adam_opts = optimizer.init("adam")
        adam_cfg = cfg["adam"]
        adam_opts.lr = adam_cfg.get("lr", 0.001)
        adam_opts.lrd = adam_cfg.get("lrd", 0.999)
        adam_opts.b1 = adam_cfg.get("b1", 0.9)
        adam_opts.b2 = adam_cfg.get("b2", 0.999)
        adam_opts.eps = adam_cfg.get("eps", 1e-8)
        adam_opts.weightDecay = adam_cfg.get("weightDecay", 0.0)
        adam_opts.amsgrad = adam_cfg.get("amsgrad", False)


def populate_message_from_toml(config: Dict[str, Any]) -> schema.InputParameters:
    """
    Populates a Cap'n Proto message from a nested config dictionary.

    Args:
        config: A dictionary loaded from a TOML configuration file.

    Returns:
        A populated InputParameters Cap'n Proto message object.
    """
    params = schema.InputParameters.new_message()

    if "cellDimensions" in config:
        dims = config["cellDimensions"]
        cell_dims_list = params.init("cellDimensions", len(dims))
        for i, dim in enumerate(dims):
            cell_dims_list[i] = dim

    if "runSelector" in config:
        cfg = config["runSelector"]
        params.runSelector.iDist = cfg.get("iDist", 0)
        params.runSelector.iRun = cfg.get("iRun", 0)
        if "sepDists" in cfg:
            dists = cfg["sepDists"]
            sep_dists_list = params.runSelector.init("sepDists", len(dists))
            for i, dist in enumerate(dists):
                sep_dists_list[i] = dist

    if "gpr" in config:
        cfg = config["gpr"]
        params.gpr.gpSigma2 = cfg.get("gpSigma2", 0.0)
        params.gpr.jitterSigma2 = cfg.get("jitterSigma2", 0.0)
        params.gpr.noiseSigma2 = cfg.get("noiseSigma2", 0.0)
        if "prior" in cfg:
            params.gpr.prior.mu = cfg["prior"].get("mu", 0.0)
            params.gpr.prior.nu = cfg["prior"].get("nu", 0.0)
            params.gpr.prior.s2 = cfg["prior"].get("s2", 0.0)

    if "early_stopping" in config:
        cfg = config["early_stopping"]
        params.earlyStopping.distMetric = cfg.get("dist_metric", "emd")
        params.earlyStopping.threshold = cfg.get("threshold", 1.2)

    if "pruning" in config:
        cfg = config["pruning"]
        params.pruning.usePrune = cfg.get("usePrune", False)
        params.pruning.pruneBegin = cfg.get("pruneBegin", 0)
        params.pruning.pruneNVals = cfg.get("pruneNVals", 0)
        params.pruning.pruneThreshold = cfg.get("pruneThreshold", 0.0)

    if "debug" in config:
        cfg = config["debug"]
        debug = params.debug
        debug.reportLevel = cfg.get("reportLevel", 0)
        debug.debugLevel = cfg.get("debugLevel", 0)
        debug.debugOutDir = cfg.get("debugOutDir", "")
        debug.debugPosFile = cfg.get("debugPosFile", "")
        debug.debugEnergyFile = cfg.get("debugEnergyFile", "")
        debug.debugGradFile = cfg.get("debugGradFile", "")
        debug.debugOutExt = cfg.get("debugOutExt", ".dat")
        debug.debugOffsetMidPoint = cfg.get("debugOffsetMidPoint", 0.0)
        debug.debugDy = cfg.get("debugDy", 0.1)
        debug.debugDz = cfg.get("debugDz", 0.1)

    if "problem" in config:
        cfg = config["problem"]
        problem = params.problem
        _populate_problem_group(config["problem"], problem)
    if "optimizer" in config:
        _populate_optimizer(config["optimizer"], params.optimizer)

    return params


def main():
    """Main execution function: parses arguments and orchestrates file conversion."""
    parser = argparse.ArgumentParser(
        description="Generates a Cap'n Proto binary buffer from a TOML config file."
    )
    parser.add_argument("input_toml", type=str, help="Path to the input TOML file.")
    parser.add_argument(
        "output_file", type=str, help="Path to write the binary output buffer."
    )
    args = parser.parse_args()

    log.info(f"Parsing TOML config file: [cyan]{args.input_toml}[/cyan]")
    try:
        with open(args.input_toml, "rb") as f:
            config_dict = tomli.load(f)
    except FileNotFoundError:
        log.critical(f"Input file not found at '{args.input_toml}'")
        sys.exit(1)
    except tomli.TOMLDecodeError as e:
        log.critical(f"Invalid TOML format in '{args.input_toml}': {e}")
        sys.exit(1)

    log.info("Populating Cap'n Proto message from TOML data...")
    try:
        capnp_message = populate_message_from_toml(config_dict)
    except Exception as e:
        log.critical(f"Failed to populate schema from TOML data: {e}", exc_info=True)
        sys.exit(1)

    log.info(f"Writing binary buffer to: [cyan]{args.output_file}[/cyan]")
    try:
        with open(args.output_file, "wb") as f:
            capnp_message.write(f)
    except IOError as e:
        log.critical(f"Failed to write to output file '{args.output_file}': {e}")
        sys.exit(1)

    log.info(f"[bold green]Successfully created '{args.output_file}'.[/bold green]")


if __name__ == "__main__":
    main()
