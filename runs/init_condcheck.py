#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "ase>=3.22",
#   "numpy>=1.26",
#   "rich>=13.0",
# ]
# ///
"""
A script to compare the initial conditions and results of two EON runs.

This script loads the starting (pos.con), displaced (displacement.con),
and saddle (saddle.con) structures, along with the initial direction vector
(direction.dat) and summary file (results.dat) for two separate EON runs
and computes their similarity.

Similarity is measured by:
1.  CSHDA RMSD for atomic structures (permutation-invariant).
2.  Cosine Similarity for the 3N-dimensional direction vectors.
3.  Key metrics from results.dat.

This script follows modern Python scripting guidelines, including inline
dependency management (PEP 723) and structured logging with rich output.

Usage:
    # This can use a PEP 723-compliant runner like uv or pipx
    uv run compare_eon_initials.py /path/to/run1 /path/to/run2
"""

import argparse
import logging
import sys
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

import ase.io
import numpy as np
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table

try:
    import ira_mod
except ImportError:
    logging.basicConfig(level="CRITICAL", handlers=[RichHandler()])
    logging.critical(
        "Error: `ira_mod.py` not found. Please install the IRA python bindings.",
    )
    sys.exit(1)


# --- Data Structures and Enums ---


class EONSaddleStatus(Enum):
    """Enumeration for EON saddle search termination reasons."""

    # See SaddleSearchJob.cpp for the ordering
    # (numerical_status, descriptive_string)
    GOOD = 0, "Success"
    INIT = 1, "Initial"  # Should never show up, before run
    BAD_NO_CONVEX = 2, "Initial displacement unable to reach convex region"
    BAD_HIGH_ENERGY = 3, "Barrier too high"
    BAD_MAX_CONCAVE_ITERATIONS = 4, "Too many iterations in concave region"
    BAD_MAX_ITERATIONS = 5, "Too many iterations"
    BAD_NOT_CONNECTED = 6, "Saddle is not connected to initial state"
    BAD_PREFACTOR = 7, "Prefactors not within window"
    BAD_HIGH_BARRIER = 8, "Energy barrier not within window"
    BAD_MINIMA = 9, "Minimizations from saddle did not converge"
    FAILED_PREFACTOR = 10, "Hessian calculation failed"
    POTENTIAL_FAILED = 11, "Potential failed"
    NONNEGATIVE_ABORT = 12, "Nonnegative initial mode, aborting"
    NONLOCAL_ABORT = 13, "Nonlocal abort"
    NEGATIVE_BARRIER = 14, "Negative barrier detected"
    BAD_MD_TRAJECTORY_TOO_SHORT = 15, "No reaction found during MD trajectory"
    BAD_NO_NEGATIVE_MODE_AT_SADDLE = (
        16,
        "Converged to stationary point with zero negative modes",
    )
    BAD_NO_BARRIER = 17, "No forward barrier was found along minimized band"
    ZEROMODE_ABORT = 18, "Zero mode abort."
    OPTIMIZER_ERROR = 19, "Optimizer error."
    UNKNOWN = -1, "Unknown status"

    def __new__(cls, value, description):
        obj = object.__new__(cls)
        obj._value_ = value
        obj.description = description
        return obj

    @classmethod
    def from_value(cls, value):
        try:
            value = int(value)
            for member in cls:
                if member.value == value:
                    return member
        except (ValueError, TypeError):
            pass
        return cls.UNKNOWN


@dataclass
class EONRunData:
    """A container for EON run files."""

    directory: Path
    start: ase.Atoms
    direction: np.ndarray
    displacement: ase.Atoms
    saddle: Optional[ase.Atoms] = None
    results: Optional[dict] = None


# --- Core Logic Functions ---


def parse_results_dat(directory: Path) -> Optional[dict]:
    """Parses the results.dat file into a dictionary."""
    results_file = directory / "results.dat"
    if not results_file.exists():
        logging.warning(
            f"'results.dat' not found in {directory}. Results comparison will be skipped."
        )
        return None

    logging.info(f"Found and loading 'results.dat' from {directory}")
    data = {}
    with open(results_file, "r") as f:
        for line in f:
            parts = line.strip().split(maxsplit=1)
            if len(parts) == 2:
                value_str, key = parts
                try:
                    # Attempt to convert value to int, then float, otherwise keep as string
                    value = int(value_str)
                except ValueError:
                    try:
                        value = float(value_str)
                    except ValueError:
                        value = value_str
                data[key] = value
    return data


def load_eon_run_data(directory: Path) -> EONRunData:
    """Loads all EON run files from a specified directory."""
    logging.info(f"Loading data from: {directory}")
    pos_file = directory / "pos.con"
    direction_file = directory / "direction.dat"
    displacement_file = directory / "displacement.con"

    if not all([f.exists() for f in [pos_file, direction_file, displacement_file]]):
        raise FileNotFoundError(
            f"Required files (pos.con, direction.dat, displacement.con) not found in {directory}"
        )

    start_atoms = ase.io.read(pos_file)
    direction_vec = np.loadtxt(direction_file)
    displacement_atoms = ase.io.read(displacement_file)

    num_atoms = len(start_atoms)
    assert direction_vec.shape[0] == num_atoms, (
        f"Mismatch in {directory}: pos.con has {num_atoms} atoms, "
        f"but direction.dat has {direction_vec.shape[0]} entries."
    )

    saddle_file = directory / "saddle.con"
    saddle_atoms = None
    if saddle_file.exists():
        logging.info(f"Found and loading 'saddle.con' from {directory}")
        saddle_atoms = ase.io.read(saddle_file)
    else:
        logging.warning(
            f"'saddle.con' not found in {directory}. Saddle comparison will be skipped."
        )

    results_data = parse_results_dat(directory)

    logging.debug(f"Successfully loaded data for {num_atoms} atoms from {directory}.")
    return EONRunData(
        directory=directory,
        start=start_atoms,
        direction=direction_vec,
        displacement=displacement_atoms,
        saddle=saddle_atoms,
        results=results_data,
    )


def calculate_cshda_rmsd(atoms1: ase.Atoms, atoms2: ase.Atoms) -> float:
    """Calculates the RMSD between two structures using the CShDA algorithm."""
    if len(atoms1) != len(atoms2) or not np.array_equal(
        atoms1.get_atomic_numbers(), atoms2.get_atomic_numbers()
    ):
        logging.warning(
            "Structures have different atom counts or symbols. CSHDA may be meaningless."
        )
        return 999.9

    ira = ira_mod.IRA()
    _, dist = ira.cshda(
        len(atoms1),
        atoms1.get_atomic_numbers(),
        atoms1.get_positions(),
        len(atoms2),
        atoms2.get_atomic_numbers(),
        atoms2.get_positions(),
    )
    return np.sqrt(np.mean(dist**2))


def calculate_cosine_similarity(vec1: np.ndarray, vec2: np.ndarray) -> float:
    """Calculates the cosine similarity between two N-dimensional vectors."""
    v1_flat, v2_flat = vec1.flatten(), vec2.flatten()
    dot_product = np.dot(v1_flat, v2_flat)
    norm_v1, norm_v2 = np.linalg.norm(v1_flat), np.linalg.norm(v2_flat)
    return dot_product / (norm_v1 * norm_v2) if norm_v1 > 0 and norm_v2 > 0 else 0.0


# --- Main Application ---


def main():
    """The main function to parse arguments, run comparisons, and print reports."""
    parser = argparse.ArgumentParser(
        description="Compare the initial conditions and results of two EON runs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("dir1", type=Path, help="Path to the first EON run directory.")
    parser.add_argument("dir2", type=Path, help="Path to the second EON run directory.")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose debug logging."
    )
    args = parser.parse_args()

    log_level = "DEBUG" if args.verbose else "INFO"
    logging.basicConfig(
        level=log_level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
    )

    console = Console()
    console.rule("[bold green]EON Run Comparison[/]")

    try:
        run1 = load_eon_run_data(args.dir1)
        run2 = load_eon_run_data(args.dir2)
    except (FileNotFoundError, AssertionError) as e:
        logging.critical(f"Failed to load data: {e}")
        sys.exit(1)

    # --- Structural and Vector Comparison ---
    logging.info("Calculating structural and vector similarity metrics...")
    rmsd_start = calculate_cshda_rmsd(run1.start, run2.start)
    rmsd_displacement = calculate_cshda_rmsd(run1.displacement, run2.displacement)
    similarity_direction = calculate_cosine_similarity(run1.direction, run2.direction)
    rmsd_saddle = (
        calculate_cshda_rmsd(run1.saddle, run2.saddle)
        if run1.saddle and run2.saddle
        else None
    )

    struct_table = Table(title="Structural & Vector Comparison")
    struct_table.add_column("Metric", justify="right", style="cyan", no_wrap=True)
    struct_table.add_column("Value", style="magenta")
    struct_table.add_row("Run 1 Directory", str(run1.directory.resolve()))
    struct_table.add_row("Run 2 Directory", str(run2.directory.resolve()))
    struct_table.add_section()
    struct_table.add_row("Start Structure CSHDA RMSD (Å)", f"{rmsd_start:.6f}")
    struct_table.add_row(
        "Displaced Structure CSHDA RMSD (Å)", f"{rmsd_displacement:.6f}"
    )
    if rmsd_saddle is not None:
        struct_table.add_row("Saddle Structure CSHDA RMSD (Å)", f"{rmsd_saddle:.6f}")
    struct_table.add_row("Direction Cosine Similarity", f"{similarity_direction:.6f}")
    console.print(struct_table)

    # --- results.dat Comparison ---
    if run1.results and run2.results:
        logging.info("Comparing results.dat metrics...")
        results_table = Table(title="results.dat Comparison")
        results_table.add_column("Metric", justify="right", style="cyan", no_wrap=True)
        results_table.add_column("Run 1", justify="right", style="magenta")
        results_table.add_column("Run 2", justify="right", style="yellow")
        results_table.add_column("Difference", justify="right", style="green")

        # Compare common keys
        common_keys = sorted(list(set(run1.results.keys()) & set(run2.results.keys())))

        for key in common_keys:
            val1, val2 = run1.results.get(key), run2.results.get(key)
            diff_str = "-"

            if key == "termination_reason":
                val1_str = EONSaddleStatus.from_value(val1).description
                val2_str = EONSaddleStatus.from_value(val2).description
            elif isinstance(val1, (int, float)) and isinstance(val2, (int, float)):
                diff = val2 - val1
                diff_str = f"{diff:+.4g}" if isinstance(diff, float) else f"{diff:+.0f}"
                val1_str = f"{val1:.4f}" if isinstance(val1, float) else str(val1)
                val2_str = f"{val2:.4f}" if isinstance(val2, float) else str(val2)
            else:
                val1_str, val2_str = str(val1), str(val2)

            results_table.add_row(key, val1_str, val2_str, diff_str)
        console.print(results_table)

    console.print("\n[bold]Notes:[/]")
    console.print(
        "- CSHDA RMSD closer to 0 means more similar structures (permutation-invariant)."
    )
    console.print(
        "- Cosine Similarity closer to 1 means more similar directions (0=orthogonal, -1=opposite)."
    )


if __name__ == "__main__":
    main()
