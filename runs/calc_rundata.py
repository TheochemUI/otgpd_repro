#!/usr/bin/env python3

# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "ase==3.23.0",
#   "click==8.1.7",
#   "polars==1.0.0",
#   "rich==13.7.1",
#   "chemparseplot>=0.1.0",
# ]
# ///

import logging
import sys
from dataclasses import asdict, dataclass
from enum import StrEnum
from pathlib import Path

import ase.io
import click
import polars as pl
from rich.logging import RichHandler
from rich.progress import track

from chemparseplot.analyze.use_ira import calculate_rmsd
from chemparseplot.basetypes import SpinID
from chemparseplot.parse.eon.saddle_search import SaddleMeasure, parse_eon_saddle

SINGLET_IDS = [f"{x:003d}" for x in range(265)]
DOUBLET_IDS = [f"{x:003d}" for x in range(235)]

GENERIC_ERROR_PATTERNS = ["error", "killed", "sigterm", "failed", "abnormally"]
TIMEOUT_PATTERN = "DUE TO TIME LIMIT"


class Spin(StrEnum):
    SINGLET = "singlet"
    DOUBLET = "doublet"


class RunStatus(StrEnum):
    SUCCESS = "Success"
    FAILED_TIMEOUT = "Timeout"
    FAILED_ERROR = "Error"


@dataclass
class CalculationResult:
    base_data: SaddleMeasure
    method: str
    rmsd_init_saddle: float | None = None


logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, markup=True)],
)
log = logging.getLogger("rich")


def process_single_run(basepath: Path, spindat: SpinID, method: str) -> dict:
    """
    Parses a single EON run and ALWAYS returns a dictionary.
    Contains full data on success or minimal data on failure.
    """
    results_dat_path = basepath / "results.dat"

    # --- Check for Timeout/Generic Error in slurm-.err file ---
    run_status_override = None
    try:
        # Find the latest slurm error file
        latest_err_file = max(
            basepath.glob("slurm-*.err"), key=lambda p: p.stat().st_mtime
        )
        content = latest_err_file.read_text(encoding="utf-8", errors="ignore")
        if TIMEOUT_PATTERN in content:
            run_status_override = RunStatus.FAILED_TIMEOUT
        elif any(p in content.lower() for p in GENERIC_ERROR_PATTERNS):
            run_status_override = RunStatus.FAILED_ERROR
    except (ValueError, FileNotFoundError):
        # Pass if no error file exists or is empty/unreadable
        pass

    if run_status_override:
        log.warning(
            f"Run {basepath.name} marked as {run_status_override.value} via .err file."
        )
        return {
            "mol_id": spindat.mol_id,
            "spin": spindat.spin,
            "method": method,
            "success": False,
            "termination_status": run_status_override.value,
        }

    # --- Failure Case 1: results.dat does not exist ---
    if not results_dat_path.exists():
        log.debug(f"No results.dat in {basepath}")
        return {
            "mol_id": spindat.mol_id,
            "spin": spindat.spin,
            "method": method,
            "success": False,
        }

    try:
        # --- Success Case ---
        base_results = parse_eon_saddle(basepath, spindat)
        if not base_results.success:
            log.info(f"Run failed for {basepath}: {base_results.termination_status}")
            return asdict(base_results) | {"method": method}
        init = ase.io.read(basepath / "pos.con")
        fin = ase.io.read(basepath / "saddle.con")
        rmsd = calculate_rmsd(init, fin, 10)

        # Flatten the data and return the dictionary directly
        return asdict(base_results) | {
            "method": method,
            "rmsd_init_saddle": rmsd,
        }
    except (KeyError, FileNotFoundError, Exception) as e:
        # --- Failure Case 2: Parsing error ---
        log.error(f"Failed to process {basepath}: {e}")
        return {
            "mol_id": spindat.mol_id,
            "spin": spindat.spin,
            "method": method,
            "success": False,
        }


@click.command()
@click.option(
    "--base_path",
    default=".",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Base path containing the 'runs' folder.",
)
@click.option(
    "--run_dir", required=True, help="Run directory, e.g., 'idimer/otgp_run'."
)
@click.option("--output", default="run_data.csv", help="Output CSV file name.")
@click.option("--method", default="OTGP", help="Method name for data collection.")
def main(base_path: Path, run_dir: str, output: str, method: str):
    """Parses EON saddle point search results and saves them to a CSV file."""
    log.info(f"Starting data collection for method: [bold cyan]{method}[/]")
    all_results = []

    # Using a list of tuples for iteration clarity.
    run_configs = [
        (Spin.SINGLET, SINGLET_IDS),
        (Spin.DOUBLET, DOUBLET_IDS),
    ]

    for spin, mol_ids in run_configs:
        description = f"Processing [bold green]{spin.value}[/]..."
        for mol_id in track(mol_ids, description=description):
            results_path = base_path / f"{run_dir}/{spin.value}/{mol_id}"
            if not results_path.exists():
                log.debug(f"Assuming: [bold cyan]{spin.value}s[/]")
                results_path = base_path / f"{run_dir}/{spin.value}s/{mol_id}"
            # for consistency with ref gprd
            spindat = SpinID(mol_id=mol_id, spin=f"{spin.value}s")
            result_dict = process_single_run(results_path, spindat, method)
            if "scf" in result_dict.keys():
                result_dict.pop("scf")
            all_results.append(result_dict)

    if not all_results:
        log.warning("No results found to save. Exiting.")
        sys.exit(0)

    df = pl.DataFrame(all_results)
    df = df.sort(["spin", "mol_id"])
    df.write_csv(output)
    log.info(f"✅ Data successfully saved to [bold magenta]{output}[/]")


if __name__ == "__main__":
    main()
