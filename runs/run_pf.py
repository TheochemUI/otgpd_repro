#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "rich",
#   "click",
#   "polars",
# ]
# ///

"""
A high-performance CLI tool to analyze and classify simulation runs.

This script recursively searches a directory for run folders, classifies each
run in parallel, extracts timing and force call information for successful
runs, and presents a formatted summary. It supports multiple output formats
(table, CSV, Markdown) and uses the Polars library for data handling.
"""

import gzip
import io
import logging
import os
import pathlib
import re
import sys
import threading
from collections import Counter
from dataclasses import dataclass
from enum import StrEnum
from queue import Queue
from concurrent.futures import ThreadPoolExecutor

import click
import polars as pl
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table

# --- Configuration & Constants ---

logging.basicConfig(
    level="INFO",
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, show_path=False)],
)

GENERIC_ERROR_PATTERNS = ["error", "killed", "sigterm", "failed", "abnormally"]
TIMEOUT_PATTERN = "DUE TO TIME LIMIT"
TIME_REGEX = re.compile(r"real\s+([\d.]+)\s+seconds")


# --- Data Structures ---


class RunStatus(StrEnum):
    SUCCESS = "✅ Success"
    FAILED_TIMEOUT = "⏰ Failed (Timeout)"
    FAILED_ERROR = "❌ Failed (Error)"
    RUNNING = "🏃 Running"
    NOT_RUN = "💤 Not-Run"


@dataclass(frozen=True)
class RunResult:
    path: pathlib.Path
    status: RunStatus
    runtime: float | None = None
    total_force_calls: int | None = None


# --- Core Logic (MODIFIED) ---


def classify_run(run_path: pathlib.Path) -> tuple[RunStatus, float | None, int | None]:
    """
    Classifies a run and extracts its runtime and force calls if successful.
    Returns a tuple of (RunStatus, runtime_in_seconds, total_force_calls).
    """
    if (run_path / "results.dat").exists():
        runtime, force_calls = None, None

        # --- Parse runtime from .log.gz file ---
        try:
            latest_log_file = max(
                run_path.glob("*.log.gz"), key=lambda p: p.stat().st_mtime
            )
            with gzip.open(
                latest_log_file, "rt", encoding="utf-8", errors="ignore"
            ) as f:
                lines = f.readlines()
                for line in reversed(lines[-20:]):
                    match = TIME_REGEX.search(line)
                    if match:
                        runtime = float(match.group(1))
                        break
        except (ValueError, FileNotFoundError):
            pass  # Fail gracefully if log not found

        # --- Parse total_force_calls from results.dat ---
        try:
            with open(run_path / "results.dat", "r") as f:
                for line in f:
                    parts = line.split()
                    if len(parts) >= 2 and parts[1] == "total_force_calls":
                        force_calls = int(parts[0])
                        break
        except (ValueError, FileNotFoundError):
            pass  # Fail gracefully if results.dat is malformed

        return RunStatus.SUCCESS, runtime, force_calls

    # --- Logic for failed/running jobs (unchanged) ---
    try:
        latest_err_file = max(
            run_path.glob("slurm-*.err"), key=lambda p: p.stat().st_mtime
        )
        content = latest_err_file.read_text(encoding="utf-8", errors="ignore")
        if TIMEOUT_PATTERN in content:
            return RunStatus.FAILED_TIMEOUT, None, None
        if any(p in content.lower() for p in GENERIC_ERROR_PATTERNS):
            return RunStatus.FAILED_ERROR, None, None
    except ValueError:
        pass
    if next(run_path.glob("*.gz"), None):
        return RunStatus.RUNNING, None, None
    return RunStatus.NOT_RUN, None, None


def find_run_directories_parallel(search_path: pathlib.Path) -> list[pathlib.Path]:
    logging.info(
        f"Recursively searching for run directories in parallel from '{search_path}'..."
    )
    dir_queue = Queue()
    dir_queue.put(search_path)
    found_dirs, lock = [], threading.Lock()

    def worker():
        while True:
            current_dir = dir_queue.get()
            if current_dir is None:
                break
            try:
                for item in current_dir.iterdir():
                    if item.is_dir():
                        if item.name.isdigit():
                            with lock:
                                found_dirs.append(item)
                        else:
                            dir_queue.put(item)
            except PermissionError:
                logging.warning(f"Permission denied: {current_dir}")
            finally:
                dir_queue.task_done()

    num_workers = (os.cpu_count() or 1) * 2
    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        for _ in range(num_workers):
            executor.submit(worker)
        dir_queue.join()
        for _ in range(num_workers):
            dir_queue.put(None)

    logging.info(f"Found {len(found_dirs)} potential run directories.")
    return found_dirs


# --- Output Formatting (MODIFIED) ---


def render_rich_tables(df: pl.DataFrame) -> tuple[Table, Table]:
    """Generates a summary and a detailed table for console output."""
    status_counts = Counter(df["status"])

    summary_table = Table(title="📊 Run Status Summary", show_header=False)
    summary_table.add_column("Status", style="bold")
    summary_table.add_column("Count", justify="right")
    for status, count in sorted(status_counts.items()):
        summary_table.add_row(status, str(count))

    success_runs = df.filter(pl.col("status") == RunStatus.SUCCESS.value)

    # --- Runtime statistics ---
    runtimes = success_runs["runtime"].drop_nulls()
    if not runtimes.is_empty():
        summary_table.add_section()
        summary_table.add_row("[bold]Successful Run Times[/bold]", "")
        avg, std, min_val, max_val = (
            runtimes.mean(),
            runtimes.std(),
            runtimes.min(),
            runtimes.max(),
        )
        unit, factor = ("min", 60) if avg > 120 else ("sec", 1)
        summary_table.add_row(
            "  Avg ± Std Dev", f"{avg / factor:.2f} ± {std / factor:.2f} {unit}"
        )
        summary_table.add_row(
            "  Min / Max",
            f"{min_val / factor:.2f} {unit} / {max_val / factor:.2f} {unit}",
        )

    # --- Force call statistics ---
    force_calls = success_runs["total_force_calls"].drop_nulls()
    if not force_calls.is_empty():
        summary_table.add_section()
        summary_table.add_row("[bold]Force Call Statistics[/bold]", "")
        avg, std, min_val, max_val = (
            force_calls.mean(),
            force_calls.std(),
            force_calls.min(),
            force_calls.max(),
        )
        summary_table.add_row("  Avg ± Std Dev", f"{avg:.1f} ± {std:.1f}")
        summary_table.add_row("  Min / Max", f"{min_val} / {max_val}")

    # --- Detailed results table ---
    results_table = Table(title="📋 Detailed Run Status")
    results_table.add_column("Path", style="cyan", no_wrap=True)
    results_table.add_column("Status", style="bold")
    status_styles = {
        RunStatus.SUCCESS: "green",
        RunStatus.FAILED_ERROR: "red",
        RunStatus.FAILED_TIMEOUT: "yellow",
        RunStatus.RUNNING: "blue",
        RunStatus.NOT_RUN: "dim",
    }
    for row in df.sort("status", "path").iter_rows(named=True):
        status = row["status"]
        style = status_styles.get(RunStatus(status), "white")
        results_table.add_row(row["path"], f"[{style}]{status}[/]")

    return summary_table, results_table


def render_markdown(df: pl.DataFrame) -> str:
    headers = (
        "| Path | Status | Runtime (s) | Total Force Calls |\n|---|---|---|---|\n"
    )
    rows = [
        f"| {row['path']} | {row['status']} | {row['runtime'] or ''} | {row['total_force_calls'] or ''} |"
        for row in df.sort("path").iter_rows(named=True)
    ]
    return headers + "\n".join(rows)


def render_csv(df: pl.DataFrame) -> str:
    buffer = io.StringIO()
    df.sort("path").write_csv(buffer)
    return buffer.getvalue()


# --- CLI Definition (MODIFIED) ---


@click.command()
@click.argument(
    "search_path",
    type=click.Path(exists=True, file_okay=False, path_type=pathlib.Path),
    default=".",
)
@click.option(
    "-f",
    "--format",
    type=click.Choice(["table", "csv", "md"], case_sensitive=False),
    default="table",
    help="Output format.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False, writable=True, path_type=pathlib.Path),
    default=None,
    help="Path to save the output file.",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Show only the summary, not the detailed list of runs.",
)
def main(search_path: pathlib.Path, format: str, output: pathlib.Path, quiet: bool):
    """Analyzes simulation runs in SEARCH_PATH and reports their status."""
    console = Console()
    run_dirs = find_run_directories_parallel(search_path)
    if not run_dirs:
        logging.warning("No run directories found to analyze.")
        sys.exit(0)

    logging.info(f"Classifying {len(run_dirs)} runs in parallel...")
    with ThreadPoolExecutor() as executor:
        results_tuples = executor.map(classify_run, run_dirs)
        results = [
            RunResult(path=d, status=s, runtime=t, total_force_calls=fc)
            for d, (s, t, fc) in zip(run_dirs, results_tuples)
        ]

    df = pl.DataFrame(
        [
            {
                "path": str(r.path),
                "status": r.status.value,
                "runtime": r.runtime,
                "total_force_calls": r.total_force_calls,
            }
            for r in results
        ]
    )

    if format == "table":
        summary_table, details_table = render_rich_tables(df)
        if output:
            with console.capture() as capture:
                console.print(summary_table)
                if not quiet:
                    console.print(details_table)
            output.write_text(capture.get(), encoding="utf-8")
            logging.info(f"Plain text table saved to '{output}'")
        else:
            console.print(summary_table)
            if not quiet:
                console.print(details_table)
        return

    output_content = ""
    if format == "csv":
        output_content = render_csv(df)
    elif format == "md":
        output_content = render_markdown(df)
    if output:
        output.write_text(output_content, encoding="utf-8")
        logging.info(f"Results written to '{output}'")
    else:
        print(output_content)


if __name__ == "__main__":
    main()
