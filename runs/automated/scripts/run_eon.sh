#!/usr/bin/env bash

set -euo pipefail

IS_HPC="${IS_HPC:-false}"

CWD=$(pwd)
# RUNDIR, SPIN, and INDEX are exported from the Snakefile
FULL_RUNDIR="$CWD/snake_runs/${RUNDIR}/${SPIN}/${INDEX}/"
BLESS_LABEL="ec_gprd_${SPIN}_${INDEX}"

# Function to run the simulation
run_simulation() {
    echo "--- Starting eonclient and NWChem ---"
    # Start nwchem in the background
    mpirun -np "$NWCHEM_MPI" "$NWCHEM_COMMAND" nwchem_socket.nwi &
    NWCHEM_PID=$!

    # Start eonclient in the background
    bless --label "$BLESS_LABEL" -- pixi r eonclient &
    EON_PID=$!

    echo "eonclient PID: $EON_PID"
    echo "NWChem PID: $NWCHEM_PID"

    # The function now only waits. The caller script will handle exit codes.
    wait $EON_PID
    wait $NWCHEM_PID
}

case "$IS_HPC" in
true | True | TRUE | 1 | yes | Yes | YES)
    echo "Running on HPC"

    # Get the job ID (example using SLURM)
    JOB_ID="${SLURM_JOB_ID:-$$}"

    # Define the full path to the scratch directory
    SCRATCH_BASE="/scratch/users/rog32"
    SCRATCH_DIR="$SCRATCH_BASE/eon_scratch_${RUNDIR}_${SPIN}_${INDEX}_${JOB_ID}"

    # Define the destination for periodic syncing and the interval in seconds (e.g., 600s = 10 minutes)
    SYNC_DEST="$FULL_RUNDIR/sync_results/"
    SYNC_INTERVAL=60

    echo "Using scratch directory: $SCRATCH_DIR"
    echo "Using CWD: $CWD"
    echo "Using FULL_RUNDIR: $FULL_RUNDIR"
    echo "Periodic sync destination: $SYNC_DEST"

    # Create the scratch and sync directories
    mkdir -p "$SCRATCH_DIR"
    mkdir -p "$SYNC_DEST"
    rsync -a "$FULL_RUNDIR/" "$SCRATCH_DIR/runner/"

    # Run eon in the scratch directory using full paths
    pushd "$SCRATCH_DIR/runner"
    ln -sf "$CWD"/{.pixi,pixi.toml,pixi.lock} .

    echo "--- Starting eonclient and NWChem ---"
    . "$HOME/spack/share/spack/setup-env.sh"
    spack load nwchem

    # Start nwchem and eonclient in the background
    mpirun -np "$NWCHEM_MPI" "$NWCHEM_COMMAND" nwchem_socket.nwi &
    NWCHEM_PID=$!
    bless --label "$BLESS_LABEL" -- pixi r eonclient &
    EON_PID=$!

    echo "eonclient PID: $EON_PID"
    echo "NWChem PID: $NWCHEM_PID"

    # Start the periodic sync process in the background
    (
        # Keep syncing as long as the main eonclient process is running
        while kill -0 $EON_PID 2>/dev/null; do
            sleep "$SYNC_INTERVAL"
            echo "--- Syncing results at $(date) ---"
            # The rsync command copies only the specified file patterns.
            # Using --include and --exclude is robust and won't fail if no files match.
            rsync -av --include='*.log' --include='*.gz' --include='*.h5' --include='*.con' --exclude='*' . "$SYNC_DEST"
        done
        echo "--- Sync process finished ---"
    ) &
    SYNC_PID=$!

    # Wait for both main simulation processes to finish and capture exit codes
    wait $EON_PID
    EON_EXIT_CODE=$?
    wait $NWCHEM_PID
    NWCHEM_EXIT_CODE=$?

    # Clean up the background sync process, if it's still sleeping
    kill "$SYNC_PID" 2>/dev/null || true

    echo "--- Simulation Finished ---"
    echo "eonclient exit code: $EON_EXIT_CODE"
    echo "NWChem exit code: $NWCHEM_EXIT_CODE"

    popd

    # Perform one final, complete sync of all results
    echo "--- Performing final sync of all results ---"
    rsync -a "$SCRATCH_DIR/runner/" "$FULL_RUNDIR/"

    # Clean up the scratch directory
    echo "Cleaning up scratch directory: $SCRATCH_DIR"
    rm -rf "$SCRATCH_DIR"

    # Exit with a non-zero code if either process failed
    if [ $EON_EXIT_CODE -ne 0 ] || [ $NWCHEM_EXIT_CODE -ne 0 ]; then
        exit 1
    fi
    ;;
false | False | FALSE | 0 | no | No | NO | "")
    echo "Running locally (IS_HPC=false, 0, no, or empty)"
    export NWCHEM_MPI="${NWCHEM_MPI:-1}"
    export NWCHEM_COMMAND="${NWCHEM_COMMAND:-nwchem}"

    echo "--- Local Execution Debug Info ---"
    echo "NWCHEM_MPI is set to: $NWCHEM_MPI"
    echo "NWCHEM_COMMAND is set to: $NWCHEM_COMMAND"
    echo "--------------------------------"
    cd "$FULL_RUNDIR"
    run_simulation
    ;;
*)
    echo "Not Running (IS_HPC set to an unexpected value: $IS_HPC)"
    exit 1
    ;;
esac
