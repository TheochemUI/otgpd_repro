#!/usr/bin/env bash

set -e

# Default values
REPO_URL="https://github.com/nwchemgit/nwchem"
GITROOT=$(git rev-parse --show-toplevel)
REPO_DIR="$GITROOT/subrepos/nwchem"

# Check if a commit hash or branch name is provided
if [ -z "$1" ]; then
    echo "Error: Git commit hash or branch name is required."
    echo "Usage: $0 <commit_hash_or_branch_name>"
    exit 1
fi

COMMIT_HASH="$1"
TAR_GZ_URL="${REPO_URL}/archive/${COMMIT_HASH}.tar.gz"

# --- Main Workflow ---

# Clone the repository if it doesn't exist
if [ -d "$REPO_DIR" ]; then
    echo "Directory '$REPO_DIR' already exists. Skipping clone."
else
    echo "Cloning repository from $REPO_URL into $REPO_DIR..."
    mkdir -p "$REPO_DIR"
    cd "$REPO_DIR"
    cd ../
    rm "$REPO_DIR"
    git clone "$REPO_URL" "$REPO_DIR"
fi

echo "Checking out commit/branch '$COMMIT_HASH'..."
cd "$REPO_DIR"
git checkout "$COMMIT_HASH"

# Set up environment variables
echo "Setting environment variables..."
export NWCHEM_TOP="$REPO_DIR"
export NWCHEM_TARGET="LINUX64"
export NWCHEM_MODULES="all"
# MPI and ARMCI network configuration.
export USE_MPI=1
export USE_MPIF=1
export USE_MPIF4=1
export ARMCI_NETWORK="MPI-PR"
# BLAS and ScaLAPACK configuration.
export BUILD_OPENBLAS=1
export BUILD_SCALAPACK=1
export BLAS_SIZE=8
export SCALAPACK_SIZE=8
export USE_SCALAPACK=1

# Configure and build NWChem
echo "--- Starting NWChem Build Process ---"
cd "$NWCHEM_TOP"/src

echo "Executing: make nwchem_config"
make nwchem_config

make -j"$(nproc)"

echo "--- NWChem Build Complete ---"
