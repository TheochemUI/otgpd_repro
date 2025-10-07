#!/bin/bash

# $1 - name of the compiler
function build_gtest {
  cmake -DCMAKE_CXX_COMPILER="$1" -DCMAKE_CXX_FLAGS="-std=c++11" ../
}

# $1 - name of the compiler
function build_gpr {
  echo ""
  echo "--- Building GPR with" $1
  make clean OBJDIR="objdir_$1" EXE_NAME="gpr_$1.out"
  make CXX="$1" OBJDIR="objdir_$1" EXE_NAME="gpr_$1.out" LIBS_PATHS=-L./googletest/googletest/build_$1/lib
}

# $1 - name of the compiler
function run_gpr {
  echo ""
  echo "--- Running unit tests with" $1
  ./gpr_$1.out
}

# Preparation step
echo ""
echo "--- Initializing submodule"
rm -rf googletest
git submodule init
git submodule update
cd googletest/googletest
patch < ../../cmake_gtest.patch
mkdir -p build_g++
mkdir -p build_icpc
cd ../../

# Build gtest with gcc
echo ""
echo "--- Building googletest lib with GCC"
module purge
module load 2019
module load foss/2018b
module load CMake/3.12.1-GCCcore-7.3.0
cd googletest/googletest/build_g++
build_gtest g++
make
cd ../../../

# Build unit tests with gcc
build_gpr g++


################################
################################
################################

# Build gtest with intel
echo ""
echo "--- Building googletest with intel"
module purge
module load 2019
module load intel/2018b
module load CMake/3.12.1-GCCcore-7.3.0
cd googletest/googletest/build_icpc
build_gtest icpc
make
cd ../../../

# Build unit tests with intel
build_gpr icpc


################################
################################
################################

# Invoke unit tests with gcc
run_gpr g++

# Invoke unit tests with intel
run_gpr icpc
