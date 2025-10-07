@0x9ffe5ee8b2bd5ed9;

# --- Reusable Enums for Type Safety ---

enum DimerOptimizer {
  lbfgs @0;
}

enum HyperparameterOptimizer {
  scg @0;
  adam @1;
}

enum DistanceMetric{
  max1DLog @0;
  rmsd @1;
  emd @2;
}

# --- Reusable Building Block Structs ---

# A generic block for defining convergence criteria.
struct ConvergenceCriteria {
  force @0 :Float64;
  angle @1 :Float64; # Used for rotations
}

# A generic block for iteration settings.
struct IterationControl {
  maxIterations @0 :UInt64;
  useGpr @1 :Bool = true; # Default to using GPR unless specified otherwise
}

struct TranslationParams {
  stepLength @0 :Float64;
  maxStepLength @1 :Float64;
  rotremThresh @2 :Float64;
}

# --- Phase-Specific Control Structs ---

# Control parameters for the initial rotation phase.
struct InitialRotationControl {
  maxIterations @0 :UInt64;    # Maps to num_iter_initrot
  useGpr @1 :Bool = true;        # Maps to initrot_nogp (inverted: false if initrot_nogp=1)
  convergenceAngle @2 :Float64;  # Maps to T_anglerot_init
}

# Control parameters for the main relaxation loop (outer/inner iterations).
struct RelaxationControl {
  maxOuterIterations @0 :UInt64; # Maps to num_bigiter
  maxInnerIterations @1 :UInt64; # Maps to num_iter
  maxRotationIterations @2 :UInt64; # Maps to num_iter_rot_gp
  convergence @3 :ConvergenceCriteria; # Holds T_anglerot_gp and the force criterion for the GP surface
  assumeManyIterations @4 :Bool; # Maps to islarge_num_iter
  useGprForInitialTranslation @5 :Bool = true; # Maps to inittrans_nogp (inverted)
}


# Parameters for the GPR model itself.
struct GprModel {
  gpSigma2 @0 :Float64;
  jitterSigma2 @1 :Float64;
  noiseSigma2 @2 :Float64; # from p.sigma2

  struct Prior {
    mu @0 :Float64;
    nu @1 :Float64;
    s2 @2 :Float64;
  }
  prior @3 :Prior;
}

# Settings for optimizing the GPR hyperparameters.
struct HyperparameterOptimization {
  algorithm @0 :HyperparameterOptimizer;
  checkDerivative @1 :Bool;
  maxIter @2 :Int32;
  tolFunc @3 :Float64;
  tolSol @4 :Float64;

  struct ScgOptions {
    lambdaLimit @0 :Float64;
    lambda @1 :Float64;
  }
  struct AdamOptions {
    lr @0 :Float64;
    lrd @1 :Float64;
    b1 @2 :Float64;
    b2 @3 :Float64;
    eps @4 :Float64;
    weightDecay @5 :Float64;
    amsgrad @6 :Bool;
  }

  union {
    scg @5 :ScgOptions;
    adam @6 :AdamOptions;
  }
}

# --- Ancillary Modules ---
struct Pruning {
  usePrune @0 :Bool;
  pruneBegin @1 :Int32;
  pruneNVals @2 :Int32;
  pruneThreshold @3 :Float64;
}

struct Debug {
  reportLevel @0 :Int32;
  debugLevel @1 :Int32;
  debugOutDir @2 :Text;
  debugPosFile @3 :Text;
  debugEnergyFile @4 :Text;
  debugGradFile @5 :Text;
  debugOutExt @6 :Text;
  debugOffsetMidPoint @7 :Float64;
  debugDy @8 :Float64;
  debugDz @9 :Float64;
}

# --- Run-Specific Setup (Not part of the core algorithm parameters) ---
struct RunSelector {
  iDist @0 :Int32;
  iRun @1 :Int32;
  sepDists @2 :List(Float64);
}

struct EarlyStopping{
  distMetric @0 :DistanceMetric;
  threshold  @1 :Float64;
}

# --- Top-Level Configuration Struct ---


struct InputParameters {
  # --- Problem Setup and Dimer Control ---
  problem :group {
    actdistFro @0 :Float64;
    dimerSep @1 :Float64;
    methodRot @2 :DimerOptimizer;
    methodTrans @3 :DimerOptimizer;

    paramTrans @4 :TranslationParams;

    # --- Iteration and Convergence Control ---
    # Grouped by calculation phase for clarity.

    # Final convergence criteria on the *true* potential energy surface.
    finalConvergenceForce @5 :Float64; # Maps to T_dimer

    # Controls for the initial setup rotations.
    initialRotations @6 :InitialRotationControl;

    # Controls for the main GPR-driven relaxation loops.
    relaxation @7 :RelaxationControl;

    # --- Other Algorithm Parameters ---
    divisorTdimerGP @8 :Int32;
    midpointMaxDisp @9 :Float64; # Maps to disp_max
    ratioAtLimit @10 :Float64;
  }

  # --- GPR Model & Optimization ---
  gpr @11 :GprModel;
  optimizer @12 :HyperparameterOptimization;

  earlyStopping @13 :EarlyStopping;

  # --- Ancillary Modules ---
  pruning @14 :Pruning;
  debug @15 :Debug;

  # --- Run-Specific Setup ---
  runSelector @16 :RunSelector;
  cellDimensions @17 :List(Float64);
}
