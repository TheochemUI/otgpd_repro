/*
 * Enums.h
 *
 *  Created on: 22 Sep 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_ENUMS_H_
#define GPR_ENUMS_H_

enum Potentials { POTENTIAL_GP, POTENTIAL_TEST1, POTENTIAL_TEST2 };

enum OptionsForGradCalculation {
    D1 = 0b00000001,
    D2 = 0b00000010,
    D12 = 0b00000100,

    D1_pt = 0b00010000,
    D2_pt = 0b00100000,
    D12_pt = 0b01000000,
};

enum Directions { x = 0, y = 1, z = 2 };

enum OptimizationAlgorithms { SCG_opt, ADAM_opt };

enum DimerAlgorithm { LBFGS_alg };

enum FrozenAndMovingAtoms {
    MOVING_ATOM = 0,
    FROZEN_ATOM = 1,
};

enum DebugLevels {
    DEBUG_L0 = 0,
    DEBUG_L1 = 1,
    DEBUG_L2 = 2,  // Write out after each outer loop
    DEBUG_L3 = 3,  // Write out after each inner relaxation loop
};

enum class DistanceMetricType { MAX_1D_LOG, RMSD, EMD };

#endif /* GPR_ENUMS_H_ */
