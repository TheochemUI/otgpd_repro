/*
 * Macros.h
 *
 *  Created on: 15 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef BACKEND_MACROS_H_
#define BACKEND_MACROS_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <cinttypes>

namespace gpr {
typedef uint32_t Index_t;
#define EigenDynamic Eigen::Dynamic
#define EigenMatrixStorage Eigen::StorageOptions::ColMajor
typedef Eigen::Matrix<double, EigenDynamic, EigenDynamic, EigenMatrixStorage>
    EigenMatrix;
}  // namespace gpr

#define EMPTY -1

#endif /* BACKEND_MACROS_H_ */
