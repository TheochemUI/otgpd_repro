/*
 * LBFGS.h
 *
 *  Created on: 31 Mar 2021
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef GPR_LBFGS_H
#define GPR_LBFGS_H

#include <Eigen/Core>

#include "../../backend/Macros.h"
#include "../../data_types/Coord.h"
#include "../../structures/Structures.h"

namespace gpr {
namespace optim {

/**
 * @brief Limited-memory BFGS algorithm.
 */
class LBFGS {
public:
    LBFGS();
    ~LBFGS();

    /**
     * @brief Executes LBFGS algorithm to find rotational direction.
     * @param F Force.
     * @param info Structure with the information from the previous external
     *             iterations.
     * @param z Final result on search direction.
     */
    void execute(const Coord& F, const LBFGSInfo& info, Eigen::VectorXd& z);

private:
    double default_scaling;
};

} /* namespace optim */
} /* namespace gpr */

#endif /* GPR_LBFGS_H */
