//
//  EAMPotentialTest.hpp
//  gpr_dimer
//
//  Created by Maxim Masterov on 02/02/2021.
//

#ifndef EAMPotentialTest_hpp
#define EAMPotentialTest_hpp

#include <gtest/gtest.h>

#include "../../../gpr/potentials/EAMPotential.h"
#include "../../../managers/io/FileManager.h"

namespace gpr {
namespace tests {

class EAMPotentialTest : public ::testing::Test {
protected:
    EAMPotentialTest();
    ~EAMPotentialTest();

    pot::EAMPotential eam_potential;
    io::FileManager io_manager;
    double threshold;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* EAMPotentialTest_hpp */
