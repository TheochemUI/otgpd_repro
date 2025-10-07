//
//  SCGTest.h
//  gpr_dimer
//
//  Created by Maxim Masterov on 25/11/2020.
//

#ifndef SCGTest_h
#define SCGTest_h

#include <gtest/gtest.h>

#include "../../../gpr/ml/SCG.h"

namespace gpr {
namespace tests {

class SCGTest : public ::testing::Test {
public:
    SCGTest();
    virtual ~SCGTest();

    funcmin::SCG scg;

    double threshold;  // Epsilon for comparison operators.
};

} /* namespace tests */
} /* namespace gpr */

#endif /* SCGTest_h */
