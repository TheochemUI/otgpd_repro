#ifndef ADAMTest_h
#define ADAMTest_h

#include <gtest/gtest.h>

#include "../../../gpr/ml/ADAM.h"

namespace gpr {
namespace tests {

class ADAMTest : public ::testing::Test {
public:
    ADAMTest();
    virtual ~ADAMTest();

    funcmin::ADAM adam;

    double threshold;  // Epsilon for comparison operators.
};

} /* namespace tests */
} /* namespace gpr */

#endif /* ADAMTest_h */
