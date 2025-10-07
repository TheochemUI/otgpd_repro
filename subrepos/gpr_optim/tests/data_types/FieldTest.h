/*
 * FieldTest.h
 *
 *  Created on: 26 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_FIELDTEST_H_
#define TESTS_FIELDTEST_H_

#include <gtest/gtest.h>

#include "../../../data_types/Field.h"

namespace gpr {
namespace tests {

class FieldTest : public ::testing::Test {
public:
    FieldTest();
    virtual ~FieldTest();
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_FIELDTEST_H_ */
