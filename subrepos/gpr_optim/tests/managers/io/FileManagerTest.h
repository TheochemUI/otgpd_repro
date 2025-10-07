/*
 * IOManagerTest.h
 *
 *  Created on: 17 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef TESTS_FILEMANAGERTEST_H_
#define TESTS_FILEMANAGERTEST_H_

#include <gtest/gtest.h>

#include "../../../managers/io/FileManager.h"

namespace gpr {
namespace tests {

class FileManagerTest : public ::testing::Test {
protected:
    FileManagerTest();
    virtual ~FileManagerTest();

    io::FileManager io_manager;
};

} /* namespace tests */
} /* namespace gpr */

#endif /* TESTS_FILEMANAGERTEST_H_ */
