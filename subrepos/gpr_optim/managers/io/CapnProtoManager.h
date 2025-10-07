#pragma once

#include <capnp/message.h>
#include <capnp/serialize.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <string>

#include "structures/Structures.h"

namespace gpr::io {

/**
 * @brief Loads parameters by reading a memory-mapped Cap'n Proto binary file.
 *
 * This function replaces the manual text parsing of FileManager::readInputFile.
 * It opens the binary file, maps it to memory, and populates the legacy
 * gpr::InputParameters struct.
 *
 * @param binary_filename The path to the Cap'n Proto binary file (e.g.,
 * "capnp_params.bin").
 * @param parameters The legacy C++ parameters struct to be filled.
 */
void loadParametersFromCapnp(const std::string& binary_filename,
                             gpr::InputParameters& parameters);
}  // namespace gpr::io
