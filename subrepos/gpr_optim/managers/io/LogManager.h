/*
 * Log.h
 *
 *  Created on: 17 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef MANAGERS_IO_LOGMANAGER_H_
#define MANAGERS_IO_LOGMANAGER_H_

#include <iostream>

namespace gpr {
namespace io {

/**
 * @brief Responsible for logging.
 */
class LogManager {
public:
    LogManager() { }
    ~LogManager() { }

    /**
     * @brief Writes data to the log.
     * @param data Data.
     * @return Reference to \e this.
     */
    template <typename T>
    inline LogManager& operator<<(const T& data)
    {
        // TODO: Implement logging into the file and/or stdout here
        std::clog << data;
        return *this;
    }
};

} /* namespace io */
} /* namespace gpr */

#endif /* MANAGERS_IO_LOGMANAGER_H_ */
