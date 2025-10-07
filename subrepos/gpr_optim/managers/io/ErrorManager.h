/*
 * ErrorManager.h
 *
 *  Created on: 17 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef MANAGERS_IO_ERRORMANAGER_H_
#define MANAGERS_IO_ERRORMANAGER_H_

#include <iostream>
#include <string>

namespace gpr {
namespace io {

/**
 * @brief Responsible for errors handling.
 */
class ErrorManager {
public:
    ErrorManager() { }
    ~ErrorManager() { }

    /**
     * @brief Writes data to the log.
     * @param data Data.
     * @return Reference to \e this.
     */
    template <typename T>
    inline ErrorManager& operator<<(const T& data)
    {
        // TODO: Implement logging into the file and/or stdout here
        std::cerr << data;
        return *this;
    }

    /**
     * @brief Prints error message and terminates the program.
     * @param expression Assert.
     * @param message Message.
     */
    //    inline void assertMsg(const bool expression,
    //            const std::string& message) const {
    //
    //        if (!expression) {
    //            std::cerr << message << "\n";
    //            exit(1);
    //        }
    //    }

    /**
     * @brief Prints error message and terminates the program.
     * @param expression Assert.
     * @param message Message.
     * @param file Filen ame.
     * @param func Function name.
     * @param line Line number.
     */
    inline static void assertMsg(const bool expression,
                                 const std::string& message, const char* file,
                                 const char* func, const int line)
    {
        if (!expression) {
            std::cerr << "[ " << file << "]::" << func << "::" << line << " "
                      << message << "\n";
            exit(1);
        }
    }
};

} /* namespace io */
} /* namespace gpr */

#define assertMsg(expression, msg)                                       \
    io::ErrorManager::assertMsg(expression, msg, __FILE__, __FUNCTION__, \
                                __LINE__)

#endif /* MANAGERS_IO_ERRORMANAGER_H_ */
