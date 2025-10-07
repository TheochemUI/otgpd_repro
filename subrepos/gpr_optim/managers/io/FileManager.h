/*
 * IOManager.h
 *
 *  Created on: 17 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef MANAGERS_IO_FILEMANAGER_H_
#define MANAGERS_IO_FILEMANAGER_H_

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>
#include <vector>

#include "../../data_types/Field.h"
#include "../../structures/Structures.h"

namespace gpr {
namespace io {
/**
 * @brief Responsible for generic IO operations.
 */
class FileManager {
public:
    FileManager();
    virtual ~FileManager();

    /**
     * @brief Read data from the input file and populates structure of \e
     * parameters. The input file may have empty lines and in-line comments that
     * should start with "%" character. The values should be assigned to
     * variable using the assignment operator "=".
     * @param file_name Name of the input file.
     * @param parameters Structure of parameters
     */
    void readInputFile(const std::string file_name,
                       gpr::InputParameters &parameters);

    /**
     * @brief Read data file.
     * The method reads data file written in ASCII format. The values of \e data
     * should point to fields that supposed to be populated. The keys in \e data
     * should have names of fields from the data file.
     * \note Template \e T should be of type \e Coord or \e Field.
     * @param file_name Name of the data file.
     * @param data Map of fields names and pointers to the fields objects.
     */
    template <typename T>
    void readDataFile(const std::string file_name,
                      std::map<std::string, T *> &data);

    /**
     * @brief Read plain ASCII data file.
     * This function reads data line-by-line from a plain text file. No memory
     * allocation and no range chacks are invoked during the execution of the
     * function. \warning \e field should be pre-allocated!
     * @param file_name File name
     * @param field 2D field of data
     * @return EXIT_SUCCESS or EXIT_FAILURE codes
     */
    template <typename T>
    bool readFromPlainFile(const std::string file_name, T &field);

    template <typename T>
    bool writePlainFile(const std::string file_name, const gpr::Field<T> &field,
                        std::ios::openmode mode);

    bool writeXYZFile(const std::string file_name, const gpr::Coord &field,
                      std::ios::openmode mode);

private:
    /**
     * @brief Split the given string using the regular expression.
     * @param str String to split.
     * @param _regex Regular expression.
     * @return Vector of strings.
     */
    inline std::vector<std::string> splitString(const std::string &str,
                                                const std::string _regex);

    /**
     * @brief Convert a string to a scalar value.
     * Method looks for a key \e parameters.key in the \e input_data map and
     * uses this key to access the corresponding value and assign it to the
     * scalar \e parameters.value. The conversion is performed only if the \e
     * parameters.key is registered in the \e input_data map.
     * @param input_data Map of strings.
     * @param parameters Pair of string key and a scalar value.
     */
    template <typename T>
    void convertStringToValue(
        std::map<std::string, std::vector<std::string> > &input_data,
        gpr::KeyValuePair<T> &parameters);

    /**
     * @brief Convert a string to a string value.
     * Method looks for a key \e parameters.key in the \e input_data map and
     * uses this key to access the corresponding value and assign it to the
     * scalar \e parameters.value. The conversion is performed only if the \e
     * parameters.key is registered in the \e input_data map.
     * @param input_data Map of strings.
     * @param parameters Pair of string key and a string value.
     */
    template <typename T>
    void convertStringToString(
        std::map<std::string, std::vector<std::string> > &input_data,
        gpr::KeyValuePair<T> &parameters);

    /**
     * @brief Convert an array of strings to an array of values.
     * Method looks for a key \e parameters.key in the \e input_data map and
     * uses this key to access the corresponding array of scalar values. The
     * found array is converted to \e parameters.value. The conversion is
     * performed only if the \e parameters.key is registered in the \e
     * input_data map.
     * @param input_data Map of strings.
     * @param parameters Pair of string key and array.
     */
    template <typename T>
    void convertStringToValueArray(
        std::map<std::string, std::vector<std::string> > &input_data,
        gpr::KeyValuePair<T> &parameters);

    /**
     * @brief Convert an array of strings to an array of values.
     * Method looks for a key \e parameters.key in the \e input_data map and
     * uses this key to access the corresponding array of scalar values. The
     * found array is converted to \e parameters.value. The conversion is
     * performed only if the \e parameters.key is registered in the \e
     * input_data map.
     * @param input_data Map of strings.
     * @param parameters Pair of string key and array of strings.
     */
    template <typename T>
    void convertStringToStringArray(
        std::map<std::string, std::vector<std::string> > &input_data,
        gpr::KeyValuePair<T> &parameters);
};

} /* namespace io */
} /* namespace gpr */

#include "FileManager.inl"

#endif /* MANAGERS_IO_FILEMANAGER_H_ */
