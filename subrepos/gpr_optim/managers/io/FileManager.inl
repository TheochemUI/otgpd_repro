/*
 * IOManager.inl
 *
 *  Created on: 17 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#ifndef MANAGERS_IO_FILEMANAGER_INL_
#define MANAGERS_IO_FILEMANAGER_INL_

namespace gpr {
namespace io {

inline std::vector<std::string> FileManager::splitString(const std::string &str,
                                                         std::string _regex)
{
    std::regex re(_regex);
    std::sregex_token_iterator re_it(str.begin(), str.end(), re, -1);
    return std::vector<std::string>(re_it, std::sregex_token_iterator());
}

template <typename T>
void FileManager::convertStringToValue(
    std::map<std::string, std::vector<std::string> > &input_data,
    gpr::KeyValuePair<T> &parameters)
{
    if (input_data.find(parameters.key) != input_data.end()) {
        if (typeid(parameters.value) == typeid(int))
            parameters.value = std::atoi(input_data[parameters.key][0].c_str());
        if (typeid(parameters.value) == typeid(double))
            parameters.value = std::atof(input_data[parameters.key][0].c_str());
    }
}

template <typename T>
void FileManager::convertStringToString(
    std::map<std::string, std::vector<std::string> > &input_data,
    gpr::KeyValuePair<T> &parameters)
{
    if (input_data.find(parameters.key) != input_data.end()) {
        parameters.value = input_data[parameters.key][0];
    }
}

template <typename T>
void FileManager::convertStringToValueArray(
    std::map<std::string, std::vector<std::string> > &input_data,
    gpr::KeyValuePair<T> &parameters)
{
    if (input_data.find(parameters.key) != input_data.end()) {
        for (size_t n = 0; n < input_data[parameters.key].size(); ++n) {
            if (typeid(parameters.value[0]) == typeid(int))
                parameters.value[n] =
                    std::atoi(input_data[parameters.key][n].c_str());
            if (typeid(parameters.value[0]) == typeid(double))
                parameters.value[n] =
                    std::atof(input_data[parameters.key][n].c_str());
        }
    }
}

template <typename T>
void FileManager::convertStringToStringArray(
    std::map<std::string, std::vector<std::string> > &input_data,
    gpr::KeyValuePair<T> &parameters)
{
    if (input_data.find(parameters.key) != input_data.end()) {
        for (size_t n = 0; n < input_data[parameters.key].size(); ++n) {
            parameters.value[n] = input_data[parameters.key][n];
        }
    }
}

template <typename T>
void FileManager::readDataFile(const std::string file_name,
                               std::map<std::string, T *> &data)
{
    std::ifstream is;
    std::string line;
    char split_char = ' ';

    is.open(file_name.c_str(), std::ios::in);

    assert(((void)("Input file cannot be opened"), is.is_open()));

    // Read the file by data sets:
    //  first - field name and dimensions
    //  second - values
    while (std::getline(is, line)) {
        if (!line.empty()) {
            std::vector<std::string> split_line;
            std::string key = "";

            split_line = splitString(line, std::string(1, split_char));

            key = split_line[0];

            if (data.find(key) != data.end()) {
                T *field = data[key];
                gpr::Index_t counter_i = 0;
                gpr::Index_t counter_j = 0;
                gpr::Index_t num_elts_i = std::atoi(split_line[1].c_str());
                gpr::Index_t num_elts_j = std::atoi(split_line[2].c_str());

                field->resize(num_elts_i, num_elts_j);

                // Get values
                std::getline(is, line);

                assert(((void)("Cannot read values from the data file"),
                        !line.empty()));

                split_line.clear();
                split_line = splitString(line, std::string(1, split_char));

                counter_i = 0;
                counter_j = 0;
                for (std::string value: split_line) {
                    field->operator()(counter_i, counter_j) =
                        std::atof(value.c_str());
                    ++counter_j;
                    if (counter_j == num_elts_j) {
                        ++counter_i;
                        counter_j = 0;
                    }
                }
            }
        }
    }

    is.close();

    // Print if needed
    //    for(std::pair<std::string, gpr::Field<T>*> elt : data) {
    //        std::cout << elt.first << "\n";
    //        for(Index_t n = 0; n < elt.second->getDimensions(); ++n)
    //            for(Index_t m = 0; m < elt.second->getSize(n); ++m)
    //                std::cout << elt.second->operator()(n, m) << " ";
    //        std::cout << "\n";
    //    }
}

template <typename T>
bool FileManager::readFromPlainFile(const std::string file_name, T &field)
{
    std::ifstream in_stream;
    in_stream.open(file_name, std::ios::in);

    if (in_stream.is_open()) {
        for (gpr::Index_t n = 0; n < field.size(); ++n) {
            in_stream >> field[n];
        }
        in_stream.close();
        return EXIT_SUCCESS;
    } else {
        io::ErrorManager err;
        err << "Error! File with the reference data not found!";
        return EXIT_FAILURE;
    }
}

template <typename T>
bool FileManager::writePlainFile(const std::string file_name,
                                 const gpr::Field<T> &field,
                                 std::ios::openmode mode)
{
    std::ofstream out_stream;

    out_stream.open(file_name, mode);

    if (out_stream.is_open()) {
        for (gpr::Index_t n = 0; n < field.getSize(); ++n) {
            out_stream << field[n] << "\n";
        }
        out_stream.close();
        return EXIT_SUCCESS;
    } else {
        io::ErrorManager err;
        err << "Error! Can't create file '" << file_name << "'!";
        return EXIT_FAILURE;
    }
}

} /* namespace io */
} /* namespace gpr */

#endif /* MANAGERS_IO_FILEMANAGER_INL_ */
