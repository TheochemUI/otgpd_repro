/*
 * FileManager.cpp
 *
 *  Created on: 17 Jun 2020
 *      Author: Maxim Masterov
 *     Company: SURFSara
 */

#include "FileManager.h"

namespace gpr {
namespace io {

FileManager::FileManager()
{
    // TODO Auto-generated constructor stub
}

FileManager::~FileManager()
{
    // TODO Auto-generated destructor stub
}

void FileManager::readInputFile(const std::string file_name,
                                gpr::InputParameters& parameters)
{
    std::ifstream is;
    std::string line;
    std::map<std::string, std::vector<std::string> > input_data;

    char comment_char = '%';     // Change these two characters to modify the
    char assignment_char = '=';  // comment symbol and the assignment operator

    char open_array_char = '[';       // Change these characters to modify the
    char close_array_char = ']';      // opening and closing brackets for arrays
    char array_separator_char = ',';  // and the separator character

    std::string empy_value =
        "none";  // Used to keep some strings empty in the file

    is.open(file_name.c_str(), std::ios::in);

    assert(((void)("Input file cannot be opened"), is.is_open()));

    // Read the file line-by-line
    while (std::getline(is, line)) {
        // Process the line only if it is not empty
        if (!line.empty()) {
            std::vector<std::string> split_line;
            size_t comment_pos = 0;
            bool array_as_value = false;

            // Remove all tabs and white spaces if there are any
            line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
            line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

            // Skip the following steps if the line is a comment
            if (line[0] == comment_char) continue;

            // Check and remove any characters after (and including) the comment
            // indicator
            comment_pos = line.find_first_of(comment_char);
            if (comment_pos != std::string::npos)
                line.erase(line.begin() + comment_pos, line.end());

            // Split the line using the assignment character
            split_line = splitString(line, std::string(1, assignment_char));

            // Check if the value is a valid array. It should have both opening
            // and closing brackets
            if (split_line[1][0] == open_array_char ||
                split_line[1][split_line[1].size() - 1] == close_array_char) {
                if (split_line[1][0] != open_array_char ||
                    split_line[1][split_line[1].size() - 1] !=
                        close_array_char) {
                    assert(((void)("The input file contains wrong format for "
                                   "the array"),
                            false));
                } else {
                    split_line[1].erase(
                        std::remove(split_line[1].begin(), split_line[1].end(),
                                    open_array_char),
                        split_line[1].end());
                    split_line[1].erase(
                        std::remove(split_line[1].begin(), split_line[1].end(),
                                    close_array_char),
                        split_line[1].end());
                    array_as_value = true;
                }
            }

            // The resulting string should contain exactly two elements,
            // otherwise something went wrong
            assert(((void)("Something went wrong when parsing the input file"),
                    split_line.size() == 2));

            // Add the result to the map
            if (!array_as_value) {
                input_data[split_line[0]].resize(1);
                if (split_line[1] == empy_value) split_line[1] = "";
                input_data[split_line[0]][0] = split_line[1];
            } else {
                std::vector<std::string> split_array;
                split_array = splitString(split_line[1],
                                          std::string(1, array_separator_char));

                input_data[split_line[0]].resize(split_array.size());

                for (size_t n = 0; n < split_array.size(); ++n) {
                    if (split_array[n] == empy_value) split_array[n] = "";
                    input_data[split_line[0]][n] = split_array[n];
                }
            }
        }
    }

    // Print the data if needed
    //    for(std::map<std::string, std::vector<std::string> >::iterator
    //            it_beg = input_data.begin(),
    //            it_end = input_data.end();
    //            it_beg != it_end; ++it_beg) {
    //        std::cout << it_beg->first << " : [";
    //        for(size_t n = 0, end = it_beg->second.size(); n < end; ++n) {
    //            std::cout << it_beg->second[n];
    //            if (n < end - 1)
    //                std::cout << " ";
    //        }
    //        std::cout << "]\n";
    //    }

    convertStringToValue(input_data, parameters.i_dist);
    convertStringToValue(input_data, parameters.i_run);
    convertStringToValue(input_data, parameters.actdist_fro);
    convertStringToValue(input_data, parameters.eval_image1);
    convertStringToValue(input_data, parameters.initrot_nogp);
    convertStringToValue(input_data, parameters.num_iter_initrot);
    convertStringToValue(input_data, parameters.num_iter_rot_gp);
    convertStringToValue(input_data, parameters.divisor_T_dimer_gp);
    convertStringToValue(input_data, parameters.T_dimer);
    convertStringToValue(input_data, parameters.T_anglerot_init);
    convertStringToValue(input_data, parameters.T_anglerot_gp);
    convertStringToValue(input_data, parameters.disp_max);
    convertStringToValue(input_data, parameters.ratio_at_limit);
    convertStringToValue(input_data, parameters.num_bigiter);
    convertStringToValue(input_data, parameters.num_iter);
    convertStringToValue(input_data, parameters.islarge_num_iter);
    convertStringToValue(input_data, parameters.dimer_sep);
    convertStringToValue(input_data, parameters.gp_sigma2);
    convertStringToValue(input_data, parameters.jitter_sigma2);
    convertStringToValue(input_data, parameters.sigma2);
    convertStringToValue(input_data, parameters.magnSigma2);
    convertStringToValue(input_data, parameters.constSigma2);
    convertStringToValue(input_data, parameters.prior_mu);
    convertStringToValue(input_data, parameters.prior_nu);
    convertStringToValue(input_data, parameters.prior_s2);
    convertStringToValue(input_data, parameters.report_level);
    convertStringToValue(input_data, parameters.max_iter);
    convertStringToValue(input_data, parameters.tolerance_func);
    convertStringToValue(input_data, parameters.tolerance_sol);
    convertStringToValue(input_data, parameters.lambda_limit);
    convertStringToValue(input_data, parameters.lambda);

    convertStringToValue(input_data, parameters.use_prune);
    convertStringToValue(input_data, parameters.start_prune_at);
    convertStringToValue(input_data, parameters.nprune_vals);
    convertStringToValue(input_data, parameters.prune_threshold);

    convertStringToString(input_data, parameters.method_rot);
    convertStringToString(input_data, parameters.method_trans);
    convertStringToString(input_data, parameters.load_file);
    convertStringToString(input_data, parameters.save_file);

    convertStringToString(input_data, parameters.debug_output_dir);
    convertStringToString(input_data, parameters.debug_output_file_R);
    convertStringToString(input_data, parameters.debug_output_file_E);
    convertStringToString(input_data, parameters.debug_output_file_G);
    convertStringToString(input_data, parameters.debug_output_file_extension);
    convertStringToValue(input_data, parameters.debug_offset_from_mid_point);
    convertStringToValue(input_data, parameters.debug_dy);
    convertStringToValue(input_data, parameters.debug_dz);
    convertStringToValue(input_data, parameters.debug_level);

    convertStringToString(input_data, parameters.optimization_alg);
    convertStringToString(input_data, parameters.check_derivative);

    convertStringToValueArray(input_data, parameters.dist_sp);
    convertStringToValueArray(input_data, parameters.param_trans);
    convertStringToValueArray(input_data, parameters.cell_dimensions);

    is.close();
}

bool FileManager::writeXYZFile(const std::string file_name,
                               const gpr::Coord& field, std::ios::openmode mode)
{
    std::ofstream out_stream;

    out_stream.open(file_name, mode);

    if (out_stream.is_open()) {
        for (gpr::Index_t n = 0; n < field.getNumPoints(); ++n) {
            out_stream << field.at(0, n) << "\n";
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
