//
// Created by John Lambert-Admin on 4/19/19.
//

#include <vector>
#include <string>
#include <sstream>

#pragma once

/*
 * Split comma separated line into individual elements (as strings).
 */
std::vector<std::string> split_comma_delimited_string(std::string string_to_split)
{
    std::vector<std::string> comma_sep_strings;
    std::stringstream ss(string_to_split);
    int i;
    while (ss.good())
    {
        std::string substr;
        getline( ss, substr, ',' );
        comma_sep_strings.push_back(substr);
    }
    return comma_sep_strings;
}

size_t str2size_t(std::string str)
{
    std::stringstream sstream(str);
    size_t result;
    sstream >> result;
    return result;
}