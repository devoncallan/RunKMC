#pragma once
#include <string>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iterator>

/**
 * @brief Functions to parse strings.
 *
 */
namespace str
{

    // ------------ Modify strings ------------

    // trim from start (in place)
    static void ltrim(std::string &s)
    {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch)
                                        { return !std::isspace(ch); }));
    }

    // trim from end (in place)
    static void rtrim(std::string &s)
    {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch)
                             { return !std::isspace(ch); })
                    .base(),
                s.end());
    }

    // trim from both ends (in place)
    static void trim(std::string &s)
    {
        ltrim(s);
        rtrim(s);
    }

    // trim a vector of strings
    static void trim(std::vector<std::string> &v)
    {
        for (auto &s : v)
            trim(s);
    }

    // ------------ Parse strings ------------

    static bool startswith(const std::string &s, const std::string &fs)
    {
        return s.rfind(fs, 0) == 0;
    }

    static std::vector<std::string> splitByDelimeter(const std::string &s, const std::string &delim)
    {
        size_t pos_start = 0, pos_end, delim_len = delim.length();
        std::string token;
        std::vector<std::string> res;

        while ((pos_end = s.find(delim, pos_start)) != std::string::npos)
        {
            token = s.substr(pos_start, pos_end - pos_start);
            pos_start = pos_end + delim_len;
            res.push_back(token);
        }

        res.push_back(s.substr(pos_start));
        return res;
    }

    static std::vector<std::string> splitByWhitespace(const std::string &input)
    {
        std::istringstream buffer(input);
        std::vector<std::string> ret((std::istream_iterator<std::string>(buffer)),
                                     std::istream_iterator<std::string>());
        return ret;
    }
};