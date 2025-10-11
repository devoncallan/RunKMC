#pragma once
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <string_view>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iterator>
#include <cmath>
#include <limits>

namespace str
{

    static std::vector<std::string_view> to_views(const std::vector<std::string> &strings)
    {
        return std::vector<std::string_view>(strings.begin(), strings.end());
    }

    // Trim string of whitespace from both ends (in place)
    static void trim(std::string &s)
    {
        // Left trim
        auto leftSpaceEnd = std::find_if(
            s.begin(), s.end(), [](unsigned char ch)
            { return !std::isspace(ch); });
        s.erase(s.begin(), leftSpaceEnd);

        // Right trim
        auto rightSpaceBegin = std::find_if(
            s.rbegin(), s.rend(), [](unsigned char ch)
            { return !std::isspace(ch); });
        s.erase(rightSpaceBegin.base(), s.end());
    }

    // Check if string starts with prefix
    static bool startswith(std::string_view s, std::string_view prefix)
    {
        return s.size() >= prefix.size() &&
               s.compare(0, prefix.size(), prefix) == 0;
    }

    // Check if string ends with suffix
    static bool endswith(std::string_view s, std::string_view suffix)
    {
        return s.size() >= suffix.size() &&
               s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    // Split string by delimiter
    static std::vector<std::string> splitByDelimeter(std::string_view s, std::string_view delim)
    {
        std::vector<std::string> res;
        size_t pos_start = 0, pos_end;
        std::string delimStr(delim); // Convert once for find()

        while ((pos_end = s.find(delimStr, pos_start)) != std::string::npos)
        {
            res.push_back(std::string(s.substr(pos_start, pos_end - pos_start)));
            pos_start = pos_end + delim.size();
        }

        res.push_back(std::string(s.substr(pos_start)));
        return res;
    }

    // Split string by whitespace
    static std::vector<std::string> splitByWhitespace(std::string_view s)
    {
        std::istringstream buffer{std::string(s)};
        std::vector<std::string> ret{
            std::istream_iterator<std::string>{buffer},
            std::istream_iterator<std::string>{}};
        return ret;
    }

    // Join a vector of strings with a delimiter
    static std::string join(const std::vector<std::string_view> &strings, std::string_view delim = ", ")
    {
        std::ostringstream os;
        for (size_t i = 0; i < strings.size(); ++i)
        {
            os << strings[i];
            if (i < strings.size() - 1)
                os << delim;
        }
        return os.str();
    }
    static std::string join(const std::vector<std::string> &strings, std::string_view delim = ", ") { return join(to_views(strings), delim); }

    template <typename T>
    size_t findInVector(const std::string &varName, const std::vector<T> &vec)
    {
        for (size_t i = 0; i < vec.size(); ++i)
        {
            if (vec[i].name == varName)
                return i;
        }
        return SIZE_T_MAX;
    }
};