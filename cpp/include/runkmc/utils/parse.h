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

#include "utils/console.h"

namespace str
{
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
    static bool startswith(const std::string &s, std::string_view prefix)
    {
        return s.size() >= prefix.size() &&
               s.compare(0, prefix.size(), prefix) == 0;
    }

    // Check if string ends with suffix
    static bool endswith(const std::string &s, std::string_view suffix)
    {
        return s.size() >= suffix.size() &&
               s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    // Split string by delimiter
    static std::vector<std::string> splitByDelimeter(const std::string &s, std::string_view delim)
    {
        std::vector<std::string> res;
        size_t pos_start = 0, pos_end;
        std::string delimStr(delim); // Convert once for find()

        while ((pos_end = s.find(delimStr, pos_start)) != std::string::npos)
        {
            res.push_back(s.substr(pos_start, pos_end - pos_start));
            pos_start = pos_end + delim.size();
        }

        res.push_back(s.substr(pos_start));
        return res;
    }

    // Split string by whitespace
    static std::vector<std::string> splitByWhitespace(const std::string &input)
    {
        std::istringstream buffer(input);
        std::vector<std::string> ret(
            (std::istream_iterator<std::string>(buffer)),
            std::istream_iterator<std::string>());
        return ret;
    }

    // Join a vector of strings with a delimiter
    static std::string join(const std::vector<std::string> &strings, const std::string &delim = ", ")
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
};

namespace input
{
    /**
     * @brief Convert string to numeric type, handling scientific notation correctly.
     *
     * This function handles all numeric formats including:
     * - Integers: "1000", "1000000000"
     * - Scientific notation: "1e9", "1.5e3", "1e-6"
     * - Decimals: "0.001", "10000.0"
     *
     * For integer types (int, uint64_t), the function:
     * 1. Parses as double to handle scientific notation
     * 2. Validates the value is a valid integer
     * 3. Converts to the target integer type
     *
     * This ensures "1e9" is correctly parsed as 1,000,000,000 rather than 1.
     */
    template <typename T>
    static T convertValue(const std::string &s)
    {
        if constexpr (std::is_same<T, double>::value)
        {
            return std::stod(s);
        }
        else if constexpr (std::is_same<T, uint64_t>::value)
        {
            // Parse as double first to handle scientific notation (e.g., "1e9")
            double d = std::stod(s);

            // Validate it's a non-negative integer value
            if (d < 0)
                console::input_error("Expected non-negative integer, got negative value: " + s);

            if (d != std::floor(d))
                console::input_error("Expected integer value, got decimal: " + s);

            // Check for overflow (uint64_t max is ~1.8e19)
            if (d > static_cast<double>(std::numeric_limits<uint64_t>::max()))
                console::input_error("Value too large for uint64_t: " + s);

            return static_cast<uint64_t>(d);
        }
        else if constexpr (std::is_same<T, int>::value)
        {
            // Parse as double first to handle scientific notation
            double d = std::stod(s);

            // Validate it's an integer value
            if (d != std::floor(d))
                console::input_error("Expected integer value, got decimal: " + s);

            // Validate it fits in int range
            if (d < static_cast<double>(std::numeric_limits<int>::min()) ||
                d > static_cast<double>(std::numeric_limits<int>::max()))
                console::input_error("Integer value out of range: " + s);

            return static_cast<int>(d);
        }
        else if constexpr (std::is_same<T, std::string>::value)
        {
            return s;
        }
        else
        {
            static_assert(sizeof(T) == 0, "Unsupported type for convertValue");
        }
    }

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
