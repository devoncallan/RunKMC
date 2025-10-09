#pragma once
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>

#include "utils/string.h"
#include "utils/console.h"

namespace input
{
    template <typename T>
    static T convertValue(const std::string &s)
    {
        if constexpr (std::is_same<T, int>::value)
            return std::stoi(s);
        else if constexpr (std::is_same<T, uint64_t>::value)
            return static_cast<uint64_t>(std::stoull(s));
        else if constexpr (std::is_same<T, double>::value)
            return std::stod(s);
        else if constexpr (std::is_same<T, std::string>::value)
            return s;
        else
            static_assert("Unsupported type for convertValue");
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
