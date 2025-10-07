#pragma once
#include <string>
#include <string_view>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iterator>

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

    // Trim a vector of strings
    static void trim(std::vector<std::string> &v)
    {
        for (auto &s : v)
            trim(s);
    }

    // Check if string starts with prefix - accepts string_view for flexibility
    static bool startswith(const std::string &s, std::string_view prefix)
    {
        return s.size() >= prefix.size() &&
               s.compare(0, prefix.size(), prefix) == 0;
    }

    // Check if string ends with suffix - accepts string_view for flexibility
    static bool endswith(const std::string &s, std::string_view suffix)
    {
        return s.size() >= suffix.size() &&
               s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    // Split string by delimiter - accepts string_view for flexibility
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

    static std::vector<std::string> splitByWhitespace(const std::string &input)
    {
        std::istringstream buffer(input);
        std::vector<std::string> ret(
            (std::istream_iterator<std::string>(buffer)),
            std::istream_iterator<std::string>());
        return ret;
    }

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

// /**
//  * @brief Functions to parse strings.
//  *
//  */
// namespace str
// {

//     // ------------ Modify strings ------------

//     // trim from start (in place)
//     static void ltrim(std::string &s)
//     {
//         s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch)
//                                         { return !std::isspace(ch); }));
//     }

//     // trim from end (in place)
//     static void rtrim(std::string &s)
//     {
//         s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch)
//                              { return !std::isspace(ch); })
//                     .base(),
//                 s.end());
//     }

//     // trim from both ends (in place)
//     static void trim(std::string &s)
//     {
//         ltrim(s);
//         rtrim(s);
//     }

//     // trim a vector of strings
//     static void trim(std::vector<std::string> &v)
//     {
//         for (auto &s : v)
//             trim(s);
//     }

//     // ------------ Parse strings ------------

//     static bool startswith(const std::string &s, const std::string &fs)
//     {
//         return s.rfind(fs, 0) == 0;
//     }

//     static bool endswith(const std::string &s, const std::string &es)
//     {
//         if (s.length() >= es.length())
//             return (0 == s.compare(s.length() - es.length(), es.length(), es));
//         else
//             return false;
//     }

//     static std::vector<std::string> splitByDelimeter(const std::string &s, const std::string &delim)
//     {
//         size_t pos_start = 0, pos_end, delim_len = delim.length();
//         std::string token;
//         std::vector<std::string> res;

//         while ((pos_end = s.find(delim, pos_start)) != std::string::npos)
//         {
//             token = s.substr(pos_start, pos_end - pos_start);
//             pos_start = pos_end + delim_len;
//             res.push_back(token);
//         }

//         res.push_back(s.substr(pos_start));
//         return res;
//     }

//     static std::vector<std::string> splitByWhitespace(const std::string &input)
//     {
//         std::istringstream buffer(input);
//         std::vector<std::string> ret((std::istream_iterator<std::string>(buffer)),
//                                      std::istream_iterator<std::string>());
//         return ret;
//     }

//     static std::string join(const std::vector<std::string> &strings, const std::string &delim = ", ")
//     {
//         std::ostringstream os;
//         for (size_t i = 0; i < strings.size(); ++i)
//         {
//             os << strings[i];
//             if (i < strings.size() - 1)
//                 os << delim;
//         }
//         return os.str();
//     }

//     static std::string join(const std::vector<const char *> &strings, const std::string &delim = ", ")
//     {
//         std::ostringstream os;
//         for (size_t i = 0; i < strings.size(); ++i)
//         {
//             os << strings[i];
//             if (i < strings.size() - 1)
//                 os << delim;
//         }
//         return os.str();
//     }
// };