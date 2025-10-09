#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

#define BLK "\x1b[0;30m"
#define RED "\x1b[0;31m"
#define GRN "\x1b[0;32m"
#define YLW "\x1b[0;33m"
#define BLU "\x1b[0;34m"
#define MAG "\x1b[0;35m"
#define CYN "\x1b[0;36m"
#define WHT "\x1b[0;37m"
#define NRM "\x1b[0m"
#define UBLU "\x1b[4;34m"

namespace console
{
    static void debug_rxn(std::string message)
    {
        std::cout << MAG << "[ DEBUG ] : ";
        std::cout << message << NRM << std::endl;
    }

    static void debug(std::string message)
    {
        std::cout << BLU << "[ DEBUG ] : ";
        std::cout << message << NRM << std::endl;
    }

    static void debug_msg(std::string message)
    {
        std::cout << UBLU << "[ DEBUG ] : ";
        std::cout << message << NRM << std::endl;
    }

    template <typename T>
    static void debug_value(std::string name, T variable)
    {
        std::cout << BLU << "[ DEBUG ] : ";
        std::cout << name << ": " << variable << NRM << std::endl;
    }

    static void log(std::string message)
    {
        std::cout << GRN << "[ LOG ] : ";
        std::cout << message << NRM << std::endl;
    }
    static void warning(std::string message)
    {
        std::cout << YLW << "[ WARNING ] : ";
        std::cout << message << NRM << std::endl;
    }
    static void input_warning(std::string message)
    {
        std::cout << YLW << "[ INPUT WARNING ] : ";
        std::cout << message << NRM << std::endl;
    }
    static void error(std::string message)
    {
        std::cout << RED << "[ ERROR ] : ";
        std::cout << message << NRM << std::endl;
        exit(EXIT_FAILURE);
    }
    static void input_error(std::string message)
    {
        std::cout << RED << "[ INPUT ERROR ] : ";
        std::cout << message << NRM << std::endl;
        exit(EXIT_FAILURE);
    }
};