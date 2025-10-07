#pragma once
#include <yaml-cpp/yaml.h>
#include "utils/console.h"

/**
 * @brief Utilities for reading variables from YAML nodes.
 * Mirrors the input::readVariable pattern from utils/parse.h but for YAML.
 *
 * Separates pure parsing (extracting values) from assignment (using values).
 */
namespace yaml
{
    /************************************************************************************************************/
    /************************************* Functions for reading YAML variables *********************************/
    /************************************************************************************************************/

    

    /**
     * @brief Read a string variable from YAML node (optional).
     *
     * @param node YAML node containing the variable
     * @param variableName name of variable to be read
     * @param variable reference to variable to be overwritten (only if found)
     */
    static void readVariable(const YAML::Node &node, const std::string &variableName, std::string &variable)
    {
        if (node[variableName])
        {
            variable = node[variableName].as<std::string>();
        }
    }

    /**
     * @brief Read a double variable from YAML node (optional).
     *
     * @param node YAML node containing the variable
     * @param variableName name of variable to be read
     * @param variable reference to variable to be overwritten (only if found)
     */
    static void readVariable(const YAML::Node &node, const std::string &variableName, double &variable)
    {
        if (node[variableName])
        {
            variable = node[variableName].as<double>();
        }
    }

    /**
     * @brief Read a uint64_t variable from YAML node (optional).
     *
     * @param node YAML node containing the variable
     * @param variableName name of variable to be read
     * @param variable reference to variable to be overwritten (only if found)
     */
    static void readVariable(const YAML::Node &node, const std::string &variableName, uint64_t &variable)
    {
        if (node[variableName])
        {
            variable = node[variableName].as<uint64_t>();
        }
    }

    /**
     * @brief Read a bool variable from YAML node (optional).
     *
     * @param node YAML node containing the variable
     * @param variableName name of variable to be read
     * @param variable reference to variable to be overwritten (only if found)
     */
    static void readVariable(const YAML::Node &node, const std::string &variableName, bool &variable)
    {
        if (node[variableName])
        {
            variable = node[variableName].as<bool>();
        }
    }

    /**
     * @brief Read a required string variable from YAML node.
     * Errors if not found.
     *
     * @param node YAML node containing the variable
     * @param variableName name of variable to be read
     * @param variable reference to variable to be overwritten
     */
    static void readVariableRequired(const YAML::Node &node, const std::string &variableName, std::string &variable)
    {
        if (!node[variableName])
        {
            console::input_error("Required variable '" + variableName + "' not found in YAML.");
        }
        variable = node[variableName].as<std::string>();
    }

    /**
     * @brief Read a required double variable from YAML node.
     * Errors if not found.
     *
     * @param node YAML node containing the variable
     * @param variableName name of variable to be read
     * @param variable reference to variable to be overwritten
     */
    static void readVariableRequired(const YAML::Node &node, const std::string &variableName, double &variable)
    {
        if (!node[variableName])
        {
            console::input_error("Required variable '" + variableName + "' not found in YAML.");
        }
        variable = node[variableName].as<double>();
    }

    /**
     * @brief Read a required uint64_t variable from YAML node.
     * Errors if not found.
     *
     * @param node YAML node containing the variable
     * @param variableName name of variable to be read
     * @param variable reference to variable to be overwritten
     */
    static void readVariableRequired(const YAML::Node &node, const std::string &variableName, uint64_t &variable)
    {
        if (!node[variableName])
        {
            console::input_error("Required variable '" + variableName + "' not found in YAML.");
        }
        variable = node[variableName].as<uint64_t>();
    }

    /**
     * @brief Read a required bool variable from YAML node.
     * Errors if not found.
     *
     * @param node YAML node containing the variable
     * @param variableName name of variable to be read
     * @param variable reference to variable to be overwritten
     */
    static void readVariableRequired(const YAML::Node &node, const std::string &variableName, bool &variable)
    {
        if (!node[variableName])
        {
            console::input_error("Required variable '" + variableName + "' not found in YAML.");
        }
        variable = node[variableName].as<bool>();
    }

    /**
     * @brief Check if a YAML node exists (for validation).
     *
     * @param node YAML node to check
     * @param variableName name of variable to check for
     * @return true if variable exists, false otherwise
     */
    static bool hasVariable(const YAML::Node &node, const std::string &variableName)
    {
        return node[variableName] && !node[variableName].IsNull();
    }
}
