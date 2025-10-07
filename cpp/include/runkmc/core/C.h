#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace C::io
{
    inline constexpr std::string_view PARAMETERS_SECTION = "parameters";
    inline constexpr std::string_view SPECIES_SECTION = "species";
    inline constexpr std::string_view RATE_CONSTANTS_SECTION = "rateconstants";
    inline constexpr std::string_view REACTIONS_SECTION = "reactions";
    inline constexpr std::string_view END_SECTION = "end";
    inline constexpr std::array<std::string_view, 4> REQUIRED_SECTIONS = {PARAMETERS_SECTION, SPECIES_SECTION, RATE_CONSTANTS_SECTION, REACTIONS_SECTION};

    // General
    inline constexpr std::string_view NAME_KEY = "name";
    inline constexpr std::string_view TYPE_KEY = "type";
    inline constexpr std::string_view VALUE_KEY = "value";

    // Parameters
    inline constexpr std::string_view NUM_UNITS_KEY = "num_units";
    inline constexpr std::string_view TERMINATION_TIME_KEY = "termination_time";
    inline constexpr std::string_view ANALYSIS_TIME_KEY = "analysis_time";

    // Species
    inline constexpr std::string_view FW_KEY = "FW";
    inline constexpr std::string_view C0_KEY = "C0";
    inline constexpr std::string_view EFFICIENCY_KEY = "f";
    inline constexpr std::string_view END_GROUP_NAMES_KEY = "end_group_units";
    inline constexpr std::string_view POLYMER_NAMES_KEY = "polymer_names";

    // Reactions
    inline constexpr std::string_view RATE_CONSTANT_KEY = "rate_constant";
    inline constexpr std::string_view REACTANTS_KEY = "reactants";
    inline constexpr std::string_view PRODUCTS_KEY = "products";
};
