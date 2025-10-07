#pragma once

#include <string>
#include <vector>

namespace C
{
    namespace io
    {
        inline constexpr const char *PARAMETERS_SECTION = "parameters";
        inline constexpr const char *SPECIES_SECTION = "species";
        inline constexpr const char *RATE_CONSTANTS_SECTION = "rateconstants";
        inline constexpr const char *REACTIONS_SECTION = "reactions";
        inline constexpr std::array<const char *, 4> REQUIRED_SECTIONS = {PARAMETERS_SECTION, SPECIES_SECTION, RATE_CONSTANTS_SECTION, REACTIONS_SECTION};

        // General
        inline constexpr const char *NAME_KEY = "name";
        inline constexpr const char *TYPE_KEY = "type";
        inline constexpr const char *VALUE_KEY = "value";

        // Parameters
        inline constexpr const char *NUM_UNITS_KEY = "num_units";
        inline constexpr const char *TERMINATION_TIME_KEY = "termination_time";
        inline constexpr const char *ANALYSIS_TIME_KEY = "analysis_time";

        // Species
        inline constexpr const char *FW_KEY = "FW";
        inline constexpr const char *C0_KEY = "C0";
        inline constexpr const char *EFFICIENCY_KEY = "f";
        inline constexpr const char *END_GROUP_NAMES_KEY = "end_group_units";
        inline constexpr const char *POLYMER_NAMES_KEY = "polymer_names";

        // Reactions
        inline constexpr const char *RATE_CONSTANT_KEY = "rate_constant";
        inline constexpr const char *REACTANTS_KEY = "reactants";
        inline constexpr const char *PRODUCTS_KEY = "products";
    }
};
