#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace C
{
    inline constexpr double NA = 6.022149e23; // Avogadro's number
};

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
    inline constexpr std::string_view C0_KEY = "[C0]";
    inline constexpr std::string_view EFFICIENCY_KEY = "f";
    inline constexpr std::string_view END_GROUP_NAMES_KEY = "end_group_units";
    inline constexpr std::string_view POLYMER_NAMES_KEY = "polymer_names";

    // Reactions
    inline constexpr std::string_view RATE_CONSTANT_KEY = "rate_constant";
    inline constexpr std::string_view REACTANTS_KEY = "reactants";
    inline constexpr std::string_view PRODUCTS_KEY = "products";
};

namespace C::state
{
    // Basic KMC state info
    inline constexpr std::string_view ITERATION_KEY = "Iteration";
    inline constexpr std::string_view KMC_TIME_KEY = "KMC Time";
    inline constexpr std::string_view KMC_STEP_KEY = "KMC Step";
    inline constexpr std::string_view NAV_KEY = "NAV";

    // Simulation time info
    inline constexpr std::string_view SIM_TIME_KEY = "Simulation Time";
    inline constexpr std::string_view SIM_TIME_PER_1E6_STEPS_KEY = "Simulation Time per 1e6 KMC Steps";

    // Species info
    inline constexpr std::string_view CONV_PREFIX = "Conv_";
    inline constexpr std::string_view COUNT_PREFIX = "Count_";
    inline constexpr std::string_view MONOMER = "Monomer";

    // Chain and molecular weight stats
    inline constexpr std::string_view NAVGCL_KEY = "nAvgCL";
    inline constexpr std::string_view WAVGCL_KEY = "wAvgCL";
    inline constexpr std::string_view DISPCL_KEY = "dispCL";
    inline constexpr std::string_view NAVGMW_KEY = "nAvgMW";
    inline constexpr std::string_view WAVGMW_KEY = "wAvgMW";
    inline constexpr std::string_view DISPMW_KEY = "dispMW";

    // Copolymerization stats
    inline constexpr std::string_view NAVGCOMP_PREFIX = "nAvgComp_";
    inline constexpr std::string_view NAVGSL_PREFIX = "nAvgSL_";
    inline constexpr std::string_view WAVGSL_PREFIX = "wAvgSL_";
    inline constexpr std::string_view DISPSL_PREFIX = "dispSL_";

    // Sequence stats
    inline constexpr std::string_view BUCKET_KEY = "Bucket";
    inline constexpr std::string_view MONCOUNT_PREFIX = "MonCount_";
    inline constexpr std::string_view SEQCOUNT_PREFIX = "SeqCount_";
    inline constexpr std::string_view SEQLEN2_PREFIX = "SeqLen2_";
}
