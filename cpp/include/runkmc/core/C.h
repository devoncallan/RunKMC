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

    inline constexpr int PRECISION = std::numeric_limits<double>::max_digits10;

    inline constexpr std::string_view PARAMETERS_SECTION = "parameters";
    inline constexpr std::string_view SPECIES_SECTION = "species";
    inline constexpr std::string_view RATE_CONSTANTS_SECTION = "rateconstants";
    inline constexpr std::string_view REACTIONS_SECTION = "reactions";
    inline constexpr std::string_view END_SECTION = "end";
    inline constexpr std::array<std::string_view, 4> REQUIRED_SECTIONS = {PARAMETERS_SECTION, SPECIES_SECTION, RATE_CONSTANTS_SECTION, REACTIONS_SECTION};

    inline constexpr std::string_view RUN_INFO_SECTION = "run_info";

    // General
    inline constexpr std::string_view NAME_KEY = "name";
    inline constexpr std::string_view TYPE_KEY = "type";
    inline constexpr std::string_view VALUE_KEY = "value";
    inline constexpr std::string_view ID_KEY = "ID";

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

    inline constexpr std::string_view UNITS_KEY = "units";
    inline constexpr std::string_view POLYMERS_KEY = "polymers";

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

namespace C::paths
{
    inline constexpr std::string_view INPUT_FILE = "input.txt";
    inline constexpr std::string_view RESULTS_FILE = "results.csv";
    inline constexpr std::string_view SEQUENCES_FILE = "sequences.csv";
    inline constexpr std::string_view POLYMERS_FILE = "polymers.dat";
    inline constexpr std::string_view METADATA_FILE = "metadata.yaml";
    inline constexpr std::string_view SPECIES_FILE = "species.yaml";
    inline constexpr std::string_view PARSED_INPUT_SUFFIX = ".kmc.yaml";
};

namespace C::io::color
{
    inline bool enabled() { return std::getenv("NO_COLOR") == nullptr; }

    inline constexpr std::string_view DEFAULT = "\x1b[0m";
    inline constexpr std::string_view ULINE = "\x1b[4m";

    inline constexpr std::string_view BLK = "\x1b[0;30m";
    inline constexpr std::string_view RED = "\x1b[0;31m";
    inline constexpr std::string_view GRN = "\x1b[0;32m";
    inline constexpr std::string_view YLW = "\x1b[0;33m";
    inline constexpr std::string_view BLU = "\x1b[0;34m";
    inline constexpr std::string_view MAG = "\x1b[0;35m";
    inline constexpr std::string_view CYN = "\x1b[0;36m";
    inline constexpr std::string_view WHT = "\x1b[0;37m";

    inline std::string_view on(std::string_view code) { return enabled() ? code : std::string_view{}; }

    inline std::string underline(std::string_view text)
    {
        if (on(ULINE).empty())
            return std::string(text);
        return std::string(ULINE) + std::string(text) + std::string(DEFAULT);
    }
};