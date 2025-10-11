"""
Needs to be kept in sync with cpp/include/runkmc/core/C.h

"""

NA = 6.022149e23  # Avogadro's number


class io:

    PARAMETERS_SECTION = "parameters"
    SPECIES_SECTION = "species"
    RATE_CONSTANTS_SECTION = "rate_constants"
    REACTIONS_SECTION = "reactions"
    END_SECTION = "end"
    REQUIRED_SECTIONS = [
        PARAMETERS_SECTION,
        SPECIES_SECTION,
        RATE_CONSTANTS_SECTION,
        REACTIONS_SECTION,
    ]

    RUN_INFO_SECTION = "run_info"

    # General
    NAME_KEY = "name"
    TYPE_KEY = "type"
    VALUE_KEY = "value"
    ID_KEY = "ID"

    # Parameters
    NUM_UNITS_KEY = "num_units"
    TERMINATION_TIME_KEY = "termination_time"
    ANALYSIS_TIME_KEY = "analysis_time"

    # Species registry keys
    SPECIES_KEY = "species"
    UNITS_KEY = "units"
    MONOMERS_KEY = "monomers"
    POLYMERS_KEY = "polymers"

    # Species
    FW_KEY = "FW"
    C0_KEY = "[C0]"
    EFFICIENCY_KEY = "f"
    END_GROUP_NAMES_KEY = "end_group_units"
    POLYMER_NAMES_KEY = "polymer_names"

    # Reactions
    RATE_CONSTANT_KEY = "rate_constant"
    REACTANTS_KEY = "reactants"
    PRODUCTS_KEY = "products"


class state:

    # Basic KMC state info
    ITERATION_KEY = "Iteration"
    KMC_TIME_KEY = "KMC Time"
    KMC_STEP_KEY = "KMC Step"
    NAV_KEY = "NAV"

    # Simulation time info
    SIM_TIME_KEY = "Simulation Time"
    SIM_TIME_PER_1E6_STEPS_KEY = "Simulation Time per 1e6 KMC Steps"

    # Species info
    CONV_PREFIX = "Conv_"
    COUNT_PREFIX = "Count_"
    MONOMER = "Monomer"

    # Chain and molecular weight stats
    NAVGCL_KEY = "nAvgCL"
    WAVGCL_KEY = "wAvgCL"
    DISPCL_KEY = "dispCL"
    NAVGMW_KEY = "nAvgMW"
    WAVGMW_KEY = "wAvgMW"
    DISPMW_KEY = "dispMW"

    # Copolymerization stats
    NAVGCOMP_PREFIX = "nAvgComp_"
    NAVGSL_PREFIX = "nAvgSL_"
    WAVGSL_PREFIX = "wAvgSL_"
    DISP_SL_PREFIX = "dispSL_"

    # Sequence stats
    BUCKET_KEY = "Bucket"
    MONCOUNT_PREFIX = "MonCount_"
    SEQCOUNT_PREFIX = "SeqCount_"
    SEQLEN2_PREFIX = "SeqLen2_"


class paths:

    INPUT_FILE = "input.txt"
    RESULTS_FILE = "results.csv"
    SEQUENCES_FILE = "sequences.csv"
    POLYMERS_FILE = "polymers.dat"
    SPECIES_FILE = "species.yaml"
    PARSED_INPUT_FILE = "parsed_input.yaml"
