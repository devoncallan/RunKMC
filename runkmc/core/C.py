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
    RUN_INFO_SECTION = "run_info"

    NAME_KEY = "name"
    TYPE_KEY = "type"
    VALUE_KEY = "value"

    NUM_UNITS_KEY = "num_units"
    TERMINATION_TIME_KEY = "termination_time"
    ANALYSIS_TIME_KEY = "analysis_time"

    FW_KEY = "FW"
    C0_KEY = "[C0]"
    EFFICIENCY_KEY = "f"
    END_GROUP_NAMES_KEY = "end_group_units"
    POLYMER_NAMES_KEY = "polymer_names"

    UNITS_KEY = "units"
    POLYMERS_KEY = "polymers"

    RATE_CONSTANT_KEY = "rate_constant"
    REACTANTS_KEY = "reactants"
    PRODUCTS_KEY = "products"


class state:

    ITERATION_KEY = "Iteration"
    KMC_TIME_KEY = "KMC Time"
    KMC_STEP_KEY = "KMC Step"
    NAV_KEY = "NAV"

    SIM_TIME_KEY = "Simulation Time"
    SIM_TIME_PER_1E6_STEPS_KEY = "Simulation Time per 1e6 KMC Steps"

    CONV_PREFIX = "Conv_"
    COUNT_PREFIX = "Count_"
    MONOMER = "Monomer"

    NAVGCL_KEY = "nAvgCL"
    WAVGCL_KEY = "wAvgCL"
    DISPCL_KEY = "dispCL"
    NAVGMW_KEY = "nAvgMW"
    WAVGMW_KEY = "wAvgMW"
    DISPMW_KEY = "dispMW"

    NAVGCOMP_PREFIX = "nAvgComp_"
    NAVGSL_PREFIX = "nAvgSL_"
    WAVGSL_PREFIX = "wAvgSL_"
    DISP_SL_PREFIX = "dispSL_"

    BUCKET_KEY = "Bucket"
    MONCOUNT_PREFIX = "MonCount_"
    SEQCOUNT_PREFIX = "SeqCount_"
    SEQLEN2_PREFIX = "SeqLen2_"
