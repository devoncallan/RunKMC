#pragma once
#include <string>
#include <vector>

static const size_t NUM_BUCKETS = 30;
typedef uint8_t SpeciesID;

class ReactionType
{
public:
    static inline const std::string ELEMENTARY = "EL";
    static inline const std::string INITIATOR_DECOMPOSITION = "ID";
    static inline const std::string INIT_DECOMP_POLY = "IDP";
    static inline const std::string INITIATION = "IN";
    static inline const std::string PROPAGATION = "PR";
    static inline const std::string DEPROPAGATION = "DP";
    static inline const std::string TERMINATION_C = "TC";
    static inline const std::string TERMINATION_D = "TD";
    static inline const std::string CHAINTRANSFER_M = "CTM";
    static inline const std::string CHAINTRANSFER_S = "CTS";
    static inline const std::string THERM_INIT_M = "TIM";

private:
    ReactionType() = delete;
    ~ReactionType() = delete;
};

class SpeciesType
{
public:
    static inline const std::string UNIT = "U";
    static inline const std::string MONOMER = "M";
    static inline const std::string INITIATOR = "I";
    static inline const std::string POLYMER = "P";
    static inline const std::string UNDEFINED = "?";
    static inline const std::string LABEL = "LABEL";

    static bool isUnitType(const std::string &type)
    {
        return type == UNIT || type == MONOMER || type == INITIATOR;
    }

    static bool isPolymerType(const std::string &type)
    {
        return type == POLYMER || type == LABEL;
    }

    static bool isValidType(const std::string &type)
    {
        for (const auto &validType : validTypes)
        {
            if (type == validType)
                return true;
        }
        return false;
    }

    static void checkValid(const std::string &type)
    {
        for (const auto &validType : validTypes)
        {
            if (type == validType)
                return;
        }
        throw std::invalid_argument("Invalid species type: " + type);
    }

private:
    SpeciesType() = delete;
    ~SpeciesType() = delete;
    static inline const std::vector<std::string> validTypes = {UNIT, MONOMER, INITIATOR, POLYMER, UNDEFINED, LABEL};
};

enum PolymerState
{
    UNINITIATED,
    ALIVE,
    TERMINATED_D,
    TERMINATED_C,
    TERMINATED_CT,
};

namespace config
{

    struct CommandLineConfig
    {
        std::string inputFilepath;
        std::string outputDir;
        bool reportPolymers = false;
        bool reportSequences = false;
    };

    struct SimulationConfig
    {
        uint64_t numParticles;
        double terminationTime;
        double analysisTime;
    };
}

namespace types
{
    template <typename T>
    struct Variable
    {
        std::string name;
        T value;
    };

    // SPECIES TYPES

    struct SpeciesRead
    {
        std::string name;
        std::string type;
    };

    struct RegisteredSpecies : SpeciesRead
    {
        SpeciesID ID;
    };

    struct SpeciesWrite : RegisteredSpecies
    {
        uint64_t count;
        uint64_t initialCount;
    };

    //

    struct UnitRead : SpeciesRead
    {
        double C0 = 0.0;
        double FW = 0.0;
        double efficiency = 1.0;
    };

    struct PolymerTypeRead : SpeciesRead
    {
        std::vector<std::string> endGroupUnitNames;
    };

    struct PolymerLabelsRead : SpeciesRead
    {
        std::vector<std::string> polymerNames;
    };

    struct SpeciesSetRead
    {
        std::vector<UnitRead> units;
        std::vector<PolymerTypeRead> polymerTypes;
        std::vector<PolymerLabelsRead> polymerLabels;
    };

    struct RateConstantRead
    {
        std::string name;
        double k;
    };

    struct ReactionRead
    {
        std::string type;
        std::string rateConstantName;
        std::vector<std::string> reactantNames;
        std::vector<std::string> productNames;
    };

    struct KMCInputRead
    {
        config::SimulationConfig config;
        SpeciesSetRead species;
        std::vector<RateConstantRead> rateConstants;
        std::vector<ReactionRead> reactions;
    };
}