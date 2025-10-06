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

    static bool isUnitType(const std::string &type)
    {
        return type == UNIT || type == MONOMER || type == INITIATOR;
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
    static inline const std::vector<std::string> validTypes = {UNIT, MONOMER, INITIATOR, POLYMER, UNDEFINED};
};

enum PolymerState
{
    UNINITIATED,
    ALIVE,
    TERMINATED_D,
    TERMINATED_C,
    TERMINATED_CT,
};