#pragma once

#include "utils/console.h"
#include "C.h"

static const size_t NUM_BUCKETS = 30;
typedef uint8_t SpeciesID;

class ReactionType
{
public:
    static constexpr std::string_view ELEMENTARY = "EL";
    static constexpr std::string_view INITIATOR_DECOMPOSITION = "ID";
    static constexpr std::string_view INIT_DECOMP_POLY = "IDP";
    static constexpr std::string_view INITIATION = "IN";
    static constexpr std::string_view PROPAGATION = "PR";
    static constexpr std::string_view DEPROPAGATION = "DP";
    static constexpr std::string_view TERMINATION_C = "TC";
    static constexpr std::string_view TERMINATION_D = "TD";
    static constexpr std::string_view CHAINTRANSFER_M = "CTM";
    static constexpr std::string_view CHAINTRANSFER_S = "CTS";
    static constexpr std::string_view THERM_INIT_M = "TIM";

private:
    ReactionType() = delete;
    ~ReactionType() = delete;
};

class SpeciesType
{
public:
    static constexpr std::string_view UNIT = "U";
    static constexpr std::string_view MONOMER = "M";
    static constexpr std::string_view INITIATOR = "I";
    static constexpr std::string_view POLYMER = "P";
    static constexpr std::string_view UNDEFINED = "?";
    static constexpr std::string_view LABEL = "LABEL";

    static bool isUnitType(std::string_view type)
    {
        return type == UNIT || type == MONOMER || type == INITIATOR;
    }

    static bool isPolymerType(std::string_view type)
    {
        return type == POLYMER || type == LABEL;
    }

    static bool isValidType(std::string_view type)
    {
        for (const auto &_validType : _validTypes)
        {
            if (type == _validType)
                return true;
        }
        return false;
    }

    static std::string invalidTypeString(std::string_view type)
    {
        std::string validTypesStr = "";
        for (const auto &_validType : _validTypes)
            validTypesStr += std::string(_validType) + " ";
        return "Invalid species type `" + std::string(type) + "`. Valid types are: " + validTypesStr + ".";
    }

    static void checkValid(std::string_view type)
    {
        if (!isValidType(type))
            console::input_error(invalidTypeString(type));
    }

private:
    SpeciesType() = delete;
    ~SpeciesType() = delete;
    static inline const std::vector<std::string_view> _validTypes = {UNIT, MONOMER, INITIATOR, POLYMER, UNDEFINED, LABEL};
};

enum PolymerState
{
    UNINITIATED,
    ALIVE,
    TERMINATED_D,
    TERMINATED_C,
    TERMINATED_CT,
};