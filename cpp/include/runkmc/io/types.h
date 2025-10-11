#pragma once
#include <string>
#include "core/types.h"

namespace io::types
{

    struct CommandLineConfig
    {
        std::string inputFilepath;
        std::string outputDir;
        bool reportPolymers = false;
        bool reportSequences = false;
        bool parseOnly = false;
        bool debug = false;
    };

    struct SimulationConfig
    {
        uint64_t numParticles;
        double terminationTime;
        double analysisTime;
    };
};

namespace io::types
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
        SimulationConfig config;
        SpeciesSetRead species;
        std::vector<RateConstantRead> rateConstants;
        std::vector<ReactionRead> reactions;
    };
};
