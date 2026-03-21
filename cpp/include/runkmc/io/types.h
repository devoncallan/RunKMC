#pragma once
#include <string>
#include "core/types.h"

namespace io::types
{

    // RunKMC program configuration options
    // Set via command line arguments
    // Instruct the program how to run

    struct CommandLineConfig
    {
        std::string inputFilepath;
        std::string outputDir;
        bool reportPolymers = false;
        bool reportChains = false;
        bool reportPositionalStats = false;
        bool reportSegmentHistogram = false;
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
        UnitRead() = default;
        UnitRead(const SpeciesRead &s) : SpeciesRead(s) {}
        double C0 = 0.0;
        double FW = 0.0;
        double efficiency = 1.0;
    };

    struct MonomerRead : UnitRead
    {
        MonomerRead() = default;
        MonomerRead(const UnitRead &u) : UnitRead(u) {}
        double rho_m = 0.0; // g/cm^3
        double rho_p = 0.0; // g/cm^3
    };

    struct PolymerTypeRead : SpeciesRead
    {
        PolymerTypeRead() = default;
        PolymerTypeRead(const SpeciesRead &s) : SpeciesRead(s) {}
        std::vector<std::string> endGroupUnitNames;
        bool report = false;
    };

    struct PolymerLabelsRead : SpeciesRead
    {
        PolymerLabelsRead() = default;
        PolymerLabelsRead(const SpeciesRead &s) : SpeciesRead(s) {}
        std::vector<std::string> polymerNames;
        bool report = false;
    };

    struct SpeciesSetRead
    {
        std::vector<MonomerRead> monomers;
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
