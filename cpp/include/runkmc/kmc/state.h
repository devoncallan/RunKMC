#pragma once
#include "common.h"
struct KMCState
{
    uint64_t iteration = 0;
    uint64_t kmcStep = 0;
    double kmcTime = 0;
    double simulationTime = 0;
    double simulationTimePer1e6Steps = 0;
    double NAV = 0;

    static std::vector<std::string> getTitles()
    {
        return {"Iteration", "KMC Step", "KMC Time", "Simulation Time", "Simulation Time per 1e6 KMC Steps", "NAV"};
    }

    /*
    Iteration, KMC Step, KMC Time, Simulation Time, Simulation Time per 1e6 KMC Steps, NAV
    */
    std::vector<std::string> getDataAsVector() const
    {
        std::vector<std::string> output;
        output.push_back(std::to_string(iteration));
        output.push_back(std::to_string(kmcStep));
        output.push_back(std::to_string(kmcTime));
        output.push_back(std::to_string(simulationTime));
        output.push_back(std::to_string(simulationTimePer1e6Steps));
        output.push_back(std::to_string(NAV));
        return output;
    }
};

struct SpeciesState
{
    std::vector<double> unitConversions;
    std::vector<uint64_t> unitCounts;
    std::vector<uint64_t> polymerCounts;
    double totalConversion;

    static std::vector<std::string> getTitles()
    {
        std::vector<std::string> names;

        auto unitNames = registry::getAllUnitNames();
        auto polymerGroupNames = registry::getNamesOf(SpeciesType::POLYMER);

        // Unit conversions
        for (const auto &name : unitNames)
            names.push_back("Conv_" + name);
        names.push_back("Conv_Total");

        // Unit counts
        for (const auto &name : unitNames)
            names.push_back("Count_" + name);

        // Polymer counts
        for (const auto &name : polymerGroupNames)
            names.push_back("Count_" + name);

        return names;
    }

    /*
    Conv_R, Conv_A, Conv_B, ..., Conv_Total,
    Count_R, Count_A, Count_B, ...,
    Count_Poly1, Count_Poly2, ...
    */
    std::vector<std::string> getDataAsVector() const
    {
        std::vector<std::string> output;

        // Unit conversions
        for (const auto &conv : unitConversions)
            output.push_back(std::to_string(conv));
        output.push_back(std::to_string(totalConversion));

        // Unit counts
        for (const auto &count : unitCounts)
            output.push_back(std::to_string(count));

        // Polymer counts
        for (const auto &count : polymerCounts)
            output.push_back(std::to_string(count));

        return output;
    }
};

struct AnalysisState
{
    double nAvgCL = 0;
    double wAvgCL = 0;
    double dispCL = 0;

    double nAvgMW = 0;
    double wAvgMW = 0;
    double dispMW = 0;

    // Sequence statistics for each monomer type
    std::vector<double> nAvgComp;
    std::vector<double> nAvgSL;
    std::vector<double> wAvgSL;
    std::vector<double> dispSL;

    AnalysisState()
    {
        nAvgComp.resize(registry::NUM_MONOMERS, 0);
        nAvgSL.resize(registry::NUM_MONOMERS, 0);
        wAvgSL.resize(registry::NUM_MONOMERS, 0);
        dispSL.resize(registry::NUM_MONOMERS, 0);
    }

    static std::vector<std::string> getTitles()
    {

        std::vector<std::string> names = {
            "nAvgCL",
            "wAvgCL",
            "dispCL",
            "nAvgMW",
            "wAvgMW",
            "dispMW",
        };

        auto monomerNames = registry::getNamesOf(SpeciesType::MONOMER);
        for (const auto &monomerName : monomerNames)
            names.push_back("nAvgComp_" + monomerName);
        for (const auto &monomerName : monomerNames)
            names.push_back("nAvgSL_" + monomerName);
        for (const auto &monomerName : monomerNames)
            names.push_back("wAvgSL_" + monomerName);
        for (const auto &monomerName : monomerNames)
            names.push_back("dispSL_" + monomerName);

        return names;
    }

    /*
    nAvgCL, wAvgCL, dispCL, nAvgMW, wAvgMW, dispMW,
    nAvgComp_A, nAvgComp_B, ..., nAvgSL_A, nAvgSL_B, ...,
    wAvgSL_A, wAvgSL_B, ..., dispSL_A, dispSL_B, ...
    */
    std::vector<std::string> getDataAsVector() const
    {
        std::vector<std::string> output;

        output.push_back(std::to_string(nAvgCL));
        output.push_back(std::to_string(wAvgCL));
        output.push_back(std::to_string(dispCL));

        output.push_back(std::to_string(nAvgMW));
        output.push_back(std::to_string(wAvgMW));
        output.push_back(std::to_string(dispMW));

        for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
            output.push_back(std::to_string(nAvgComp[i]));
        for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
            output.push_back(std::to_string(nAvgSL[i]));
        for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
            output.push_back(std::to_string(wAvgSL[i]));
        for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
            output.push_back(std::to_string(dispSL[i]));

        return output;
    }
};

struct SequenceState
{
    KMCState kmcState;
    std::vector<analysis::SequenceStats> stats;

    static std::vector<std::string> getTitles()
    {
        std::vector<std::string> names = {"Iteration", "KMC Time", "Bucket"};
        auto monomerNames = registry::getNamesOf(SpeciesType::MONOMER);
        for (const auto &monomerName : monomerNames)
            names.push_back("monCount_" + monomerName);
        for (const auto &monomerName : monomerNames)
            names.push_back("seqCount_" + monomerName);
        for (const auto &monomerName : monomerNames)
            names.push_back("seqLengths2_" + monomerName);

        return names;
    }

    /*
    Iteration, KMC Time, Bucket, monCount_A, monCount_B, ...,
    seqCount_A, seqCount_B, ..., seqLengths2_A, seqLengths2_B, ...
    */
    std::vector<std::string> getDataAsVector(size_t bucket) const
    {
        std::vector<std::string> output;
        output.push_back(std::to_string(kmcState.iteration));
        output.push_back(std::to_string(kmcState.kmcTime));
        output.push_back(std::to_string(bucket));

        for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
            output.push_back(std::to_string(stats[bucket].monCounts[i]));
        for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
            output.push_back(std::to_string(stats[bucket].seqCounts[i]));
        for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
            output.push_back(std::to_string(stats[bucket].seqLengths2[i]));

        return output;
    }
};

struct SystemState
{
    KMCState kmc;
    SpeciesState species;
    AnalysisState analysis;
    SequenceState sequence;
};