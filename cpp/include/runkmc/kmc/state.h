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
        return {
            std::string(C::state::ITERATION_KEY),
            std::string(C::state::KMC_STEP_KEY),
            std::string(C::state::KMC_TIME_KEY),
            std::string(C::state::SIM_TIME_KEY),
            std::string(C::state::SIM_TIME_PER_1E6_STEPS_KEY),
            std::string(C::state::NAV_KEY)};
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
    double monomerConversion = 0;

    static std::vector<std::string> getTitles()
    {
        std::vector<std::string> names;

        auto unitNames = registry::getAllUnitNames();
        auto polymerGroupNames = registry::getPolymerNames();

        // Unit conversions
        for (const auto &name : unitNames)
            names.push_back(std::string(C::state::CONV_PREFIX) + name);
        names.push_back(std::string(C::state::CONV_PREFIX) + std::string(C::state::MONOMER));

        // Unit counts
        for (const auto &name : unitNames)
            names.push_back(std::string(C::state::COUNT_PREFIX) + name);

        // Polymer counts
        for (const auto &name : polymerGroupNames)
            names.push_back(std::string(C::state::COUNT_PREFIX) + name);

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
        output.push_back(std::to_string(monomerConversion));

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
        auto numMonomers = registry::getNumMonomers();
        if (numMonomers <= 1)
            return;

        nAvgComp.resize(numMonomers, 0);
        nAvgSL.resize(numMonomers, 0);
        wAvgSL.resize(numMonomers, 0);
        dispSL.resize(numMonomers, 0);
    }

    static std::vector<std::string> getTitles()
    {

        // Chain length statistics and molecular weight statistics
        std::vector<std::string> names = {
            std::string(C::state::NAVGCL_KEY),
            std::string(C::state::WAVGCL_KEY),
            std::string(C::state::DISPCL_KEY),
            std::string(C::state::NAVGMW_KEY),
            std::string(C::state::WAVGMW_KEY),
            std::string(C::state::DISPMW_KEY)};

        // If no monomer or homopolymer, don't add copolymer stats
        if (registry::getNumMonomers() <= 1)
            return names;
        auto monomerNames = registry::getMonomerNames();

        // Monomer composition
        for (const auto &monomerName : monomerNames)
            names.push_back(std::string(C::state::NAVGCOMP_PREFIX) + monomerName);

        // Number-averaged sequence lengths
        for (const auto &monomerName : monomerNames)
            names.push_back(std::string(C::state::NAVGSL_PREFIX) + monomerName);

        // Weight-averaged sequence lengths
        for (const auto &monomerName : monomerNames)
            names.push_back(std::string(C::state::WAVGSL_PREFIX) + monomerName);

        // Dispersity of sequence lengths
        for (const auto &monomerName : monomerNames)
            names.push_back(std::string(C::state::DISPSL_PREFIX) + monomerName);

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

        // Chain length and molecular weight statistics
        output.push_back(std::to_string(nAvgCL));
        output.push_back(std::to_string(wAvgCL));
        output.push_back(std::to_string(dispCL));
        output.push_back(std::to_string(nAvgMW));
        output.push_back(std::to_string(wAvgMW));
        output.push_back(std::to_string(dispMW));

        // If no monomer or homopolymer, don't add copolymer stats
        auto numMonomers = registry::getNumMonomers();
        if (numMonomers <= 1)
            return output;

        // Monomer composition
        for (size_t i = 0; i < numMonomers; ++i)
            output.push_back(std::to_string(nAvgComp[i]));

        // Number-averaged sequence lengths
        for (size_t i = 0; i < numMonomers; ++i)
            output.push_back(std::to_string(nAvgSL[i]));

        // Weight-averaged sequence lengths
        for (size_t i = 0; i < numMonomers; ++i)
            output.push_back(std::to_string(wAvgSL[i]));

        // Dispersity of sequence lengths
        for (size_t i = 0; i < numMonomers; ++i)
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
        // If no monomer or homopolymer, return empty vector
        if (registry::getNumMonomers() <= 1)
            return {};
        auto monomerNames = registry::getMonomerNames();

        // Basic KMC state info and bucket index
        std::vector<std::string> names = {
            std::string(C::state::ITERATION_KEY),
            std::string(C::state::KMC_TIME_KEY),
            std::string(C::state::BUCKET_KEY)};

        // Monomer counts
        for (const auto &monomerName : monomerNames)
            names.push_back(std::string(C::state::MONCOUNT_PREFIX) + monomerName);

        // Sequence counts
        for (const auto &monomerName : monomerNames)
            names.push_back(std::string(C::state::SEQCOUNT_PREFIX) + monomerName);

        // Sum of squared sequence lengths
        for (const auto &monomerName : monomerNames)
            names.push_back(std::string(C::state::SEQLEN2_PREFIX) + monomerName);

        return names;
    }

    /*
    Iteration, KMC Time, Bucket, monCount_A, monCount_B, ...,
    seqCount_A, seqCount_B, ..., seqLengths2_A, seqLengths2_B, ...
    */
    std::vector<std::string> getDataAsVector(size_t bucket) const
    {
        // If no monomer or homopolymer, return empty vector
        auto numMonomers = registry::getNumMonomers();
        if (numMonomers <= 1)
            return {};

        // Basic KMC state info and bucket index
        std::vector<std::string> output;
        output.push_back(std::to_string(kmcState.iteration));
        output.push_back(std::to_string(kmcState.kmcTime));
        output.push_back(std::to_string(bucket));

        // Monomer counts
        for (size_t i = 0; i < numMonomers; ++i)
            output.push_back(std::to_string(stats[bucket].monCounts[i]));

        // Sequence counts
        for (size_t i = 0; i < numMonomers; ++i)
            output.push_back(std::to_string(stats[bucket].seqCounts[i]));

        // Sum of squared sequence lengths
        for (size_t i = 0; i < numMonomers; ++i)
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