#pragma once
#include "common.h"
#include "kmc/analysis/types.h"
#include "kmc/analysis/analysis.h"

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

struct ChainStatsState
{
    analysis::ChainStats stats;

    static std::vector<std::string> getTitles()
    {
        std::vector<std::string> names = {
            std::string(C::state::NAVGCL_KEY),
            std::string(C::state::WAVGCL_KEY),
            std::string(C::state::DISPCL_KEY),
            std::string(C::state::NAVGMW_KEY),
            std::string(C::state::WAVGMW_KEY),
            std::string(C::state::DISPMW_KEY)};

        if (registry::getNumMonomers() <= 1)
            return names;

        const auto monomerNames = registry::getMonomerNames();
        for (const auto &name : monomerNames)
            names.push_back(std::string(C::state::NAVGCOMP_PREFIX) + name);
        for (const auto &name : monomerNames)
            names.push_back(std::string(C::state::NAVGSL_PREFIX) + name);
        for (const auto &name : monomerNames)
            names.push_back(std::string(C::state::WAVGSL_PREFIX) + name);
        for (const auto &name : monomerNames)
            names.push_back(std::string(C::state::DISPSL_PREFIX) + name);

        return names;
    }

    std::vector<std::string> getDataAsVector() const
    {
        using analysis::safeDivide;
        std::vector<std::string> output;
        output.push_back(std::to_string(stats.chainLength.nAvg()));
        output.push_back(std::to_string(stats.chainLength.wAvg()));
        output.push_back(std::to_string(stats.chainLength.disp()));

        const bool hasMW = stats.chainMW.count > 0;
        output.push_back(std::to_string(hasMW ? stats.chainMW.nAvg() : stats.chainLength.nAvg()));
        output.push_back(std::to_string(hasMW ? stats.chainMW.wAvg() : stats.chainLength.wAvg()));
        output.push_back(std::to_string(hasMW ? stats.chainMW.disp() : stats.chainLength.disp()));

        const auto numMonomers = stats.sequenceLengths.size();
        if (numMonomers <= 1)
            return output;

        double totalMonCount = 0;
        for (const auto &seq : stats.sequenceLengths)
            totalMonCount += seq.sum;

        for (size_t i = 0; i < numMonomers; ++i)
            output.push_back(std::to_string(safeDivide(stats.sequenceLengths[i].sum, totalMonCount)));
        for (size_t i = 0; i < numMonomers; ++i)
            output.push_back(std::to_string(stats.sequenceLengths[i].nAvg()));
        for (size_t i = 0; i < numMonomers; ++i)
            output.push_back(std::to_string(stats.sequenceLengths[i].wAvg()));
        for (size_t i = 0; i < numMonomers; ++i)
            output.push_back(std::to_string(stats.sequenceLengths[i].disp()));

        return output;
    }
};

struct SystemState
{
    KMCState kmc;
    SpeciesState species;
    ChainStatsState chainStats;
    SequenceState sequence;
};
