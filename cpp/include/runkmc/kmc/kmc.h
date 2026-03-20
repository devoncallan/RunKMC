#pragma once
#include <map>

#include "common.h"
#include "kmc/state.h"
#include "kmc/reactions/reaction_set.h"
#include "kmc/species/species_set.h"
#include "kmc/analysis/analysis.h"
#include "kmc/plugin.h"
#include "results/state.h"
#include "results/polymers.h"

/**
 * @brief Kinetic Monte Carlo simulation class
 *
 */
class KMC
{
public:
    KMC(SpeciesSet &species, ReactionSet &reactions, io::types::CommandLineConfig config_, io::types::SimulationConfig options_)
        : speciesSet(std::move(species)), reactionSet(std::move(reactions)), config(config_), options(options_)
    {
        paths = SimulationPaths(config_);
        state = SystemState();

        state.kmc.NAV = speciesSet.getNAV();

        speciesSet.updatePolymerContainers();

        reactionSet.updateReactionProbabilities(state.kmc.NAV);

        // Initialize per-container accumulators for reporting containers
        for (const auto *container : speciesSet.getReportingContainers())
        {
            containerChainStats[container->name] = analysis::ChainStats{};
            containerPosStats[container->name] = analysis::PositionalSequenceStats{};
        }

        output::writeStateHeaders(paths, config);

        state.species = speciesSet.getStateData();
    }

    void run()
    {
        for (auto &plugin : plugins)
            plugin->onSimulationStart(*this);

        if (reactionSet.cantProceed())
            console::error("No reactions can occur with the initial species set. Stopping simulation.");

        startTime = std::chrono::steady_clock::now();

        // Print initial state
        output::writeResults(state, paths);

        // Main simulation loop
        while (state.kmc.kmcTime < options.terminationTime)
        {
            state.kmc.iteration += 1;

            auto targetTime = state.kmc.kmcTime + options.analysisTime;
            bool success = runToTime(targetTime);

            if (!success)
            {
                console::warning(
                    "Could not reach termination time (" +
                    std::to_string(options.terminationTime) +
                    ") - no more reactions can occur. Stopping simulation at " +
                    std::to_string(state.kmc.kmcTime) + ".");
                break;
            }

            // Snapshot species counts before sinks are cleared
            snapshotSpeciesState();

            // Accumulate chain stats, write per-container files, clear sinks
            writeReportingContainers();

            // Copy accumulated chain stats into state then write
            updateSystemState();

            output::writeResults(state, paths);
        }

        // Final flush of reporting containers
        writeReportingContainers();

        if (config.reportPolymers)
            output::writePolymers(paths, speciesSet);

        for (auto &plugin : plugins)
            plugin->onSimulationEnd(*this);
    }

    double getNAV() const { return state.kmc.NAV; }

    void setNAV(double newNAV)
    {
        if (newNAV <= 0)
            console::error("Attempted to set NAV to non-positive value: " + std::to_string(newNAV) + ".");
        state.kmc.NAV = newNAV;
        reactionSet.setNAV(newNAV);
    }

    void scaleNAV(double factor)
    {
        if (factor <= 0)
            console::error("Attempted to scale NAV by non-positive factor: " + std::to_string(factor) + ".");
        setNAV(state.kmc.NAV * factor);
    }

    void registerPlugin(SimulationPluginPtr plugin)
    {
        if (!plugin)
            return;

        plugin->onRegistered(*this);
        plugins.emplace_back(std::move(plugin));
    }

    const io::types::CommandLineConfig &getConfig() const { return config; };
    const io::types::SimulationConfig &getOptions() const { return options; };
    const SimulationPaths &getPaths() const { return paths; };
    const SystemState &getState() const { return state; };
    const SpeciesSet &getSpeciesSet() const { return speciesSet; };
    const ReactionSet &getReactionSet() const { return reactionSet; };

    // Returns cumulative chain stats for a reporting container by name (or empty stats if not found)
    const analysis::ChainStats &getContainerChainStats(const std::string &name) const
    {
        static analysis::ChainStats empty{};
        auto it = containerChainStats.find(name);
        return (it != containerChainStats.end()) ? it->second : empty;
    }

    const std::map<std::string, analysis::ChainStats> &getAllContainerChainStats() const
    {
        return containerChainStats;
    }

private:
    // ********** Simulation functions **********

    bool runToTime(double time)
    {
        while (state.kmc.kmcTime < time)
        {
            if (reactionSet.cantProceed())
                return false;

            step();
        }
        return true;
    }

    // Core Kinetic Monte Carlo Simulation Step
    void step()
    {
        size_t reactionIndex = reactionSet.chooseRandomReactionIndex();

        Reaction *reaction = reactionSet.getReaction(reactionIndex);

        reaction->react();

        speciesSet.updatePolymerContainers();

        for (auto &plugin : plugins)
            plugin->afterReaction(*this, *reaction, reactionIndex);

        for (auto &plugin : plugins)
            plugin->beforePropensityUpdate(*this, reactionSet);

        reactionSet.updateReactionProbabilities(state.kmc.NAV);

        if (reactionSet.cantProceed())
            return;

        // Update time
        double rn = rng::rand();
        state.kmc.kmcTime -= log(rn) / reactionSet.getTotalReactionRate();
        state.kmc.kmcStep += 1;
    }

    // ********** Output functions **********

    // Process all reporting containers: accumulate stats, write records, clear sinks
    void writeReportingContainers()
    {
        const auto &FWs = registry::getMonomerFWs();

        for (auto *container : speciesSet.getReportingContainers())
        {
            const auto &name = container->name;
            const auto &polymers = container->getPolymers();

            if (polymers.empty())
                continue;

            // Accumulate into persistent chain stats
            auto &ccs = containerChainStats[name];
            for (const auto *polymer : polymers)
                ccs += polymer->toChainStats(FWs);

            // Accumulate into per-interval positional stats (copolymer + multi-bucket only)
            if (config.reportPositionalStats && NUM_BUCKETS > 1)
            {
                auto &pcs = containerPosStats[name];
                for (const auto *polymer : polymers)
                    pcs += polymer->getPositionalStats();
            }

            // Write per-chain records
            if (config.reportChains)
                output::writeChainRecords(*container, state.kmc, paths);

            // Write positional stats block
            if (config.reportPositionalStats && NUM_BUCKETS > 1)
            {
                output::writePosChainStats(containerPosStats[name], name, state.kmc, paths);
                containerPosStats[name] = analysis::PositionalSequenceStats{};
            }

            // Write segment histogram
            if (config.reportSegmentHistogram)
                output::writeSegmentHistogram(*container, state.kmc, paths);

            // Clear if this container is a sink
            if (container->isSink)
                container->clearPolymers();
        }
    }

    // ********** State functions **********

    // Snapshot species counts and timing into state (call before clearing sinks)
    void snapshotSpeciesState()
    {
        auto currentTime = std::chrono::steady_clock::now();
        state.kmc.simulationTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.;
        if (state.kmc.kmcStep > 0)
            state.kmc.simulationTimePer1e6Steps = state.kmc.simulationTime / (state.kmc.kmcStep / 1e6);
        state.species = speciesSet.getStateData();
    }

    // Aggregate all reporting container chain stats into state.chainStats (call after writeReportingContainers)
    void updateSystemState()
    {
        analysis::ChainStats aggregate{};
        for (const auto &[name, ccs] : containerChainStats)
            aggregate += ccs;
        state.chainStats.stats = aggregate;
    }

    // Simulation inputs
    io::types::CommandLineConfig config;
    io::types::SimulationConfig options;

    // Managing outputs
    SimulationPaths paths;
    SystemState state;

    // Core simulation objects
    ReactionSet reactionSet;
    SpeciesSet speciesSet;

    std::vector<SimulationPluginPtr> plugins;

    // Per-container accumulators
    std::map<std::string, analysis::ChainStats> containerChainStats;         // persistent across flushes
    std::map<std::string, analysis::PositionalSequenceStats> containerPosStats; // reset each interval

    // Simulation start time
    std::chrono::steady_clock::time_point startTime;
};
