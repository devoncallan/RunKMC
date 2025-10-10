#pragma once
#include "common.h"
#include "kmc/state.h"
#include "kmc/reactions/reaction_set.h"
#include "kmc/species/species_set.h"
#include "kmc/analysis/analysis.h"
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

        output::writeStateHeaders(paths, config);

        state.species = speciesSet.getStateData();
    }

    void run()
    {
        if (reactionSet.cantProceed())
            console::error("No reactions can occur with the initial species set. Stopping simulation.");

        startTime = std::chrono::steady_clock::now();

        // Print initial state
        output::writeState(state, paths, config);

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

            // Analyze current state
            updateSystemState();

            output::writeState(state, paths, config);
        }

        if (config.reportPolymers)
            output::writePolymers(paths, speciesSet);
    }

    const io::types::CommandLineConfig &getConfig() const { return config; };
    const io::types::SimulationConfig &getOptions() const { return options; };
    const SimulationPaths &getPaths() const { return paths; };
    const SystemState &getState() const { return state; };
    const SpeciesSet &getSpeciesSet() const { return speciesSet; };
    const ReactionSet &getReactionSet() const { return reactionSet; };

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

        reactionSet.updateReactionProbabilities(state.kmc.NAV);

        if (reactionSet.cantProceed())
            return;

        // Update time
        double rn = rng_utils::dis(rng_utils::rng) + 1e-40;
        state.kmc.kmcTime -= log(rn) / reactionSet.getTotalReactionRate();
        state.kmc.kmcStep += 1;
    }

    // ********** State functions **********

    void updateSystemState()
    {
        auto currentTime = std::chrono::steady_clock::now();
        state.kmc.simulationTime = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime - startTime).count() / 1000.;
        if (state.kmc.kmcStep > 0)
            state.kmc.simulationTimePer1e6Steps = state.kmc.simulationTime / (state.kmc.kmcStep / 1e6);

        state.species = speciesSet.getStateData();

        analysis::analyze(speciesSet, state);
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

    // Simulation start time
    std::chrono::steady_clock::time_point startTime;
};