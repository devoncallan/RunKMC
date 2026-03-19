#pragma once

#include <algorithm>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include "common.h"
#include "kmc/kmc.h"
#include "kmc/plugin.h"
#include "kmc/species/types.h"

/**
 * @brief Volume management plugin.
 *
 * Consumes simple density data to estimate how each reaction changes
 * system volume. The plugin keeps the stochastic scaling (NAV) in sync
 * with the evolving volume while leaving the base simulation code
 * unchanged.
 */
class VolumePlugin : public SimulationPlugin
{
public:
    struct SpeciesConfig
    {
        std::string name;
        double monomerDensity = 0.0; // g/cm^3
        double polymerDensity = 0.0; // g/cm^3
    };

    struct Config
    {
        double initialVolumeL = 1.0;
        std::vector<SpeciesConfig> species;
    };

    explicit VolumePlugin(Config config_) : config(std::move(config_)) {}

    std::string_view getName() const override { return pluginName; }

    void onRegistered(KMC &kmc) override
    {
        cacheUnitCoefficients(kmc.getSpeciesSet());
        precomputeReactionDeltas(kmc.getReactionSet());

        initialVolume = config.initialVolumeL;
        if (initialVolume <= 0.0)
        {
            console::warning("VolumePlugin: initial volume must be positive. Resetting to 1 L.");
            initialVolume = 1.0;
        }
        initialNAV = kmc.getNAV();
        currentVolume = initialVolume;
    }

    void onSimulationStart(KMC &kmc) override
    {
        initialNAV = kmc.getNAV();
        currentVolume = initialVolume;
    }

    void afterReaction(KMC &kmc, const Reaction &, size_t reactionIndex) override
    {
        if (reactionIndex >= volumeDeltas.size())
            return;

        const double deltaVolume = volumeDeltas[reactionIndex];
        if (deltaVolume == 0.0)
            return;

        currentVolume += deltaVolume;
        if (currentVolume <= 0.0)
        {
            console::warning("VolumePlugin: volume became non-positive (" + std::to_string(currentVolume) + "). Clamping to epsilon.");
            currentVolume = std::numeric_limits<double>::min();
        }

        const double newNAV = initialNAV * (initialVolume / currentVolume);
        kmc.setNAV(newNAV);
    }

    double getCurrentVolume() const { return currentVolume; }
    double getInitialVolume() const { return initialVolume; }

private:
    void cacheUnitCoefficients(const SpeciesSet &speciesSet)
    {
        unitVolumeCoeff.clear();
        const auto &units = speciesSet.getUnits();

        for (const auto &entry : config.species)
        {
            auto it = std::find_if(
                units.begin(),
                units.end(),
                [&](const Unit &unit) { return unit.name == entry.name; });

            if (it == units.end())
            {
                console::warning("VolumePlugin: unit '" + entry.name + "' not found. Skipping.");
                continue;
            }

            if (entry.monomerDensity <= 0 || entry.polymerDensity <= 0)
            {
                console::warning("VolumePlugin: densities for '" + entry.name + "' must be positive. Skipping.");
                continue;
            }

            const double coeff = it->FW * (1.0 / entry.polymerDensity - 1.0 / entry.monomerDensity);
            unitVolumeCoeff[it->ID] = coeff;
        }
    }

    void precomputeReactionDeltas(const ReactionSet &reactionSet)
    {
        const size_t numReactions = reactionSet.getNumReactions();
        volumeDeltas.assign(numReactions, 0.0);

        for (size_t i = 0; i < numReactions; ++i)
        {
            const Reaction *reaction = reactionSet.getReaction(i);
            if (!reaction)
                continue;

            const auto &species = reaction->getSpecies();
            double delta = 0.0;

            for (const auto *reactant : species.reactants)
            {
                if (!reactant || !SpeciesType::isUnitType(reactant->type))
                    continue;
                const auto *unit = static_cast<const Unit *>(reactant);
                auto coeffIt = unitVolumeCoeff.find(unit->ID);
                if (coeffIt != unitVolumeCoeff.end())
                    delta -= coeffIt->second;
            }

            for (const auto *product : species.products)
            {
                if (!product || !SpeciesType::isUnitType(product->type))
                    continue;
                const auto *unit = static_cast<const Unit *>(product);
                auto coeffIt = unitVolumeCoeff.find(unit->ID);
                if (coeffIt != unitVolumeCoeff.end())
                    delta += coeffIt->second;
            }

            volumeDeltas[i] = delta;
        }
    }

    static inline constexpr std::string_view pluginName = "VolumePlugin";

    Config config;
    std::unordered_map<SpeciesID, double> unitVolumeCoeff;
    std::vector<double> volumeDeltas;
    double initialVolume = 1.0;
    double currentVolume = 1.0;
    double initialNAV = 1.0;
};
