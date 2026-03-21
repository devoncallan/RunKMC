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
        double rho_m = 0.0; // g/cm^3
        double rho_p = 0.0; // g/cm^3
    };

    struct Config
    {
        std::vector<SpeciesConfig> species;
    };

    explicit VolumePlugin(Config config_) : config(std::move(config_)) {}

    std::string_view getName() const override { return pluginName; }

    void onRegistered(KMC &kmc) override
    {
        initialNAV = kmc.getNAV();
        currentNAV = initialNAV;
        buildVolumeDeltas(kmc.getSpeciesSet(), kmc.getReactionSet());
    }

    void onSimulationStart(KMC &kmc) override
    {
        initialNAV = kmc.getNAV();
        currentNAV = initialNAV;
    }

    void afterReaction(KMC &kmc, const Reaction &, size_t reactionIndex) override
    {
        if (reactionIndex >= deltaNAVs.size())
            return;

        const double deltaNAV = deltaNAVs[reactionIndex];
        if (deltaNAV == 0.0)
            return;

        currentNAV += deltaNAV;
        if (currentNAV <= 0.0)
        {
            console::warning("VolumePlugin: volume became non-positive (" + std::to_string(currentNAV * C::NA) + "). Clamping to epsilon.");
            currentNAV = std::numeric_limits<double>::min();
        }
        kmc.setNAV(currentNAV);
    }

    double getCurrentVolume() const { return currentNAV / C::NA; }
    double getInitialVolume() const { return initialNAV / C::NA; }

private:
    // Look up a monomer's volume delta from a species list; returns 0 if not found
    static double getMonomerDeltaNAV(const std::vector<Species *> &speciesList,
                                     const std::unordered_map<SpeciesID, double> &coeffs)
    {
        for (const auto *s : speciesList)
        {
            if (!s || s->type != SpeciesType::MONOMER)
                continue;
            auto it = coeffs.find(static_cast<const Unit *>(s)->ID);
            return it != coeffs.end() ? it->second : 0.0;
        }
        return 0.0;
    }

    void buildVolumeDeltas(const SpeciesSet &speciesSet, const ReactionSet &reactionSet)
    {
        // Build name -> coeff map from config
        std::unordered_map<SpeciesID, double> dNAV;
        const auto &units = speciesSet.getUnits();
        for (const auto &entry : config.species)
        {
            auto it = std::find_if(units.begin(), units.end(),
                                   [&](const Unit &u)
                                   { return u.name == entry.name; });
            if (it == units.end())
            {
                console::warning("VolumePlugin: monomer '" + entry.name + "' not found. Skipping.");
                continue;
            }
            // ΔNAV = Change in NA * volume (L * molecules/mol)
            dNAV[it->ID] = it->FW * (1.0 / entry.rho_p - 1.0 / entry.rho_m) / 1000.0;
        }

        const size_t n = reactionSet.getNumReactions();
        deltaNAVs.assign(n, 0.0);

        for (size_t i = 0; i < n; ++i)
        {
            const Reaction *reaction = reactionSet.getReaction(i);
            if (!reaction)
                continue;

            const auto type = reaction->getType();
            const auto &sp = reaction->getSpecies();

            if (type == ReactionType::PROPAGATION ||
                type == ReactionType::INITIATION ||
                type == ReactionType::CHAINTRANSFER_M)
                deltaNAVs[i] = getMonomerDeltaNAV(sp.reactants, dNAV);
            else if (type == ReactionType::DEPROPAGATION)
                deltaNAVs[i] = -getMonomerDeltaNAV(sp.products, dNAV);
        }
    }

    static inline constexpr std::string_view pluginName = "VolumePlugin";

    Config config;
    std::vector<double> deltaNAVs;
    double currentNAV = 1.0;
    double initialNAV = 1.0;
};
