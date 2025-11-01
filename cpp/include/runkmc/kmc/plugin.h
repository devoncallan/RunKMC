#pragma once

#include <cstddef>
#include <memory>
#include <string_view>

class KMC;
class ReactionSet;
class Reaction;

/**
 * @brief Minimal interface for optional simulation extensions.
 *
 * A SimulationPlugin instance is owned by the KMC driver for the
 * duration of a run. Hooks are intentionally narrow and provide
 * sensible no-op defaults so callers can override only what they
 * need. Implementations should keep hook bodies lightweight – most
 * callbacks fire on the hot reaction path.
 */
class SimulationPlugin
{
public:
    virtual ~SimulationPlugin() = default;

    /**
     * @brief Identifier used for logging/debugging.
     */
    virtual std::string_view getName() const = 0;

    /**
     * @brief Called immediately after the plugin is registered.
     * Hook runs before the simulation starts so the plugin can
     * cache references or precompute data.
     */
    virtual void onRegistered(KMC &kmc) {}

    /**
     * @brief Fired once right before the simulation loop begins.
     */
    virtual void onSimulationStart(KMC &kmc) {}

    /**
     * @brief Invoked on every step after a reaction has occurred.
     *
     * The reaction index corresponds to the entry in the ReactionSet
     * selected for the event that just fired.
     */
    virtual void afterReaction(KMC &kmc, const Reaction &reaction, size_t reactionIndex) {}

    /**
     * @brief Called just before reaction propensities are recomputed.
     * Plugins can adjust global scalars (e.g., NAV, volume) or update
     * cached coefficients before the ReactionSet recalculates rates.
     */
    virtual void beforePropensityUpdate(KMC &kmc, ReactionSet &reactions) {}

    /**
     * @brief Fired once after the simulation loop terminates.
     */
    virtual void onSimulationEnd(KMC &kmc) {}
};

using SimulationPluginPtr = std::unique_ptr<SimulationPlugin>;
