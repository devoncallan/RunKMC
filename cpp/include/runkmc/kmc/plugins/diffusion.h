#pragma once

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "common.h"
#include "kmc/kmc.h"
#include "kmc/plugin.h"
#include "kmc/reactions/reaction_set.h"
#include "kmc/species/types.h"

/**
 * @brief Diffusion-controlled kinetics plugin.
 *
 * Applies free-volume-based diffusion corrections to kp, kt, and f
 * at a configurable step interval. Modeled after the FRP1_Diff equations
 * from the sparks-hybrid codebase.
 *
 * Parameters are loaded from a separate YAML data file specified in the
 * main input under the plugins: section.
 */
class DiffusionPlugin : public SimulationPlugin
{
public:
    struct MonomerDiffProps
    {
        std::string name;
        double Tg_m = 0.0;    // Glass transition temp of monomer (K)
        double Tg_p = 0.0;    // Glass transition temp of polymer (K)
        double Vf_m = 0.025;  // Free volume of pure monomer
        double Vf_p = 0.025;  // Free volume of pure polymer
        double alpha_m = 0.001; // Thermal expansion coefficient monomer (1/K)
        double alpha_p = 0.00048; // Thermal expansion coefficient polymer (1/K)
        double Vf_c = 0.05;   // Critical free volume for kp/kfm diffusion control
        double B = 1.0;       // Translational diffusion decay for kp
        double Vf_i = 0.05;   // Critical free volume for f diffusion control
        double C = 1.0;       // Initiation diffusion decay
        double delta = 0.001; // Segmental diffusion reaction radius (L/g)
        double m = 0.5;       // Gel-effect model exponent (Mw)
        double n = 1.75;      // Post-gel translational diffusion exponent
        double A = 1.11;      // Gel-effect model parameter
        double K3 = 1e10;     // Gel-point onset criterion
        double n_s = 0.9;     // Avg monomer units per chain segment
        double l0 = 2.52e-8;  // Monomer unit length (cm)
        double rho_m = 1.0;   // Density of monomer (g/cm^3)
        double rho_p = 1.0;   // Density of polymer (g/cm^3)
    };

    struct Config
    {
        std::vector<MonomerDiffProps> monomers;
        uint64_t updateInterval = 1;
        double temperature = 60.0; // °C
    };

    explicit DiffusionPlugin(Config config_) : config(std::move(config_)) {}

    std::string_view getName() const override { return pluginName; }

    void onRegistered(KMC &kmc) override
    {
        const auto &units = kmc.getSpeciesSet().getUnits();
        const auto &reactionSet = kmc.getReactionSet();
        const size_t n = reactionSet.getNumReactions();

        monomerIndices.clear();
        for (const auto &props : config.monomers)
        {
            auto it = std::find_if(units.begin(), units.end(),
                                   [&](const Unit &u) { return u.name == props.name; });
            if (it == units.end())
            {
                console::warning("DiffusionPlugin: monomer '" + props.name + "' not found in species set.");
                monomerIndices.push_back(SIZE_MAX);
            }
            else
            {
                monomerIndices.push_back(static_cast<size_t>(it - units.begin()));
            }
        }

        // Classify reactions into buckets for targeted multiplier application
        propagationReactions.clear();
        terminationDReactions.clear();
        terminationCReactions.clear();
        initiationReactions.clear();
        initiatorDecompReactions.clear();

        for (size_t i = 0; i < n; ++i)
        {
            const Reaction *r = reactionSet.getReaction(i);
            if (!r) continue;
            const auto type = r->getType();
            if (type == ReactionType::PROPAGATION || type == ReactionType::CHAINTRANSFER_M)
                propagationReactions.push_back(i);
            else if (type == ReactionType::TERMINATION_D)
                terminationDReactions.push_back(i);
            else if (type == ReactionType::TERMINATION_C)
                terminationCReactions.push_back(i);
            else if (type == ReactionType::INITIATION)
                initiationReactions.push_back(i);
            else if (type == ReactionType::INITIATOR_DECOMPOSITION ||
                     type == ReactionType::INIT_DECOMP_POLY)
                initiatorDecompReactions.push_back(i);
        }
    }

    void onSimulationStart(KMC &kmc) override
    {
        gelPointReached = false;
        Mw_gel = 0.0;
        Vf_gel = 0.0;
        kt_gel_factor = 1.0;
    }

    void beforePropensityUpdate(KMC &kmc, ReactionSet &reactions) override
    {
        if (kmc.getState().kmc.kmcStep % config.updateInterval != 0)
            return;
        computeAndApply(kmc, reactions);
    }

private:
    void computeAndApply(KMC &kmc, ReactionSet &reactions)
    {
        if (config.monomers.empty() || monomerIndices.empty())
            return;

        const double T_K = config.temperature + 273.15;
        const double NAV = kmc.getNAV();
        const auto &units = kmc.getSpeciesSet().getUnits();

        // Use first monomer's properties (for homopolymer; extend to copolymer later)
        const size_t mIdx = monomerIndices[0];
        if (mIdx == SIZE_MAX) return;
        const MonomerDiffProps &props = config.monomers[0];
        const Unit &monomer = units[mIdx];

        // V [L] = NAV / NA  (since NAV = NA * V)
        // rho [g/cm^3] * 1000 = [g/L]; vol [L] = moles * FW [g/mol] / (rho * 1000 [g/L])
        const double V = NAV / C::NA;
        const double nM_mol = static_cast<double>(monomer.count) / C::NA;
        const double consumed = static_cast<double>(monomer.getInitCount()) - static_cast<double>(monomer.count);
        const double nPM_mol = consumed / C::NA;
        const double vol_m_L = nM_mol * monomer.FW / (props.rho_m * 1000.0);
        const double vol_p_L = nPM_mol * monomer.FW / (props.rho_p * 1000.0);

        if (V <= 0.0) return;

        // Free volume
        const double Vf_m_contrib = (props.Vf_m + props.alpha_m * (T_K - props.Tg_m)) * vol_m_L / V;
        const double Vf_p_contrib = (props.Vf_p + props.alpha_p * (T_K - props.Tg_p)) * vol_p_L / V;
        const double Vf = Vf_m_contrib + Vf_p_contrib;

        if (Vf <= 0.0) return;

        // --- kp/kfm diffusion factor ---
        double kp_factor = 1.0;
        if (Vf < props.Vf_c)
            kp_factor = std::exp(-props.B * (1.0 / Vf - 1.0 / props.Vf_c));

        // --- f (initiator efficiency) diffusion factor ---
        double f_factor = 1.0;
        if (Vf < props.Vf_i)
            f_factor = std::exp(-props.C * (1.0 / Vf - 1.0 / props.Vf_i));

        // --- kt diffusion factor ---
        const double c_poly = nPM_mol * monomer.FW / V; // [g/L]
        double kt_factor = 1.0;
        const double Mw = kmc.getState().chainStats.stats.chainMW.wAvg();

        if (!gelPointReached)
        {
            kt_factor = 1.0 + props.delta * c_poly;

            const double K3_test = std::pow(Mw + 1.0, props.m) * std::exp(props.A / Vf);
            if (K3_test >= props.K3)
            {
                gelPointReached = true;
                Mw_gel = Mw;
                Vf_gel = Vf;
                kt_gel_factor = kt_factor;
                console::log("DiffusionPlugin: gel point reached at Mw=" + std::to_string(Mw) +
                             ", Vf=" + std::to_string(Vf));
            }
        }
        else
        {
            // Post-gel: translational diffusion
            const double kt_trans = kt_gel_factor *
                                    std::pow((Mw_gel + 1.0) / (Mw + 1.0), props.n) *
                                    std::exp(-props.A * (1.0 / Vf - 1.0 / Vf_gel));
            // Reaction-diffusion termination
            const double cM = nM_mol / V; // [mol/L]
            const double D = props.n_s * (props.l0 * props.l0) * cM / 6.0;
            const double delta_rd = std::pow(6.0 * (vol_m_L / (static_cast<double>(monomer.count) + 1.0)) / (M_PI * C::NA), 1.0 / 3.0);
            const double kt_rd = 8.0 * M_PI * C::NA * D * delta_rd / 1000.0;
            kt_factor = kt_trans + kt_rd;
        }

        // Apply multipliers
        for (size_t i : propagationReactions)
            reactions.setRateMultiplier(i, kp_factor);
        for (size_t i : terminationDReactions)
            reactions.setRateMultiplier(i, kt_factor);
        for (size_t i : terminationCReactions)
            reactions.setRateMultiplier(i, kt_factor);
        for (size_t i : initiatorDecompReactions)
            reactions.setRateMultiplier(i, f_factor);
    }

    static inline constexpr std::string_view pluginName = "DiffusionPlugin";

    Config config;
    std::vector<size_t> monomerIndices;

    std::vector<size_t> propagationReactions;
    std::vector<size_t> terminationDReactions;
    std::vector<size_t> terminationCReactions;
    std::vector<size_t> initiationReactions;
    std::vector<size_t> initiatorDecompReactions;

    bool gelPointReached = false;
    double Mw_gel = 0.0;
    double Vf_gel = 0.0;
    double kt_gel_factor = 1.0;
};

// +--------------------------
// | Data file loader
// +--------------------------
namespace diffusion
{
    inline DiffusionPlugin::MonomerDiffProps readMonomerProps(const YAML::Node &node)
    {
        DiffusionPlugin::MonomerDiffProps p;
        if (node["name"])    p.name    = node["name"].as<std::string>();
        if (node["Tg_m"])    p.Tg_m    = node["Tg_m"].as<double>();
        if (node["Tg_p"])    p.Tg_p    = node["Tg_p"].as<double>();
        if (node["Vf_m"])    p.Vf_m    = node["Vf_m"].as<double>();
        if (node["Vf_p"])    p.Vf_p    = node["Vf_p"].as<double>();
        if (node["alpha_m"]) p.alpha_m = node["alpha_m"].as<double>();
        if (node["alpha_p"]) p.alpha_p = node["alpha_p"].as<double>();
        if (node["Vf_c"])    p.Vf_c    = node["Vf_c"].as<double>();
        if (node["B"])       p.B       = node["B"].as<double>();
        if (node["Vf_i"])    p.Vf_i    = node["Vf_i"].as<double>();
        if (node["C"])       p.C       = node["C"].as<double>();
        if (node["delta"])   p.delta   = node["delta"].as<double>();
        if (node["m"])       p.m       = node["m"].as<double>();
        if (node["n"])       p.n       = node["n"].as<double>();
        if (node["A"])       p.A       = node["A"].as<double>();
        if (node["K3"])      p.K3      = node["K3"].as<double>();
        if (node["n_s"])     p.n_s     = node["n_s"].as<double>();
        if (node["l0"])      p.l0      = node["l0"].as<double>();
        if (node["rho_m"])   p.rho_m   = node["rho_m"].as<double>();
        if (node["rho_p"])   p.rho_p   = node["rho_p"].as<double>();
        return p;
    }

    inline DiffusionPlugin::Config loadDataFile(const std::string &dataFilePath,
                                                 uint64_t updateInterval,
                                                 double temperature)
    {
        if (!std::filesystem::exists(dataFilePath))
            console::input_error("DiffusionPlugin: data file not found: " + dataFilePath);

        YAML::Node root = YAML::LoadFile(dataFilePath);
        DiffusionPlugin::Config cfg;
        cfg.updateInterval = updateInterval;
        cfg.temperature = temperature;

        if (!root["monomers"] || !root["monomers"].IsSequence())
            console::input_error("DiffusionPlugin: data file '" + dataFilePath + "' must have a 'monomers' list.");

        for (const auto &node : root["monomers"])
            cfg.monomers.push_back(readMonomerProps(node));

        return cfg;
    }
}
