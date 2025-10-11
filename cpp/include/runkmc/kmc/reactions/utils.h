#pragma once
#include "common.h"
#include "kmc/species/polymer_type.h"

namespace rxn_print
{

    std::string speciesToString(std::string_view name, const uint64_t count, bool with_counts = true)
    {
        console::log("Getting species string for " + std::string(name));
        console::log("Count: " + std::to_string(count) + ", with_counts: " + (with_counts ? "true" : "false"));
        if (with_counts)
            return std::string(name) + " (" + std::to_string(count) + ")";
        console::debug("Returning species name without count.");
        return std::string(name);
    }

    std::vector<std::string> getReactantStrings(const std::vector<Unit *> &unitReactants, const std::vector<PolymerContainer *> &polyReactants, bool with_counts = true)
    {
        console::debug("Getting reactant strings...");
        std::vector<std::string> reactantStrings;
        for (const auto &polyReactant : polyReactants)
            reactantStrings.push_back(speciesToString(polyReactant->name, polyReactant->count, with_counts));
        console::debug("Got polymer reactant strings.");
        for (const auto &unitReactant : unitReactants)
            reactantStrings.push_back(speciesToString(unitReactant->name, unitReactant->count, with_counts));
        console::debug("Got unit reactant strings.");
        return reactantStrings;
    }

    std::vector<std::string> getProductStrings(const std::vector<Unit *> &unitProducts, const std::vector<PolymerContainer *> &polyProducts, bool with_counts = true)
    {
        std::vector<std::string> productStrings;
        console::log("Getting product strings...");
        for (const auto &polyProduct : polyProducts)
            productStrings.push_back(speciesToString(polyProduct->name, polyProduct->count, with_counts));
        console::log("Got polymer product strings.");
        for (const auto &unitProduct : unitProducts)
            productStrings.push_back(speciesToString(unitProduct->name, unitProduct->count, with_counts));
        console::log("Got unit product strings.");
        return productStrings;
    }

    std::string reactionToString(
        const std::vector<Unit *> &unitReactants,
        const std::vector<PolymerContainer *> &polyReactants,
        const std::vector<Unit *> &unitProducts,
        const std::vector<PolymerContainer *> &polyProducts,
        bool with_counts = false)
    {

        std::string reactionString = "";

        console::log("Generating reaction string...");

        auto reactantStrings = getReactantStrings(unitReactants, polyReactants, with_counts);
        console::log("Reactant strings:" + str::join(reactantStrings, ", "));
        auto productStrings = getProductStrings(unitProducts, polyProducts, with_counts);

        console::logVector(reactantStrings);
        console::logVector(productStrings);

        // Add reactants
        if (!reactantStrings.empty())
        {
            for (size_t i = 0; i < reactantStrings.size() - 1; ++i)
                reactionString += reactantStrings[i] + " + ";
            reactionString += reactantStrings[reactantStrings.size() - 1];
        }

        reactionString += " -> ";

        // Add products
        if (!productStrings.empty())
        {
            for (size_t i = 0; i < productStrings.size() - 1; ++i)
                reactionString += productStrings[i] + " + ";
            reactionString += productStrings[productStrings.size() - 1];
        }

        return reactionString;
    }
}
