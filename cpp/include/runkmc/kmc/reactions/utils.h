#pragma once
#include "common.h"
#include "kmc/species/polymer_type.h"

namespace rxn_print
{

    std::string speciesToString(const std::string name, const uint64_t count, bool with_counts = true)
    {
        if (with_counts)
            return name + " (" + std::to_string(count) + ")";
        return name;
    }

    std::vector<std::string> getReactantStrings(const std::vector<Unit *> &unitReactants, const std::vector<PolymerContainerPtr> &polyReactants, bool with_counts = true)
    {
        std::vector<std::string> reactantStrings;
        for (const auto &polyReactant : polyReactants)
            reactantStrings.push_back(speciesToString(polyReactant->name, polyReactant->count, with_counts));
        for (const auto &unitReactant : unitReactants)
            reactantStrings.push_back(speciesToString(unitReactant->name, unitReactant->count, with_counts));
        return reactantStrings;
    }

    std::vector<std::string> getProductStrings(const std::vector<Unit *> &unitProducts, const std::vector<PolymerContainerPtr> &polyProducts, bool with_counts = true)
    {
        std::vector<std::string> productStrings;
        for (const auto &polyProduct : polyProducts)
            productStrings.push_back(speciesToString(polyProduct->name, polyProduct->count, with_counts));
        for (const auto &unitProduct : unitProducts)
            productStrings.push_back(speciesToString(unitProduct->name, unitProduct->count, with_counts));
        return productStrings;
    }

    std::string reactionToString(
        const std::vector<Unit *> &unitReactants,
        const std::vector<PolymerContainerPtr> &polyReactants,
        const std::vector<Unit *> &unitProducts,
        const std::vector<PolymerContainerPtr> &polyProducts,
        bool with_counts = false)
    {

        std::string reactionString = "";

        auto reactantStrings = getReactantStrings(unitReactants, polyReactants, with_counts);
        auto productStrings = getProductStrings(unitProducts, polyProducts, with_counts);

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
