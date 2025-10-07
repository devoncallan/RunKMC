#pragma once

#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>

#include "common.h"
#include "kmc/kmc.h"
#include "io/text_parsers.h"

namespace builder
{

    KMC fromFile(config::CommandLineConfig config)
    {
        auto sections = io::parsers::text::parseFile(config.inputFilepath);
        config::SimulationConfig simConfig = io::parsers::text::parseSimulationConfig(sections.parameters);

        auto speciesSetRead = io::parsers::text::parseSpeciesSet(sections.species);
        SpeciesSet speciesSet = buildSpeciesSet(speciesSetRead);

        registry::builder.build();

        auto rateConstants = io::parsers::text::parseRateConstants(sections.rateConstants);
        auto reactionsRead = io::parsers::text::parseReactions(sections.reactions);
        ReactionSet reactionSet = buildReactionSet(reactionsRead, speciesSet);
    };

    SpeciesSet buildSpeciesSet(const types::SpeciesSetRead &data)
    {
        size_t numUnits = data.units.size();
        size_t numPolyTypes = data.polymerTypes.size();
        size_t numPolyLabels = data.polymerLabels.size();

        std::vector<types::UnitRead> sortedUnits = data.units;
        std::stable_sort(sortedUnits.begin(), sortedUnits.end(), [](const types::UnitRead &a, const types::UnitRead &b)
                         {
            auto priority = [](const std::string &type) {
                if (type == SpeciesType::MONOMER) return 0;
                if (type == SpeciesType::INITIATOR) return 1;
                if (type == SpeciesType::UNIT) return 2;
                if (type == SpeciesType::POLYMER) return 3;
                return 4;
            };
            return priority(a.type) < priority(b.type); });

        for (const auto &unit : sortedUnits)
        {
            registry::builder.registerNewSpecies(unit.name, unit.type);
        }

        for (const auto &polyType : data.polymerTypes)
        {
            registry::builder.registerNewSpecies(polyType.name, polyType.type);
            for (const auto &unitName : polyType.endGroupUnitNames)
            {
                if (!registry::builder.isRegistered(unitName))
                    console::input_error("End group unit " + unitName + " for polymer " + polyType.name + " is not registered. Exiting.");
            }
        }

        // for (const auto &label : data.polymerLabels)
        // {
        //     registry::registerNewSpecies(label.name, label.type);
        //     for (const auto &polyName : label.polymerNames)
        //     {
        //         if (!registry::isRegistered(polyName))
        //             console::input_error("Polymer " + polyName + " for label " + label.name + " is not registered. Exiting.");
        //     }
        // }
    }

    std::vector<RateConstant> buildRateConstants(const std::vector<types::RateConstantRead> &data)
    {
        std::vector<RateConstant> rateConstants;
        rateConstants.reserve(data.size());

        for (const auto &rc : data)
        {
            RateConstant rateConstant(rc.name, rc.k);
            rateConstants.push_back(rateConstant);
        }
        return rateConstants;
    }

    ReactionSet buildReactionSet(const types::ReactionSetRead &data, SpeciesSet &speciesSet)
    {
        std::vector<Reaction *> reactions;

        std::vector<RateConstant> rateConstants = buildRateConstants(data.rateConstants);

        for (const auto &reactionData : data.reactions)
        {

            std::vector<Unit *> unitReactants;
            std::vector<Unit *> unitProducts;
            std::vector<PolymerTypeGroupPtr> polyReactants;
            std::vector<PolymerTypeGroupPtr> polyProducts;
            unitReactants.reserve(3);
            unitProducts.reserve(3);
            polyReactants.reserve(3);
            polyProducts.reserve(3);

            // Find rate constant
            size_t index = input::findInVector(reactionData.rateConstantName, rateConstants);
            // if (index <

            // Find reactants
            for (const auto &reactantName : reactionData.reactantNames)
            {
                size_t index = input::findInVector(reactantName, speciesSet.getUnits());
                if (index < speciesSet.getUnits().size())
                    unitReactants.push_back(&speciesSet.getUnits()[index]);
                else
                {
                    index = input::findInVector(reactantName, speciesSet.getPolymerTypeGroups());
                    if (index < speciesSet.getPolymerTypeGroups().size())
                        polyReactants.push_back(&speciesSet.getPolymerTypeGroups()[index]);
                    else
                        console::input_error("Species " + reactantName + " not found. Exiting.");
                }
            }

            // Find products
            for (const auto &productName : reactionData.productNames)
            {
                size_t index = input::findInVector(productName, speciesSet.getUnits());
                if (index < speciesSet.getUnits().size())
                    unitProducts.push_back(&speciesSet.getUnits()[index]);
                else
                {
                    index = input::findInVector(productName, speciesSet.getPolymerTypeGroups());
                    if (index < speciesSet.getPolymerTypeGroups().size())
                        polyProducts.push_back(&speciesSet.getPolymerTypeGroups()[index]);
                    else
                        console::input_error("Species " + productName + " not found. Exiting.");
                }
            }
            uint8_t sameReactant = 0;

            if (reactionData.type == Elementary::TYPE)
                reactions.push_back(new Elementary(rateConstant, unitReactants, unitProducts));
            else if (reactionData.type == InitiatorDecomposition::TYPE)
                reactions.push_back(new InitiatorDecomposition(rateConstant, unitReactants[0], unitProducts[0], unitProducts[1], unitReactants[0]->efficiency));
        }

        ReactionSet reactionSet(reactions, rateConstants);
        return reactionSet;
    };
};