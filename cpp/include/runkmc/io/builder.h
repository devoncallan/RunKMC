#pragma once

#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>

#include "common.h"
#include "kmc/kmc.h"
#include "io/text_parsers.h"
#include "io/yaml_parsers.h"

namespace builder
{

    static KMC buildModel(const config::CommandLineConfig &config, const types::KMCInputRead &data)
    {
        SpeciesSet speciesSet = buildSpeciesSet(data.species);

        registry::builder.build();

        ReactionSet reactionSet = buildReactionSet(data.reactions, data.rateConstants, speciesSet);

        // KMC kmc(simConfig, speciesSet, reactionSet, config.outputDir, config.reportPolymers, config.reportSequences);
        // return kmc;
    }

    SpeciesSet buildSpeciesSet(const types::SpeciesSetRead &data)
    {
        size_t numUnits = data.units.size();
        size_t numPolyTypes = data.polymerTypes.size();
        size_t numPolyLabels = data.polymerLabels.size();
        size_t numPolyCollections = numPolyTypes + numPolyLabels;

        std::vector<types::UnitRead> unitsRead = data.units;
        std::stable_sort(
            unitsRead.begin(),
            unitsRead.end(),
            [](const types::UnitRead &a, const types::UnitRead &b)
            {
            auto priority = [](const std::string &type) {
                if (type == SpeciesType::MONOMER) return 0;
                if (type == SpeciesType::INITIATOR) return 1;
                if (type == SpeciesType::UNIT) return 2;
                if (type == SpeciesType::POLYMER) return 3;
                return 4;
            };
            return priority(a.type) < priority(b.type); });

        std::vector<Unit> units;
        std::vector<PolymerType> polymerTypes;
        std::vector<PolymerContainerMap> polymerContainers;

        // Register and create unit species
        for (const auto &unitRead : unitsRead)
        {
            SpeciesID id = registry::builder.registerNewSpecies(unitRead.name, unitRead.type);
            units.push_back(Unit(unitRead.type, unitRead.name, id, unitRead.C0, unitRead.FW, unitRead.efficiency));
        }

        // Register and create polymer types
        for (const auto &polyType : data.polymerTypes)
        {
            std::vector<SpeciesID> endGroupIDs;
            endGroupIDs.reserve(polyType.endGroupUnitNames.size());

            for (const auto &unitName : polyType.endGroupUnitNames)
            {
                if (!registry::builder.isRegistered(unitName))
                    console::input_error("End group unit " + unitName + " for polymer " + polyType.name + " is not registered. Exiting.");
                endGroupIDs.push_back(registry::builder.getSpeciesID(unitName));
            }

            SpeciesID id = registry::builder.registerNewSpecies(polyType.name, polyType.type);
            polymerTypes.push_back(PolymerType(id, polyType.name, endGroupIDs));
            polymerContainers.push_back(PolymerContainerMap(id, polyType.name));
        }

        // Create polymer collections
        for (const auto &label : data.polymerLabels)
        {
            std::vector<SpeciesID> labelPolyIDs;
            labelPolyIDs.reserve(label.polymerNames.size());

            for (const auto &polyName : label.polymerNames)
            {
                if (!registry::builder.isRegistered(polyName))
                    console::input_error("Polymer " + polyName + " for label " + label.name + " is not registered. Exiting.");
                labelPolyIDs.push_back(registry::builder.getSpeciesID(polyName));
            }

            SpeciesID id = registry::builder.registerNewSpecies(label.name, label.type);
            polymerContainers.push_back(PolymerContainerMap(id, label.name, labelPolyIDs));
        }
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

    ReactionSet buildReactionSet(const std::vector<types::ReactionRead> &reactionsRead, const std::vector<types::RateConstantRead> &rateConstantsRead, SpeciesSet &speciesSet)
    {
        std::vector<Reaction *> reactions;

        const std::vector<RateConstant> rateConstants = buildRateConstants(rateConstantsRead);

        for (const auto &reactionData : reactionsRead)
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
            if (index < rateConstants.size())
                rateConstant = rateConstants[index];
            else
                console::input_error("Rate constant " + reactionData.rateConstantName + " not found. Exiting.");
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

namespace builder::utils
{

};

// namespace build
// {
//     // Build from file
//     KMC fromFile(config::CommandLineConfig config);
//     KMC _fromModelFile(config::CommandLineConfig config);
//     KMC _fromYamlFile(config::CommandLineConfig config);

//     // Species
//     static SpeciesSet buildSpeciesSet(const types::SpeciesSetRead &data);

//     // Rate constants
//     static std::vector<RateConstant> buildRateConstants(const std::vector<types::RateConstantRead> &data);

//     // Reactions
//     static ReactionSet buildReactionSet(const types::ReactionSetRead &data, SpeciesSet &speciesSet);

// };