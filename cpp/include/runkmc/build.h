#pragma once

#include "common.h"
#include "io/text.h"
#include "io/yaml.h"
#include "io/cli.h"
#include "utils/parse.h"
#include "kmc/kmc.h"

namespace builder
{
    static SpeciesSet buildSpeciesSet(const types::SpeciesSetRead &data, const config::SimulationConfig &config);
    static std::vector<RateConstant> buildRateConstants(const std::vector<types::RateConstantRead> &data);
    static ReactionSet buildReactionSet(const std::vector<types::ReactionRead> &reactionsRead, const std::vector<types::RateConstantRead> &rateConstantsRead, SpeciesSet &speciesSet);

    static config::CommandLineConfig parseArguments(int arc, char **argv)
    {
        return io::parse::cli::parseArguments(arc, argv);
    }

    static types::KMCInputRead parseKMCInput(const std::string &filepath)
    {
        if (str::endswith(filepath, ".yaml") || str::endswith(filepath, ".yml"))
            return io::parse::parseKMCInputFromYaml(filepath);
        else
            return io::parse::text::parseKMCInput(filepath);
    }

    static KMC buildModel(const config::CommandLineConfig &config, const types::KMCInputRead &data)
    {
        SpeciesSet speciesSet = buildSpeciesSet(data.species, data.config);

        ReactionSet reactionSet = buildReactionSet(data.reactions, data.rateConstants, speciesSet);

        KMC kmc(speciesSet, reactionSet, config, data.config);

        console::debug("Built KMC object successfully.");

        return kmc;
    }
};
namespace builder
{
    static SpeciesSet buildSpeciesSet(const types::SpeciesSetRead &data, const config::SimulationConfig &config)
    {

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

        // Register and create unit species
        std::vector<Unit> units;
        size_t numUnits = data.units.size();
        units.reserve(numUnits);

        for (const auto &unitRead : unitsRead)
        {
            SpeciesID id = registry::builder.registerNewSpecies(unitRead.name, unitRead.type);
            units.push_back(Unit(unitRead.type, unitRead.name, id, unitRead.C0, unitRead.FW, unitRead.efficiency));
        }

        // Register and create polymer types
        std::vector<PolymerType> polymerTypes;
        size_t numPolyTypes = data.polymerTypes.size();
        polymerTypes.reserve(numPolyTypes);

        std::vector<PolymerContainerMap> polymerContainerMap;
        size_t numPolyContainers = data.polymerLabels.size() + numPolyTypes;
        polymerContainerMap.reserve(numPolyContainers);

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
            polymerContainerMap.push_back(PolymerContainerMap(id, polyType.name, {polymerTypes.size() - 1}));
        }

        for (const auto &label : data.polymerLabels)
        {
            std::vector<size_t> labelPolyIndices;
            labelPolyIndices.reserve(label.polymerNames.size());

            for (const auto &polyName : label.polymerNames)
            {
                if (!registry::builder.isRegistered(polyName))
                    console::input_error("Polymer " + polyName + " for label " + label.name + " is not registered. Exiting.");

                // Find index of this polymer type
                size_t index = input::findInVector(polyName, polymerTypes);
                if (index == SIZE_T_MAX)
                    console::input_error("Polymer " + polyName + " for label " + label.name + " is not registered. Exiting.");

                labelPolyIndices.push_back(index);
            }

            SpeciesID id = registry::builder.registerNewSpecies(label.name, label.type);
            polymerContainerMap.push_back(PolymerContainerMap(id, label.name, labelPolyIndices));
        }

        // Finalize registry
        registry::initialize();

        // Create species set
        SpeciesSet speciesSet(std::move(units), std::move(polymerTypes), std::move(polymerContainerMap), config.numParticles);
        return speciesSet;
    }

    static std::vector<RateConstant> buildRateConstants(const std::vector<types::RateConstantRead> &data)
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

    static ReactionSet buildReactionSet(const std::vector<types::ReactionRead> &reactionsRead, const std::vector<types::RateConstantRead> &rateConstantsRead, SpeciesSet &speciesSet)
    {
        std::vector<Reaction *> reactions;

        const std::vector<RateConstant> rateConstants = buildRateConstants(rateConstantsRead);

        const std::vector<Unit> units = speciesSet.getUnits();
        const std::vector<PolymerContainer> &polymerContainers = speciesSet.getPolymerContainers();

        for (const auto &reactionData : reactionsRead)
        {
            std::vector<Unit *> unitReactants;
            std::vector<Unit *> unitProducts;
            std::vector<PolymerContainerPtr> polyReactants;
            std::vector<PolymerContainerPtr> polyProducts;
            unitReactants.reserve(3);
            unitProducts.reserve(3);
            polyReactants.reserve(3);
            polyProducts.reserve(3);

            // Find rate constant
            size_t index = input::findInVector(reactionData.rateConstantName, rateConstants);
            if (index == SIZE_T_MAX)
                console::input_error("Rate constant " + reactionData.rateConstantName + " not found. Exiting.");

            RateConstant rateConstant = rateConstants[index];

            // Find reactants
            for (const auto &reactantName : reactionData.reactantNames)
            {
                auto species = registry::getSpecies(reactantName);

                if (SpeciesType::isUnitType(species.type))
                    unitReactants.push_back(&speciesSet.getUnits()[registry::getUnitIndex(species.ID)]);
                else if (SpeciesType::isPolymerType(species.type))
                    polyReactants.push_back(&speciesSet.getPolymerContainers()[registry::getPolymerIndex(species.ID)]);
                else
                    console::input_error("Species type " + species.type + " for species " + species.name + " not recognized. Exiting.");
            }

            for (const auto &productName : reactionData.productNames)
            {
                auto species = registry::getSpecies(productName);

                if (SpeciesType::isUnitType(species.type))
                    unitProducts.push_back(&speciesSet.getUnits()[registry::getUnitIndex(species.ID)]);
                else if (SpeciesType::isPolymerType(species.type))
                    polyProducts.push_back(&speciesSet.getPolymerContainers()[registry::getPolymerIndex(species.ID)]);
                else
                    console::input_error("Species type " + species.type + " for species " + species.name + " not recognized. Exiting.");
            }

            uint8_t sameReactant = 0;

            if (reactionData.type == Elementary::TYPE)
                reactions.push_back(new Elementary(rateConstant, unitReactants, unitProducts));
            else if (reactionData.type == InitiatorDecomposition::TYPE)
                reactions.push_back(new InitiatorDecomposition(rateConstant, unitReactants[0], unitProducts[0], unitProducts[1], unitReactants[0]->efficiency));
            else if (reactionData.type == InitiatorDecompositionPolymer::TYPE)
                reactions.push_back(new InitiatorDecompositionPolymer(rateConstant, unitReactants[0], polyProducts[0], polyProducts[1], unitReactants[0]->efficiency));
            else if (reactionData.type == Initiation::TYPE)
                reactions.push_back(new Initiation(rateConstant, unitReactants[0], unitReactants[1], polyProducts[0]));
            else if (reactionData.type == Propagation::TYPE)
                reactions.push_back(new Propagation(rateConstant, polyReactants[0], unitReactants[0], polyProducts[0]));
            else if (reactionData.type == Depropagation::TYPE)
                reactions.push_back(new Depropagation(rateConstant, polyReactants[0], polyProducts[0], unitProducts[0]));
            else if (reactionData.type == TerminationCombination::TYPE)
            {
                if (polyReactants[0]->name == polyReactants[1]->name)
                    sameReactant = 1;
                reactions.push_back(new TerminationCombination(rateConstant, polyReactants[0], polyReactants[1], polyProducts[0], sameReactant));
            }
            else if (reactionData.type == TerminationDisproportionation::TYPE)
            {
                if (polyReactants[0]->name == polyReactants[1]->name)
                    sameReactant = 1;
                reactions.push_back(new TerminationDisproportionation(rateConstant, polyReactants[0], polyReactants[1], polyProducts[0], polyProducts[1], sameReactant));
            }
            else if (reactionData.type == ChainTransferToMonomer::TYPE)
                reactions.push_back(new ChainTransferToMonomer(rateConstant, polyReactants[0], unitReactants[0], polyProducts[0], polyProducts[1]));
            else if (reactionData.type == ThermalInitiationMonomer::TYPE)
                reactions.push_back(new ThermalInitiationMonomer(rateConstant, unitReactants[0], unitReactants[1], unitReactants[2], polyProducts[0], polyProducts[1]));
            else
                console::input_error(reactionData.type + " is not a valid reaction type.");
        }

        ReactionSet reactionSet(reactions, rateConstants);

        reactionSet.printSummary();

        return reactionSet;
    };
};