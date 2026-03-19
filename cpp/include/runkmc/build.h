#pragma once

#include <unordered_set>

#include "common.h"
#include "io/text.h"
#include "io/yaml.h"
#include "io/cli.h"
#include "kmc/kmc.h"

namespace build
{
    inline constexpr console::LogContext logger("KMC::Build");

    // Forward declarations of builder functions
    static SpeciesSet buildSpeciesSet(const io::types::SpeciesSetRead &data, const io::types::SimulationConfig &config);
    static std::vector<RateConstant> buildRateConstants(const std::vector<io::types::RateConstantRead> &data);
    static ReactionSet buildReactionSet(const std::vector<io::types::ReactionRead> &reactionsRead, const std::vector<io::types::RateConstantRead> &rateConstantsRead, SpeciesSet &speciesSet);

    static io::types::CommandLineConfig parseArguments(int arc, char **argv) { return io::cli::parseArguments(arc, argv); }

    static io::types::KMCInputRead parseModelFile(const std::string &filepath)
    {
        if (str::endswith(filepath, ".yaml") || str::endswith(filepath, ".yml"))
            return io::parseYamlModelFile(filepath);
        else if (str::endswith(filepath, ".txt"))
            return io::parseTextModelFile(filepath);

        console::input_error("Unrecognized file extension for model file: " + filepath + ".");
    }

    static KMC buildModel(const io::types::CommandLineConfig &config, const io::types::KMCInputRead &data)
    {
        SpeciesSet speciesSet = buildSpeciesSet(data.species, data.config);

        ReactionSet reactionSet = buildReactionSet(data.reactions, data.rateConstants, speciesSet);

        // Sink detection: collect all SpeciesIDs that appear as reactants
        std::unordered_set<SpeciesID> reactantIDs;
        for (const auto &reactionRead : data.reactions)
            for (const auto &name : reactionRead.reactantNames)
                if (registry::isRegistered(name))
                    reactantIDs.insert(registry::getSpecies(name).ID);

        // Mark polymer containers as sinks if they never appear as reactants
        for (auto &container : speciesSet.getPolymerContainers())
            if (reactantIDs.find(container.ID) == reactantIDs.end())
                container.isSink = true;

        KMC kmc(speciesSet, reactionSet, config, data.config);

        console::debug("Built KMC object successfully.");

        return kmc;
    }
};
namespace build
{
    static SpeciesSet buildSpeciesSet(const io::types::SpeciesSetRead &data, const io::types::SimulationConfig &config)
    {

        std::vector<io::types::UnitRead> unitsRead = data.units;
        std::stable_sort(
            unitsRead.begin(),
            unitsRead.end(),
            [](const io::types::UnitRead &a, const io::types::UnitRead &b)
            {
            auto priority = [](std::string_view type) {
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
            units.push_back(Unit(id, unitRead.name, unitRead.type, unitRead.C0, unitRead.FW, unitRead.efficiency));
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
                    logger.error("End group unit " + unitName + " for polymer " + polyType.name + " is not registered. Exiting.");
                // console::input_error("End group unit " + unitName + " for polymer " + polyType.name + " is not registered. Exiting.");
                endGroupIDs.push_back(registry::builder.getSpeciesID(unitName));
            }

            SpeciesID id = registry::builder.registerNewSpecies(polyType.name, polyType.type);
            PolymerType type(id, polyType.name, endGroupIDs);
            polymerTypes.push_back(type);
            polymerContainerMap.push_back(PolymerContainerMap(id, polyType.name, {polymerTypes.size() - 1}, polyType.report));
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
                size_t index = str::findInVector(polyName, polymerTypes);
                if (index == SIZE_MAX)
                    console::input_error("Polymer " + polyName + " for label " + label.name + " is not registered. Exiting.");

                labelPolyIndices.push_back(index);
            }

            SpeciesID id = registry::builder.registerNewSpecies(label.name, label.type);
            polymerContainerMap.push_back(PolymerContainerMap(id, label.name, labelPolyIndices, label.report));
        }

        // Finalize registry
        registry::initialize();

        auto chainType = ChainType::Sequence;
        if (registry::getNumMonomers() == 1)
            chainType = ChainType::Homopolymer;
        for (auto &polyType : polymerTypes)
            polyType.setChainType(chainType);

        // Create species set
        SpeciesSet speciesSet(std::move(units), std::move(polymerTypes), std::move(polymerContainerMap), config.numParticles);
        return speciesSet;
    }

    static std::vector<RateConstant> buildRateConstants(const std::vector<io::types::RateConstantRead> &data)
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

    static ReactionSpecies buildReactionSpecies(const io::types::ReactionRead &reactionData, SpeciesSet &speciesSet)
    {
        ReactionSpecies species;
        species.reactants.reserve(reactionData.reactantNames.size());
        species.products.reserve(reactionData.productNames.size());

        for (const auto &reactantName : reactionData.reactantNames)
        {
            if (!registry::isRegistered(reactantName))
                console::input_error("Reactant species " + reactantName + " not registered. Exiting.");

            auto speciesInfo = registry::getSpecies(reactantName);
            if (SpeciesType::isUnitType(speciesInfo.type))
                species.reactants.push_back(&speciesSet.getUnits()[registry::getUnitIndex(speciesInfo.ID)]);
            else if (SpeciesType::isPolymerType(speciesInfo.type))
                species.reactants.push_back(&speciesSet.getPolymerContainers()[registry::getPolymerIndex(speciesInfo.ID)]);
            else
                console::input_error("Species type " + std::string(speciesInfo.type) + " for species " + speciesInfo.name + " not recognized. Exiting.");
        }

        for (const auto &productName : reactionData.productNames)
        {
            if (!registry::isRegistered(productName))
                console::input_error("Product species " + productName + " not registered. Exiting.");

            auto speciesInfo = registry::getSpecies(productName);
            if (SpeciesType::isUnitType(speciesInfo.type))
                species.products.push_back(&speciesSet.getUnits()[registry::getUnitIndex(speciesInfo.ID)]);
            else if (SpeciesType::isPolymerType(speciesInfo.type))
                species.products.push_back(&speciesSet.getPolymerContainers()[registry::getPolymerIndex(speciesInfo.ID)]);
            else
                console::input_error("Species type " + std::string(speciesInfo.type) + " for species " + speciesInfo.name + " not recognized. Exiting.");
        }

        return species;
    }

    static ReactionSet buildReactionSet(const std::vector<io::types::ReactionRead> &reactionsRead, const std::vector<io::types::RateConstantRead> &rateConstantsRead, SpeciesSet &speciesSet)
    {
        std::vector<Reaction *> reactions;

        const std::vector<RateConstant> rateConstants = buildRateConstants(rateConstantsRead);

        for (const auto &reactionData : reactionsRead)
        {
            // 1. Find rate constant
            size_t index = str::findInVector(reactionData.rateConstantName, rateConstants);
            if (index == SIZE_MAX)
                console::input_error("Rate constant " + reactionData.rateConstantName + " not found. Exiting.");
            RateConstant rateConstant = rateConstants[index];

            // 2. Build ReactionSpecies
            ReactionSpecies species = buildReactionSpecies(reactionData, speciesSet);

            // 3. Create new Reaction based on type
            if (reactionData.type == ReactionType::ELEMENTARY)
                reactions.push_back(new Elementary(rateConstant, species));
            else if (reactionData.type == ReactionType::INITIATOR_DECOMPOSITION)
                reactions.push_back(new InitiatorDecomposition(rateConstant, species));
            else if (reactionData.type == ReactionType::INIT_DECOMP_POLY)
                reactions.push_back(new InitiatorDecompositionPolymer(rateConstant, species));
            else if (reactionData.type == ReactionType::INITIATION)
                reactions.push_back(new Initiation(rateConstant, species));
            else if (reactionData.type == ReactionType::PROPAGATION)
                reactions.push_back(new Propagation(rateConstant, species));
            else if (reactionData.type == ReactionType::DEPROPAGATION)
                reactions.push_back(new Depropagation(rateConstant, species));
            else if (reactionData.type == ReactionType::TERMINATION_D)
                reactions.push_back(new TerminationDisproportionation(rateConstant, species));
            else if (reactionData.type == ReactionType::TERMINATION_C)
                reactions.push_back(new TerminationCombination(rateConstant, species));
            else if (reactionData.type == ReactionType::CHAINTRANSFER_M)
                reactions.push_back(new ChainTransferToMonomer(rateConstant, species));
            else if (reactionData.type == ReactionType::THERM_INIT_M)
                reactions.push_back(new ThermalInitiationMonomer(rateConstant, species));
            else
                console::input_error("Reaction type " + std::string(reactionData.type) + " not recognized. Exiting.");
        }

        ReactionSet reactionSet(reactions, rateConstants);
        reactionSet.printSummary();
        return reactionSet;
    };
};
