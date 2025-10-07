#pragma once
#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>

#include "common.h"
#include "kmc/kmc.h"
#include "utils/yaml.h"
#include "types.h"
#include "core/C.h"

// namespace io::parsers
// {
//     class yaml
//     {
//     public:
//         types::SpeciesRead parseSpecies(const YAML::Node &node)
//         {
//             types::SpeciesRead species;
//             readVar(node, "name", species.name, true);
//             readVar(node, "type", species.type, true);
//             return species;
//         };

//         types::UnitRead parseUnit(const YAML::Node &node)
//         {
//             types::SpeciesRead species = parseSpecies(node);

//             types::UnitRead unit;
//             unit.name = species.name;
//             unit.type = species.type;

//             readVar(node, C::C0_VAR, unit.C0);
//             readVar(node, C::FW_VAR, unit.FW);
//             readVar(node, C::EFFICIENCY_VAR, unit.efficiency); // optional, default 1

//             return unit;
//         }

//         types::PolymerTypeRead parsePolymerType(const YAML::Node &node)
//         {
//             types::SpeciesRead species = parseSpecies(node);
//             types::PolymerTypeRead polymer;
//             polymer.name = species.name;
//             polymer.type = species.type;

//             readVar(node, "end_group_units", polymer.endGroupUnitNames, true); // required

//             return polymer;
//         }

//         types::PolymerLabelsRead parsePolymerLabels(const YAML::Node &node)
//         {
//             types::SpeciesRead species = parseSpecies(node);

//             types::PolymerLabelsRead labels;
//             labels.name = species.name;
//             labels.type = species.type;

//             readVar(node, "polymer_names", labels.polymerNames, true); // required

//             return labels;
//         }

//     private:
//         template <typename T>
//         static void readVar(const YAML::Node &node, const std::string &varName, T &var, bool required = false)
//         {
//             if (node[varName])
//                 var = node[varName].as<T>();
//             else if (required)
//                 throw std::runtime_error("Missing required variable: " + varName);
//         }
//     };
// };

namespace io::parse::yaml
{

    struct SectionNodes
    {
        YAML::Node parameters;
        YAML::Node species;
        YAML::Node rateConstants;
        YAML::Node reactions;
    };

    static SectionNodes parseFile(const std::string &filepath)
    {
        YAML::Node root = YAML::LoadFile(filepath);

        if (!root["parameters"] || !root["species"] || !root["rate_constants"] || !root["reactions"])
            console::input_error("YAML input file is missing one or more required sections: 'parameters', 'species', 'rate_constants', 'reactions'.");

        SectionNodes sections;
        sections.parameters = root["parameters"];
        sections.species = root["species"];
        sections.rateConstants = root["rate_constants"];
        sections.reactions = root["reactions"];

        return sections;
    }

    // Parameters
    static config::SimulationConfig parseSimulationConfig(const YAML::Node &node)
    {
        config::SimulationConfig config;
        utils::readVar(node, C::io::NUM_UNITS_KEY, config.numParticles, true);
        utils::readVar(node, C::io::TERMINATION_TIME_KEY, config.terminationTime, true);
        utils::readVar(node, C::io::ANALYSIS_TIME_KEY, config.analysisTime, true);
        return config;
    }

    // Species
    static types::SpeciesRead parseBaseSpecies(const YAML::Node &node)
    {
        types::SpeciesRead species;
        utils::readVar(node, C::io::NAME_KEY, species.name, true);
        utils::readVar(node, C::io::TYPE_KEY, species.type, true);

        if (!SpeciesType::isValidType(species.type))
            console::input_error("Species type " + species.type + " is not valid.");

        return species;
    }

    static types::UnitRead parseUnit(const YAML::Node &node)
    {
        types::SpeciesRead species = parseBaseSpecies(node);

        types::UnitRead unit;
        unit.name = species.name;
        unit.type = species.type;

        utils::readVar(node, C::io::C0_KEY, unit.C0);
        utils::readVar(node, C::io::FW_KEY, unit.FW);

        if (unit.type == SpeciesType::INITIATOR)
            utils::readVar(node, C::io::EFFICIENCY_KEY, unit.efficiency, true); // required for initiators
        else
            utils::readVar(node, C::io::EFFICIENCY_KEY, unit.efficiency);

        return unit;
    }
    static types::PolymerTypeRead parsePolymerType(const YAML::Node &node)
    {
        types::SpeciesRead species = parseBaseSpecies(node);
        types::PolymerTypeRead polymer;
        polymer.name = species.name;
        polymer.type = species.type;

        utils::readVar(node, C::io::END_GROUP_NAMES_KEY, polymer.endGroupUnitNames, true); // required

        return polymer;
    }

    static types::PolymerLabelsRead parsePolymerLabels(const YAML::Node &node)
    {
        types::SpeciesRead species = parseBaseSpecies(node);

        types::PolymerLabelsRead labels;
        labels.name = species.name;
        labels.type = species.type;

        utils::readVar(node, C::io::POLYMER_NAMES_KEY, labels.polymerNames, true); // required

        return labels;
    }

    static types::SpeciesSetRead parseSpecies(const YAML::Node &node)
    {
        // ===== PASS 1: Separate nodes by species type =====
        std::vector<YAML::Node> unitNodes;
        std::vector<YAML::Node> polymerNodes;
        std::vector<YAML::Node> labelNodes;

        for (const auto &specNode : node)
        {
            types::SpeciesRead species = parseBaseSpecies(specNode);
            if (SpeciesType::isUnitType(species.type))
                unitNodes.push_back(specNode);
            else if (species.type == SpeciesType::POLYMER)
                polymerNodes.push_back(specNode);
            else if (species.type == SpeciesType::LABEL)
                labelNodes.push_back(specNode);
            else
                console::input_error("Unknown species type: " + species.type);
        }

        // ===== PASS 2: Parse in order: Units, Polymers, Labels =====
        std::vector<types::UnitRead> units;
        std::vector<types::PolymerTypeRead> polymerTypes;
        std::vector<types::PolymerLabelsRead> labels;
        units.reserve(unitNodes.size());
        polymerTypes.reserve(polymerNodes.size());
        labels.reserve(labelNodes.size());

        for (const auto &specNode : unitNodes)
            units.push_back(parseUnit(specNode));

        for (const auto &specNode : polymerNodes)
            polymerTypes.push_back(parsePolymerType(specNode));

        for (const auto &specNode : labelNodes)
            labels.push_back(parsePolymerLabels(specNode));

        types::SpeciesSetRead speciesSet;
        speciesSet.units = std::move(units);
        speciesSet.polymerTypes = std::move(polymerTypes);
        speciesSet.polymerLabels = std::move(labels);

        return speciesSet;
    }

    // Rate constants
    static std::vector<types::RateConstantRead> parseRateConstants(const YAML::Node &node)
    {
        std::vector<types::RateConstantRead> rateConstants;
        rateConstants.reserve(node.size());
        for (const auto &rcNode : node)
        {
            types::RateConstantRead rateConstant;
            utils::readVar(rcNode, C::io::NAME_KEY, rateConstant.name, true);
            utils::readVar(rcNode, C::io::VALUE_KEY, rateConstant.k, true);
            rateConstants.push_back(rateConstant);
        }
        return rateConstants;
    }

    // Reactions
    static types::ReactionRead parseReaction(const YAML::Node &node)
    {
        types::ReactionRead reaction;
        utils::readVar(node, C::io::TYPE_KEY, reaction.type, true);
        utils::readVar(node, C::io::RATE_CONSTANT_KEY, reaction.rateConstantName, true);
        utils::readVar(node, C::io::REACTANTS_KEY, reaction.reactantNames, true);
        utils::readVar(node, C::io::PRODUCTS_KEY, reaction.productNames, true);
        return reaction;
    }
    static std::vector<types::ReactionRead> parseReactions(const YAML::Node &node)
    {
        std::vector<types::ReactionRead> reactions;
        reactions.reserve(node.size());
        for (const auto &rxnNode : node)
            reactions.push_back(parseReaction(rxnNode));

        return reactions;
    }
};

namespace io::parse::yaml::utils
{
    template <typename T>
    static void readVar(const YAML::Node &node, const std::string &varName, T &var, bool required = false)
    {
        if (node[varName])
            var = node[varName].as<T>();
        else if (required)
            throw std::runtime_error("Missing required variable: " + varName);
    }
};