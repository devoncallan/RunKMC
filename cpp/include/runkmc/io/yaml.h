#pragma once
#include <yaml-cpp/yaml.h>

#include "common.h"
#include "io/types.h"
namespace io::yaml
{
    template <typename T>
    struct Parser
    {
        static T read(const YAML::Node &node);
        static YAML::Node write(const T &obj);
    };

    template <typename T>
    struct Parser<std::vector<T>>
    {
        static std::vector<T> read(const YAML::Node &node)
        {

            if (!node.IsSequence())
                throw std::runtime_error("Expected a sequence node.");

            std::vector<T> data;
            data.reserve(node.size());
            for (const auto &item : node)
                data.push_back(Parser<T>::read(item));
            return data;
        }

        static YAML::Node write(const std::vector<T> &data)
        {
            YAML::Node node;
            for (const auto &item : data)
                node.push_back(Parser<T>::write(item));
            return node;
        }
    };

    static bool hasKey(const YAML::Node &node, const std::string_view &key)
    {
        if (node[std::string(key)])
            return true;
        return false;
    }

    template <typename T>
    static void readVar(const YAML::Node &node, const std::string_view &key, T &value, bool required = false)
    {
        auto keyStr = std::string(key);
        if (hasKey(node, key))
            value = node[keyStr].as<T>();
        else if (required)
            console::input_error("Missing required key `" + std::string(key) + "` in node `" + node.Tag() + "`.");
    }

    template <typename T>
    static void readVarRequired(const YAML::Node &node, const std::string_view &key, T &value) { readVar(node, key, value, true); }

    static YAML::Node getRequiredNode(const YAML::Node &root, const std::string_view &sectionName)
    {
        if (!hasKey(root, sectionName))
            console::input_error("YAML input file is missing required section: " + std::string(sectionName) + ".");
        return root[std::string(sectionName)];
    }

    // +--------------------------
    // | Parameters - Config
    // +--------------------------
    template <>
    struct Parser<types::SimulationConfig>
    {
        static types::SimulationConfig read(const YAML::Node &node)
        {
            types::SimulationConfig data;
            readVarRequired(node, C::io::NUM_UNITS_KEY, data.numParticles);
            readVarRequired(node, C::io::TERMINATION_TIME_KEY, data.terminationTime);
            readVarRequired(node, C::io::ANALYSIS_TIME_KEY, data.analysisTime);
            return data;
        }

        static YAML::Node write(const types::SimulationConfig &data)
        {
            YAML::Node node;
            node[C::io::NUM_UNITS_KEY] = data.numParticles;
            node[C::io::TERMINATION_TIME_KEY] = data.terminationTime;
            node[C::io::ANALYSIS_TIME_KEY] = data.analysisTime;
            return node;
        }
    };

    // +--------------------------
    // | Species - Base
    // +--------------------------
    template <>
    struct Parser<types::SpeciesRead>
    {
        static types::SpeciesRead read(const YAML::Node &node)
        {
            types::SpeciesRead data;
            readVarRequired(node, C::io::NAME_KEY, data.name);
            readVarRequired(node, C::io::TYPE_KEY, data.type);

            SpeciesType::checkValid(data.type);

            return data;
        }

        static YAML::Node write(const types::SpeciesRead &data)
        {
            YAML::Node node;
            node[C::io::NAME_KEY] = data.name;
            node[C::io::TYPE_KEY] = data.type;
            return node;
        }
    };

    // +--------------------------
    // | Species - Registered
    // +--------------------------
    template <>
    struct Parser<RegisteredSpecies>
    {
        static RegisteredSpecies read(const YAML::Node &node)
        {
            RegisteredSpecies data;
            readVarRequired(node, C::io::NAME_KEY, data.name);
            readVarRequired(node, C::io::TYPE_KEY, data.type);
            readVarRequired(node, C::io::ID_KEY, data.ID);

            SpeciesType::checkValid(data.type);

            return data;
        }

        static YAML::Node write(const RegisteredSpecies &data)
        {
            YAML::Node node;
            node[C::io::NAME_KEY] = data.name;
            node[C::io::TYPE_KEY] = data.type;
            node[C::io::ID_KEY] = int(data.ID);
            return node;
        }
    };

    // +--------------------------
    // | Species - Unit
    // +--------------------------
    template <>
    struct Parser<types::UnitRead> : public Parser<types::SpeciesRead>
    {
        static types::UnitRead read(const YAML::Node &node)
        {
            types::SpeciesRead species = Parser<types::SpeciesRead>::read(node);
            types::UnitRead data;
            data.name = species.name;
            data.type = species.type;
            readVar(node, C::io::C0_KEY, data.C0);
            readVar(node, C::io::FW_KEY, data.FW);
            readVar(node, C::io::EFFICIENCY_KEY, data.efficiency);
            return data;
        }

        static YAML::Node write(const types::UnitRead &data)
        {
            YAML::Node node = Parser<types::SpeciesRead>::write(data);
            node[C::io::C0_KEY] = data.C0;
            node[C::io::FW_KEY] = data.FW;
            node[C::io::EFFICIENCY_KEY] = data.efficiency;
            return node;
        }
    };

    // +--------------------------
    // | Species - Polymer Type
    // +--------------------------
    template <>
    struct Parser<types::PolymerTypeRead> : public Parser<types::SpeciesRead>
    {
        static types::PolymerTypeRead read(const YAML::Node &node)
        {
            types::SpeciesRead species = Parser<types::SpeciesRead>::read(node);
            types::PolymerTypeRead data;
            data.name = species.name;
            data.type = species.type;
            readVarRequired(node, C::io::END_GROUP_NAMES_KEY, data.endGroupUnitNames);
            return data;
        }

        static YAML::Node write(const types::PolymerTypeRead &data)
        {
            YAML::Node node = Parser<types::SpeciesRead>::write(data);
            node[C::io::END_GROUP_NAMES_KEY] = data.endGroupUnitNames;
            return node;
        }
    };

    // +--------------------------
    // | Species - Polymer Labels
    // +--------------------------
    template <>
    struct Parser<types::PolymerLabelsRead> : public Parser<types::SpeciesRead>
    {
        static types::PolymerLabelsRead read(const YAML::Node &node)
        {
            types::SpeciesRead species = Parser<types::SpeciesRead>::read(node);
            types::PolymerLabelsRead data;
            data.name = species.name;
            data.type = species.type;
            readVarRequired(node, C::io::POLYMER_NAMES_KEY, data.polymerNames);
            return data;
        }

        static YAML::Node write(const types::PolymerLabelsRead &data)
        {
            YAML::Node node = Parser<types::SpeciesRead>::write(data);
            node[C::io::POLYMER_NAMES_KEY] = data.polymerNames;
            return node;
        }
    };

    // +--------------------------
    // | Species - Section
    // +--------------------------
    template <>
    struct Parser<types::SpeciesSetRead>
    {
        static types::SpeciesSetRead read(const YAML::Node &node)
        {
            // ===== PASS 1: Separate nodes by species type =====
            std::vector<YAML::Node> unitNodes;
            std::vector<YAML::Node> polymerNodes;
            std::vector<YAML::Node> labelNodes;

            for (const auto &specNode : node)
            {
                types::SpeciesRead species = Parser<types::SpeciesRead>::read(specNode);

                if (SpeciesType::isUnitType(species.type))
                    unitNodes.push_back(specNode);
                else if (species.type == SpeciesType::POLYMER)
                    polymerNodes.push_back(specNode);
                else if (species.type == SpeciesType::LABEL)
                    labelNodes.push_back(specNode);

                SpeciesType::checkValid(species.type);
            }

            // ===== PASS 2: Parse in order: Units, Polymers, Labels =====
            std::vector<types::UnitRead> units;
            std::vector<types::PolymerTypeRead> polymerTypes;
            std::vector<types::PolymerLabelsRead> labels;
            units.reserve(unitNodes.size());
            polymerTypes.reserve(polymerNodes.size());
            labels.reserve(labelNodes.size());

            for (const auto &specNode : unitNodes)
                units.push_back(Parser<types::UnitRead>::read(specNode));

            for (const auto &specNode : polymerNodes)
                polymerTypes.push_back(Parser<types::PolymerTypeRead>::read(specNode));

            for (const auto &specNode : labelNodes)
                labels.push_back(Parser<types::PolymerLabelsRead>::read(specNode));

            types::SpeciesSetRead speciesSet;
            speciesSet.units = std::move(units);
            speciesSet.polymerTypes = std::move(polymerTypes);
            speciesSet.polymerLabels = std::move(labels);

            return speciesSet;
        }

        static YAML::Node write(const types::SpeciesSetRead &data)
        {
            YAML::Node node;
            for (const auto &unit : data.units)
                node.push_back(Parser<types::UnitRead>::write(unit));
            for (const auto &polymer : data.polymerTypes)
                node.push_back(Parser<types::PolymerTypeRead>::write(polymer));
            for (const auto &label : data.polymerLabels)
                node.push_back(Parser<types::PolymerLabelsRead>::write(label));
            return node;
        }
    };

    // +--------------------------
    // | Reactions - Rate Constant
    // +--------------------------

    template <>
    struct Parser<types::RateConstantRead>
    {
        static types::RateConstantRead read(const YAML::Node &node)
        {
            types::RateConstantRead data;
            readVarRequired(node, C::io::NAME_KEY, data.name);
            readVarRequired(node, C::io::VALUE_KEY, data.k);
            return data;
        }

        static YAML::Node write(const types::RateConstantRead &data)
        {
            YAML::Node node;
            node[C::io::NAME_KEY] = data.name;
            node[C::io::VALUE_KEY] = data.k;
            return node;
        }
    };

    // +--------------------------
    // | Reactions - Base
    // +--------------------------

    template <>
    struct Parser<types::ReactionRead>
    {
        static types::ReactionRead read(const YAML::Node &node)
        {
            types::ReactionRead data;
            readVarRequired(node, C::io::TYPE_KEY, data.type);
            readVarRequired(node, C::io::RATE_CONSTANT_KEY, data.rateConstantName);
            readVarRequired(node, C::io::REACTANTS_KEY, data.reactantNames);
            readVarRequired(node, C::io::PRODUCTS_KEY, data.productNames);
            return data;
        }

        static YAML::Node write(const types::ReactionRead &data)
        {
            YAML::Node node;
            node[C::io::TYPE_KEY] = data.type;
            node[C::io::RATE_CONSTANT_KEY] = data.rateConstantName;
            node[C::io::REACTANTS_KEY] = data.reactantNames;
            node[C::io::PRODUCTS_KEY] = data.productNames;
            return node;
        }
    };

    template <>
    struct Parser<types::KMCInputRead>
    {
        static types::KMCInputRead read(const YAML::Node &node)
        {
            auto parameters = getRequiredNode(node, C::io::PARAMETERS_SECTION);
            auto species = getRequiredNode(node, C::io::SPECIES_SECTION);
            auto rateConstants = getRequiredNode(node, C::io::RATE_CONSTANTS_SECTION);
            auto reactions = getRequiredNode(node, C::io::REACTIONS_SECTION);

            types::KMCInputRead data;
            data.config = Parser<types::SimulationConfig>::read(parameters);
            data.species = Parser<types::SpeciesSetRead>::read(species);
            data.rateConstants = Parser<std::vector<types::RateConstantRead>>::read(rateConstants);
            data.reactions = Parser<std::vector<types::ReactionRead>>::read(reactions);

            return data;
        }

        static YAML::Node write(const types::KMCInputRead &data)
        {
            YAML::Node node;
            node[C::io::PARAMETERS_SECTION] = Parser<types::SimulationConfig>::write(data.config);
            node[C::io::SPECIES_SECTION] = Parser<types::SpeciesSetRead>::write(data.species);
            node[C::io::RATE_CONSTANTS_SECTION] = Parser<std::vector<types::RateConstantRead>>::write(data.rateConstants);
            node[C::io::REACTIONS_SECTION] = Parser<std::vector<types::ReactionRead>>::write(data.reactions);
            return node;
        }
    };

    // +--------------------------
    // | File I/O
    // +--------------------------
    static void writeYamlToFile(const std::filesystem::path &filepath, const YAML::Node &node)
    {
        std::ofstream file(filepath);

        if (!file.is_open())
        {
            throw std::runtime_error("Failed to open file for writing: " + filepath.string());
        }

        file << std::setprecision(C::io::PRECISION);
        file << node;
    }
};
namespace io
{
    static types::KMCInputRead parseYamlModelFile(const std::string &filepath)
    {
        YAML::Node root = YAML::LoadFile(filepath);
        types::KMCInputRead input = io::yaml::Parser<types::KMCInputRead>::read(root);
        return input;
    }
}
