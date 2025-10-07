#pragma once
#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>

#include "common.h"
#include "kmc/kmc.h"
#include "utils/yaml.h"
#include "types.h"
#include "core/C.h"

namespace io::parse::text
{
    struct SectionLines
    {
        std::vector<std::string> parameters;
        std::vector<std::string> species;
        std::vector<std::string> rateConstants;
        std::vector<std::string> reactions;
    };

    static SectionLines parseFile(const std::string &filepath)
    {
        std::ifstream modelFile(filepath);
        if (!modelFile.is_open())
            console::input_error("Cannot open model file: " + filepath);

        SectionLines sections;
        std::string line;
        while (std::getline(modelFile, line))
        {
            str::trim(line);
            if (input::canIgnoreLine(line))
                continue;

            if (sections.parameters.empty() && str::startswith(line, C::io::PARAMETERS_SECTION))
                sections.parameters = input::parseSection(C::io::PARAMETERS_SECTION, modelFile);

            else if (sections.species.empty() && str::startswith(line, C::io::SPECIES_SECTION))
                sections.species = input::parseSection(C::io::SPECIES_SECTION, modelFile);

            else if (sections.rateConstants.empty() && str::startswith(line, C::io::RATE_CONSTANTS_SECTION))
                sections.rateConstants = input::parseSection(C::io::RATE_CONSTANTS_SECTION, modelFile);

            else if (sections.reactions.empty() && str::startswith(line, C::io::REACTIONS_SECTION))
                sections.reactions = input::parseSection(C::io::REACTIONS_SECTION, modelFile);
        }

        if (sections.parameters.empty() || sections.species.empty() ||
            sections.rateConstants.empty() || sections.reactions.empty())
            console::input_error("Missing required sections in model file. Required sections:"); //+ str::join(C::REQUIRED_SECTIONS));

        return sections;
    }

    // Parameters
    static config::SimulationConfig parseSimulationConfig(const std::vector<std::string> &paramLines)
    {
        config::SimulationConfig config;
        utils::readVar(paramLines, C::io::NUM_UNITS_KEY, config.numParticles, true);
        utils::readVar(paramLines, C::io::TERMINATION_TIME_KEY, config.terminationTime, true);
        utils::readVar(paramLines, C::io::ANALYSIS_TIME_KEY, config.analysisTime, true);
        return config;
        // + more when I think of them
    }

    // Species
    static types::SpeciesRead parseBaseSpecies(const std::vector<std::string> &args)
    {
        if (args.size() < 2)
            console::input_error("Species definition requires at least 2 arguments: <type> <name>");

        types::SpeciesRead species;
        species.type = args[0];
        species.name = args[1];

        if (!SpeciesType::isValidType(species.type))
            console::input_error("Species type " + species.type + " is not valid.");

        return species;
    }

    static types::UnitRead parseUnit(const std::vector<std::string> &args)
    {
        types::SpeciesRead species = parseBaseSpecies(args);

        types::UnitRead unit;
        unit.name = species.name;
        unit.type = species.type;

        std::vector<std::string> vars = {args.begin() + 2, args.end()};

        utils::readVar(vars, C::io::C0_KEY, unit.C0);
        utils::readVar(vars, C::io::FW_KEY, unit.FW);

        if (unit.type == SpeciesType::INITIATOR)
            utils::readVar(vars, C::io::EFFICIENCY_KEY, unit.efficiency, true); // required for initiators
        else
            utils::readVar(vars, C::io::EFFICIENCY_KEY, unit.efficiency);

        return unit;
    }

    static types::PolymerTypeRead parsePolymerType(const std::vector<std::string> &args)
    {
        types::SpeciesRead species = parseBaseSpecies(args);

        types::PolymerTypeRead polymer;
        polymer.name = species.name;
        polymer.type = species.type;

        // Check if the polymer name has brackets (e.g., "P[A.A]")
        size_t seqStart = species.name.find("[");
        size_t seqEnd = species.name.find("]");
        bool hasSequence = (seqStart != std::string::npos && seqEnd != std::string::npos && seqEnd > seqStart);
        std::vector<std::string> endGroupUnitNames;

        if (!hasSequence)
        {
            polymer.endGroupUnitNames = {};
            return polymer;
        }

        // Get string inside the brackets and split by period
        // e.g., "P[A.A]" -> "A.A" -> {"A","A"}
        std::string sequenceStr = species.name.substr(seqStart + 1, seqEnd - seqStart - 1);
        std::vector<std::string> endGroupNames = str::splitByDelimeter(sequenceStr, ".");
        polymer.endGroupUnitNames = endGroupNames;
        return polymer;
    }

    static types::PolymerLabelsRead parsePolymerLabels(const std::vector<std::string> &args)
    {
        types::SpeciesRead species = parseBaseSpecies(args);

        types::PolymerLabelsRead labels;
        labels.name = species.name;
        labels.type = species.type;

        std::vector<std::string> vars = {args.begin() + 2, args.end()};
        if (vars.size() != 1)
            console::input_error("Polymer label definition requires polymer names separated by '|', e.g., LABEL MyLabel P[A.A]|P[B.B]");

        std::vector<std::string> polymerNames = str::splitByDelimeter(vars[0], "|");
        labels.polymerNames = polymerNames;

        return labels;
    }

    static types::SpeciesSetRead parseSpecies(const std::vector<std::string> &lines)
    {
        // ===== PASS 1: Separate lines by species type =====
        std::vector<std::string> unitLines;
        std::vector<std::string> polymerLines;
        std::vector<std::string> labelLines;

        for (const auto &line : lines)
        {
            std::vector<std::string> args = str::splitByWhitespace(line);
            if (args.empty())
                continue;

            types::SpeciesRead species = parseBaseSpecies(args);
            if (SpeciesType::isUnitType(species.type))
                unitLines.push_back(line);
            else if (species.type == SpeciesType::POLYMER)
                polymerLines.push_back(line);
            else if (species.type == SpeciesType::LABEL)
                labelLines.push_back(line);
            else
                console::input_error("Unknown species type: " + species.type);
        }

        // ===== PASS 2: Parse in order: Units, Polymers, Labels =====
        std::vector<types::UnitRead> units;
        std::vector<types::PolymerTypeRead> polymerTypes;
        std::vector<types::PolymerLabelsRead> labels;
        units.reserve(unitLines.size());
        polymerTypes.reserve(polymerLines.size());
        labels.reserve(labelLines.size());

        for (const auto &line : unitLines)
            units.push_back(parseUnit(str::splitByWhitespace(line)));

        for (const auto &line : polymerLines)
            polymerTypes.push_back(parsePolymerType(str::splitByWhitespace(line)));

        for (const auto &line : labelLines)
            labels.push_back(parsePolymerLabels(str::splitByWhitespace(line)));

        types::SpeciesSetRead speciesSet;
        speciesSet.units = std::move(units);
        speciesSet.polymerTypes = std::move(polymerTypes);
        speciesSet.polymerLabels = std::move(labels);

        return speciesSet;
    }

    // Rate constants
    static std::vector<types::RateConstantRead> parseRateConstants(const std::vector<std::string> &lines)
    {
        std::vector<types::RateConstantRead> rateConstants;
        rateConstants.reserve(lines.size());
        for (const auto &line : lines)
        {
            types::RateConstantRead rateConstant;
            auto var = utils::parseVariable<double>(line);
            rateConstant.name = var.name;
            rateConstant.k = var.value;
            rateConstants.push_back(rateConstant);
        }
        return rateConstants;
    }

    // Reactions
    static types::ReactionRead parseReaction(const std::string &line)
    {
        std::vector<std::string> args = str::splitByWhitespace(line);

        types::ReactionRead reaction;
        reaction.type = args[0];

        if (args.size() < 4)
            console::input_error("Reaction definition requires at least 4 arguments: <type> <reactants> -<rateConstantName>-> <products>. Provided: " + line);

        // Loop through args (ignoring first arg)
        bool isReactants = true;
        for (size_t i = 1; i < args.size(); ++i)
        {
            std::string reactionArg = args[i];

            // Check if string is "-name->"
            if (str::startswith(reactionArg, "-") && str::endswith(reactionArg, "->"))
            {
                std::string rateConstantName = reactionArg.substr(1, reactionArg.size() - 3); // remove "-" and "->"
                reaction.rateConstantName = rateConstantName;
                isReactants = false;
                continue;
            }

            if (reactionArg == "+")
                continue;

            if (isReactants)
                reaction.reactantNames.push_back(reactionArg);
            else
                reaction.productNames.push_back(reactionArg);
        }

        return reaction;
    }

    static std::vector<types::ReactionRead> parseReactions(const std::vector<std::string> &lines)
    {
        std::vector<types::ReactionRead> reactions;
        reactions.reserve(lines.size());
        for (const auto &line : lines)
            reactions.push_back(parseReaction(line));

        return reactions;
    }

};

namespace io::parse::text::utils
{
    template <typename T>
    static void readVar(const std::vector<std::string> &vars, const std::string &varName, T &value, bool required = false)
    {
        bool found = false;
        for (const auto &varString : vars)
        {
            if (str::startsWith(varString, varName))
            {
                types::Variable<T> var = parseVariable<T>(varString);
                value = var.value;
                found = true;
                break;
            }
        }
        if (required && !found)
            console::input_error("Required variable '" + varName + "' not found in text input.");
    }

    template <typename T>
    static types::Variable<T> parseVariable(const std::string &s)
    {
        std::vector<std::string> tokens = str::splitByDelimeter(s, "=");
        if (tokens.size() != 2)
            console::input_error("Variable '" + s + "' is not in 'name=value' format.");

        str::trim(tokens[0]);
        str::trim(tokens[1]);

        return types::Variable<T>{tokens[0], input::convertValue<T>(tokens[1])};
    }
};

// namespace io::parsers
// {
//     struct SectionLines
//     {
//         std::vector<std::string> parameters;
//         std::vector<std::string> species;
//         std::vector<std::string> rateConstants;
//         std::vector<std::string> reactions;
//     };

//     class text
//     {
//     public:
//         static SectionLines parseFile(const std::string &filepath)
//         {
//             std::ifstream modelFile(filepath);
//             if (!modelFile.is_open())
//                 console::input_error("Cannot open model file: " + filepath);

//             SectionLines sections;
//             std::string line;
//             while (std::getline(modelFile, line))
//             {
//                 str::trim(line);
//                 if (input::canIgnoreLine(line))
//                     continue;
//                 if (sections.parameters.empty() && str::startswith(line, C::PARAMETERS_SECTION))
//                 {
//                     sections.parameters = input::parseSection("parameters", modelFile);
//                 }
//                 else if (sections.species.empty() && str::startswith(line, "species"))
//                 {
//                     sections.species = input::parseSection("species", modelFile);
//                     // sections.species = sortSpeciesLinesByType(sections.species);
//                 }
//                 else if (sections.rateConstants.empty() && str::startswith(line, "rateconstants"))
//                 {
//                     sections.rateConstants = input::parseSection("rateconstants", modelFile);
//                 }
//                 else if (sections.reactions.empty() && str::startswith(line, "reactions"))
//                 {
//                     sections.reactions = input::parseSection("reactions", modelFile);
//                 }
//             }
//             return sections;
//         }

//         static config::SimulationConfig parseSimulationConfig(const std::vector<std::string> &lines)
//         {
//             config::SimulationConfig config;
//             readVar(lines, "num_units", config.numParticles, true);
//             readVar(lines, "termination_time", config.terminationTime, true);
//             readVar(lines, "analysis_time", config.analysisTime, true);
//             return config;
//             // + more when I think of them
//         }

//         static std::vector<types::RateConstantRead> parseRateConstants(const std::vector<std::string> &lines)
//         {
//             std::vector<types::RateConstantRead> rateConstants;
//             rateConstants.reserve(lines.size());
//             for (const auto &line : lines)
//             {
//                 types::RateConstantRead rateConstant;
//                 auto var = parseVariable<double>(line);
//                 rateConstant.name = var.name;
//                 rateConstant.k = var.value;
//                 rateConstants.push_back(rateConstant);
//             }
//             return rateConstants;
//         }

//         static types::ReactionSetRead parseReactions(const std::vector<std::string> &lines)
//         {
//             std::vector<types::ReactionRead> reactions;
//         }

//         static types::ReactionRead parseReaction(const std::string &line)
//         {
//             std::vector<std::string> args = str::splitByWhitespace(line);

//             types::ReactionRead reaction;
//             reaction.type = args[0];

//             if (args.size() < 4)
//                 console::input_error("Reaction definition requires at least 4 arguments: <type> <reactants> -<rateConstantName>-> <products>. Provided: " + line);

//             // Loop through args (ignoring first arg)
//             bool isReactants = true;
//             for (size_t i = 1; i < args.size(); ++i)
//             {
//                 std::string reactionArg = args[i];

//                 // Check if string is "-name->"
//                 if (str::startswith(reactionArg, "-") && str::endswith(reactionArg, "->"))
//                 {
//                     std::string rateConstantName = reactionArg.substr(1, reactionArg.size() - 3); // remove "-" and "->"
//                     reaction.rateConstantName = rateConstantName;
//                     isReactants = false;
//                     continue;
//                 }

//                 if (reactionArg == "+")
//                     continue;

//                 if (isReactants)
//                     reaction.reactantNames.push_back(reactionArg);
//                 else
//                     reaction.productNames.push_back(reactionArg);
//             }
//         }

//         static types::SpeciesRead parseSpecies(const std::vector<std::string> &args)
//         {

//             if (args.size() < 2)
//                 console::input_error("Species definition requires at least 2 arguments: <type> <name>");

//             types::SpeciesRead species;
//             species.type = args[0];
//             species.name = args[1];

//             if (!SpeciesType::isValidType(species.type))
//                 console::input_error("Species type " + species.type + " is not valid. Exiting.");

//             return species;
//         }

//         static types::UnitRead parseUnit(const std::vector<std::string> &args)
//         {

//             types::SpeciesRead species = parseSpecies(args);

//             types::UnitRead unit;
//             unit.name = species.name;
//             unit.type = species.type;

//             std::vector<std::string> vars = {args.begin() + 2, args.end()};

//             readVar(vars, "C0", unit.C0);
//             readVar(vars, "FW", unit.FW);

//             if (unit.type == SpeciesType::INITIATOR)
//                 readVar(vars, "f", unit.efficiency, true); // required for initiators
//             else
//                 readVar(vars, "efficiency", unit.efficiency);

//             return unit;
//         }

//         static types::PolymerTypeRead parsePolymerType(const std::vector<std::string> &args)
//         {

//             types::SpeciesRead species = parseSpecies(args);

//             types::PolymerTypeRead polymer;
//             polymer.name = species.name;
//             polymer.type = species.type;

//             // Check if the polymer name has brackets (e.g., "P[A.A]")
//             size_t seqStart = species.name.find("[");
//             size_t seqEnd = species.name.find("]");
//             bool hasSequence = (seqStart != std::string::npos && seqEnd != std::string::npos && seqEnd > seqStart);
//             std::vector<std::string> endGroupUnitNames;

//             if (!hasSequence)
//             {
//                 polymer.endGroupUnitNames = {};
//                 return polymer;
//             }

//             // Get string inside the brackets and split by period
//             // e.g., "P[A.A]" -> "A.A" -> {"A","A"}
//             std::string sequenceStr = species.name.substr(seqStart + 1, seqEnd - seqStart - 1);
//             std::vector<std::string> endGroupNames = str::splitByDelimeter(sequenceStr, ".");
//             polymer.endGroupUnitNames = endGroupNames;
//             return polymer;
//         };

//         static types::PolymerLabelsRead parsePolymerLabels(const std::vector<std::string> &args)
//         {
//             types::SpeciesRead species = parseSpecies(args);

//             types::PolymerLabelsRead labels;
//             labels.name = species.name;
//             labels.type = species.type;

//             std::vector<std::string> vars = {args.begin() + 2, args.end()};
//             if (vars.size() != 1)
//                 console::input_error("Polymer label definition requires polymer names separated by '|', e.g., LABEL MyLabel P[A.A]|P[B.B]");

//             std::vector<std::string> polymerNames = str::splitByDelimeter(vars[0], "|");
//             labels.polymerNames = polymerNames;

//             return labels;
//         }

//         static types::SpeciesSetRead parseSpeciesSet(const std::vector<std::string> &lines)
//         {

//             // ===== PASS 1: Separate lines by species type =====
//             std::vector<std::string> unitLines;
//             std::vector<std::string> polymerLines;
//             std::vector<std::string> labelLines;

//             for (const auto &line : lines)
//             {
//                 std::vector<std::string> args = str::splitByWhitespace(line);
//                 if (args.empty())
//                     continue;

//                 types::SpeciesRead species = parseSpecies(args);
//                 if (SpeciesType::isUnitType(species.type))
//                     unitLines.push_back(line);
//                 else if (species.type == SpeciesType::POLYMER)
//                     polymerLines.push_back(line);
//                 else if (species.type == SpeciesType::LABEL)
//                     labelLines.push_back(line);
//                 else
//                     console::input_error("Unknown species type: " + species.type);
//             }

//             // ===== PASS 2: Parse in order: Units, Polymers, Labels =====
//             std::vector<types::UnitRead> units(unitLines.size());
//             std::vector<types::PolymerTypeRead> polymerTypes(polymerLines.size());
//             std::vector<types::PolymerLabelsRead> labels(labelLines.size());

//             for (const auto &line : unitLines)
//                 units.push_back(parseUnit(str::splitByWhitespace(line)));

//             for (const auto &line : polymerLines)
//                 polymerTypes.push_back(parsePolymerType(str::splitByWhitespace(line)));

//             for (const auto &line : labelLines)
//                 labels.push_back(parsePolymerLabels(str::splitByWhitespace(line)));

//             types::SpeciesSetRead speciesSet;
//             speciesSet.units = std::move(units);
//             speciesSet.polymerTypes = std::move(polymerTypes);
//             speciesSet.polymerLabels = std::move(labels);

//             return speciesSet;
//         }

//     private:
//         template <typename T>
//         static void readVar(const std::vector<std::string> &vars, const std::string &varName, T &value, bool required = false)
//         {
//             bool found = false;
//             for (const auto &varString : vars)
//             {
//                 if (str::startsWith(varString, varName))
//                 {
//                     types::Variable<T> var = parseVariable<T>(varString);
//                     value = var.value;
//                     found = true;
//                     break;
//                 }
//             }
//             if (required && !found)
//                 console::input_error("Required variable '" + varName + "' not found in text input.");
//         }

//         template <typename T>
//         static types::Variable<T> parseVariable(const std::string &s)
//         {
//             std::vector<std::string> tokens = str::splitByDelimeter(s, "=");
//             if (tokens.size() != 2)
//                 console::input_error("Variable '" + s + "' is not in 'name=value' format.");

//             str::trim(tokens[0]);
//             str::trim(tokens[1]);

//             return types::Variable<T>{tokens[0], input::convertValue<T>(tokens[1])};
//         }
//     };
// };

// namespace io::parse::text
// {
//     template <typename T>
//     static void readVar(const std::vector<std::string> &vars, const std::string &varName, T &value, bool required = false)
//     {
//         bool found = false;
//         for (const auto &varString : vars)
//         {
//             if (str::startsWith(varString, varName))
//             {
//                 types::Variable<T> var = parseVariable<T>(varString);
//                 value = var.value;
//                 found = true;
//                 break;
//             }
//         }
//         if (required && !found)
//             console::input_error("Required variable '" + varName + "' not found in text input.");
//     }

//     template <typename T>
//     static types::Variable<T> parseVariable(const std::string &s)
//     {
//         std::vector<std::string> tokens = str::splitByDelimeter(s, "=");
//         if (tokens.size() != 2)
//             console::input_error("Variable '" + s + "' is not in 'name=value' format.");

//         str::trim(tokens[0]);
//         str::trim(tokens[1]);

//         return types::Variable<T>{tokens[0], input::convertValue<T>(tokens[1])};
//     }
// };