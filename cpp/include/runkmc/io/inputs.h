#pragma once
#include <string>
#include <vector>
#include <fstream>

#include "core/types.h"
#include "core/C.h"
#include "utils/parse.h"
#include "utils/console.h"

#include "io/text.h"

namespace io::inputs
{
    static types::KMCInputRead parseTextModelFile(const std::string &filepath)
    {

        std::ifstream modelFile(filepath);
        if (!modelFile.is_open())
            console::input_error("Cannot open model file: " + filepath);

        std::vector<std::string> parameters;
        std::vector<std::string> species;
        std::vector<std::string> rateConstants;
        std::vector<std::string> reactions;

        std::string line;
        while (std::getline(modelFile, line))
        {
            str::trim(line);
            if (utils::canIgnoreLine(line))
                continue;

            if (parameters.empty() && str::startswith(line, C::io::PARAMETERS_SECTION))
                parameters = utils::parseSection(modelFile, C::io::PARAMETERS_SECTION);

            else if (species.empty() && str::startswith(line, C::io::SPECIES_SECTION))
                species = utils::parseSection(modelFile, C::io::SPECIES_SECTION);

            else if (rateConstants.empty() && str::startswith(line, C::io::RATE_CONSTANTS_SECTION))
                rateConstants = utils::parseSection(modelFile, C::io::RATE_CONSTANTS_SECTION);

            else if (reactions.empty() && str::startswith(line, C::io::REACTIONS_SECTION))
                reactions = utils::parseSection(modelFile, C::io::REACTIONS_SECTION);
        }

        if (parameters.empty() || species.empty() || rateConstants.empty() || reactions.empty())
            console::input_error("Missing required sections in model file.");

        types::KMCInputRead input;
        input.config = parseSimulationConfig(parameters);
        input.species = parseSpecies(species);
        input.rateConstants = parseRateConstants(rateConstants);
        input.reactions = parseReactions(reactions);

        return input;
    }

    static types::KMCInputRead parseYamlModelFile(const std::string &filepath)
    {
        YAML::Node root = YAML::LoadFile(filepath);

        types::KMCInputRead input = yaml::Parser<types::KMCInputRead>::read(root);

        return input;
    }

    

}