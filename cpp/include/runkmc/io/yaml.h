#pragma once
#include <string>
#include <vector>

#include "common.h"
#include "utils/yaml.h"
#include "core/serialization.h"

namespace io::parse
{
    static types::KMCInputRead parseKMCInput(const std::string &filepath)
    {
        YAML::Node root = YAML::LoadFile(filepath);

        types::KMCInputRead input;

        if (!yaml::hasKey(root, C::io::PARAMETERS_SECTION))
            console::input_error("YAML input file is missing required section: " + std::string(C::io::PARAMETERS_SECTION) + ".");
        if (!yaml::hasKey(root, C::io::SPECIES_SECTION))
            console::input_error("YAML input file is missing required section: " + std::string(C::io::SPECIES_SECTION) + ".");
        if (!yaml::hasKey(root, C::io::RATE_CONSTANTS_SECTION))
            console::input_error("YAML input file is missing required section: " + std::string(C::io::RATE_CONSTANTS_SECTION) + ".");
        if (!yaml::hasKey(root, C::io::REACTIONS_SECTION))
            console::input_error("YAML input file is missing required section: " + std::string(C::io::REACTIONS_SECTION) + ".");

        input.config = yaml::Parser<config::SimulationConfig>::read(root[C::io::PARAMETERS_SECTION]);
        input.species = yaml::Parser<types::SpeciesSetRead>::read(root[C::io::SPECIES_SECTION]);
        input.rateConstants = yaml::Parser<std::vector<types::RateConstantRead>>::read(root[C::io::RATE_CONSTANTS_SECTION]);
        input.reactions = yaml::Parser<std::vector<types::ReactionRead>>::read(root[C::io::REACTIONS_SECTION]);

        return input;
    }
};