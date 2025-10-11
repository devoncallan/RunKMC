#pragma once
#include <yaml-cpp/yaml.h>

#include "common.h"
#include "kmc/kmc.h"

namespace output
{
    using namespace io;

    void writeSpeciesRegistry(const std::filesystem::path &filepath)
    {
        std::vector<RegisteredSpecies> species = registry::getAllSpecies();

        YAML::Node root;
        root[C::io::SPECIES_KEY] = yaml::Parser<std::vector<RegisteredSpecies>>::write(species);
        root[C::io::UNITS_KEY] = registry::getAllUnitNames();
        root[C::io::MONOMERS_KEY] = registry::getMonomerNames();
        root[C::io::POLYMERS_KEY] = registry::getPolymerNames();

        io::yaml::writeYamlToFile(filepath, root);
    }

    void writeInputFile(const types::KMCInputRead &data, const std::filesystem::path &filepath)
    {
        YAML::Node root;
        root[C::io::PARAMETERS_SECTION] = yaml::Parser<io::types::SimulationConfig>::write(data.config);
        root[C::io::SPECIES_SECTION] = yaml::Parser<types::SpeciesSetRead>::write(data.species);
        root[C::io::RATE_CONSTANTS_SECTION] = yaml::Parser<std::vector<types::RateConstantRead>>::write(data.rateConstants);
        root[C::io::REACTIONS_SECTION] = yaml::Parser<std::vector<types::ReactionRead>>::write(data.reactions);

        io::yaml::writeYamlToFile(filepath, root);
    }
};