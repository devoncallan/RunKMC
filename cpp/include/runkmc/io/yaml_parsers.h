#pragma once
#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>

#include "common.h"
#include "kmc/kmc.h"
#include "utils/yaml.h"
#include "types.h"

namespace io::parsers
{
    class yaml
    {
    public:
        types::SpeciesRead parseSpecies(const YAML::Node &node)
        {
            types::SpeciesRead species;
            readVar(node, "name", species.name, true);
            readVar(node, "type", species.type, true);
            return species;
        };

        types::UnitRead parseUnit(const YAML::Node &node)
        {
            types::SpeciesRead species = parseSpecies(node);

            types::UnitRead unit;
            unit.name = species.name;
            unit.type = species.type;

            readVar(node, "C0", unit.C0);
            readVar(node, "FW", unit.FW);
            readVar(node, "efficiency", unit.efficiency); // optional, default 1

            return unit;
        }

        types::PolymerTypeRead parsePolymerType(const YAML::Node &node)
        {
            types::SpeciesRead species = parseSpecies(node);
            types::PolymerTypeRead polymer;
            polymer.name = species.name;
            polymer.type = species.type;

            readVar(node, "end_group_units", polymer.endGroupUnitNames, true); // required

            return polymer;
        }

        types::PolymerLabelsRead parsePolymerLabels(const YAML::Node &node)
        {
            types::SpeciesRead species = parseSpecies(node);

            types::PolymerLabelsRead labels;
            labels.name = species.name;
            labels.type = species.type;

            readVar(node, "polymer_names", labels.polymerNames, true); // required

            return labels;
        }

    private:
        template <typename T>
        static void readVar(const YAML::Node &node, const std::string &varName, T &var, bool required = false)
        {
            if (node[varName])
                var = node[varName].as<T>();
            else if (required)
                throw std::runtime_error("Missing required variable: " + varName);
        }
    };
};