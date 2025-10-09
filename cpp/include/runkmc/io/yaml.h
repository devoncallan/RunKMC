#pragma once
#include <string>
#include <vector>

#include "common.h"
#include "core/serialization.h"

namespace io::parse
{
    static types::KMCInputRead parseKMCInputFromYaml(const std::string &filepath)
    {
        YAML::Node root = YAML::LoadFile(filepath);

        types::KMCInputRead input = yaml::Parser<types::KMCInputRead>::read(root);

        return input;
    }
};