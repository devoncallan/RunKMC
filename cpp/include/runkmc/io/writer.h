#pragma once
#include <yaml-cpp/yaml.h>
#include <string>
#include <vector>

#include "common.h"
#include "kmc/kmc.h"
#include "utils/yaml.h"
#include "types.h"

namespace io::writers
{
    class yaml
    {
    public:
        yaml() = delete;

        static YAML::Node writeUnit(const Unit &unit);
        static YAML::Node writePolymerType(const PolymerType &polyType);
        static YAML::Node writePolymerGroup(const PolymerTypeGroup &polyGroup);
        static YAML::Node writeSpeciesSet(const SpeciesSet &speciesSet);

        static YAML::Node writeRateConstant(const RateConstant &rateConstant);
        static YAML::Node writeReaction(const Reaction &reaction);
        static YAML::Node writeReactionSet(const ReactionSet &reactionSet);
    };
}