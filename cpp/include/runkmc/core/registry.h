#pragma once
#include <unordered_map>
#include <species/unit.h>
#include <vector>
#include <algorithm>

#include "core/types.h"

#define AVOGADROS 6.022149e+23;

typedef std::string SpeciesNameStr;
typedef std::string SpeciesTypeStr;

struct RegisteredSpecies
{
    SpeciesNameStr name;
    SpeciesTypeStr type;
    SpeciesID ID;
};

namespace registry
{
    static std::vector<RegisteredSpecies> REGISTERED_SPECIES;
    static std::unordered_map<SpeciesTypeStr, std::vector<SpeciesID>> SPECIES_IDS;
    static std::unordered_map<SpeciesTypeStr, std::vector<SpeciesNameStr>> SPECIES_NAMES;

    static size_t NUM_MONOMERS;
    static std::vector<SpeciesID> MONOMER_IDS;

    RegisteredSpecies getByID(SpeciesID id)
    {
        auto it = std::find_if(REGISTERED_SPECIES.begin(), REGISTERED_SPECIES.end(),
                               [id](const RegisteredSpecies &species)
                               { return species.ID == id; });
        if (it != REGISTERED_SPECIES.end())
            return *it;
        throw std::invalid_argument("Species with ID " + std::to_string(id) + " not found");
    }

    RegisteredSpecies getByName(SpeciesNameStr name)
    {
        auto it = std::find_if(REGISTERED_SPECIES.begin(), REGISTERED_SPECIES.end(),
                               [name](const RegisteredSpecies &species)
                               { return species.name == name; });
        if (it != REGISTERED_SPECIES.end())
            return *it;
        throw std::invalid_argument("Species with name " + name + " not found");
    }

    static std::vector<SpeciesNameStr> getNamesOf(SpeciesTypeStr type)
    {
        SpeciesType::checkValid(type);
        if (SPECIES_NAMES.find(type) != SPECIES_NAMES.end())
            return SPECIES_NAMES[type];
        return {};
    }

    static std::vector<SpeciesID> getIDsOf(SpeciesTypeStr type)
    {
        SpeciesType::checkValid(type);
        if (SPECIES_IDS.find(type) != SPECIES_IDS.end())
            return SPECIES_IDS[type];
        return {};
    }

    static size_t getNumOf(SpeciesTypeStr type) { return getIDsOf(type).size(); }

    static std::vector<SpeciesID> getAllUnitIDs()
    {
        std::vector<SpeciesID> unitIDs;
        auto unitTypes = {SpeciesType::UNIT, SpeciesType::MONOMER, SpeciesType::INITIATOR};
        for (const auto &type : unitTypes)
        {
            auto ids = getIDsOf(type);
            unitIDs.insert(unitIDs.end(), ids.begin(), ids.end());
        }
        std::sort(unitIDs.begin(), unitIDs.end());

        return unitIDs;
    }

    static std::vector<SpeciesNameStr> getAllUnitNames()
    {
        std::vector<SpeciesNameStr> unitNames;
        auto unitIDs = getAllUnitIDs();
        for (const auto &id : unitIDs)
            unitNames.push_back(getByID(id).name);

        return unitNames;
    }

    size_t getIndex(SpeciesID id, SpeciesTypeStr type)
    {
        SpeciesType::checkValid(type);
        auto ids = getIDsOf(type);
        auto it = std::find(ids.begin(), ids.end(), id);
        if (it != ids.end())
            return std::distance(ids.begin(), it);
        return SIZE_MAX;
    }

    size_t getIndex(SpeciesNameStr name, SpeciesTypeStr type)
    {
        SpeciesType::checkValid(type);
        auto names = getNamesOf(type);
        auto it = std::find(names.begin(), names.end(), name);
        if (it != names.end())
            return std::distance(names.begin(), it);
        return SIZE_MAX;
    }

    bool isType(SpeciesID id, SpeciesTypeStr type) { return registry::getIndex(id, type) != SIZE_MAX; };
    bool isType(SpeciesNameStr name, SpeciesTypeStr type) { return registry::getIndex(name, type) != SIZE_MAX; };

    static SpeciesID registerNewSpecies(SpeciesNameStr name, SpeciesTypeStr type)
    {
        SpeciesType::checkValid(type);

        if (registry::isType(name, type))
            throw std::invalid_argument("Species with name " + name + " and type " + type + " already registered");

        SpeciesID newID = REGISTERED_SPECIES.size() + 1;
        REGISTERED_SPECIES.push_back(RegisteredSpecies{name, type, newID});
        SPECIES_IDS[type].push_back(newID);
        SPECIES_NAMES[type].push_back(name);
        return newID;
    }

    static void finalizeRegistry()
    {
        if (REGISTERED_SPECIES.empty())
            console::error("No species registered.");
        assert(SPECIES_IDS.size() == SPECIES_NAMES.size());
        for (const auto &[type, ids] : SPECIES_IDS)
        {
            if (SPECIES_NAMES.find(type) == SPECIES_NAMES.end())
                console::error("Species type " + type + " has IDs but no names registered.");
            auto names = SPECIES_NAMES[type];
            if (ids.size() != names.size())
                console::error("Species type " + type + " has mismatched number of IDs and names registered.");
        }

        NUM_MONOMERS = getNumOf(SpeciesType::MONOMER);
        MONOMER_IDS = getIDsOf(SpeciesType::MONOMER);
    }

    static void printRegisteredSpecies()
    {
        console::log("Registered Species:");
        for (const auto &species : REGISTERED_SPECIES)
            console::log("\tID: " + std::to_string(species.ID) + ", Name: " + species.name + ", Type: " + species.type);

        console::debug("Registered Species by Type:");
        for (const auto &[type, ids] : SPECIES_IDS)
        {
            console::debug("Type: " + type);
            for (const auto &id : ids)
            {
                RegisteredSpecies regSpecies = getByID(id);
                console::debug("\tID: " + std::to_string(id) + ", Name: " + regSpecies.name);
            }
        }
    }
}