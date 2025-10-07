

#pragma once
#include <string>
#include <vector>
#include "types.h"

class Registry
{
public:
    // Immutable, optimized, cached lookups
    const std::vector<SpeciesID> &getAllUnitIDs() const { return _allUnitIDs; }
    const std::vector<std::string> &getAllUnitNames() const { return _allUnitNames; }
    const std::vector<std::string> &getNamesOf(const std::string &type) const { return _typeToNames.at(type); }
    const std::vector<SpeciesID> getMonomerIDs() const { return _monomerIDs; }
    size_t getNumMonomers() const { return _numMonomers; }

    bool isMonomer(SpeciesID id) const { return _idToSpecies.at(id).type == SpeciesType::MONOMER; }
    size_t getMonomerIndex(SpeciesID id) const { return _monomerIDtoIndices.at(id); }

private:
    friend class RegistryBuilder;

    Registry(
        const std::vector<types::RegisteredSpecies> &species)
    {

        _species = species;
        for (const auto &s : species)
        {
            _nameToSpecies[s.name] = s;
            _idToSpecies[s.ID] = s;

            size_t index = _typeToNames[s.type].size();
            _typeToNames[s.type].push_back(s.name);

            if (SpeciesType::isUnitType(s.type))
            {
                _allUnitNames.push_back(s.name);
                _allUnitIDs.push_back(s.ID);
            }
            if (s.type == SpeciesType::MONOMER)
            {
                size_t monomerIndex = _monomerIDs.size();
                _monomerIDtoIndices[s.ID] = monomerIndex;
                _monomerIDs.push_back(s.ID);
            }
        }
        _numMonomers = _monomerIDs.size();
    };

    // getAllUnitIDs
    // MONOMER_IDS
    // NUM_MONOMERS
    // getAllUnitNames
    // getIndex
    // isType

    std::vector<types::RegisteredSpecies> _species;
    std::unordered_map<SpeciesID, types::RegisteredSpecies> _idToSpecies;
    std::unordered_map<std::string, types::RegisteredSpecies> _nameToSpecies;
    std::unordered_map<std::string, std::vector<std::string>> _typeToNames;

    size_t _numMonomers;
    std::vector<SpeciesID> _monomerIDs;
    std::unordered_map<SpeciesID, size_t> _monomerIDtoIndices;

    std::vector<SpeciesID> _allUnitIDs;
    std::vector<std::string> _allUnitNames;
};

class RegistryBuilder
{

public:
    bool isRegistered(const std::string &name) const
    {
        for (const auto &species : registered_species)
        {
            if (species.name == name)
                return true;
        }
        return false;
    }

    SpeciesID registerNewSpecies(const std::string &name, const std::string &type)
    {
        if (finalized)
            throw std::runtime_error("Cannot register new species after registry has been finalized.");

        SpeciesType::checkValid(type);

        if (isRegistered(name))
            console::input_error("Species with name " + name + " already registered.");

        SpeciesID newID = registered_species.size() + 1;
        registered_species.push_back({name, type, newID});
        return newID;
    }

    SpeciesID getSpeciesID(const std::string &name) const
    {
        for (const auto &species : registered_species)
        {
            if (species.name == name)
                return species.ID;
        }
        throw std::runtime_error("Species with name " + name + " is not registered.");
    }

    Registry build()
    {
        finalized = true;
        return Registry(registered_species);
    }

private:
    bool finalized = false;
    std::vector<types::RegisteredSpecies> registered_species;
};

namespace registry
{
    static RegistryBuilder builder;
    static Registry _instance;

    static void initialize() { _instance = builder.build(); }
    inline std::vector<SpeciesID> getAllUnitIDs() { return _instance.getAllUnitIDs(); }
    inline std::vector<std::string> getAllUnitNames() { return _instance.getAllUnitNames(); }
    inline size_t getNumMonomers() { return _instance.getNumMonomers(); }
}