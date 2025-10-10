#pragma once
#include <string>
#include <vector>

#include "utils/parse.h"
#include "utils/console.h"
#include "types.h"
#include "io/types.h"

/**
 * @brief Immutable, optimized, cached lookups for species information.
 *
 */
class Registry
{
public:
    Registry() = default;
    ~Registry() = default;

    // Helper functions
    const std::vector<io::types::RegisteredSpecies> &getAllSpecies() const { return _species; }
    const io::types::RegisteredSpecies &getSpecies(const SpeciesID &id) const { return _idToSpecies.at(id); }
    const io::types::RegisteredSpecies &getSpecies(const std::string &name) const { return _nameToSpecies.at(name); }

    // Unit helper functions
    const std::vector<SpeciesID> &getAllUnitIDs() const { return _allUnitIDs; }
    const std::vector<std::string> &getAllUnitNames() const { return _allUnitNames; }
    size_t getUnitIndex(const SpeciesID &id) const { return _unitIDtoIndices.at(id); }

    // Monomer helper functions
    const std::vector<SpeciesID> getMonomerIDs() const { return _monomerIDs; }
    const std::vector<std::string> getMonomerNames() const { return _monomerNames; }
    size_t getNumMonomers() const { return _numMonomers; }
    bool isMonomer(SpeciesID id) const { return _idToSpecies.at(id).type == SpeciesType::MONOMER; }
    size_t getMonomerIndex(SpeciesID id) const { return _monomerIDtoIndices.at(id); }

    // Polymer helper functions
    const std::vector<std::string> getPolymerNames() const { return _polymerContainerNames; }
    size_t getPolymerIndex(SpeciesID id) const { return _polymerContainerIDtoIndices.at(id); }

private:
    friend class RegistryBuilder;

    Registry(
        const std::vector<io::types::RegisteredSpecies> &species)
    {
        _species = species; // in order of registration
        for (const auto &s : _species)
        {
            _nameToSpecies[s.name] = s;
            _idToSpecies[s.ID] = s;
            _typeToNames[s.type].push_back(s.name);

            if (SpeciesType::isUnitType(s.type))
            {
                _allUnitNames.push_back(s.name);
                _allUnitIDs.push_back(s.ID);
                _unitIDtoIndices[s.ID] = _allUnitNames.size() - 1;
            }
            if (s.type == SpeciesType::MONOMER)
            {
                _monomerIDs.push_back(s.ID);
                _monomerNames.push_back(s.name);
                _monomerIDtoIndices[s.ID] = _monomerIDs.size() - 1;
            }
            if (s.type == SpeciesType::POLYMER)
            {
                _polymerTypeNames.push_back(s.name);
                _polymerTypeIDtoIndices[s.ID] = _polymerTypeNames.size() - 1;
            }
            if (s.type == SpeciesType::POLYMER || s.type == SpeciesType::LABEL)
            {
                _polymerContainerNames.push_back(s.name);
                _polymerContainerIDtoIndices[s.ID] = _polymerContainerNames.size() - 1;
            }
        }
        _numMonomers = _monomerIDs.size();
    };

    // All species
    std::vector<io::types::RegisteredSpecies> _species;

    std::unordered_map<std::string, io::types::RegisteredSpecies> _nameToSpecies;
    std::unordered_map<SpeciesID, io::types::RegisteredSpecies> _idToSpecies;
    std::unordered_map<std::string, std::vector<std::string>> _typeToNames;

    // Unit data (cached)
    std::vector<SpeciesID> _allUnitIDs;
    std::vector<std::string> _allUnitNames;
    std::unordered_map<SpeciesID, size_t> _unitIDtoIndices;

    // Monomer data (cached)
    size_t _numMonomers;
    std::vector<SpeciesID> _monomerIDs;
    std::vector<std::string> _monomerNames;
    std::unordered_map<SpeciesID, size_t> _monomerIDtoIndices;

    // Polymer data (cached)
    std::vector<std::string> _polymerTypeNames;
    std::unordered_map<SpeciesID, size_t> _polymerTypeIDtoIndices;

    std::vector<std::string> _polymerContainerNames;
    std::unordered_map<SpeciesID, size_t> _polymerContainerIDtoIndices;
};

class RegistryBuilder
{
public:
    size_t findSpecies(const std::string &name) const { return input::findInVector(name, registered_species); }
    bool isRegistered(const std::string &name) const { return findSpecies(name) != SIZE_T_MAX; }

    io::types::RegisteredSpecies getSpecies(const std::string &name) const
    {
        size_t index = findSpecies(name);
        if (index != SIZE_T_MAX)
            return registered_species[index];
        console::input_error("Species with name " + name + " is not registered.");
    }

    SpeciesID getSpeciesID(const std::string &name) const { return getSpecies(name).ID; }

    SpeciesID registerNewSpecies(const std::string &name, const std::string &type)
    {
        if (finalized)
            console::error("Cannot register new species after registry has been finalized.");

        SpeciesType::checkValid(type);

        if (isRegistered(name))
            console::input_error("Species with name " + name + " already registered.");

        // Assign IDs sequentially
        SpeciesID newID = registered_species.size() + 1;
        registered_species.push_back({name, type, newID});
        return newID;
    }

    Registry build()
    {
        finalized = true;
        return Registry(registered_species);
    }

private:
    bool finalized = false;
    std::vector<io::types::RegisteredSpecies> registered_species;
};

namespace registry
{
    static RegistryBuilder builder;
    static Registry _instance;

    // Initialize registry
    static void initialize() { _instance = builder.build(); }

    // Helper functions
    static inline const std::vector<io::types::RegisteredSpecies> &getAllSpecies() { return _instance.getAllSpecies(); }
    static inline const io::types::RegisteredSpecies &getSpecies(const SpeciesID &id) { return _instance.getSpecies(id); }
    static inline const io::types::RegisteredSpecies &getSpecies(const std::string &name) { return _instance.getSpecies(name); }

    // Unit helpers
    static inline std::vector<SpeciesID> getAllUnitIDs() { return _instance.getAllUnitIDs(); }
    static inline std::vector<std::string> getAllUnitNames() { return _instance.getAllUnitNames(); }
    static inline size_t getUnitIndex(const SpeciesID &id) { return _instance.getUnitIndex(id); }

    // Monomer helpers
    static inline std::vector<std::string> getMonomerNames() { return _instance.getMonomerNames(); }
    static inline std::vector<SpeciesID> getMonomerIDs() { return _instance.getMonomerIDs(); }
    static inline size_t getNumMonomers() { return _instance.getNumMonomers(); }
    static inline size_t getMonomerIndex(const SpeciesID &id) { return _instance.getMonomerIndex(id); }
    static inline bool isMonomer(const SpeciesID &id) { return _instance.isMonomer(id); }

    // Polymer helpers
    static inline std::vector<std::string> getPolymerNames() { return _instance.getPolymerNames(); }
    static inline size_t getPolymerIndex(const SpeciesID &id) { return _instance.getPolymerIndex(id); }
}