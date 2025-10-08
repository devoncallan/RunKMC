#pragma once
#include "common.h"
#include "polymer.h"
#include "analysis/types.h"

/**
 * @brief Stores pointers to polymer objects of a specific type.
 *
 */
class PolymerType
{
public:
    SpeciesID ID;
    std::string name;
    uint64_t count;

    PolymerType(const SpeciesID &ID_, const std::string &name_, const std::vector<SpeciesID> &endGroup_) : ID(ID_), name(name_), endGroup(endGroup_), count(0) {};

    ~PolymerType() {};

    void insertPolymer(Polymer *polymer)
    {
        ++count;
        polymers.push_back(polymer);
    }

    Polymer *removeRandomPolymer()
    {
        --count;
        size_t randomIndex = int(rng_utils::dis(rng_utils::rng) * polymers.size());
        Polymer *polymer = polymers[randomIndex];           // get random polymer
        polymers[randomIndex] = std::move(polymers.back()); // swap
        polymers.pop_back();                                // and pop!
        return polymer;
    }

    const std::vector<Polymer *> &getPolymers() const { return polymers; }

    const std::vector<SpeciesID> &getEndGroup() const { return endGroup; }

private:
    std::vector<Polymer *> polymers;
    std::vector<SpeciesID> endGroup; // endGroup to identify the terminal units on the chain end.
};

typedef PolymerType *PolymerTypePtr;

/**
 * @brief Stores a collection of PolymerType pointers.
 *
 */
class PolymerContainer
{
public:
    SpeciesID ID;
    std::string name;
    uint64_t count = 0;

    PolymerContainer(const SpeciesID &ID_, const std::string &name_, const std::vector<PolymerTypePtr> &polymerTypePtrs_)
        : ID(ID_), name(name_), polymerTypePtrs(polymerTypePtrs_)
    {
        polymerTypeCounts.resize(polymerTypePtrs_.size());
    };
    ~PolymerContainer() {}

    Polymer *removeRandomPolymer()
    {
        if (polymerTypePtrs.size() == 1)
        {
            --count;
            --polymerTypeCounts[0];
            return polymerTypePtrs[0]->removeRandomPolymer();
        }

        std::discrete_distribution<size_t> discrete_dis(polymerTypeCounts.begin(), polymerTypeCounts.end());
        size_t typeIndex = discrete_dis(rng_utils::rng);
        --count;
        --polymerTypeCounts[typeIndex];
        return polymerTypePtrs[typeIndex]->removeRandomPolymer();
    }

    void insertPolymer(Polymer *polymer)
    {
        // No classification needed. Directly store the polymer.
        if (polymerTypePtrs.size() == 1)
        {
            ++count;
            ++polymerTypeCounts[0];
            polymerTypePtrs[0]->insertPolymer(polymer);
            return;
        }

        // Classify the polymer based on its end group.
        for (int i = 0; i < polymerTypePtrs.size(); ++i)
        {
            if (polymer->endGroupIs(polymerTypePtrs[i]->getEndGroup()))
            {
                ++count;
                ++polymerTypeCounts[i];
                polymerTypePtrs[i]->insertPolymer(polymer);
                return;
            }
        }
        console::error("End sequence for inserted polymer does not match. Exiting.....");
    }

    void updatePolymerCounts()
    {
        uint64_t totalCount = 0;
        for (size_t i = 0; i < polymerTypePtrs.size(); ++i)
        {
            auto typeCount = polymerTypePtrs[i]->count;
            polymerTypeCounts[i] = typeCount;
            totalCount += typeCount;
        }
        count = totalCount;
    }

    const std::string toString() const
    {
        return name + ": " + std::to_string(count);
    }

    const std::vector<PolymerTypePtr> &getPolymerTypes() const { return polymerTypePtrs; }

private:
    std::vector<PolymerTypePtr> polymerTypePtrs;
    std::vector<uint64_t> polymerTypeCounts;
};

typedef PolymerContainer *PolymerContainerPtr;

struct PolymerContainerMap
{
    SpeciesID ID;
    std::string name;
    std::vector<size_t> polymerTypeIndices;

    PolymerContainerMap(SpeciesID ID_, std::string name_, std::vector<size_t> polymerTypeIndices_)
        : ID(ID_), name(name_), polymerTypeIndices(polymerTypeIndices_) {};
};