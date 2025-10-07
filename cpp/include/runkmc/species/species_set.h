#pragma once
#include "common.h"
#include "species/polymer_type.h"
#include "kmc/state.h"

class SpeciesSet
{
public:
    SpeciesSet() {};

    SpeciesSet(
        std::vector<Unit> &&units_,
        std::vector<PolymerType> &&polymerTypes_,
        std::vector<PolymerCollectionMap> &&polymerCollections_,
        size_t numParticles_

    )
    {
        // Calculate NAV
        double totalC0 = 0;
        auto unitIDs = registry::getAllUnitIDs();
        for (const auto &unitID : unitIDs)
            totalC0 += units_[unitID].C0;
        NAV = numParticles_ / totalC0;

        // Set initial counts
        for (const auto &unitID : unitIDs)
        {
            double initAmount = this->units[unitID].C0 * NAV;
            uint64_t initCount = static_cast<uint64_t>(initAmount), uint64_t(1);

            if (initAmount < 1)
            {
                console::input_warning("Initial amount of " + this->units[unitID].name + " is less than 1 (" + std::to_string(initAmount) + "). Setting initial count to 1.");
                initCount = 1;
            }

            double roundingError = std::abs((initAmount - double(initCount)) / initAmount);
            if (roundingError > 0.10)
                console::input_error("Initial amount of " + this->units[unitID].name + " has a abs rounding error of " + std::to_string(roundingError * 100) + "%. Consider increasing num_units to reduce this error. Exiting.....");

            this->units[unitID].setInitialCount(initCount);
        }

        polymerGroups.reserve(polymerCollections_.size());
        polymerGroupPtrs.reserve(polymerCollections_.size());

        // Creating polymer groups
        for (const auto &polymerGroup : polymerCollections_)
        {
            auto indices = polymerGroup.polymerTypeIndices;
            std::vector<PolymerTypePtr> polymerSubTypePtrs;
            polymerSubTypePtrs.reserve(indices.size());
            for (const auto &index : indices)
                polymerSubTypePtrs.push_back(&polymerTypes[index]);

            polymerGroups.push_back(PolymerTypeGroup(polymerGroup.name, polymerSubTypePtrs));
            polymerGroupPtrs.push_back(&polymerGroups.back());
        }
    }

    SpeciesSet(
        std::vector<PolymerType> &&polymerTypes_,
        std::vector<PolymerGroupStruct> &&PolymerGroupStructs_,
        std::vector<Unit> &&units_,
        size_t numParticles_) : polymerTypes(std::move(polymerTypes_)), units(std::move(units_)), numParticles(numParticles_)
    {
        // Calculate NAV
        double totalC0 = 0;
        auto unitIDs = registry::getAllUnitIDs();
        for (const auto &unitID : unitIDs)
            totalC0 += this->units[unitID].C0;
        NAV = numParticles / totalC0;

        // Set initial counts
        for (const auto &unitID : unitIDs)
        {
            double initAmount = this->units[unitID].C0 * NAV;
            uint64_t initCount = static_cast<uint64_t>(initAmount), uint64_t(1);

            if (initAmount < 1)
            {
                console::input_warning("Initial amount of " + this->units[unitID].name + " is less than 1 (" + std::to_string(initAmount) + "). Setting initial count to 1.");
                initCount = 1;
            }

            double roundingError = std::abs((initAmount - double(initCount)) / initAmount);
            if (roundingError > 0.10)
                console::input_error("Initial amount of " + this->units[unitID].name + " has a abs rounding error of " + std::to_string(roundingError * 100) + "%. Consider increasing num_units to reduce this error. Exiting.....");

            this->units[unitID].setInitialCount(initCount);
        }

        polymerGroups.reserve(PolymerGroupStructs_.size());
        polymerGroupPtrs.reserve(PolymerGroupStructs_.size());

        // Creating polymer groups
        for (const auto &polymerGroup : PolymerGroupStructs_)
        {
            auto indices = polymerGroup.polymerTypeIndices;
            std::vector<PolymerTypePtr> polymerSubTypePtrs;
            polymerSubTypePtrs.reserve(indices.size());
            for (const auto &index : indices)
                polymerSubTypePtrs.push_back(&polymerTypes[index]);

            polymerGroups.push_back(PolymerTypeGroup(polymerGroup.name, polymerSubTypePtrs));
            polymerGroupPtrs.push_back(&polymerGroups.back());
        }
    };

    void updatePolyTypeGroups()
    {
        for (const auto &polymerGroupPtr : polymerGroupPtrs)
            polymerGroupPtr->updatePolymerCounts();
    };

    SpeciesState getStateData() const
    {
        SpeciesState data;

        // Unit counts / conversions
        auto unitIDs = registry::getAllUnitIDs();
        for (auto &id : unitIDs)
        {
            data.unitCounts.push_back(units[id].count);
            data.unitConversions.push_back(units[id].calculateConversion());
        }

        // Polymer counts
        for (const auto &polymerGroup : polymerGroups)
            data.polymerCounts.push_back(polymerGroup.count);

        // Total conversion
        data.totalConversion = calculateConversion();

        return data;
    }

    double calculateConversion() const
    {
        double numerator = 0;
        double denominator = 0;
        uint64_t initialCount;

        auto monomerIDs = registry::MONOMER_IDS;
        for (const auto &id : monomerIDs)
        {
            initialCount = units[id].getInitialCount();
            numerator += initialCount - units[id].count;
            denominator += initialCount;
        }
        if (denominator == 0)
            return 0;

        return numerator / denominator;
    };

    std::vector<Polymer *> getPolymers() const
    {

        std::vector<Polymer *> polymers;

        // Reserve space for all polymers
        uint64_t numPolymers = 0;
        for (const auto &polymerType : polymerTypes)
            numPolymers += polymerType.count;
        polymers.reserve(numPolymers);

        // Add all polymer pointers to the reserved space
        for (const auto &polymerType : polymerTypes)
        {
            const auto &typePolymers = polymerType.getPolymers();
            polymers.insert(polymers.end(), typePolymers.begin(), typePolymers.end());
        }
        return polymers;
    }

    analysis::RawSequenceData getRawSequenceData() const
    {
        const auto &polymers = getPolymers();
        auto sequenceData = analysis::RawSequenceData(polymers.size());

        for (const auto *polymer : polymers)
        {
            if (!polymer->isCompressed())
                sequenceData.sequences.push_back(polymer->getSequence());
            else
                sequenceData.precomputedStats.push_back(polymer->getPositionalStats());
        }

        return sequenceData;
    };

    void printSummary() const
    {
        console::log("Units:");
        for (const auto &unit : units)
            console::log("\t" + unit.toString());

        // console::log("Polymer Types:");
        // for (const auto &polyType : polymerTypes)
        //     console::log("\t" + polyType.toString());

        console::log("Polymer Groups:");
        for (const auto &polyGroup : polymerGroups)
            console::log("\t" + polyGroup.toString());
    }

    std::vector<double> getMonomerFWs() const
    {
        std::vector<double> monomerFWs;
        monomerFWs.reserve(registry::NUM_MONOMERS);
        for (const auto &id : registry::MONOMER_IDS)
            monomerFWs.push_back(units[id].FW);
        return monomerFWs;
    }

    std::vector<Unit> &getUnits() { return units; }
    const std::vector<Unit> &getUnits() const { return units; }
    std::vector<PolymerTypeGroup> &getPolyTypeGroups() { return polymerGroups; }
    const std::vector<PolymerTypeGroup> &getPolyTypeGroups() const { return polymerGroups; }
    std::vector<PolymerTypeGroupPtr> &getPolymerGroupPtrs() { return polymerGroupPtrs; }
    const std::vector<PolymerTypeGroupPtr> &getPolymerGroupPtrs() const { return polymerGroupPtrs; }
    double getNAV() const { return NAV; }

private:
    std::vector<PolymerType> polymerTypes;
    std::vector<PolymerTypeGroup> polymerGroups;
    std::vector<PolymerTypeGroupPtr> polymerGroupPtrs;

    std::vector<Unit> units;
    size_t numParticles;
    double NAV;
};