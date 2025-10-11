#pragma once
#include "common.h"
#include "kmc/state.h"
#include "kmc/species/polymer_type.h"
#include "kmc/analysis/analysis.h"

namespace species
{
    constexpr console::LogContext logger("KMC::Species");
}
class SpeciesSet
{
public:
    SpeciesSet() {};

    SpeciesSet(
        std::vector<Unit> &&units_,
        std::vector<PolymerType> &&polymerTypes_,
        std::vector<PolymerContainerMap> &&polymerContainerMaps_,
        size_t numParticles_)
    {
        units = std::move(units_);
        polymerTypes = std::move(polymerTypes_);
        numParticles = numParticles_;

        // Calculate NAV
        double totalC0 = 0;
        for (const auto &unit : units)
            totalC0 += unit.C0;
        NAV = numParticles / totalC0;

        // Set initial counts

        for (auto &unit : units)
        {
            double initAmount = unit.C0 * NAV;
            uint64_t initCount = static_cast<uint64_t>(initAmount);
            if (initAmount > 0)
            {
                double roundingError = std::abs((initAmount - double(initCount)) / initAmount);
                if (initAmount < 1)
                {
                    console::input_warning("Initial amount of " + unit.name + " is less than 1 (" + std::to_string(initAmount) + "). Setting initial count to 1.");
                    initCount = 1;
                }
                if (roundingError > 0.10)
                    console::input_error("Initial amount of " + unit.name + " has a abs rounding error of " + std::to_string(roundingError * 100) + "%. Consider increasing num_units to reduce this error. Exiting.....");
            }
            unit.setInitCount(initCount);
        }

        polymerContainers.reserve(polymerContainerMaps_.size());
        polymerContainerPtrs.reserve(polymerContainerMaps_.size());

        // Create polymer containers
        for (const auto &containerMap : polymerContainerMaps_)
        {
            // Get pointers to the polymer types in this container
            auto indices = containerMap.polymerTypeIndices;
            std::vector<PolymerType *> polymerTypePtrs;
            polymerTypePtrs.reserve(indices.size());
            for (const auto &index : indices)
                polymerTypePtrs.push_back(&polymerTypes[index]);

            polymerContainers.push_back(PolymerContainer(containerMap.ID, containerMap.name, polymerTypePtrs));
            polymerContainerPtrs.push_back(&polymerContainers.back());
        }

        printSummary();
    }

    void updatePolymerContainers()
    {
        for (const auto &polymerContainerPtr : polymerContainerPtrs)
            polymerContainerPtr->updatePolymerCounts();
    }

    SpeciesState getStateData() const
    {
        SpeciesState data;

        // Unit counts / conversions
        auto unitIDs = registry::getAllUnitIDs();
        for (auto &id : unitIDs)
        {
            auto idx = registry::getUnitIndex(id);
            data.unitCounts.push_back(units[idx].count);
            data.unitConversions.push_back(units[idx].calculateConversion());
        }

        // Monomer conversion
        data.monomerConversion = calculateMonomerConversion();

        // Polymer counts
        for (const auto &container : polymerContainers)
            data.polymerCounts.push_back(container.count);

        return data;
    }

    double calculateMonomerConversion() const
    {
        double numerator = 0;
        double denominator = 0;
        uint64_t initialCount;

        auto monomerIDs = registry::getMonomerIDs();
        for (const auto &id : monomerIDs)
        {
            auto monIdx = registry::getMonomerIndex(id);
            initialCount = units[monIdx].getInitCount();
            numerator += initialCount - units[monIdx].count;
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

    void analyze(SystemState &systemState)
    {
        auto sequenceData = getRawSequenceData();
        auto summary = analysis::calculateSequenceSummary(sequenceData);

        AnalysisState analysisState;
        analysis::analyzeChainLengthDist(summary.sequenceStatsMatrix, getMonomerFWs(), analysisState);
        systemState.analysis = analysisState;

        if (registry::getNumMonomers() <= 1)
            return;

        SequenceState sequenceState = SequenceState{systemState.kmc, summary.positionalStats};
        analysis::analyzeSequenceLengthDist(summary.sequenceStatsMatrix, analysisState);
        systemState.analysis = analysisState;
        systemState.sequence = sequenceState;
    }

    void printSummary() const
    {
        species::logger.info("Units:");
        for (const auto &unit : units)
            species::logger.info("\t" + unit.printSummary());

        species::logger.info("Polymer Containers:");
        for (const auto &container : polymerContainers)
            species::logger.info("\t" + container.printSummary());
    }

    std::vector<double> getMonomerFWs() const
    {
        std::vector<double> monomerFWs;
        monomerFWs.reserve(registry::getNumMonomers());

        for (const auto &id : registry::getMonomerIDs())
            monomerFWs.push_back(units[registry::getMonomerIndex(id)].FW);

        return monomerFWs;
    }

    std::vector<Unit> &getUnits() { return units; }
    const std::vector<Unit> &getUnits() const { return units; }

    std::vector<PolymerContainer> &getPolymerContainers() { return polymerContainers; }
    const std::vector<PolymerContainer> &getPolymerContainers() const { return polymerContainers; }
    const std::vector<PolymerContainer *> &getPolymerContainerPtrs() const { return polymerContainerPtrs; }

    double getNAV() const { return NAV; }

private:
    std::vector<PolymerType> polymerTypes;
    std::vector<PolymerContainer> polymerContainers;
    std::vector<PolymerContainer *> polymerContainerPtrs;

    std::vector<Unit> units;
    size_t numParticles;
    double NAV;
};