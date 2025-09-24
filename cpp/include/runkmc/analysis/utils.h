#pragma once
#include "common.h"
#include "analysis/types.h"

namespace analysis
{

    static size_t getBucketIndex(size_t position, size_t chainLength, size_t numBuckets)
    {
        if (chainLength <= 1)
            return 0;

        double normalizedPos = static_cast<double>(position) / (chainLength);

        size_t bucket = static_cast<size_t>(normalizedPos * numBuckets);
        return (bucket == numBuckets) ? numBuckets - 1 : bucket;
    }

    // Calculate sequence statistics for a single polymer sequence, divided into buckets
    std::vector<SequenceStats> calculatePositionalSequenceStats(const std::vector<SpeciesID> &sequence, const size_t &numBuckets)
    {
        std::vector<SequenceStats> stats(numBuckets);
        if (sequence.empty())
            return stats;

        SpeciesID currentMonomer = 0;
        size_t currentSequenceLength = 0;
        bool inSequence = false;

        for (size_t i = 0; i < sequence.size(); ++i)
        {
            size_t bucket = getBucketIndex(i, sequence.size(), numBuckets);
            SpeciesID id = sequence[i];

            // Skip non-monomer units
            if (!registry::isType(id, SpeciesType::MONOMER))
                continue;

            size_t monomerIndex = registry::getIndex(id, SpeciesType::MONOMER);
            stats[bucket].monCounts[monomerIndex]++;

            if (id == currentMonomer)
            {
                currentSequenceLength++;
            }
            else if (inSequence)
            {
                size_t prevIndex = registry::getIndex(currentMonomer, SpeciesType::MONOMER);
                stats[bucket].seqCounts[prevIndex] += 1;
                stats[bucket].seqLengths2[prevIndex] += currentSequenceLength * currentSequenceLength;
                currentSequenceLength = 1;
            }
            else
            {
                inSequence = true;
                currentSequenceLength = 1;
            }

            currentMonomer = id;
        }

        // Add the stats for the last sequence
        size_t bucket = getBucketIndex(sequence.size() - 1, sequence.size(), numBuckets);
        size_t lastMonomerIdx = registry::getIndex(currentMonomer, SpeciesType::MONOMER);
        stats[bucket].seqCounts[lastMonomerIdx] += 1;
        stats[bucket].seqLengths2[lastMonomerIdx] += currentSequenceLength * currentSequenceLength;

        return stats;
    }

    template <typename Func>
    void forEachStats(const RawSequenceData &sequenceData, size_t numBuckets, Func callback)
    {
        // Process sequences on-the-fly
        auto numSequences = sequenceData.sequences.size();
        for (size_t i = 0; i < numSequences; ++i)
        {
            auto stats = calculatePositionalSequenceStats(sequenceData.sequences[i], numBuckets);
            callback(i, stats);
        }

        // Process precomputed
        auto numPrecomputed = sequenceData.precomputedStats.size();
        for (size_t i = 0; i < numPrecomputed; ++i)
        {
            callback(numSequences + i, sequenceData.precomputedStats[i]);
        }
    }

    /*
     */
    SequenceSummary calculateSequenceSummary(const analysis::RawSequenceData &sequenceData)
    {
        // Calculate sequence stats matrix (polymers x (monomers*fields)) -> Summed across all buckets
        // Calculate positional average stats (buckets x (monomers*fields)) -> Summed across all polymers
        Eigen::MatrixXd sequenceStatsMatrix = Eigen::MatrixXd::Zero(sequenceData.length, SequenceStats::SIZE());
        std::vector<SequenceStats> positionalStats(NUM_BUCKETS);

        forEachStats(
            sequenceData,
            NUM_BUCKETS,
            [&](size_t index, const std::vector<SequenceStats> &allStats)
            {
                for (size_t bucket = 0; bucket < NUM_BUCKETS; ++bucket)
                {
                    const auto &stats = allStats[bucket];
                    sequenceStatsMatrix.row(index) += stats.toEigen();
                    positionalStats[bucket] += stats;
                }
            });

        return SequenceSummary{sequenceStatsMatrix, positionalStats};
    }
}
