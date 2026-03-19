#pragma once
#include "common.h"
#include "kmc/analysis/types.h"

namespace analysis
{

    namespace utils
    {
        static size_t getBucketIndex(size_t position, size_t chainLength, size_t numBuckets);
    }

    // Calculate sequence statistics for a single polymer sequence, divided into buckets
    std::vector<SequenceStats> calculatePositionalSequenceStats(const std::vector<SpeciesID> &sequence, const size_t &numBuckets)
    {
        std::vector<SequenceStats> stats(numBuckets);
        if (sequence.empty())
            return stats;

        SpeciesID currentMonomerID = 0;
        size_t currentSequenceLength = 0;

        for (size_t i = 0; i < sequence.size(); ++i)
        {
            size_t bucket = utils::getBucketIndex(i, sequence.size(), numBuckets);
            SpeciesID id = sequence[i];

            if (!registry::isMonomer(id))
                continue;

            if (id == currentMonomerID)
            {
                currentSequenceLength++;
                continue;
            }

            if (currentSequenceLength > 0)
                stats[bucket].addSequence(currentMonomerID, currentSequenceLength);

            currentMonomerID = id;
            currentSequenceLength = 1;
        }

        size_t bucket = utils::getBucketIndex(sequence.size() - 1, sequence.size(), numBuckets);
        if (currentSequenceLength > 0)
            stats[bucket].addSequence(currentMonomerID, currentSequenceLength);

        return stats;
    }

} // namespace analysis

namespace analysis::utils
{
    static size_t getBucketIndex(size_t position, size_t chainLength, size_t numBuckets)
    {
        if (chainLength <= 1)
            return 0;

        double normalizedPos = static_cast<double>(position) / (chainLength);
        size_t bucket = static_cast<size_t>(normalizedPos * numBuckets);
        return (bucket == numBuckets) ? numBuckets - 1 : bucket;
    }
}
