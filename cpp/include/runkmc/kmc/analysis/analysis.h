#pragma once
#include "common.h"
#include "kmc/species/types.h"
#include "kmc/analysis/types.h"

namespace analysis
{
    inline CompressedCopolymerBuffer compressCopolymer(const CopolymerBuffer &buf)
    {
        const size_t numBuckets = NUM_BUCKETS;
        const size_t numMonomers = registry::getNumMonomers();
        const auto &units = buf.units;

        CompressedCopolymerBuffer result;
        result.length = static_cast<uint32_t>(units.size());
        result.posStats.resize(numBuckets);

        if (units.empty())
            return result;

        SpeciesID currentMonomerID = 0;
        size_t currentSequenceLength = 0;

        auto flushRun = [&](size_t bucket)
        {
            if (currentSequenceLength == 0)
                return;
            result.posStats[bucket].addSequence(currentMonomerID, currentSequenceLength);
            size_t monIdx = registry::getMonomerIndex(currentMonomerID);
            auto &counts = result.segHist.data[static_cast<uint32_t>(currentSequenceLength)];
            if (counts.size() <= monIdx)
                counts.resize(numMonomers, 0);
            counts[monIdx]++;
        };

        for (size_t i = 0; i < units.size(); ++i)
        {
            size_t bucket = utils::getBucketIndex(i, units.size(), numBuckets);
            SpeciesID id = units[i];

            if (!registry::isMonomer(id))
                continue;

            if (id == currentMonomerID)
            {
                currentSequenceLength++;
                continue;
            }

            flushRun(bucket);

            currentMonomerID = id;
            currentSequenceLength = 1;
        }

        size_t lastBucket = utils::getBucketIndex(units.size() - 1, units.size(), numBuckets);
        flushRun(lastBucket);

        return result;
    }
}
