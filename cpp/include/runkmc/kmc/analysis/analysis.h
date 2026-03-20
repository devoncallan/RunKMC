#pragma once
#include "common.h"
#include "kmc/species/types.h"
#include "kmc/analysis/types.h"

namespace analysis
{
    CompressedCopolymerBuffer compressCopolymer(const CopolymerBuffer &buf)
    {
        const size_t numBuckets = NUM_BUCKETS;
        const size_t numMonomers = registry::getNumMonomers();
        const auto &units = buf.units;

        CompressedCopolymerBuffer result;
        result.length = static_cast<uint32_t>(units.size());

        if (units.empty())
            return result;

        SpeciesID currentMonomerID = 0;
        size_t currentSequenceLength = 0;

        auto flushRun = [&](size_t bucket)
        {
            if (currentSequenceLength == 0)
                return;
            result.posStats.buckets[bucket].addSequence(currentMonomerID, currentSequenceLength);
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

    ChainStats toChainStats(const HomopolymerBuffer &_chain, const std::vector<double> &FWs)
    {
        ChainStats cs;
        const double L = static_cast<double>(_chain.length);
        cs.chainLength.add(L);
        if (_chain.monomer != INVALID_SPECIES_ID && _chain.length > 0)
        {
            size_t monIdx = registry::getMonomerIndex(_chain.monomer);
            cs.sequenceLengths.stats[monIdx].add(L);
            cs.chainMW.add(L * FWs[monIdx]);
        }
        return cs;
    }

    ChainStats toChainStats(const CompressedCopolymerBuffer &_chain, const std::vector<double> &FWs)
    {
        ChainStats cs;
        cs.chainLength.add(static_cast<double>(_chain.length));

        double mw = 0;
        for (const auto &s : _chain.posStats.buckets)
        {
            cs.sequenceLengths += s;
            for (size_t i = 0; i < FWs.size(); ++i)
                mw += s.stats[i].sum * FWs[i];
        }
        cs.chainMW.add(mw);

        return cs;
    }
}
