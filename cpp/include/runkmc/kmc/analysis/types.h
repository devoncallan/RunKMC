#pragma once
#include <map>
#include <vector>

#include "common.h"

namespace analysis
{
    namespace utils
    {
        double safeDivide(double num, double denom)
        {
            if (denom == 0)
                return 0;
            return num / denom;
        };

        static size_t getBucketIndex(size_t position, size_t chainLength, size_t numBuckets)
        {
            if (chainLength <= 1)
                return 0;

            double normalizedPos = static_cast<double>(position) / (chainLength);
            size_t bucket = static_cast<size_t>(normalizedPos * numBuckets);
            return (bucket == numBuckets) ? numBuckets - 1 : bucket;
        };
    }

    struct MomentAccumulator
    {
        double count = 0;
        double sum = 0;
        double sumSq = 0;

        void add(double value)
        {
            count += 1;
            sum += value;
            sumSq += value * value;
        }

        MomentAccumulator &operator+=(const MomentAccumulator &other)
        {
            count += other.count;
            sum += other.sum;
            sumSq += other.sumSq;
            return *this;
        }

        double nAvg() const { return utils::safeDivide(sum, count); }

        double wAvg() const { return utils::safeDivide(sumSq, sum); }

        double disp() const { return utils::safeDivide(wAvg(), nAvg()); }
    };

    struct SequenceStats
    {
        std::vector<MomentAccumulator> sequences; // one per monomer

        SequenceStats() { sequences.resize(registry::getNumMonomers()); }

        SequenceStats &operator+=(const SequenceStats &other)
        {
            for (size_t i = 0; i < sequences.size(); ++i)
                sequences[i] += other.sequences[i];
            return *this;
        }

        void addSequence(SpeciesID id, size_t length)
        {
            sequences[registry::getMonomerIndex(id)].add(static_cast<double>(length));
        }
    };

    // Segment length histogram: maps segment length → per-monomer counts.
    struct SegmentHistogram
    {
        std::map<uint32_t, std::vector<uint64_t>> data;

        bool empty() const { return data.empty(); }
        size_t size() const { return data.size(); }
        auto begin() const { return data.begin(); }
        auto end() const { return data.end(); }

        SegmentHistogram &operator+=(const SegmentHistogram &other)
        {
            for (const auto &[length, counts] : other.data)
            {
                auto &dest = data[length];
                if (dest.size() < counts.size())
                    dest.resize(counts.size(), 0);
                for (size_t i = 0; i < counts.size(); ++i)
                    dest[i] += counts[i];
            }
            return *this;
        }
    };

    struct ChainStats
    {
        size_t numChains = 0;
        MomentAccumulator chainLength;
        MomentAccumulator chainMW;
        std::vector<MomentAccumulator> sequenceLengths; // one per monomer: run-length moments

        ChainStats() { sequenceLengths.resize(registry::getNumMonomers()); }

        double totalMonomerSum() const {
            double total = 0;
            for (const auto &s : sequenceLengths)
                total += s.sum;
            return total;
        }
        double nAvgComp(size_t monomerIdx) const { return utils::safeDivide(sequenceLengths[monomerIdx].sum, totalMonomerSum()); }

        ChainStats &operator+=(const ChainStats &other)
        {
            numChains += other.numChains;
            chainLength += other.chainLength;
            chainMW += other.chainMW;
            for (size_t i = 0; i < sequenceLengths.size(); ++i)
                sequenceLengths[i] += other.sequenceLengths[i];
            return *this;
        }

        void add(size_t length, const std::vector<SequenceStats> &posStats, const std::vector<double> &FWs = {})
        {
            const double L = static_cast<double>(length);
            chainLength.add(L);
            numChains++;

            if (!FWs.empty())
            {
                double mw = 0;
                if (posStats.empty())
                    mw = L * FWs[0]; // homopolymer
                else
                    for (const auto &stats : posStats)
                        for (size_t i = 0; i < FWs.size(); ++i)
                            mw += stats.sequences[i].sum * FWs[i];
                chainMW.add(mw);
            }

            for (const auto &stats : posStats)
                for (size_t i = 0; i < sequenceLengths.size(); ++i)
                    sequenceLengths[i] += stats.sequences[i];
        }
    };

    // Per-interval positional statistics, structured as a vector of ChainStats — one per bucket.
    // Accumulates sequence statistics broken down by normalized position along the chain.
    // Reset after each reporting interval.
    struct PositionalChainStats
    {
        std::vector<ChainStats> buckets; // one ChainStats per positional bucket

        PositionalChainStats(size_t numBuckets = NUM_BUCKETS)
        {
            buckets.resize(numBuckets);
        }

        void reset()
        {
            for (auto &b : buckets)
                b = ChainStats{};
        }

        // Accumulate posStats[b] into buckets[b] for each bucket b.
        void add(const std::vector<SequenceStats> &posStats)
        {
            for (size_t b = 0; b < posStats.size() && b < buckets.size(); ++b)
            {
                const auto &s = posStats[b];
                for (size_t m = 0; m < buckets[b].sequenceLengths.size(); ++m)
                    buckets[b].sequenceLengths[m] += s.sequences[m];
                // chainLength tracks number of chains contributing to this bucket
                buckets[b].chainLength.add(1.0);
            }
        }
    };
}
