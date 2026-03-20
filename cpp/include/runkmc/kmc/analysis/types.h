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
        std::vector<MomentAccumulator> stats; // one per monomer

        SequenceStats() { stats.resize(registry::getNumMonomers()); }

        SequenceStats &operator+=(const SequenceStats &other)
        {
            for (size_t i = 0; i < stats.size(); ++i)
                stats[i] += other.stats[i];
            return *this;
        }

        void addSequence(SpeciesID id, size_t length)
        {
            stats[registry::getMonomerIndex(id)].add(static_cast<double>(length));
        }

        double totalMonomerSum() const
        {
            double total = 0;
            for (const auto &s : stats)
                total += s.sum;
            return total;
        }

        double nAvgComp(size_t monomerIdx) const { return utils::safeDivide(stats[monomerIdx].sum, totalMonomerSum()); }
    };

    struct PositionalSequenceStats
    {
        std::vector<SequenceStats> buckets;

        PositionalSequenceStats() { buckets.resize(NUM_BUCKETS); }

        PositionalSequenceStats &operator+=(const PositionalSequenceStats &other)
        {
            for (size_t b = 0; b < buckets.size() && b < other.buckets.size(); ++b)
                buckets[b] += other.buckets[b];
            return *this;
        }

        SequenceStats collapse() const
        {
            SequenceStats result;
            for (const auto &b : buckets)
                result += b;
            return result;
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
        MomentAccumulator chainLength;
        MomentAccumulator chainMW;
        SequenceStats sequenceLengths;

        ChainStats() {}

        size_t numChains() const { return static_cast<size_t>(chainLength.count); }

        ChainStats &operator+=(const ChainStats &other)
        {
            chainLength += other.chainLength;
            chainMW += other.chainMW;
            sequenceLengths += other.sequenceLengths;
            return *this;
        }
    };

}
