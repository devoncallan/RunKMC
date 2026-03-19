#pragma once
#include <map>
#include <vector>

#include "common.h"

namespace analysis
{

    double safeDivide(double num, double denom)
    {
        if (denom == 0)
            return 0;
        return num / denom;
    };

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

        double nAvg() const { return safeDivide(sum, count); }

        double wAvg() const { return safeDivide(sumSq, sum); }

        double disp() const { return safeDivide(wAvg(), nAvg()); }
    };

    struct SequenceStats_
    {
        std::vector<MomentAccumulator> sequences;

        SequenceStats_() { sequences.resize(registry::getNumMonomers()); }

        SequenceStats_ &operator+=(const SequenceStats_ &other)
        {
            for (size_t i = 0; i < registry::getNumMonomers(); ++i)
                sequences[i] += other.sequences[i];
            return *this;
        }

        void addSequence(SpeciesID id, size_t length)
        {
            size_t monIdx = registry::getMonomerIndex(id);
            sequences[monIdx].add(length);
        }
    };

    struct SequenceStats
    {
        std::vector<uint64_t> monCounts;
        std::vector<uint64_t> seqCounts;
        std::vector<uint64_t> seqLengths2;
        std::vector<std::map<uint32_t, uint64_t>> segmentHist;

        const static size_t NUM_METRICS = 3;

        SequenceStats()
        {
            monCounts.resize(registry::getNumMonomers(), 0);
            seqCounts.resize(registry::getNumMonomers(), 0);
            seqLengths2.resize(registry::getNumMonomers(), 0);
            segmentHist.resize(registry::getNumMonomers());
        }

        static size_t SIZE() { return registry::getNumMonomers() * NUM_METRICS; }

        SequenceStats &operator+=(const SequenceStats &other)
        {
            for (size_t i = 0; i < registry::getNumMonomers(); ++i)
            {
                monCounts[i] += other.monCounts[i];
                seqCounts[i] += other.seqCounts[i];
                seqLengths2[i] += other.seqLengths2[i];
                auto &destHist = segmentHist[i];
                const auto &srcHist = other.segmentHist[i];
                for (const auto &[length, count] : srcHist)
                    destHist[length] += count;
            }
            return *this;
        }

        void addSequence(SpeciesID id, size_t length)
        {
            size_t monIdx = registry::getMonomerIndex(id);
            monCounts[monIdx] += length;
            seqCounts[monIdx] += 1;
            seqLengths2[monIdx] += length * length;
            segmentHist[monIdx][static_cast<uint32_t>(length)] += 1;
        }
    };

    struct ChainStats
    {
        size_t numChains = 0;
        MomentAccumulator chainLength;
        MomentAccumulator chainMW;
        std::vector<MomentAccumulator> sequenceLengths; // one per monomer: run-length moments

        ChainStats() { sequenceLengths.resize(registry::getNumMonomers()); }

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
                            mw += static_cast<double>(stats.monCounts[i]) * FWs[i];
                chainMW.add(mw);
            }

            for (const auto &stats : posStats)
                for (size_t i = 0; i < sequenceLengths.size(); ++i)
                    sequenceLengths[i] += MomentAccumulator{
                        static_cast<double>(stats.seqCounts[i]),
                        static_cast<double>(stats.monCounts[i]),
                        static_cast<double>(stats.seqLengths2[i])};
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
        // posStats is indexed by bucket, each element has seqCounts, monCounts, seqLengths2 per monomer.
        void add(const std::vector<SequenceStats> &posStats)
        {
            for (size_t b = 0; b < posStats.size() && b < buckets.size(); ++b)
            {
                const auto &s = posStats[b];
                for (size_t m = 0; m < buckets[b].sequenceLengths.size(); ++m)
                    buckets[b].sequenceLengths[m] += MomentAccumulator{
                        static_cast<double>(s.seqCounts[m]),
                        static_cast<double>(s.monCounts[m]),
                        static_cast<double>(s.seqLengths2[m])};
                // chainLength tracks number of chains contributing to this bucket
                buckets[b].chainLength.add(1.0);
            }
        }
    };
}
