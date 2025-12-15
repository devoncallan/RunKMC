#pragma once
#include <Eigen/Core>
#include <functional>
#include <unordered_map>
#include <vector>

#include "common.h"

namespace analysis
{
    struct SequenceStats
    {
        std::vector<uint64_t> monCounts;
        std::vector<uint64_t> seqCounts;
        std::vector<uint64_t> seqLengths2;

        const static size_t NUM_METRICS = 3;

        SequenceStats()
        {
            monCounts.resize(registry::getNumMonomers(), 0);
            seqCounts.resize(registry::getNumMonomers(), 0);
            seqLengths2.resize(registry::getNumMonomers(), 0);
        }

        static size_t SIZE() { return registry::getNumMonomers() * NUM_METRICS; }

        SequenceStats &operator+=(const SequenceStats &other)
        {
            for (size_t i = 0; i < registry::getNumMonomers(); ++i)
            {
                monCounts[i] += other.monCounts[i];
                seqCounts[i] += other.seqCounts[i];
                seqLengths2[i] += other.seqLengths2[i];
            }
            return *this;
        }

        /*
        MonCounts_A, MonCounts_B, ...,
        SeqCounts_A, SeqCounts_B, ...,
        SeqLengths2_A, SeqLengths2_B, ...
        */
        Eigen::VectorXd toEigen() const
        {
            Eigen::VectorXd result(SIZE());
            for (size_t i = 0; i < registry::getNumMonomers(); ++i)
                result(0 * registry::getNumMonomers() + i) = static_cast<double>(monCounts[i]);
            for (size_t i = 0; i < registry::getNumMonomers(); ++i)
                result(1 * registry::getNumMonomers() + i) = static_cast<double>(seqCounts[i]);
            for (size_t i = 0; i < registry::getNumMonomers(); ++i)
                result(2 * registry::getNumMonomers() + i) = static_cast<double>(seqLengths2[i]);
            return result;
        }

        void addSequence(SpeciesID id, size_t length)
        {
            size_t monIdx = registry::getMonomerIndex(id);
            monCounts[monIdx] += length;
            seqCounts[monIdx] += 1;
            seqLengths2[monIdx] += length * length;
        }
    };

    struct SequenceSummary
    {
        Eigen::MatrixXd sequenceStatsMatrix;        // (polymers x (SequenceStats))
        std::vector<SequenceStats> positionalStats; // (buckets x (monomers*fields))
    };

    struct RawSequenceData
    {
        std::vector<std::vector<SpeciesID>> sequences;
        std::vector<std::vector<SequenceStats>> precomputedStats;
        size_t length;

        RawSequenceData(size_t n)
        {
            sequences.reserve(n);
            precomputedStats.reserve(n);
            length = n;
        }
    };

    struct MonomerCountKey
    {
        std::vector<uint64_t> counts;

        MonomerCountKey() = default;
        explicit MonomerCountKey(std::vector<uint64_t> values) : counts(std::move(values)) {}

        bool operator==(const MonomerCountKey &other) const noexcept
        {
            return counts == other.counts;
        }
    };

    struct MonomerCountKeyHasher
    {
        std::size_t operator()(const MonomerCountKey &key) const noexcept
        {
            std::size_t seed = key.counts.size();
            for (auto value : key.counts)
            {
                seed ^= std::hash<uint64_t>{}(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    struct ChainHistogramBin
    {
        uint64_t chainCount = 0;
        SequenceStats aggregatedStats;
    };

    using ChainHistogramMap = std::unordered_map<MonomerCountKey, ChainHistogramBin, MonomerCountKeyHasher>;

    struct ChainHistogram
    {
        ChainHistogramMap bins;

        void clear() { bins.clear(); }

        void addChain(const std::vector<uint64_t> &counts, const SequenceStats &stats)
        {
            MonomerCountKey key{counts};
            auto &bin = bins[key];
            if (bin.chainCount == 0)
            {
                bin.aggregatedStats = SequenceStats();
            }

            ++bin.chainCount;
            auto &aggregated = bin.aggregatedStats;
            const auto numMonomers = registry::getNumMonomers();
            for (size_t i = 0; i < numMonomers; ++i)
            {
                aggregated.seqCounts[i] += stats.seqCounts[i];
                aggregated.seqLengths2[i] += stats.seqLengths2[i];
            }
        }
    };
}
