#pragma once
#include <Eigen/Core>

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
            monCounts.resize(registry::NUM_MONOMERS, 0);
            seqCounts.resize(registry::NUM_MONOMERS, 0);
            seqLengths2.resize(registry::NUM_MONOMERS, 0);
        }

        static size_t SIZE() { return registry::NUM_MONOMERS * NUM_METRICS; }

        SequenceStats &operator+=(const SequenceStats &other)
        {
            for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
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
            for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
                result(0 * registry::NUM_MONOMERS + i) = static_cast<double>(monCounts[i]);
            for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
                result(1 * registry::NUM_MONOMERS + i) = static_cast<double>(seqCounts[i]);
            for (size_t i = 0; i < registry::NUM_MONOMERS; ++i)
                result(2 * registry::NUM_MONOMERS + i) = static_cast<double>(seqLengths2[i]);
            return result;
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
}