#pragma once
#include <Eigen/Core>
#include <functional>
#include <map>
#include <unordered_map>
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
            segmentHist[monIdx][static_cast<uint32_t>(length)] += 1;
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
