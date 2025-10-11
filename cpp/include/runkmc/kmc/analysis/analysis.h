#pragma once
#include "common.h"
#include "kmc/state.h"

namespace analysis
{

    namespace utils
    {
        static size_t getBucketIndex(size_t position, size_t chainLength, size_t numBuckets);
        template <typename Func>
        void forEachStats(const RawSequenceData &sequenceData, size_t numBuckets, Func callback);
    }

    SequenceSummary calculateSequenceSummary(const analysis::RawSequenceData &sequenceData)
    {
        // Calculate sequence stats matrix (polymers x (monomers*fields)) -> Summed across all buckets
        // Calculate positional average stats (buckets x (monomers*fields)) -> Summed across all polymers
        Eigen::MatrixXd sequenceStatsMatrix = Eigen::MatrixXd::Zero(sequenceData.length, SequenceStats::SIZE());
        std::vector<SequenceStats> positionalStats(NUM_BUCKETS);

        utils::forEachStats(
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

            // Skip non-monomer units
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

        // Add the stats for the last sequence
        size_t bucket = utils::getBucketIndex(sequence.size() - 1, sequence.size(), numBuckets);
        if (currentSequenceLength > 0)
            stats[bucket].addSequence(currentMonomerID, currentSequenceLength);

        return stats;
    }

    void analyzeChainLengthDist(Eigen::MatrixXd &sequenceStatsMatrix, const std::vector<double> &monomerFWs, AnalysisState &state)
    {
        if (sequenceStatsMatrix.rows() == 0 || sequenceStatsMatrix.cols() == 0)
            return;

        // SequenceStatsMatrix: (numPolymers x (A Count, B Count, ..., A SeqCount, B SeqCount, ..., A SeqLen2, B SeqLen2, ...))
        // Extract monomer count distribution (shape: numPolymers x numMonomers)
        Eigen::MatrixXd monomerCountDist = sequenceStatsMatrix.leftCols(registry::getNumMonomers());

        // Chain length calculations (Get chain lengths by summing monomer counts)
        Eigen::VectorXd chainLengths = monomerCountDist.rowwise().sum();
        state.nAvgCL = chainLengths.mean();
        if (state.nAvgCL != 0.0)
        {
            state.wAvgCL = chainLengths.array().square().mean() / state.nAvgCL;
            state.dispCL = state.wAvgCL / state.nAvgCL;
        }

        // If any monomer has FW of 0, skip molecular weight calculations
        // and set molecular weight averages to chain length averages
        Eigen::VectorXd FWs = Eigen::Map<const Eigen::VectorXd>(monomerFWs.data(), monomerFWs.size());
        if ((FWs.array() == 0.0).any())
        {
            state.nAvgMW = state.nAvgCL;
            state.wAvgMW = state.wAvgCL;
            state.dispMW = state.dispCL;
            return;
        }

        // Molecular weight calculations
        Eigen::MatrixXd weightedMonomerCountDist = monomerCountDist.array().rowwise() * FWs.transpose().array();
        Eigen::VectorXd molecularWeights = weightedMonomerCountDist.rowwise().sum();
        state.nAvgMW = molecularWeights.mean();
        if (state.nAvgMW != 0.0)
        {
            state.wAvgMW = molecularWeights.array().square().mean() / state.nAvgMW;
            state.dispMW = state.wAvgMW / state.nAvgMW;
        }
    }

    // SequenceStatsMatrix: (numPolymers x (A Count, B Count, ..., A SeqCount, B SeqCount, ..., A SeqLen2, B SeqLen2, ...))
    void analyzeSequenceLengthDist(Eigen::MatrixXd &sequenceStatsMatrix, AnalysisState &state)
    {

        if (sequenceStatsMatrix.rows() == 0 || sequenceStatsMatrix.cols() < SequenceStats::SIZE())
            return;

        // Sum stats over all polymers -> (1 x (Total A Count, B Count, ..., A SeqCount, B SeqCount, ..., A SeqLen2, B SeqLen2, ...))
        auto totalStats = sequenceStatsMatrix.colwise().sum();
        auto totalMonomerCounts = totalStats.leftCols(registry::getNumMonomers()).sum(); // Total chain length (i.e. total monomer count across all types)

        auto numMonomers = registry::getNumMonomers();
        for (size_t i = 0; i < numMonomers; ++i)
        {
            auto monomerCounts = totalStats(0 * numMonomers + i);    // Total count of monomer type i across all polymers
            auto sequenceCounts = totalStats(1 * numMonomers + i);   // Total number of sequences of monomer type i across all polymers
            auto sequenceLengths2 = totalStats(2 * numMonomers + i); // Sum of squared sequence lengths of monomer type i across all polymers

            if (sequenceCounts > 0 && monomerCounts > 0)
            {
                state.nAvgComp[i] = monomerCounts / totalMonomerCounts; // Total count of monomer i / total count of all monomers
                state.nAvgSL[i] = monomerCounts / sequenceCounts;       // Total count of monomer i / total number of sequences of monomer i
                state.wAvgSL[i] = sequenceLengths2 / monomerCounts;     // Sum of squared sequence lengths of monomer i / total count of monomer i
                state.dispSL[i] = state.wAvgSL[i] / state.nAvgSL[i];    // Sequence length dispersity of monomer i
            }
        }
    }

};

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
}