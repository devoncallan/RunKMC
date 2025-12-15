#pragma once

#include <Eigen/Dense>

#include "core/species.h"
#include "kmc/analysis/types.h"

namespace analysis
{
ChainHistogram buildHistogramFromSequenceStats(const Eigen::MatrixXd &sequenceStatsMatrix, std::size_t numMonomers);
ChainHistogram buildHistogramFromChainLengths(const std::vector<uint64_t> &chainLengths, SpeciesID monomerID = INVALID_SPECIES_ID);

} // namespace analysis


#include <cmath>

namespace analysis
{

inline ChainHistogram buildHistogramFromSequenceStats(const Eigen::MatrixXd &sequenceStatsMatrix, std::size_t numMonomers)
{
    ChainHistogram histogram;

    for (Eigen::Index row = 0; row < sequenceStatsMatrix.rows(); ++row)
    {
        std::vector<uint64_t> counts(numMonomers, 0);
        SequenceStats stats;

        for (std::size_t col = 0; col < numMonomers; ++col)
        {
            const auto monCount = static_cast<uint64_t>(std::llround(sequenceStatsMatrix(row, static_cast<Eigen::Index>(col))));
            const auto seqCount = static_cast<uint64_t>(std::llround(sequenceStatsMatrix(row, static_cast<Eigen::Index>(numMonomers + col))));
            const auto seqLen2 = static_cast<uint64_t>(std::llround(sequenceStatsMatrix(row, static_cast<Eigen::Index>(2 * numMonomers + col))));

            counts[col] = monCount;
            stats.seqCounts[col] = seqCount;
            stats.seqLengths2[col] = seqLen2;
        }

        histogram.addChain(counts, stats);
    }

    return histogram;
}

inline ChainHistogram buildHistogramFromChainLengths(const std::vector<uint64_t> &chainLengths, SpeciesID monomerID)
{
    ChainHistogram histogram;

    const bool hasValidMonomer = (monomerID != INVALID_SPECIES_ID) && (registry::getNumMonomers() > 0);
    const std::size_t targetIndex = hasValidMonomer ? registry::getMonomerIndex(monomerID) : 0;

    for (const auto length : chainLengths)
    {
        if (hasValidMonomer)
        {
            std::vector<uint64_t> counts(registry::getNumMonomers(), 0);
            counts[targetIndex] = length;

            SequenceStats stats;
            stats.seqCounts[targetIndex] = length > 0 ? 1 : 0;
            stats.seqLengths2[targetIndex] = length * length;

            histogram.addChain(counts, stats);
        }
        else
        {
            histogram.addChain({length}, SequenceStats());
        }
    }

    return histogram;
}

} // namespace analysis
