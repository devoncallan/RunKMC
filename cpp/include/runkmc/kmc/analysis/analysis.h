#pragma once
#include "common.h"
#include "kmc/state.h"
#include "kmc/analysis/utils.h"
#include "kmc/species/species_set.h"

namespace analysis
{
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

    void analyze(const SpeciesSet &speciesSet, SystemState &systemState)
    {
        auto sequenceData = speciesSet.getRawSequenceData();
        auto summary = analysis::calculateSequenceSummary(sequenceData);

        AnalysisState analysisState;
        analysis::analyzeChainLengthDist(summary.sequenceStatsMatrix, speciesSet.getMonomerFWs(), analysisState);
        systemState.analysis = analysisState;

        if (registry::getNumMonomers() <= 1)
            return;

        SequenceState sequenceState = SequenceState{systemState.kmc, summary.positionalStats};
        analysis::analyzeSequenceLengthDist(summary.sequenceStatsMatrix, analysisState);
        systemState.analysis = analysisState;
        systemState.sequence = sequenceState;
    }
}
