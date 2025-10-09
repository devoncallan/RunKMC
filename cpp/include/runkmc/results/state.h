#pragma once
#include "common.h"
#include "core/state.h"
#include "analysis/types.h"
#include "results/paths.h"

namespace output
{
    class ResultsWriter
    {
    public:
        ResultsWriter(const KMCState &kmc, const SpeciesState &species, const AnalysisState &analysis)
            : kmcState(kmc), speciesState(species), analysisState(analysis) {}

        static void writeHeader(std::ostream &out)
        {
            auto kmcHeaders = KMCState::getTitles();
            auto speciesHeaders = SpeciesState::getTitles();
            auto analysisHeaders = AnalysisState::getTitles();

            std::string headerRow = "";
            for (const auto &header : kmcHeaders)
                headerRow += header + ",";
            for (const auto &header : speciesHeaders)
                headerRow += header + ",";
            for (const auto &header : analysisHeaders)
                headerRow += header + ",";

            headerRow.pop_back(); // remove last comma
            out << headerRow << std::endl;
        }

        void writeState(std::ostream &out) const
        {
            auto kmcData = kmcState.getDataAsVector();
            auto speciesData = speciesState.getDataAsVector();
            auto analysisData = analysisState.getDataAsVector();

            std::string row = "";
            for (const auto &data : kmcData)
                row += data + ",";
            for (const auto &data : speciesData)
                row += data + ",";
            for (const auto &data : analysisData)
                row += data + ",";

            row.pop_back(); // remove last comma
            out << row << std::endl;
        }

    private:
        const KMCState &kmcState;
        const SpeciesState &speciesState;
        const AnalysisState &analysisState;
    };

    class SequenceWriter
    {
    public:
        SequenceWriter(const SequenceState &seq) : sequenceState(seq) {}

        void writeState(std::ostream &out) const
        {
            for (size_t bucket = 0; bucket < sequenceState.stats.size(); ++bucket)
            {
                auto data = sequenceState.getDataAsVector(bucket);

                std::string row = "";
                for (const auto &d : data)
                    row += d + ",";

                row.pop_back(); // remove last comma
                out << row << std::endl;
            }
        }

        static void writeHeader(std::ostream &out)
        {
            auto headers = SequenceState::getTitles();

            std::string headerRow = "";
            for (const auto &header : headers)
                headerRow += header + ",";

            headerRow.pop_back(); // remove last comma
            out << headerRow << std::endl;
        }

    private:
        const SequenceState &sequenceState;
    };

    void writeStateHeaders(const SimulationPaths &paths, const config::CommandLineConfig &config)
    {

        console::debug("Preparing to write state headers to " + paths.resultsFile().string());

        auto resultsFile = std::ofstream(paths.resultsFile());

        console::debug("Writing state headers to " + paths.resultsFile().string());

        ResultsWriter::writeHeader(resultsFile);

        console::debug("Finished writing state headers.");

        if (config.reportSequences)
        {
            auto sequenceFile = std::ofstream(paths.sequencesFile());
            SequenceWriter::writeHeader(sequenceFile);
        }

        console::debug("Finished writing sequence headers.");
    }

    void writeState(const SystemState &state, const SimulationPaths &paths, const config::CommandLineConfig &config)
    {
        auto resultsFile = std::ofstream(paths.resultsFile(), std::ios::app);
        ResultsWriter writer(state.kmc, state.species, state.analysis);
        writer.writeState(resultsFile);

        if (config.reportSequences)
        {
            auto sequenceFile = std::ofstream(paths.sequencesFile(), std::ios::app);
            SequenceWriter seqWriter(state.sequence);
            seqWriter.writeState(sequenceFile);
        }
    }
};
