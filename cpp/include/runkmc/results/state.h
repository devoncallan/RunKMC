#pragma once
#include "common.h"
#include "kmc/state.h"
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
            out << str::join(KMCState::getTitles(), ",", true);
            out << str::join(SpeciesState::getTitles(), ",", true);
            out << str::join(AnalysisState::getTitles(), ",") << std::endl;
        }

        void writeState(std::ostream &out) const
        {
            out << str::join(kmcState.getDataAsVector(), ",", true);
            out << str::join(speciesState.getDataAsVector(), ",", true);
            out << str::join(analysisState.getDataAsVector(), ",") << std::endl;
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

        static void writeHeader(std::ostream &out)
        {
            out << str::join(SequenceState::getTitles(), ",") << std::endl;
        }

        void writeState(std::ostream &out) const
        {
            for (size_t bucket = 0; bucket < sequenceState.stats.size(); ++bucket)
                out << str::join(sequenceState.getDataAsVector(bucket), ",") << std::endl;
        }

    private:
        const SequenceState &sequenceState;
    };

    void writeStateHeaders(const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        console::debug("Writing results to " + paths.resultsFile().string());

        auto resultsFile = std::ofstream(paths.resultsFile());
        ResultsWriter::writeHeader(resultsFile);

        if (config.reportSequences)
        {
            auto sequenceFile = std::ofstream(paths.sequencesFile());
            SequenceWriter::writeHeader(sequenceFile);
        }
    }

    void writeState(const SystemState &state, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
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
