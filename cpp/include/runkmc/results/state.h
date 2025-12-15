#pragma once
#include <algorithm>

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

    class ChainWriter
    {
    public:
        ChainWriter(const ChainState &chain) : chainState(chain) {}

        static void writeHeader(std::ostream &out) {}

        void writeState(std::ostream &out) const
        {
            std::vector<std::pair<const analysis::MonomerCountKey *, const analysis::ChainHistogramBin *>> entries;
            entries.reserve(chainState.histogram.bins.size());
            for (const auto &entry : chainState.histogram.bins)
                entries.emplace_back(&entry.first, &entry.second);

            auto bins = static_cast<uint64_t>(entries.size());
            if (bins == 0)
                return;

            auto comparator = [](const auto &lhs, const auto &rhs)
            {
                const auto &a = lhs.first->counts;
                const auto &b = rhs.first->counts;
                if (a.size() != b.size())
                    return a.size() < b.size();
                return a < b;
            };
            std::sort(entries.begin(), entries.end(), comparator);

            out << '#' << std::string(C::state::ITERATION_KEY) << '='
                << std::to_string(chainState.kmcState.iteration)
                << ','
                << std::string(C::state::KMC_TIME_KEY)
                << '='
                << std::to_string(chainState.kmcState.kmcTime)
                << ','
                << std::string(C::state::BINS_KEY)
                << '='
                << std::to_string(static_cast<uint64_t>(entries.size()))
                << std::endl;

            const auto titles = ChainState::getTitles();
            if (!titles.empty())
                out << str::join(titles, ",") << std::endl;

            for (const auto &entry : entries)
                out << str::join(chainState.buildRow(*entry.first, *entry.second), ",") << std::endl;

            out << std::endl;
        }

    private:
        const ChainState &chainState;
    };

    void writeStateHeaders(const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        console::debug("Writing results to " + paths.resultsFile().string());

        auto resultsFile = std::ofstream(paths.resultsFile());
        ResultsWriter::writeHeader(resultsFile);

        if (config.reportChains)
        {
            console::debug("Writing chain stats to " + paths.chainStatsFile().string());
            auto chainFile = std::ofstream(paths.chainStatsFile());
            ChainWriter::writeHeader(chainFile);
        }
        if (config.reportSequences)
        {
            console::debug("Writing sequences to " + paths.sequencesFile().string());
            auto sequenceFile = std::ofstream(paths.sequencesFile());
            SequenceWriter::writeHeader(sequenceFile);
        }
    }

    

    void writeResults(const SystemState &state, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        auto resultsFile = std::ofstream(paths.resultsFile(), std::ios::app);
        ResultsWriter writer(state.kmc, state.species, state.analysis);
        writer.writeState(resultsFile);
    }

    void writeSequences(const SystemState &state, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        auto sequenceFile = std::ofstream(paths.sequencesFile(), std::ios::app);
        SequenceWriter seqWriter(state.sequence);
        seqWriter.writeState(sequenceFile);
    }

    void writeChainStats(const SystemState &state, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        auto chainFile = std::ofstream(paths.chainStatsFile(), std::ios::app);
        ChainWriter chainWriter(state.chains);
        chainWriter.writeState(chainFile);
    }

    void writeState(const SystemState &state, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        writeResults(state, paths, config);
        if (config.reportChains)
            writeChainStats(state, paths, config);
        if (config.reportSequences)
            writeSequences(state, paths, config);
    }
};
