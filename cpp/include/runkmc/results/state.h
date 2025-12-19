#pragma once
#include <algorithm>
#include <map>

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

    class SegmentHistogramWriter
    {
    public:
        SegmentHistogramWriter(const ChainState &chain) : chainState(chain) {}

        void writeState(std::ostream &out) const
        {
            const auto monomerNames = registry::getMonomerNames();
            if (monomerNames.empty())
                return;

            std::map<uint32_t, std::vector<uint64_t>> histogram;
            const auto numMonomers = monomerNames.size();

            for (const auto &entry : chainState.histogram.bins)
            {
                const auto &segmentHist = entry.second.segmentHist;
                if (segmentHist.empty())
                    continue;

                for (size_t idx = 0; idx < numMonomers && idx < segmentHist.size(); ++idx)
                {
                    for (const auto &[length, count] : segmentHist[idx])
                    {
                        auto &row = histogram[length];
                        if (row.size() < numMonomers)
                            row.resize(numMonomers, 0);
                        row[idx] += count;
                    }
                }
            }

            if (histogram.empty())
                return;

            out << '#' << std::string(C::state::ITERATION_KEY) << '='
                << std::to_string(chainState.kmcState.iteration)
                << ','
                << std::string(C::state::KMC_TIME_KEY)
                << '='
                << std::to_string(chainState.kmcState.kmcTime)
                << ','
                << std::string(C::state::BINS_KEY)
                << '='
                << std::to_string(static_cast<uint64_t>(histogram.size()))
                << std::endl;

            std::vector<std::string> titles;
            titles.emplace_back(std::string(C::state::SEGMENT_LENGTH_KEY));
            for (const auto &name : monomerNames)
                titles.emplace_back(std::string(C::state::SEGMENT_COUNT_PREFIX) + name);
            out << str::join(titles, ",") << std::endl;

            for (const auto &[length, counts] : histogram)
            {
                std::vector<std::string> row;
                row.reserve(numMonomers + 1);
                row.push_back(std::to_string(length));
                for (size_t idx = 0; idx < numMonomers; ++idx)
                {
                    uint64_t value = (idx < counts.size()) ? counts[idx] : 0;
                    row.push_back(std::to_string(value));
                }
                out << str::join(row, ",") << std::endl;
            }

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
        if (config.reportSegmentHistogram)
        {
            console::debug("Writing segment histograms to " + paths.segmentHistFile().string());
            auto segFile = std::ofstream(paths.segmentHistFile());
            (void)segFile;
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

    void writeSegmentHistogram(const SystemState &state, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        auto segFile = std::ofstream(paths.segmentHistFile(), std::ios::app);
        SegmentHistogramWriter writer(state.chains);
        writer.writeState(segFile);
    }

    void writeState(const SystemState &state, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        writeResults(state, paths, config);
        if (config.reportChains)
            writeChainStats(state, paths, config);
        if (config.reportSegmentHistogram)
            writeSegmentHistogram(state, paths, config);
        if (config.reportSequences)
            writeSequences(state, paths, config);
    }
};
