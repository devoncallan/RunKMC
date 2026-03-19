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
        ResultsWriter(const KMCState &kmc, const SpeciesState &species, const ChainStatsState &chainStats)
            : kmcState(kmc), speciesState(species), chainStatsState(chainStats) {}

        static void writeHeader(std::ostream &out)
        {
            out << str::join(KMCState::getTitles(), ",", true);
            out << str::join(SpeciesState::getTitles(), ",", true);
            out << str::join(ChainStatsState::getTitles(), ",") << std::endl;
        }

        void writeState(std::ostream &out) const
        {
            out << str::join(kmcState.getDataAsVector(), ",", true);
            out << str::join(speciesState.getDataAsVector(), ",", true);
            out << str::join(chainStatsState.getDataAsVector(), ",") << std::endl;
        }

    private:
        const KMCState &kmcState;
        const SpeciesState &speciesState;
        const ChainStatsState &chainStatsState;
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

    class ChainRecordWriter
    {
    public:
        ChainRecordWriter(const std::vector<Polymer *> &polymers)
            : polymers(polymers) {}

        static void writeHeader(std::ostream &out)
        {
            std::vector<std::string> cols;
            auto monomerNames = registry::getMonomerNames();

            // Monomer counts only
            for (const auto &name : monomerNames)
                cols.push_back(std::string(C::state::MONCOUNT_PREFIX) + name);

            // Sequence stats (copolymer only)
            if (monomerNames.size() > 1)
            {
                for (const auto &name : monomerNames)
                    cols.push_back(std::string(C::state::SEQCOUNT_PREFIX) + name);
                for (const auto &name : monomerNames)
                    cols.push_back(std::string(C::state::SEQLEN2_PREFIX) + name);
            }

            out << str::join(cols, ",") << std::endl;
        }

        void writeRecords(std::ostream &out) const
        {
            const auto numMonomers = registry::getNumMonomers();

            if (numMonomers == 1)
            {
                for (const auto *polymer : polymers)
                    out << polymer->getDegreeOfPolymerization() << "\n";
                return;
            }

            for (const auto *polymer : polymers)
            {
                std::vector<std::string> row;

                // Extract monomer counts from compressed stats
                const auto &posStats = polymer->getPositionalStats();
                std::vector<uint64_t> monCounts(numMonomers, 0);
                for (const auto &stats : posStats)
                {
                    for (size_t i = 0; i < numMonomers; ++i)
                        monCounts[i] += stats.monCounts[i];
                }
                for (auto count : monCounts)
                    row.push_back(std::to_string(count));

                // Sequence stats (copolymer only)
                analysis::SequenceStats aggregated;
                for (const auto &stats : posStats)
                    aggregated += stats;

                for (size_t i = 0; i < numMonomers; ++i)
                    row.push_back(std::to_string(aggregated.seqCounts[i]));
                for (size_t i = 0; i < numMonomers; ++i)
                    row.push_back(std::to_string(aggregated.seqLengths2[i]));

                out << str::join(row, ",") << std::endl;
            }
        }

    private:
        const std::vector<Polymer *> &polymers;
    };

    class SegmentHistogramWriter
    {
    public:
        SegmentHistogramWriter(const PolymerContainer &container, const KMCState &kmcState)
            : container(container), kmcState(kmcState) {}

        void writeState(std::ostream &out) const
        {
            const auto monomerNames = registry::getMonomerNames();
            if (monomerNames.empty())
                return;

            std::map<uint32_t, std::vector<uint64_t>> histogram;
            const auto numMonomers = monomerNames.size();

            for (const auto *polymer : container.getPolymers())
            {
                for (const auto &stats : polymer->getPositionalStats())
                {
                    for (size_t idx = 0; idx < numMonomers && idx < stats.segmentHist.size(); ++idx)
                    {
                        for (const auto &[length, count] : stats.segmentHist[idx])
                        {
                            auto &row = histogram[length];
                            if (row.size() < numMonomers)
                                row.resize(numMonomers, 0);
                            row[idx] += count;
                        }
                    }
                }
            }

            if (histogram.empty())
                return;

            out << '#' << std::string(C::state::ITERATION_KEY) << '='
                << std::to_string(kmcState.iteration)
                << ','
                << std::string(C::state::KMC_TIME_KEY)
                << '='
                << std::to_string(kmcState.kmcTime)
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
        const PolymerContainer &container;
        const KMCState &kmcState;
    };

    void writeStateHeaders(const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        console::debug("Writing results to " + paths.resultsFile().string());

        auto resultsFile = std::ofstream(paths.resultsFile());
        ResultsWriter::writeHeader(resultsFile);

        if (config.reportChains)
        {
            console::debug("Writing chain records to " + paths.chainRecordsFile().string());
            auto chainFile = std::ofstream(paths.chainRecordsFile());
            ChainRecordWriter::writeHeader(chainFile);
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
        ResultsWriter writer(state.kmc, state.species, state.chainStats);
        writer.writeState(resultsFile);
    }

    void writeSequences(const SystemState &state, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        auto sequenceFile = std::ofstream(paths.sequencesFile(), std::ios::app);
        SequenceWriter seqWriter(state.sequence);
        seqWriter.writeState(sequenceFile);
    }

    void writeChainRecords(const PolymerContainer &container, const SimulationPaths &paths)
    {
        auto file = paths.chainRecordsFile();
        bool fileExists = std::filesystem::exists(file);
        std::ofstream out(file, std::ios::app);

        if (!fileExists)
            ChainRecordWriter::writeHeader(out);

        ChainRecordWriter writer(container.getPolymers());
        writer.writeRecords(out);
    }

    void writeSegmentHistogram(const PolymerContainer &container, const KMCState &kmc, const SimulationPaths &paths)
    {
        auto segFile = std::ofstream(paths.segmentHistFile(), std::ios::app);
        SegmentHistogramWriter writer(container, kmc);
        writer.writeState(segFile);
    }

    void writeDeadPolymers(SpeciesSet &speciesSet, analysis::ChainStats &acc, const KMCState &kmc, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        const auto &container = speciesSet.getDeadPolymerContainer();
        if (container.getPolymers().empty())
            return;

        const auto FWs = speciesSet.getMonomerFWs();
        for (const auto *polymer : container.getPolymers())
            acc.add(polymer->getDegreeOfPolymerization(), polymer->getPositionalStats(), FWs);

        if (config.reportChains)
            writeChainRecords(container, paths);
        if (config.reportSegmentHistogram)
            writeSegmentHistogram(container, kmc, paths);

        speciesSet.getDeadPolymerContainer().clearPolymers();
    }

    void writeState(const SystemState &state, const SimulationPaths &paths, const io::types::CommandLineConfig &config)
    {
        writeResults(state, paths, config);
        if (config.reportSequences)
            writeSequences(state, paths, config);
    }
};
