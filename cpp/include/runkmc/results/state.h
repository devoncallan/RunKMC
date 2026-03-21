#pragma once
#include <algorithm>
#include <map>

#include "common.h"
#include "kmc/state.h"
#include "results/paths.h"

namespace output
{
    // Shared block header writer: #Iteration=N,KMC Time=T,Count=R
    static void writeBlockHeader(std::ostream &out, const KMCState &kmc, size_t count)
    {
        out << '#'
            << std::string(C::state::ITERATION_KEY) << '=' << std::to_string(kmc.iteration)
            << ',' << std::string(C::state::KMC_TIME_KEY) << '=' << std::to_string(kmc.kmcTime)
            << ',' << std::string(C::state::BINS_KEY) << '=' << std::to_string(count)
            << std::endl;
    }

    class ResultsWriter
    {
    public:
        ResultsWriter(const KMCState &kmc, const SpeciesState &species,
                      const ChainStatsState &chainStats)
            : kmcState(kmc), speciesState(species), chainStats(chainStats) {}

        static void writeHeader(std::ostream &out)
        {
            out << str::join(KMCState::getTitles(), ",", true);
            out << str::join(SpeciesState::getTitles(), ",", true);
            out << str::join(ChainStatsState::getTitles(), ",") << std::endl;
        }

        void writeData(std::ostream &out) const
        {
            out << str::join(kmcState.getDataAsVector(), ",", true);
            out << str::join(speciesState.getDataAsVector(), ",", true);
            out << str::join(chainStats.getDataAsVector(), ",") << std::endl;
        }

    private:
        const KMCState &kmcState;
        const SpeciesState &speciesState;
        const ChainStatsState &chainStats;
    };

    // Writes positional chain statistics for one container.
    // Format per interval: block header line, column header, one row per bucket.
    class PositionalStatsWriter
    {
    public:
        PositionalStatsWriter(const analysis::PositionalSequenceStats &posStats,
                              const KMCState &kmcState)
            : posStats(posStats), kmcState(kmcState) {}

        static void writeHeader(std::ostream &out)
        {
            std::vector<std::string> cols = {std::string(C::state::BUCKET_KEY)};
            const auto monomerNames = registry::getMonomerNames();
            for (const auto &name : monomerNames)
                cols.push_back(std::string(C::state::NAVGCOMP_PREFIX) + name);
            for (const auto &name : monomerNames)
                cols.push_back(std::string(C::state::NAVGSL_PREFIX) + name);
            for (const auto &name : monomerNames)
                cols.push_back(std::string(C::state::WAVGSL_PREFIX) + name);
            for (const auto &name : monomerNames)
                cols.push_back(std::string(C::state::DISPSL_PREFIX) + name);
            out << str::join(cols, ",") << std::endl;
        }

        void writeData(std::ostream &out) const
        {
            const auto numMonomers = registry::getNumMonomers();
            for (size_t b = 0; b < posStats.buckets.size(); ++b)
            {
                std::vector<std::string> row;

                const auto &cs = posStats.buckets[b];
                row.push_back(std::to_string(b));

                for (size_t m = 0; m < numMonomers; ++m)
                    row.push_back(std::to_string(cs.nAvgComp(m)));
                for (size_t m = 0; m < numMonomers; ++m)
                    row.push_back(std::to_string(cs.stats[m].nAvg()));
                for (size_t m = 0; m < numMonomers; ++m)
                    row.push_back(std::to_string(cs.stats[m].wAvg()));
                for (size_t m = 0; m < numMonomers; ++m)
                    row.push_back(std::to_string(cs.stats[m].disp()));

                out << str::join(row, ",") << std::endl;
            }
        }

        void writeBlock(std::ostream &out) const
        {
            if (registry::getNumMonomers() <= 1 || posStats.buckets.empty())
                return;

            writeBlockHeader(out, kmcState, posStats.buckets.size());
            writeHeader(out);
            writeData(out);
        }

    private:
        const analysis::PositionalSequenceStats &posStats;
        const KMCState &kmcState;
    };

    class ChainRecordWriter
    {
    public:
        ChainRecordWriter(const std::vector<Polymer *> &polymers)
            : polymers(polymers) {}

        static void writeHeader(std::ostream &out)
        {
            std::vector<std::string> cols;
            const auto monomerNames = registry::getMonomerNames();

            for (const auto &name : monomerNames)
                cols.push_back(std::string(C::state::MONCOUNT_PREFIX) + name);

            if (monomerNames.size() > 1)
            {
                for (const auto &name : monomerNames)
                    cols.push_back(std::string(C::state::SEQCOUNT_PREFIX) + name);
                for (const auto &name : monomerNames)
                    cols.push_back(std::string(C::state::SEQLEN2_PREFIX) + name);
            }
            out << str::join(cols, ",") << std::endl;
        }

        void writeData(std::ostream &out) const
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
                const auto record = polymer->getChainRecord();
                std::vector<std::string> row;
                for (size_t i = 0; i < numMonomers; ++i)
                    row.push_back(std::to_string(static_cast<uint64_t>(record.stats[i].sum)));
                for (size_t i = 0; i < numMonomers; ++i)
                    row.push_back(std::to_string(static_cast<uint64_t>(record.stats[i].count)));
                for (size_t i = 0; i < numMonomers; ++i)
                    row.push_back(std::to_string(static_cast<uint64_t>(record.stats[i].sumSq)));
                out << str::join(row, ",") << "\n";
            }
        }

        void writeBlock(std::ostream &out, const KMCState &kmc) const
        {
            writeBlockHeader(out, kmc, polymers.size());
            writeHeader(out);
            writeData(out);
        }

        const std::vector<Polymer *> &polymers;
    };

    class SegmentHistogramWriter
    {
    public:
        SegmentHistogramWriter(const analysis::SegmentHistogram &histogram, const KMCState &kmcState)
            : histogram(histogram), kmcState(kmcState) {}

        static void writeHeader(std::ostream &out)
        {
            const auto monomerNames = registry::getMonomerNames();
            std::vector<std::string> cols = {std::string(C::state::SEGMENT_LENGTH_KEY)};
            for (const auto &name : monomerNames)
                cols.emplace_back(std::string(C::state::SEGMENT_COUNT_PREFIX) + name);
            out << str::join(cols, ",") << std::endl;
        }

        void writeData(std::ostream &out) const
        {
            const auto numMonomers = registry::getMonomerNames().size();
            for (const auto &[length, counts] : histogram)
            {
                std::vector<std::string> row;
                row.reserve(numMonomers + 1);
                row.push_back(std::to_string(length));
                for (size_t idx = 0; idx < numMonomers; ++idx)
                    row.push_back(std::to_string(idx < counts.size() ? counts[idx] : 0));
                out << str::join(row, ",") << std::endl;
            }
        }

        void writeBlock(std::ostream &out) const
        {
            if (histogram.empty())
                return;
            writeBlockHeader(out, kmcState, histogram.size());
            writeHeader(out);
            writeData(out);
        }

    private:
        const analysis::SegmentHistogram &histogram;
        const KMCState &kmcState;
    };

    // Write results.csv header only (chain/segment/positional files need no init — all blocks are self-contained).
    void writeStateHeaders(const SimulationPaths &paths,
                           const io::types::CommandLineConfig &config)
    {
        console::debug("Writing results to " + paths.resultsFile().string());
        auto resultsFile = std::ofstream(paths.resultsFile());
        ResultsWriter::writeHeader(resultsFile);
    }

    void writeResults(const SystemState &state, const SimulationPaths &paths)
    {
        auto resultsFile = std::ofstream(paths.resultsFile(), std::ios::app);
        ResultsWriter writer(state.kmc, state.species, state.chainStats);
        writer.writeData(resultsFile);
    }

    void writeChainRecords(const PolymerContainer &container, const KMCState &kmc,
                           const SimulationPaths &paths)
    {
        auto file = paths.chainsFile(container.name);
        std::ofstream out(file, std::ios::app);
        ChainRecordWriter writer(container.getPolymers());
        writer.writeBlock(out, kmc);
    }

    void writeSegmentHistogram(const PolymerContainer &container, const KMCState &kmc,
                               const SimulationPaths &paths)
    {
        analysis::SegmentHistogram histogram;
        for (const auto *polymer : container.getPolymers())
            histogram += polymer->getSegmentHistogram();

        auto segFile = std::ofstream(paths.segmentHistFile(container.name), std::ios::app);
        SegmentHistogramWriter writer(histogram, kmc);
        writer.writeBlock(segFile);
    }

    void writePosChainStats(const analysis::PositionalSequenceStats &posStats,
                            const std::string &containerName,
                            const KMCState &kmc,
                            const SimulationPaths &paths)
    {
        auto file = std::ofstream(paths.posChainFile(containerName), std::ios::app);
        PositionalStatsWriter writer(posStats, kmc);
        writer.writeBlock(file);
    }

};
