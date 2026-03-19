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

        void writeState(std::ostream &out) const
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
        PositionalStatsWriter(const analysis::PositionalChainStats &posStats,
                              const KMCState &kmcState)
            : posStats(posStats), kmcState(kmcState) {}

        static std::vector<std::string> getColumnNames()
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
            return cols;
        }

        void writeBlock(std::ostream &out) const
        {
            if (registry::getNumMonomers() <= 1)
                return;

            const auto numMonomers = registry::getNumMonomers();
            const auto &buckets = posStats.buckets;
            if (buckets.empty())
                return;

            writeBlockHeader(out, kmcState, buckets.size());
            out << str::join(getColumnNames(), ",") << std::endl;

            // One row per bucket
            for (size_t b = 0; b < buckets.size(); ++b)
            {
                const auto &cs = buckets[b];
                std::vector<std::string> row;
                row.push_back(std::to_string(b));

                double totalMonSum = 0;
                for (size_t m = 0; m < numMonomers; ++m)
                    totalMonSum += cs.sequenceLengths[m].sum;

                for (size_t m = 0; m < numMonomers; ++m)
                    row.push_back(std::to_string(analysis::safeDivide(cs.sequenceLengths[m].sum, totalMonSum)));
                for (size_t m = 0; m < numMonomers; ++m)
                    row.push_back(std::to_string(cs.sequenceLengths[m].nAvg()));
                for (size_t m = 0; m < numMonomers; ++m)
                    row.push_back(std::to_string(cs.sequenceLengths[m].wAvg()));
                for (size_t m = 0; m < numMonomers; ++m)
                    row.push_back(std::to_string(cs.sequenceLengths[m].disp()));

                out << str::join(row, ",") << std::endl;
            }
        }

    private:
        const analysis::PositionalChainStats &posStats;
        const KMCState &kmcState;
    };

    class ChainRecordWriter
    {
    public:
        ChainRecordWriter(const std::vector<Polymer *> &polymers)
            : polymers(polymers) {}

        static std::vector<std::string> getColumnNames()
        {
            std::vector<std::string> cols;
            auto monomerNames = registry::getMonomerNames();

            if (monomerNames.size() == 1)
            {
                cols.push_back("ChainLength");
            }
            else
            {
                for (const auto &name : monomerNames)
                    cols.push_back(std::string(C::state::MONCOUNT_PREFIX) + name);
                for (const auto &name : monomerNames)
                    cols.push_back(std::string(C::state::SEQCOUNT_PREFIX) + name);
                for (const auto &name : monomerNames)
                    cols.push_back(std::string(C::state::SEQLEN2_PREFIX) + name);
            }
            return cols;
        }

        // Write a complete block: header + column names + records
        void writeBlock(std::ostream &out, const KMCState &kmc) const
        {
            writeBlockHeader(out, kmc, polymers.size());
            out << str::join(getColumnNames(), ",") << std::endl;

            writeRecords(out);
        }

    private:
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

                const auto &posStats = polymer->getPositionalStats();
                std::vector<uint64_t> monCounts(numMonomers, 0);
                for (const auto &stats : posStats)
                {
                    for (size_t i = 0; i < numMonomers; ++i)
                        monCounts[i] += stats.monCounts[i];
                }
                for (auto count : monCounts)
                    row.push_back(std::to_string(count));

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

        const std::vector<Polymer *> &polymers;
    };

    class SegmentHistogramWriter
    {
    public:
        SegmentHistogramWriter(const PolymerContainer &container, const KMCState &kmcState)
            : container(container), kmcState(kmcState) {}

        void writeBlock(std::ostream &out) const
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

            writeBlockHeader(out, kmcState, histogram.size());

            // Column header
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
        }

    private:
        const PolymerContainer &container;
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
        writer.writeState(resultsFile);
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
        auto segFile = std::ofstream(paths.segmentHistFile(container.name), std::ios::app);
        SegmentHistogramWriter writer(container, kmc);
        writer.writeBlock(segFile);
    }

    void writePosChainStats(const analysis::PositionalChainStats &posStats,
                            const std::string &containerName,
                            const KMCState &kmc,
                            const SimulationPaths &paths)
    {
        auto file = std::ofstream(paths.posChainFile(containerName), std::ios::app);
        PositionalStatsWriter writer(posStats, kmc);
        writer.writeBlock(file);
    }

};
