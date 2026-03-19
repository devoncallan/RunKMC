#pragma once
#include <filesystem>

#include "common.h"
#include "io/types.h"

class SimulationPaths
{
    std::filesystem::path _baseDir;
    std::filesystem::path _inputFilepath;

public:
    SimulationPaths() = default;

    SimulationPaths(const io::types::CommandLineConfig &config)
    {
        _baseDir = std::filesystem::path(config.outputDir);
        _inputFilepath = std::filesystem::path(config.inputFilepath);

        if (!std::filesystem::exists(_inputFilepath))
            console::error("Input file does not exist: " + _inputFilepath.string());

        if (!std::filesystem::exists(_baseDir))
            std::filesystem::create_directories(_baseDir);

        if (config.reportPolymers && !std::filesystem::exists(polymerFile()))
            std::filesystem::create_directories(polymerFile().parent_path());

        // Copy input file to output directory for record-keeping
        if (!std::filesystem::exists(sourceInputFile()))
            std::filesystem::copy(config.inputFilepath, localInputFile());
    }

    std::filesystem::path baseDirectory() const { return _baseDir; }

    // Return the original input file path the user provided (may be outside output dir)
    std::filesystem::path sourceInputFile() const { return _inputFilepath; }
    std::string localInputFile() const { return _baseDir / _inputFilepath.filename().string(); }

    // Processed/output files inside the base output directory
    std::string parsedInputFile() const { return _baseDir / C::paths::PARSED_INPUT_FILE; }
    std::filesystem::path speciesFile() const { return _baseDir / C::paths::SPECIES_FILE; }
    std::filesystem::path resultsFile() const { return _baseDir / C::paths::RESULTS_FILE; }
    std::filesystem::path polymerFile() const { return _baseDir / C::paths::POLYMERS_FILE; }

    // Per-container output files
    std::filesystem::path chainsFile(const std::string &name) const
    {
        return _baseDir / (std::string(C::paths::CHAINS_PREFIX) + name + ".dat");
    }
    std::filesystem::path posChainFile(const std::string &name) const
    {
        return _baseDir / (std::string(C::paths::POS_CHAIN_PREFIX) + name + ".dat");
    }
    std::filesystem::path segmentHistFile(const std::string &name) const
    {
        return _baseDir / (std::string(C::paths::SEGMENT_HIST_PREFIX) + name + ".dat");
    }
};
