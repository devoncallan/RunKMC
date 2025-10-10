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

        if (config.reportSequences && !std::filesystem::exists(sequencesFile()))
            std::filesystem::create_directories(sequencesFile().parent_path());

        // Copy input file to output directory for record-keeping
        if (!std::filesystem::exists(sourceInputFile()))
            std::filesystem::copy(config.inputFilepath, localInputFile());
    }

    std::filesystem::path baseDirectory() const { return _baseDir; }

    // Return the original input file path the user provided (may be outside output dir)
    std::filesystem::path sourceInputFile() const { return _inputFilepath; }
    std::string localInputFile() const { return _baseDir / _inputFilepath.filename().string(); }
    std::string parsedInputFile() const { return _baseDir / (_inputFilepath.stem().string() + std::string(C::paths::PARSED_INPUT_SUFFIX)); }

    // Processed/output files inside the base output directory
    std::filesystem::path speciesFile() const { return _baseDir / C::paths::SPECIES_FILE; }
    std::filesystem::path resultsFile() const { return _baseDir / C::paths::RESULTS_FILE; }
    std::filesystem::path polymerFile() const { return _baseDir / C::paths::POLYMERS_FILE; }
    std::filesystem::path sequencesFile() const { return _baseDir / C::paths::SEQUENCES_FILE; }

    std::filesystem::path metadataFile() const { return _baseDir / "metadata.yaml"; }
    // std::filesystem::path registryFile() const { return _baseDir / "registry.yaml"; }
    // std::filesystem::path inputFileYaml() const { return _baseDir / "input.yaml"; }
};