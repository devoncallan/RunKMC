#pragma once
#include "common.h"
#include "io/types.h"

namespace io::cli
{
    // Forward declarations of helper functions
    static bool validateInputFile(const std::string &filepath);
    static bool prepareOutputDir(const std::string &dirPath);

    static types::CommandLineConfig parseArguments(int argc, char **argv)
    {
        if (argc < 3)
        {
            std::cerr
                << "Usage: " << argv[0]
                << " <inputFilePath> <outputDirectory>"
                << " [--report-polymers] [--report-sequences]\n";
            exit(EXIT_FAILURE);
        }

        types::CommandLineConfig config;
        config.inputFilepath = argv[1];
        config.outputDir = argv[2];

        // Parse optional flags
        for (int i = 3; i < argc; ++i)
        {
            std::string arg = argv[i];

            if (arg == "--report-polymers")
                config.reportPolymers = true;
            else if (arg == "--report-sequences")
                config.reportSequences = true;
            else
            {
                std::cerr << "Unknown argument: " << arg << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        if (!validateInputFile(config.inputFilepath))
            exit(EXIT_FAILURE);

        if (!prepareOutputDir(config.outputDir))
            exit(EXIT_FAILURE);

        return config;
    }

    static bool validateInputFile(const std::string &filepath)
    {
        std::ifstream file(filepath);
        if (!file.is_open())
        {
            std::cerr << "Cannot open input file: " << filepath << std::endl;
            return false;
        }
        file.close();
        return true;
    }

    static bool prepareOutputDir(const std::string &dirPath)
    {
        try
        {
            std::filesystem::create_directories(dirPath);
            return true;
        }
        catch (const std::filesystem::filesystem_error &e)
        {
            std::cerr << "Failed to create output directory: " << dirPath << "\n"
                      << e.what() << std::endl;
            return false;
        }
    }
};
