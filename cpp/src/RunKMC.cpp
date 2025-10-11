#include "build.h"
#include "results/metadata.h"

int main(int argc, char **argv)
{
    auto config = build::parseArguments(argc, argv);
    auto input = build::parseModelFile(config.inputFilepath);

    SimulationPaths paths(config);
    output::writeInputFile(input, paths.parsedInputFile());

    auto model = build::buildModel(config, input);

    output::writeSpeciesRegistry(paths.speciesFile());

    output::writeMetadata(model);

    if (config.parseOnly)
        return EXIT_SUCCESS;

    model.run();

    return EXIT_SUCCESS;
}