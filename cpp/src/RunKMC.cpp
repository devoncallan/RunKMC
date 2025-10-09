#include "build.h"
#include "results/metadata.h"

int main(int argc, char **argv)
{
    auto config = builder::parseArguments(argc, argv);
    auto input = builder::parseKMCInput(config.inputFilepath);

    auto model = builder::buildModel(config, input);

    output::writeMetadata(model);
    output::writeRegistry(model);
    output::writeInputFile(model, input);

    model.run();

    return EXIT_SUCCESS;
}