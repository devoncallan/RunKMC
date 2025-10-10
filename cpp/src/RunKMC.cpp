#include "build.h"
#include "results/metadata.h"

int main(int argc, char **argv)
{
    auto config = build::parseArguments(argc, argv);
    auto input = build::parseModelFile(config.inputFilepath);
    auto model = build::buildModel(config, input);

    output::writeInputFile(model, input);
    output::writeRegistry(model);
    output::writeMetadata(model);

    if (config.parseOnly)
        return EXIT_SUCCESS;

    model.run();

    return EXIT_SUCCESS;
}