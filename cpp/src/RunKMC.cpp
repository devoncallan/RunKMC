#include "build.h"
#include "results/metadata.h"

int main(int argc, char **argv)
{
    auto config = build::parseArguments(argc, argv);
    auto input = build::parseModelFile(config.inputFilepath);
    auto model = build::buildModel(config, input);

    output::writeMetadata(model);
    output::writeRegistry(model);
    output::writeInputFile(model, input);

    model.run();

    return EXIT_SUCCESS;
}