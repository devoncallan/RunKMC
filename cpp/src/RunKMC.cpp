// #include "kmc/builder.h"
#include "kmc/builder.h"
#include "outputs/metadata.h"

int main(int argc, char **argv)
{
    auto config = builder::parseArguments(argc, argv);

    console::log("Building model from file: " + config.inputFilepath);
    console::log("Output directory: " + config.outputDir);

    auto model = builder::fromModelFile(config);

    output::writeMetadata(model);

    return EXIT_SUCCESS;

    model.run();

    return EXIT_SUCCESS;
}