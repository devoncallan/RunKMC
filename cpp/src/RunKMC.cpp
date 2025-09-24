#include "kmc/builder.h"
#include "outputs/metadata.h"

int main(int argc, char **argv)
{
    auto config = KMCBuilder::parseArguments(argc, argv);

    auto model = KMCBuilder::fromFile(config);

    output::writeMetadata(model);

    model.run();

    return EXIT_SUCCESS;
}