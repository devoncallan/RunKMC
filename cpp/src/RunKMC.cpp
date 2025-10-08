// #include "kmc/builder.h"
#include "kmc/builder.h"
#include "outputs/metadata.h"

int main(int argc, char **argv)
{
    auto config = builder::parseArguments(argc, argv);

    auto model = builder::fromModelFile(config);

    output::writeMetadata(model);
    // output::

    model.run();

    return EXIT_SUCCESS;
}