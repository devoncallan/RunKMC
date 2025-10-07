#pragma once
#include "common.h"
#include "utils/parse.h"
#include "kmc.h"

/**
 * Parses the model file to build a valid KMC simulator.
 */
class KMCBuilder
{

public:
    static KMC fromFile(config::CommandLineConfig config)
    {
        std::string line;
        std::ifstream modelFile(config.inputFilepath);
        if (!modelFile.is_open())
            console::input_error("Cannot open model file: " + config.inputFilepath);

        std::vector<std::string> parameterLines, speciesLines, rateConstantLines, reactionLines;

        // Parse model file for required sections
        while (std::getline(modelFile, line))
        {
            str::trim(line);
            if (input::canIgnoreLine(line))
                continue;

            if (parameterLines.empty() && str::startswith(line, "parameters"))
            {
                parameterLines = input::parseSection("parameters", modelFile);
            }
            else if (speciesLines.empty() && str::startswith(line, "species"))
            {
                speciesLines = input::parseSection("species", modelFile);
                speciesLines = KMCBuilder::sortSpeciesLinesByType(speciesLines);
            }
            else if (rateConstantLines.empty() && str::startswith(line, "rateconstants"))
            {
                rateConstantLines = input::parseSection("rateconstants", modelFile);
            }
            else if (reactionLines.empty() && str::startswith(line, "reactions"))
            {
                reactionLines = input::parseSection("reactions", modelFile);
            }
        }

        // Validate that all sections were read
        if (parameterLines.empty() || speciesLines.empty() ||
            rateConstantLines.empty() || reactionLines.empty())
        {
            console::input_error("Missing required sections in model file");
        }

        // Build inputs for the KMC model from the parsed sections
        auto simConfig = buildSimulationConfig(parameterLines);
        auto speciesSet = buildSpeciesSet(speciesLines, simConfig);
        auto rateConstants = buildRateConstants(rateConstantLines);
        auto reactionSet = buildReactionSet(reactionLines, speciesSet, rateConstants);

        registry::finalizeRegistry();

        KMC model(speciesSet, reactionSet, config, simConfig);

        return model;
    }

    static config::CommandLineConfig parseArguments(int argc, char **argv)
    {
        if (argc < 3)
        {
            std::cerr
                << "Usage: " << argv[0]
                << " <inputFilePath> <outputDirectory>"
                << " [--report-polymers] [--report-sequences]\n";
            exit(EXIT_FAILURE);
        }

        config::CommandLineConfig config;
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

private:
    static config::SimulationConfig buildSimulationConfig(const std::vector<std::string> &parameterLines)
    {
        config::SimulationConfig config;
        input::readVariableRequired(parameterLines, "num_units", config.numParticles);
        input::readVariableRequired(parameterLines, "termination_time", config.terminationTime);
        input::readVariableRequired(parameterLines, "analysis_time", config.analysisTime);
        return config;
        // + more when I think of them
    }

    static SpeciesSet buildSpeciesSet(const std::vector<std::string> &speciesLines, const config::SimulationConfig &config)
    {
        size_t totalSpecies = speciesLines.size() + 1;
        std::vector<Unit> unitSpecies = {UNIT_UNDEF}; // Initialize with placeholder (IDs start at 1)

        std::vector<PolymerType> polymerTypes;
        std::vector<PolymerGroupStruct> polymerGroupStructs;

        unitSpecies.reserve(totalSpecies);
        polymerTypes.reserve(totalSpecies);
        polymerGroupStructs.reserve(totalSpecies);

        for (const auto &line : speciesLines)
        {
            std::vector<std::string> args = str::splitByWhitespace(line);
            std::string speciesType = args[0];
            std::string speciesName = args[1];
            std::vector<std::string> vars = {args.begin() + 2, args.end()};

            // Default values
            double FW = 0.;
            double C0 = 0.;

            if (speciesType == SpeciesType::UNIT)
            {
                input::readVariable(vars, "FW", FW);
                input::readVariable(vars, "[C0]", C0);

                SpeciesID id = registry::registerNewSpecies(speciesName, speciesType);
                Unit unit = Unit(speciesType, speciesName, id, C0, FW);
                unitSpecies.push_back(unit);
            }
            else if (speciesType == SpeciesType::MONOMER)
            {
                input::readVariable(vars, "FW", FW);
                input::readVariable(vars, "[C0]", C0);

                SpeciesID id = registry::registerNewSpecies(speciesName, speciesType);
                Unit monomer = Unit(speciesType, speciesName, id, C0, FW);
                unitSpecies.push_back(monomer);
            }
            else if (speciesType == SpeciesType::INITIATOR)
            {
                double efficiency;
                input::readVariable(vars, "FW", FW);
                input::readVariable(vars, "[C0]", C0);
                input::readVariableRequired(vars, "f", efficiency);

                SpeciesID id = registry::registerNewSpecies(speciesName, speciesType);
                Unit initiator = Unit(speciesType, speciesName, id, C0, FW, efficiency);
                unitSpecies.push_back(initiator);
            }
            else if (speciesType == SpeciesType::POLYMER)
            {
                std::vector<size_t> polyTypeIndices;

                SpeciesID id = registry::registerNewSpecies(speciesName, speciesType);

                // Is this doing anything?
                for (const auto &polymerGroupStruct : polymerGroupStructs)
                    if (polymerGroupStruct.name == speciesName)
                        continue;

                // Check if this polymer has arguments (subtype)
                if (!vars.empty())
                {
                    assert(vars.size() == 1);
                    // Split "P[I,A]|P[A,A]|P[B,A]" into {"P[I,A]", "P[A,A]", "P[B,A]"}
                    std::vector<std::string> subPolymerTypeNames = str::splitByDelimeter(vars[0], "|");
                    std::vector<size_t> subPolymerTypeIndices;
                    subPolymerTypeIndices.reserve(subPolymerTypeNames.size());

                    for (const auto &subPolymerTypeName : subPolymerTypeNames)
                    {
                        size_t index = input::findInVector(subPolymerTypeName, polymerTypes);
                        if (index < polymerTypes.size())
                            subPolymerTypeIndices.push_back(index);
                        else
                            console::input_error("No PolymerType associated with " + subPolymerTypeName + ". Exiting.");
                    }
                    polymerGroupStructs.push_back(PolymerGroupStruct(speciesName, subPolymerTypeIndices));
                    continue;
                }

                // Check if this polymer has sequence information.
                size_t seqStart = speciesName.find("[");
                size_t seqEnd = speciesName.find("]");
                std::vector<SpeciesID> endSequence;

                // Check if the polymer name has brackets (e.g., "P[A,A]" or "P[B,B]")
                if (seqStart != std::string::npos && seqEnd != std::string::npos && seqStart < seqEnd)
                {
                    // Get string inside the brackets (e.g., "A,A")
                    std::string sequenceString = speciesName.substr(seqStart + 1, seqEnd - 2);
                    // Split sequenceString by a comma (e.g., {"A","A"})
                    std::vector<std::string> unitStrings = str::splitByDelimeter(sequenceString, ".");
                    // Find unit string
                    for (const auto &unitString : unitStrings)
                    {
                        size_t index = input::findInVector(unitString, unitSpecies);
                        if (index < unitSpecies.size())
                            endSequence.push_back(unitSpecies[index].ID);
                        else
                        {
                            console::input_error("Unit corresponding to " + unitString + " not found. Exiting.");
                        }
                    }
                }
                polymerTypes.push_back(PolymerType(speciesName, endSequence));
                polymerGroupStructs.push_back(PolymerGroupStruct(speciesName, {polymerTypes.size() - 1}));
            }
            else
            {
                console::input_error("Species type " + speciesType + " not recognized. Exiting.");
            }
        }

        return SpeciesSet(std::move(polymerTypes), std::move(polymerGroupStructs), std::move(unitSpecies), config.numParticles);
    }

    static std::vector<RateConstant> buildRateConstants(const std::vector<std::string> &rateConstantLines)
    {
        std::vector<RateConstant> rateConstants;
        rateConstants.reserve(rateConstantLines.size());

        for (const auto &line : rateConstantLines)
        {
            std::vector<std::string> var = input::parseVariable(line);
            RateConstant rateConstant(var[0], std::stod(var[1]));
            rateConstants.push_back(rateConstant);
        }
        return rateConstants;
    }

    static ReactionSet buildReactionSet(const std::vector<std::string> &reactionLines, SpeciesSet &speciesSet, const std::vector<RateConstant> &rateConstants)
    {
        std::vector<PolymerTypeGroup> polyTypeGroups = speciesSet.getPolyTypeGroups();
        std::vector<PolymerTypeGroupPtr> polyGroupPtrs = speciesSet.getPolymerGroupPtrs();

        std::vector<Unit> &units = speciesSet.getUnits();
        std::vector<Reaction *> reactions;
        reactions.reserve(reactionLines.size());

        for (const auto &line : reactionLines)
        {
            std::vector<Unit *> unitReactants;
            unitReactants.reserve(3);
            std::vector<Unit *> unitProducts;
            unitProducts.reserve(3);
            std::vector<PolymerTypeGroupPtr> polyReactants;
            polyReactants.reserve(3);
            std::vector<PolymerTypeGroupPtr> polyProducts;
            polyProducts.reserve(3);
            RateConstant rateConstant;
            size_t index;

            // split reaction string into {"PR","P[-,A]","+","A","-kAA->","P[A,A]"}
            std::vector<std::string> splitReactionString = str::splitByWhitespace(line);

            std::string reactionType = splitReactionString[0]; // (e.g., "PR" for propagation)

            bool isReactants = true;

            // Loop through all arguments in reaction string (ignoring first argument, type)
            for (size_t i = 1; i < splitReactionString.size(); ++i)
            {

                std::string reactionArg = splitReactionString[i];
                if (str::startswith(reactionArg, "-")) // Check if string is "-rateconstant->"
                {
                    std::string rateConstantName = reactionArg.substr(1, reactionArg.size() - 3); // remove "-" and "->"
                    size_t index = input::findInVector(rateConstantName, rateConstants);
                    if (index < rateConstants.size())
                        rateConstant = rateConstants[index];
                    else
                        console::input_error("Rate constant " + rateConstantName + " not found. Exiting.");
                    isReactants = false; // switch from reactants to products
                    continue;
                }

                if (reactionArg == "+")
                    continue;

                // Check if the reaction argument is a polymer
                bool foundPoly = false;
                index = input::findInVector(reactionArg, polyTypeGroups);
                if (index < polyTypeGroups.size())
                {
                    foundPoly = true;
                    if (isReactants)
                        polyReactants.push_back(polyGroupPtrs[index]);
                    else
                        polyProducts.push_back(polyGroupPtrs[index]);
                }

                // Check if the reaction argument is a unit
                bool foundUnit = false;
                index = input::findInVector(reactionArg, units);
                if (index < units.size())
                {
                    foundUnit = true;
                    if (isReactants)
                        unitReactants.push_back(&units[index]);
                    else
                        unitProducts.push_back(&units[index]);
                }

                if (!(foundUnit || foundPoly))
                {
                    auto unitIDs = registry::getAllUnitIDs();
                    for (const auto &id : unitIDs)
                    {
                        std::cout << units[id].name << std::endl;
                    }
                    console::input_error("Species " + reactionArg + " not found. Exiting.");
                }
            }
            uint8_t sameReactant = 0;

            if (reactionType == Elementary::TYPE)
                reactions.push_back(new Elementary(rateConstant, unitReactants, unitProducts));
            else if (reactionType == InitiatorDecomposition::TYPE)
                reactions.push_back(new InitiatorDecomposition(rateConstant, unitReactants[0], unitProducts[0], unitProducts[1], unitReactants[0]->efficiency));
            else if (reactionType == InitiatorDecompositionPolymer::TYPE)
                reactions.push_back(new InitiatorDecompositionPolymer(rateConstant, unitReactants[0], polyProducts[0], polyProducts[1], unitReactants[0]->efficiency));
            else if (reactionType == Initiation::TYPE)
                reactions.push_back(new Initiation(rateConstant, unitReactants[0], unitReactants[1], polyProducts[0]));
            else if (reactionType == Propagation::TYPE)
                reactions.push_back(new Propagation(rateConstant, polyReactants[0], unitReactants[0], polyProducts[0]));
            else if (reactionType == Depropagation::TYPE)
                reactions.push_back(new Depropagation(rateConstant, polyReactants[0], polyProducts[0], unitProducts[0]));
            else if (reactionType == TerminationCombination::TYPE)
            {
                if (polyReactants[0]->name == polyReactants[1]->name)
                    sameReactant = 1;
                reactions.push_back(new TerminationCombination(rateConstant, polyReactants[0], polyReactants[1], polyProducts[0], sameReactant));
            }
            else if (reactionType == TerminationDisproportionation::TYPE)
            {
                if (polyReactants[0]->name == polyReactants[1]->name)
                    sameReactant = 1;
                reactions.push_back(new TerminationDisproportionation(rateConstant, polyReactants[0], polyReactants[1], polyProducts[0], polyProducts[1], sameReactant));
            }
            else if (reactionType == ChainTransferToMonomer::TYPE)
                reactions.push_back(new ChainTransferToMonomer(rateConstant, polyReactants[0], unitReactants[0], polyProducts[0], polyProducts[1]));
            else if (reactionType == ThermalInitiationMonomer::TYPE)
                reactions.push_back(new ThermalInitiationMonomer(rateConstant, unitReactants[0], unitReactants[1], unitReactants[2], polyProducts[0], polyProducts[1]));
            else
                console::input_error(reactionType + " is not a valid reaction type.");
        }
        ReactionSet reactionSet(reactions, rateConstants);

        return reactionSet;
    }

    static std::vector<std::string> sortSpeciesLinesByType(std::vector<std::string> &speciesLines)
    {
        std::vector<std::string> sortedLines = speciesLines;
        // Sort species lines by type priority: MONOMER, INITIATOR, UNIT, POLYMER
        // Preserves order within each type
        std::stable_sort(sortedLines.begin(), sortedLines.end(), [](const std::string &a, const std::string &b)
                         {
            auto getTypePriority = [](const std::string &line) -> int {
                std::string type = str::splitByWhitespace(line)[0];
                if (type == SpeciesType::MONOMER) return 0;
                if (type == SpeciesType::INITIATOR) return 1;
                if (type == SpeciesType::UNIT) return 2;
                if (type == SpeciesType::POLYMER) return 3;
                return 4; // unknown types last
            };
            return getTypePriority(a) < getTypePriority(b); });
        return sortedLines;
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

class KMCYamlBuilder
{
public:
    static KMC fromFile(config::CommandLineConfig config)
    {
        console::error("YAML input not yet supported. Exiting.");
        // return KMC();
    }

private:
};