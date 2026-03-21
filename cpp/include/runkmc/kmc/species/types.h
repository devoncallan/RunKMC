#pragma once
#include "common.h"
#include "kmc/analysis/types.h"

/**
 * Minimal unit class for KMC simulation.
 * Can include any molecule/non-distributed species.
 * Example: initiator, monomer, dyads, etc.
 */
class Unit : public Species
{
public:
    double C0;
    double FW;
    double efficiency; // for initiators

    Unit(SpeciesID ID, std::string name, std::string_view type, double C0_, double FW_, double efficiency_ = 1.0)
        : Species(ID, name, type), C0(C0_), FW(FW_), efficiency(efficiency_) {}
};


enum class ChainType
{
    Homopolymer,
    Copolymer,
};

struct CopolymerBuffer
{
    std::vector<SpeciesID> units;
};

struct HomopolymerBuffer
{
    SpeciesID monomer = INVALID_SPECIES_ID;
    uint32_t length = 0;
};

struct CompressedCopolymerBuffer
{
    uint32_t length = 0;
    analysis::PositionalSequenceStats posStats;
    analysis::SegmentHistogram segHist;
};