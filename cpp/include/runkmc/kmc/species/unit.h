#pragma once
#include "common.h"

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

	std::string toString() const
	{
		return name + " (" + std::to_string(ID) + "): " + std::to_string(count) + " / " + std::to_string(getInitialCount());
	}
};

typedef Unit *UnitPtr;

Unit UNIT_UNDEF = Unit(0, "UNDEFINED", SpeciesType::UNDEFINED, 0.0, 0.0, 0.0);
