#pragma once
#include "common.h"

/**
 * Minimal unit class for KMC simulation.
 * Can include any molecule/non-distributed species.
 * Example: initiator, monomer, dyads, etc.
 */
class Unit
{
public:
	std::string type;
	std::string name;
	SpeciesID ID;

	double C0;
	double FW;
	double efficiency; // for initiators

	uint64_t count_0;
	uint64_t count;

	Unit(std::string type_, std::string name_, SpeciesID ID_, double C0_, double FW_, double efficiency_ = 1.0)
		: type(type_), name(name_), ID(ID_), C0(C0_), FW(FW_), efficiency(efficiency_) {}

	void setInitialCount(uint64_t initialCount)
	{
		count_0 = initialCount;
		count = count_0;
	}

	uint64_t getInitialCount() const
	{
		return count_0;
	};

	double calculateConversion() const
	{
		double initialCount = double(getInitialCount());
		if (initialCount == 0)
			return 0;
		return (initialCount - count) / initialCount;
	}

	std::string toString() const
	{
		return name + " (" + std::to_string(ID) + "): " + std::to_string(count) + " / " + std::to_string(getInitialCount());
	}
};

typedef Unit *UnitPtr;

Unit UNIT_UNDEF = Unit(SpeciesType::UNDEFINED, "UNDEFINED", 0, 0.0, 0.0);
