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
};

class Monomer : public Unit
{
public:
	double monomerDensity; // g/cm^3
	double polymerDensity; // g/cm^3

	Monomer(SpeciesID ID, std::string name, std::string_view type, double C0_, double FW_, double monomerDensity_, double polymerDensity_, double efficiency_ = 1.0)
		: Unit(ID, name, type, C0_, FW_, efficiency_), monomerDensity(monomerDensity_), polymerDensity(polymerDensity_) {}
};