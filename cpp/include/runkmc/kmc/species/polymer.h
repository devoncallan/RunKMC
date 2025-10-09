#pragma once
#include "common.h"
#include "kmc/species/unit.h"
#include "kmc/analysis/utils.h"

class Polymer
{
private:
	PolymerState state;
	std::vector<analysis::SequenceStats> posStats;
	std::vector<SpeciesID> sequence;
	SpeciesID initiator;

public:
	Polymer(uint32_t maxDOP = 1000)
	{
		sequence.reserve(maxDOP);
		state = ALIVE;
		posStats.reserve(NUM_BUCKETS);
	};

	~Polymer() = default;

	/***************** Modify functions ****************/

	void updateState(PolymerState ps) { state = ps; }

	void addUnitToEnd(const SpeciesID unit)
	{
		sequence.push_back(unit);
	}

	void initiate(const SpeciesID unit)
	{
		initiator = unit;
		state = ALIVE;
	}

	void removeUnitFromEnd()
	{
		if (sequence.empty())
			console::error("Trying to remove unit from empty polymer.");
		if (getDegreeOfPolymerization() <= 1)
			console::error("Trying to remove last unit from polymer.");

		sequence.pop_back();
	}

	void clearSequence()
	{
		sequence.clear();
		std::vector<SpeciesID>().swap(sequence);
	}

	/***************** State functions *****************/
	bool isUninitiated() const { return state == PolymerState::UNINITIATED; }

	bool isAlive() const { return state == PolymerState::ALIVE; }

	size_t getDegreeOfPolymerization() const { return sequence.size(); }

	bool endGroupIs(const std::vector<SpeciesID> &endGroup) const
	{
		if (!isAlive() || endGroup.size() > getDegreeOfPolymerization() + 1)
			return false;
		return equal(sequence.end() - endGroup.size(), sequence.end(), endGroup.begin());
	}

	bool isCompressed() const
	{
		if (sequence.empty() && !posStats.empty())
			return true;
		return false;
	}

	std::string getSequenceString() const
	{
		if (isCompressed())
			return "";

		std::string sequenceString;
		for (const auto &id : sequence)
			sequenceString += std::to_string(id) + " ";
		return sequenceString;
	}

	PolymerState getState() const { return state; }

	const std::vector<SpeciesID> &getSequence() const { return sequence; }

	const std::vector<analysis::SequenceStats> &getPositionalStats() const { return posStats; }

	/***************** Reaction functions *****************/

	void terminate()
	{
		posStats = analysis::calculatePositionalSequenceStats(sequence, NUM_BUCKETS);
		clearSequence();
	}

	void terminateByChainTransfer()
	{
		state = PolymerState::TERMINATED_CT;
		terminate();
	}

	void terminateByDisproportionation()
	{
		state = PolymerState::TERMINATED_D;
		terminate();
	}

	void terminateByCombination(Polymer *&polymer)
	{
		sequence.insert(sequence.end(), polymer->sequence.rbegin(), polymer->sequence.rend());
		state = PolymerState::TERMINATED_C;
		terminate();
	}
};
