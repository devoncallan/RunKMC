#pragma once

#include <type_traits>
#include <variant>

#include "common.h"
#include "kmc/species/unit.h"
#include "kmc/analysis/analysis.h"

enum class ChainType
{
    Sequence,
    Homopolymer
};

struct SequenceBuffer
{
    std::vector<SpeciesID> units;
    std::vector<analysis::SequenceStats> posStats;

    void reserve(size_t maxDOP)
    {
        units.reserve(maxDOP);
        posStats.reserve(NUM_BUCKETS);
    }

    void clear()
    {
        units.clear();
        std::vector<SpeciesID>().swap(units);
    }
};

struct HomopolymerBuffer
{
    SpeciesID monomer = INVALID_SPECIES_ID;
    uint32_t length = 0;

    void push(SpeciesID unit)
    {
        if (monomer == INVALID_SPECIES_ID)
        {
            monomer = unit;
        }
        else if (monomer != unit)
        {
            console::error("Homopolymer buffer received mismatched monomer ID; mixed compositions are unsupported.");
        }
        ++length;
    }

    void pop()
    {
        if (length == 0)
            console::error("Trying to remove unit from empty polymer.");
        --length;
        if (length == 0)
            monomer = INVALID_SPECIES_ID;
    }

    void reset()
    {
        monomer = INVALID_SPECIES_ID;
        length = 0;
    }
};

class Polymer
{
private:
    using Chain = std::variant<SequenceBuffer, HomopolymerBuffer>;

    PolymerState state;
    ChainType chainType;
    Chain chain;
    SpeciesID initiator;

    template <typename FuncSequence, typename FuncHomo>
    decltype(auto) visitChain(FuncSequence &&seqFn, FuncHomo &&homoFn)
    {
        if (auto *seq = std::get_if<SequenceBuffer>(&chain))
            return std::forward<FuncSequence>(seqFn)(*seq);
        auto *homo = std::get_if<HomopolymerBuffer>(&chain);
        if (!homo)
            console::error("Invalid polymer representation state.");
        return std::forward<FuncHomo>(homoFn)(*homo);
    }

    template <typename FuncSequence, typename FuncHomo>
    decltype(auto) visitChain(FuncSequence &&seqFn, FuncHomo &&homoFn) const
    {
        if (auto *seq = std::get_if<SequenceBuffer>(&chain))
            return std::forward<FuncSequence>(seqFn)(*seq);
        auto *homo = std::get_if<HomopolymerBuffer>(&chain);
        if (!homo)
            console::error("Invalid polymer representation state.");
        return std::forward<FuncHomo>(homoFn)(*homo);
    }

public:
    static ChainType defaultChainType()
    {
        return registry::getNumMonomers() <= 1 ? ChainType::Homopolymer : ChainType::Sequence;
    }

    Polymer(ChainType type = ChainType::Sequence, uint32_t maxDOP = 1000)
        : state(ALIVE), chainType(type), chain(), initiator(INVALID_SPECIES_ID)
    {
        if (chainType == ChainType::Sequence)
        {
            chain.emplace<SequenceBuffer>();
            std::get<SequenceBuffer>(chain).reserve(maxDOP);
        }
        else
        {
            chain.emplace<HomopolymerBuffer>();
        }
    };

    ~Polymer() = default;

    /***************** Modify functions ****************/

    void updateState(PolymerState ps) { state = ps; }

    ChainType kind() const { return chainType; }

    void addUnitToEnd(const SpeciesID unit)
    {
        visitChain(
            [&](SequenceBuffer &seq)
            {
                seq.units.push_back(unit);
            },
            [&](HomopolymerBuffer &homo)
            {
                homo.push(unit);
            });
    }

    void initiate(const SpeciesID unit)
    {
        initiator = unit;
        state = ALIVE;
    }

    void removeUnitFromEnd()
    {
        visitChain(
            [&](SequenceBuffer &seq)
            {
                if (seq.units.empty())
                    console::error("Trying to remove unit from empty polymer.");
                if (seq.units.size() <= 1)
                    console::error("Trying to remove last unit from polymer.");
                seq.units.pop_back();
            },
            [&](HomopolymerBuffer &homo)
            {
                if (homo.length == 0)
                    console::error("Trying to remove unit from empty polymer.");
                if (homo.length <= 1)
                    console::error("Trying to remove last unit from polymer.");
                homo.pop();
            });
    }

    void clearSequence()
    {
        visitChain(
            [&](SequenceBuffer &seq)
            {
                seq.clear();
            },
            [&](HomopolymerBuffer &homo)
            {
                homo.reset();
            });
    }

    /***************** State functions *****************/
    bool isUninitiated() const { return state == PolymerState::UNINITIATED; }

    bool isAlive() const { return state == PolymerState::ALIVE; }

    size_t getDegreeOfPolymerization() const
    {
        return visitChain(
            [](const SequenceBuffer &seq)
            { return seq.units.size(); },
            [](const HomopolymerBuffer &homo)
            { return static_cast<size_t>(homo.length); });
    }

    bool endGroupIs(const std::vector<SpeciesID> &endGroup) const
    {
        if (!isAlive())
            return false;

        return visitChain(
            [&](const SequenceBuffer &seq)
            {
                if (endGroup.size() > seq.units.size() + 1)
                    return false;
                return equal(seq.units.end() - endGroup.size(), seq.units.end(), endGroup.begin());
            },
            [&](const HomopolymerBuffer &)
            {
                return true;
            });
    }

    bool isCompressed() const
    {
        return visitChain(
            [](const SequenceBuffer &seq)
            { return seq.units.empty() && !seq.posStats.empty(); },
            [](const HomopolymerBuffer &)
            { return true; });
    }

    std::string getSequenceString() const
    {
        return visitChain(
            [](const SequenceBuffer &seq)
            {
                if (seq.units.empty())
                    return std::string();
                std::string out;
                out.reserve(seq.units.size() * 2);
                for (const auto id : seq.units)
                {
                    out += std::to_string(id);
                    out.push_back(' ');
                }
                return out;
            },
            [](const HomopolymerBuffer &homo)
            {
                if (homo.length == 0 || homo.monomer == INVALID_SPECIES_ID)
                    return std::string();
                std::string out;
                out.reserve(static_cast<size_t>(homo.length) * 2);
                for (uint32_t i = 0; i < homo.length; ++i)
                {
                    out += std::to_string(homo.monomer);
                    out.push_back(' ');
                }
                return out;
            });
    }

    PolymerState getState() const { return state; }

    const std::vector<SpeciesID> &getSequence() const
    {
        return visitChain(
            [](const SequenceBuffer &seq) -> const std::vector<SpeciesID> &
            { return seq.units; },
            [](const HomopolymerBuffer &) -> const std::vector<SpeciesID> &
            {
                static const std::vector<SpeciesID> empty;
                return empty;
            });
    }

    const std::vector<analysis::SequenceStats> &getPositionalStats() const
    {
        return visitChain(
            [](const SequenceBuffer &seq) -> const std::vector<analysis::SequenceStats> &
            { return seq.posStats; },
            [](const HomopolymerBuffer &) -> const std::vector<analysis::SequenceStats> &
            {
                static const std::vector<analysis::SequenceStats> empty;
                return empty;
            });
    }

    /***************** Reaction functions *****************/

    void terminate()
    {
        visitChain(
            [&](SequenceBuffer &seq)
            {
                seq.posStats = analysis::calculatePositionalSequenceStats(seq.units, NUM_BUCKETS);
                seq.clear();
            },
            [&](HomopolymerBuffer &)
            {
                // Preserve length for homopolymer analysis.
            });
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
        if (chainType == ChainType::Sequence)
        {
            if (polymer->chainType != ChainType::Sequence)
                console::error("Combination between sequence and homopolymer chains is not supported.");

            auto &selfSeq = std::get<SequenceBuffer>(chain);
            const auto &otherSeq = std::get<SequenceBuffer>(polymer->chain);
            selfSeq.units.insert(selfSeq.units.end(), otherSeq.units.rbegin(), otherSeq.units.rend());
        }
        else
        {
            if (polymer->chainType != ChainType::Homopolymer)
                console::error("Combination between homopolymer and sequence chains is not supported.");

            auto &selfHomo = std::get<HomopolymerBuffer>(chain);
            const auto &otherHomo = std::get<HomopolymerBuffer>(polymer->chain);

            if (selfHomo.monomer == INVALID_SPECIES_ID)
                selfHomo.monomer = otherHomo.monomer;

            if (otherHomo.monomer != INVALID_SPECIES_ID && selfHomo.monomer != otherHomo.monomer)
                console::error("Cannot combine homopolymers with different monomer identities.");

            selfHomo.length += otherHomo.length;
        }

        state = PolymerState::TERMINATED_C;
        terminate();
    }
};
