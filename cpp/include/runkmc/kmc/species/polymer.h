#pragma once

// #include <type_traits>
#include <variant>

#include "common.h"
#include "kmc/species/types.h"
#include "kmc/analysis/analysis.h"

class Polymer
{
private:
    using Chain = std::variant<HomopolymerBuffer, CopolymerBuffer, CompressedCopolymerBuffer>;

    PolymerState state;
    ChainType chainType;
    Chain chain;
    SpeciesID initiator;

    template <typename FnHomo, typename FnCopoly, typename FnCompressed>
    decltype(auto) visitChain(FnHomo &&homoFn, FnCopoly &&copolyFn, FnCompressed &&compressedFn)
    {
        return std::visit([&](auto &c) -> decltype(auto)
                          {
            using T = std::decay_t<decltype(c)>;
            if constexpr (std::is_same_v<T, HomopolymerBuffer>)
                return std::forward<FnHomo>(homoFn)(c);
            else if constexpr (std::is_same_v<T, CopolymerBuffer>)
                return std::forward<FnCopoly>(copolyFn)(c);
            else
                return std::forward<FnCompressed>(compressedFn)(c); }, chain);
    }

    template <typename FnHomo, typename FnCopoly, typename FnCompressed>
    decltype(auto) visitChain(FnHomo &&homoFn, FnCopoly &&copolyFn, FnCompressed &&compressedFn) const
    {
        return std::visit([&](const auto &c) -> decltype(auto)
                          {
            using T = std::decay_t<decltype(c)>;
            if constexpr (std::is_same_v<T, HomopolymerBuffer>)
                return std::forward<FnHomo>(homoFn)(c);
            else if constexpr (std::is_same_v<T, CopolymerBuffer>)
                return std::forward<FnCopoly>(copolyFn)(c);
            else
                return std::forward<FnCompressed>(compressedFn)(c); }, chain);
    }

public:
    static ChainType defaultChainType()
    {
        return registry::getNumMonomers() <= 1 ? ChainType::Homopolymer : ChainType::Copolymer;
    }

    Polymer(ChainType type = ChainType::Copolymer, uint32_t maxDOP = 1000)
        : state(ALIVE), chainType(type), chain(), initiator(INVALID_SPECIES_ID)
    {
        if (chainType == ChainType::Homopolymer)
        {
            chain.emplace<HomopolymerBuffer>();
        }
        else if (chainType == ChainType::Copolymer)
        {
            auto &buf = chain.emplace<CopolymerBuffer>();
            buf.units.reserve(maxDOP);
        }
        else
            console::error("Invalid polymer chain type.");
    };

    ~Polymer() = default;

    /***************** Modify functions ****************/

    void updateState(PolymerState ps) { state = ps; }

    ChainType kind() const { return chainType; }

    void addUnitToEnd(const SpeciesID unit)
    {
        visitChain(
            [&](HomopolymerBuffer &_chain)
            {
                if (_chain.monomer == INVALID_SPECIES_ID)
                    _chain.monomer = unit;
                else if (_chain.monomer != unit)
                    console::error("Homopolymer buffer received mismatched monomer ID; mixed compositions are unsupported.");
                ++_chain.length;
            },
            [&](CopolymerBuffer &_chain)
            {
                _chain.units.push_back(unit);
            },
            [](CompressedCopolymerBuffer &)
            {
                console::error("Cannot add unit to compressed polymer.");
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
            [](HomopolymerBuffer &_chain)
            {
                if (_chain.length == 0)
                    console::error("Trying to remove unit from empty polymer.");
                --_chain.length;
                if (_chain.length == 0)
                    _chain.monomer = INVALID_SPECIES_ID;
            },
            [](CopolymerBuffer &_chain)
            {
                if (_chain.units.empty())
                    console::error("Trying to remove unit from empty polymer.");
                if (_chain.units.size() <= 1)
                    console::error("Trying to remove last unit from polymer.");
                _chain.units.pop_back();
            },
            [](CompressedCopolymerBuffer &)
            { console::error("Cannot remove unit from compressed polymer."); });
    }

    void clearSequence()
    {
        visitChain(
            [](HomopolymerBuffer &_chain)
            {
                _chain.monomer = INVALID_SPECIES_ID;
                _chain.length = 0;
            },
            [](CopolymerBuffer &_chain)
            {
                _chain.units.clear();
            },
            [](CompressedCopolymerBuffer &) {});
    }

    /***************** State functions *****************/
    bool isUninitiated() const { return state == PolymerState::UNINITIATED; }

    bool isAlive() const { return state == PolymerState::ALIVE; }

    size_t getDegreeOfPolymerization() const
    {
        return visitChain(
            [](const HomopolymerBuffer &_chain)
            { return static_cast<size_t>(_chain.length); },
            [](const CopolymerBuffer &_chain)
            { return _chain.units.size(); },
            [](const CompressedCopolymerBuffer &_chain)
            { return static_cast<size_t>(_chain.length); });
    }

    bool endGroupIs(const std::vector<SpeciesID> &endGroup) const
    {
        if (!isAlive())
            return false;

        return visitChain(
            [](const HomopolymerBuffer &)
            {
                return true;
            },
            [&](const CopolymerBuffer &_chain)
            {
                if (endGroup.size() > _chain.units.size() + 1)
                    return false;
                return equal(_chain.units.end() - endGroup.size(), _chain.units.end(), endGroup.begin());
            },
            [](const CompressedCopolymerBuffer &)
            {
                return false;
            });
    }

    std::string getSequenceString() const
    {
        return visitChain(
            [](const HomopolymerBuffer &_chain)
            {
                if (_chain.length == 0 || _chain.monomer == INVALID_SPECIES_ID)
                    return std::string();
                std::string out;
                out.reserve(static_cast<size_t>(_chain.length) * 2);
                for (uint32_t i = 0; i < _chain.length; ++i)
                {
                    out += std::to_string(_chain.monomer);
                    out.push_back(' ');
                }
                return out;
            },
            [](const CopolymerBuffer &_chain)
            {
                if (_chain.units.empty())
                    return std::string();
                std::string out;
                out.reserve(_chain.units.size() * 2);
                for (const auto id : _chain.units)
                {
                    out += std::to_string(id);
                    out.push_back(' ');
                }
                return out;
            },
            [](const CompressedCopolymerBuffer &)
            { return std::string(); });
    }

    PolymerState getState() const { return state; }

    const std::vector<SpeciesID> &getSequence() const
    {
        static const std::vector<SpeciesID> empty;
        return visitChain(
            [&](const HomopolymerBuffer &) -> const std::vector<SpeciesID> &
            { return empty; },
            [](const CopolymerBuffer &_chain) -> const std::vector<SpeciesID> &
            { return _chain.units; },
            [&](const CompressedCopolymerBuffer &) -> const std::vector<SpeciesID> &
            { return empty; });
    }

    const analysis::PositionalSequenceStats &getPositionalStats() const
    {
        static const analysis::PositionalSequenceStats empty;
        return visitChain(
            [&](const HomopolymerBuffer &) -> const analysis::PositionalSequenceStats &
            { return empty; },
            [&](const CopolymerBuffer &) -> const analysis::PositionalSequenceStats &
            { return empty; },
            [](const CompressedCopolymerBuffer &_chain) -> const analysis::PositionalSequenceStats &
            { return _chain.posStats; });
    }

    analysis::SequenceStats getChainRecord() const
    {
        return visitChain(
            [](const HomopolymerBuffer &) { return analysis::SequenceStats{}; },
            [](const CopolymerBuffer &) { return analysis::SequenceStats{}; },
            [](const CompressedCopolymerBuffer &_chain) {
                return _chain.posStats.collapse();
            });
    }

    const analysis::SegmentHistogram &getSegmentHistogram() const
    {
        static const analysis::SegmentHistogram empty;
        return visitChain(
            [&](const HomopolymerBuffer &) -> const analysis::SegmentHistogram &
            { return empty; },
            [&](const CopolymerBuffer &) -> const analysis::SegmentHistogram &
            { return empty; },
            [](const CompressedCopolymerBuffer &_chain) -> const analysis::SegmentHistogram &
            { return _chain.segHist; });
    }

    analysis::ChainStats toChainStats(const std::vector<double> &FWs) const
    {
        return visitChain(
            [&](const HomopolymerBuffer &_chain)
            { return analysis::toChainStats(_chain, FWs); },
            [](const CopolymerBuffer &)
            { return analysis::ChainStats{}; },
            [&](const CompressedCopolymerBuffer &_chain)
            { return analysis::toChainStats(_chain, FWs); });
    }

    /***************** Reaction functions *****************/

    void compress()
    {
        visitChain(
            [](HomopolymerBuffer &) {},
            [&](CopolymerBuffer &_chain)
            {
                chain.emplace<CompressedCopolymerBuffer>(analysis::compressCopolymer(_chain));
            },
            [](CompressedCopolymerBuffer &) {});
    }

    void terminateByChainTransfer()
    {
        state = PolymerState::TERMINATED_CT;
    }

    void terminateByDisproportionation()
    {
        state = PolymerState::TERMINATED_D;
    }

    void terminateByCombination(Polymer *&polymer)
    {
        if (chainType == ChainType::Copolymer)
        {
            if (polymer->chainType != ChainType::Copolymer)
                console::error("Combination between copolymer and homopolymer chains is not supported.");

            auto &selfSeq = std::get<CopolymerBuffer>(chain);
            const auto &otherSeq = std::get<CopolymerBuffer>(polymer->chain);
            selfSeq.units.insert(selfSeq.units.end(), otherSeq.units.rbegin(), otherSeq.units.rend());
        }
        else
        {
            if (polymer->chainType != ChainType::Homopolymer)
                console::error("Combination between homopolymer and copolymer chains is not supported.");

            auto &selfHomo = std::get<HomopolymerBuffer>(chain);
            const auto &otherHomo = std::get<HomopolymerBuffer>(polymer->chain);

            if (selfHomo.monomer == INVALID_SPECIES_ID)
                selfHomo.monomer = otherHomo.monomer;

            if (otherHomo.monomer != INVALID_SPECIES_ID && selfHomo.monomer != otherHomo.monomer)
                console::error("Cannot combine homopolymers with different monomer identities.");

            selfHomo.length += otherHomo.length;
        }

        state = PolymerState::TERMINATED_C;
    }
};
