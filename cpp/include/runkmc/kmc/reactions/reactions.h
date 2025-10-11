#pragma once
#include "common.h"
#include "kmc/species/polymer_type.h"
#include "kmc/reactions/utils.h"

struct RateConstant
{
    std::string name;
    double value;
    RateConstant(const std::string name_ = "undefined", const double value_ = 0.0)
        : name(name_), value(value_) {};
};

// Container for reactant and product species in a reaction
struct ReactionSpecies
{
    std::vector<Species *> reactants;
    std::vector<Species *> products;

    // ----- Helper functions to get species as specific types -----

    // Get a reactant species by index as a PolymerContainer pointer
    template <size_t Idx>
    PolymerContainer *r_poly() const { return static_cast<PolymerContainer *>(reactants[Idx]); }

    // Get a product species by index as a PolymerContainer pointer
    template <size_t Idx>
    PolymerContainer *p_poly() const { return static_cast<PolymerContainer *>(products[Idx]); }

    // Get a reactant species by index as a Unit pointer
    template <size_t Idx>
    Unit *r_unit() const { return static_cast<Unit *>(reactants[Idx]); }

    // Get a product species by index as a Unit pointer
    template <size_t Idx>
    Unit *p_unit() const { return static_cast<Unit *>(products[Idx]); }
};

struct ReactionSchema
{
    std::string_view type;
    std::vector<std::string_view> reactantTypes;
    std::vector<std::string_view> productTypes;

    std::string toString() const
    {
        return std::string(type) + ": " + str::join(reactantTypes, " + ") + " --> " + str::join(productTypes, " + ");
    }

    constexpr void validateSchema() const
    {
        for (const auto &type : reactantTypes)
            SpeciesType::checkValid(type);
        for (const auto &type : productTypes)
            SpeciesType::checkValid(type);
    }

    void validate(const ReactionSpecies &species) const
    {
        // Elementary reactions have no type / size constraints
        if (type == ReactionType::ELEMENTARY)
            return;

        if (species.reactants.size() != reactantTypes.size())
            console::input_error("Reaction " + std::string(type) + " expects " + std::to_string(reactantTypes.size()) + " reactants, got " + std::to_string(species.reactants.size()) + ".");
        if (species.products.size() != productTypes.size())
            console::input_error("Reaction " + std::string(type) + " expects " + std::to_string(productTypes.size()) + " products, got " + std::to_string(species.products.size()) + ".");

        for (size_t i = 0; i < species.reactants.size(); ++i)
        {
            if (!species.reactants[i])
                console::input_error("Reaction " + std::string(type) + " has null reactant at position " + std::to_string(i) + ". Expected type: " + std::string(reactantTypes[i]) + ".");
            if (!typesMatch(reactantTypes[i], species.reactants[i]->type))
                console::input_error("Reaction " + std::string(type) + " has reactant type mismatch at position " + std::to_string(i) + ". Expected type: " + std::string(reactantTypes[i]) + ", got type: " + std::string(species.reactants[i]->type) + ".");
        }

        for (size_t i = 0; i < species.products.size(); ++i)
        {
            if (!species.products[i])
                console::input_error("Reaction " + std::string(type) + " has null product at position " + std::to_string(i) + ". Expected type: " + std::string(productTypes[i]) + ".");
            if (!typesMatch(productTypes[i], species.products[i]->type))
                console::input_error("Reaction " + std::string(type) + " has product type mismatch at position " + std::to_string(i) + ". Expected type: " + std::string(productTypes[i]) + ", got type: " + std::string(species.products[i]->type) + ".");
        }
    }

private:
    static bool typesMatch(std::string_view expected, std::string_view actual)
    {
        if (expected == actual)
            return true;
        if (expected == SpeciesType::UNIT && SpeciesType::isUnitType(actual))
            return true;
        if (expected == SpeciesType::POLYMER && SpeciesType::isPolymerType(actual))
            return true;
        return false;
    }
};

class Reaction
{
public:
    Reaction(RateConstant rc, const ReactionSchema &schema_, const ReactionSpecies &species_)
        : rateConstant(rc), schema(schema_), species(species_)
    {
        schema.validateSchema();
        schema.validate(species);
    }

    virtual void react() = 0;
    virtual double calculateRate(double NAV) const = 0;
    virtual std::string_view getType() const { return schema.type; }

    std::string toStringWithCounts() const { return ""; }
    // { return rxn_print::reactionToString()

protected:
    RateConstant rateConstant;
    ReactionSchema schema;
    ReactionSpecies species;
};

/**
 * @brief Elementary reaction (e.g., A + B ––> C)
 * Arbitrary number of unit reactants forming arbitrary number of unit products.
 */
class Elementary : public Reaction
{
public:
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::ELEMENTARY,
        {},
        {}};

    Elementary(RateConstant rateConstant, const ReactionSpecies &species)
        : Reaction(rateConstant, SCHEMA, species) {}

    void react() override
    {
        for (size_t i = 0; i < species.reactants.size(); ++i)
            --species.reactants[i]->count;
        for (size_t i = 0; i < species.products.size(); ++i)
            ++species.products[i]->count;
    }

    double calculateRate(double NAV) const override
    {
        double rate = rateConstant.value;
        for (size_t i = 0; i < species.reactants.size(); ++i)
            rate *= species.reactants[i]->count;
        return rate;
    }
};

/**
 * @brief Initiator decomposition reaction (e.g., AIBN ––> I + I)
 * Decomposition of initiator molecule to form two active primary radicals.
 */
class InitiatorDecomposition : public Reaction
{
public:
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::INITIATOR_DECOMPOSITION,
        {SpeciesType::INITIATOR},
        {SpeciesType::UNIT, SpeciesType::UNIT}};

    InitiatorDecomposition(RateConstant rateConstant, const ReactionSpecies &species)
        : Reaction(rateConstant, SCHEMA, species) {}

    void react() override
    {
        auto init = species.r_unit<0>();
        --init->count;

        if (rng::rand() <= init->efficiency)
            ++species.p_unit<0>()->count;
        if (rng::rand() <= init->efficiency)
            ++species.p_unit<1>()->count;
    }

    double calculateRate(double NAV) const override
    {
        return rateConstant.value * species.r_unit<0>()->count;
    }
};

/**
 * @brief Initiator decomposition reaction to polymer (e.g., AIBN ––> R + R)
 * Decomoposition of initiator molecule to form two active primary radicals (polymer)
 */
class InitiatorDecompositionPolymer : public Reaction
{
public:
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::INIT_DECOMP_POLY,
        {SpeciesType::INITIATOR},
        {SpeciesType::POLYMER, SpeciesType::POLYMER}};

    InitiatorDecompositionPolymer(RateConstant rateConstant, ReactionSpecies species)
        : Reaction(rateConstant, SCHEMA, species) {}

    void react() override
    {
        auto init = species.r_unit<0>();
        --init->count;

        if (rng::rand() <= init->efficiency)
        {
            Polymer *polymer = new Polymer();
            polymer->initiate(init->ID);
            species.p_poly<0>()->insertPolymer(polymer);
        }
        if (rng::rand() <= init->efficiency)
        {
            Polymer *polymer = new Polymer();
            polymer->initiate(init->ID);
            species.p_poly<1>()->insertPolymer(polymer);
        }
    }

    double calculateRate(double NAV) const override
    {
        return rateConstant.value * species.r_unit<0>()->count;
    }
};

/**
 * @brief Initiation reaction (e.g., I + A ––> IA)
 * Reaction between a radical molecule and monomer to create a polymer.
 */
class Initiation : public Reaction
{
public:
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::INITIATION,
        {SpeciesType::UNIT, SpeciesType::UNIT},
        {SpeciesType::POLYMER}};

    Initiation(RateConstant rateConstant, const ReactionSpecies &species)
        : Reaction(rateConstant, SCHEMA, species) {}

    void react() override
    {
        auto init = species.r_unit<0>();
        auto mon = species.r_unit<1>();
        auto prod = species.p_poly<0>();

        --init->count;
        --mon->count;

        Polymer *polymer = new Polymer();
        polymer->initiate(init->ID);
        polymer->addUnitToEnd(mon->ID);
        prod->insertPolymer(polymer);
    }

    double calculateRate(double NAV) const override
    {
        return rateConstant.value * species.r_unit<0>()->count * species.r_unit<1>()->count / NAV;
    }
};

class Propagation : public Reaction
{
public:
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::PROPAGATION,
        {SpeciesType::POLYMER, SpeciesType::UNIT},
        {SpeciesType::POLYMER}};

    Propagation(RateConstant rateConstant, const ReactionSpecies &species)
        : Reaction(rateConstant, SCHEMA, species) {}

    void react() override
    {
        auto poly = species.r_poly<0>();
        auto mon = species.r_unit<1>();
        auto prod = species.p_poly<0>();

        --mon->count;
        Polymer *polymer = poly->removeRandomPolymer();
        polymer->addUnitToEnd(mon->ID);
        prod->insertPolymer(polymer);
    }

    double calculateRate(double NAV) const override
    {
        return rateConstant.value * species.r_poly<0>()->count * species.r_unit<1>()->count / NAV;
    }
};

/**
 * @brief Depropagation reaction (e.g., P[A,A] ––> P[?,A] + A).
 * Removes the terminal chain end unit.
 */
class Depropagation : public Reaction
{
public:
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::DEPROPAGATION,
        {SpeciesType::POLYMER},
        {SpeciesType::POLYMER, SpeciesType::UNIT}};

    Depropagation(RateConstant rateConstant, const ReactionSpecies &species)
        : Reaction(rateConstant, SCHEMA, species) {}

    void react() override
    {

        auto r_poly = species.r_poly<0>();
        auto p_poly = species.p_poly<0>();
        auto p_unit = species.p_unit<0>();

        Polymer *polymer = r_poly->removeRandomPolymer();
        polymer->removeUnitFromEnd();
        p_poly->insertPolymer(polymer);
        ++p_unit->count;
    }

    double calculateRate(double NAV) const override
    {
        return rateConstant.value * species.r_poly<0>()->count;
    }
};

/**
 * @brief Termination by disproportionation (e.g., P[A,A] + P[B,A] ––> D + D).
 *
 */
class TerminationDisproportionation : public Reaction
{
public:
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::TERMINATION_D,
        {SpeciesType::POLYMER, SpeciesType::POLYMER},
        {SpeciesType::POLYMER, SpeciesType::POLYMER}};

    TerminationDisproportionation(RateConstant rateConstant, const ReactionSpecies &species)
        : Reaction(rateConstant, SCHEMA, species)
    {
        sameReactant = (species.r_poly<0>()->name == species.r_poly<1>()->name) ? 1 : 0;
    }

    void react() override
    {
        Polymer *poly1 = species.r_poly<0>()->removeRandomPolymer();
        Polymer *poly2 = species.r_poly<1>()->removeRandomPolymer();
        poly1->terminateByDisproportionation();
        poly2->terminateByDisproportionation();
        species.p_poly<0>()->insertPolymer(poly1);
        species.p_poly<1>()->insertPolymer(poly2);
    }

    double calculateRate(double NAV) const override
    {
        return rateConstant.value * species.r_poly<0>()->count * (species.r_poly<1>()->count - sameReactant) / NAV;
    }

private:
    uint8_t sameReactant; // True = 1, False = 0
};

/**
 * @brief Termination by combination (e.g., P[A,A] + P[B,A] ––> D).
 *
 */
class TerminationCombination : public Reaction
{
public:
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::TERMINATION_C,
        {SpeciesType::POLYMER, SpeciesType::POLYMER},
        {SpeciesType::POLYMER}};

    TerminationCombination(RateConstant rateConstant, const ReactionSpecies &species)
        : Reaction(rateConstant, SCHEMA, species)
    {
        sameReactant = (species.r_poly<0>()->name == species.r_poly<1>()->name) ? 1 : 0;
    }

    void react() override
    {
        Polymer *poly1 = species.r_poly<0>()->removeRandomPolymer();
        Polymer *poly2 = species.r_poly<1>()->removeRandomPolymer();
        poly1->terminateByCombination(poly2);
        species.p_poly<0>()->insertPolymer(poly1);
    }

    double calculateRate(double NAV) const override
    {
        return rateConstant.value * species.r_poly<0>()->count * (species.r_poly<1>()->count - sameReactant) / NAV;
    }

private:
    uint8_t sameReactant; // True = 1, False = 0
};

class ChainTransferToMonomer : public Reaction
{
public:
    // P + M --> D + R
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::CHAINTRANSFER_M,
        {SpeciesType::POLYMER, SpeciesType::UNIT},
        {SpeciesType::POLYMER, SpeciesType::POLYMER}};

    ChainTransferToMonomer(RateConstant rateConstant, const ReactionSpecies &species)
        : Reaction(rateConstant, SCHEMA, species) {}

    void react() override
    {
        auto mon = species.r_unit<1>();

        // Terminate a polymer
        Polymer *poly = species.r_poly<0>()->removeRandomPolymer();
        poly->terminateByChainTransfer();
        species.p_poly<0>()->insertPolymer(poly);
        --mon->count;

        // Create a new monomer radical
        Polymer *newRadical = new Polymer();
        newRadical->initiate(mon->ID);
        newRadical->addUnitToEnd(mon->ID);
        species.p_poly<1>()->insertPolymer(newRadical);
    }

    double calculateRate(double NAV) const override
    {
        return rateConstant.value * species.r_poly<0>()->count * species.r_unit<1>()->count / NAV;
    }
};

class ThermalInitiationMonomer : public Reaction
{
public:
    // M + M + M --> P + P
    static inline const ReactionSchema &SCHEMA = {
        ReactionType::THERM_INIT_M,
        {SpeciesType::UNIT, SpeciesType::UNIT, SpeciesType::UNIT},
        {SpeciesType::POLYMER, SpeciesType::POLYMER}};
    ThermalInitiationMonomer(RateConstant rateConstant, const ReactionSpecies &species)
        : Reaction(rateConstant, SCHEMA, species) {}

    void react() override
    {
        --species.r_unit<0>()->count;
        --species.r_unit<1>()->count;
        --species.r_unit<2>()->count;

        Polymer *polymer1 = new Polymer();
        polymer1->addUnitToEnd((species.r_unit<0>())->ID);
        species.p_poly<0>()->insertPolymer(polymer1);

        Polymer *polymer2 = new Polymer();
        polymer2->addUnitToEnd((species.r_unit<1>())->ID);
        species.p_poly<1>()->insertPolymer(polymer2);
    }

    double calculateRate(double NAV) const override
    {
        return rateConstant.value * pow(species.r_unit<0>()->count, 3) / NAV;
    }
};
