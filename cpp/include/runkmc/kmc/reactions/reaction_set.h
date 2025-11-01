#pragma once
#include "common.h"
#include "kmc/reactions/reactions.h"
#include <algorithm>

namespace reactions
{
    constexpr console::LogContext logger("KMC::Reactions");
}

/**
 * @brief Stores set of all reactions and can calculate cumulative properties such as
 * reaction rates and probabilities.
 *
 */
class ReactionSet
{
public:
    ReactionSet(const std::vector<Reaction *> &reactions_, const std::vector<RateConstant> &rateConstants_) : reactions(reactions_), rateConstants(rateConstants_)
    {
        numReactions = reactions.size();
        reactionRates.resize(numReactions);
        reactionProbabilities.resize(numReactions);
        reactionCumulativeProbabilities.resize(numReactions);
        rateMultipliers.assign(numReactions, 1.0);
    };

    ReactionSet() {};
    ~ReactionSet() {};

    /**
     * @brief Calculate and update reaction probabilities. First updates the reaction rates,
     * then calculates the probability and cumulative probability vectors.
     */
    void updateReactionProbabilities(double NAV_)
    {
        NAV = NAV_;
        updateReactionRates();

        if (totalReactionRate <= 0)
        {
            std::fill(reactionProbabilities.begin(), reactionProbabilities.end(), 0.0);
            std::fill(reactionCumulativeProbabilities.begin(), reactionCumulativeProbabilities.end(), 0.0);
            return;
        }

        reactionProbabilities[0] = reactionRates[0] / totalReactionRate;
        reactionCumulativeProbabilities[0] = reactionProbabilities[0];
        for (size_t i = 1; i < numReactions; ++i)
        {
            reactionProbabilities[i] = reactionRates[i] / totalReactionRate;
            reactionCumulativeProbabilities[i] = reactionProbabilities[i] + reactionCumulativeProbabilities[i - 1];
        }
    }

    size_t chooseRandomReactionIndex() const
    {
        double rn = rng::rand();
        for (size_t reactionIndex = 0; reactionIndex < numReactions; ++reactionIndex)
        {
            if (rn <= reactionCumulativeProbabilities[reactionIndex])
                return reactionIndex;
        }
        console::error("Uh oh! No reaction was chosen - something is wrong with the cumulative probability vector. Exiting.");
        return 0; // Not reached as console::error will exit
    }

    void printSummary() const
    {
        reactions::logger.info("Reaction Set (" + std::to_string(numReactions) + " reactions):");
        for (size_t i = 0; i < numReactions; ++i)
            reactions::logger.info("\t" + std::to_string(i + 1) + ": " + reactions[i]->toString());
    }

    Reaction *getReaction(size_t reactionIndex) const { return reactions[reactionIndex]; }
    size_t getNumReactions() const { return numReactions; }
    const std::vector<RateConstant> &getRateConstants() const { return rateConstants; }
    double getTotalReactionRate() const { return totalReactionRate; }
    bool cantProceed() const { return totalReactionRate == 0; }
    void setNAV(double NAV) { this->NAV = NAV; }
    double getNAV() const { return NAV; }
    const std::vector<double> &getRateMultipliers() const { return rateMultipliers; }

    void setRateMultiplier(size_t reactionIndex, double value)
    {
        if (reactionIndex >= numReactions)
            console::error("Rate multiplier index out of bounds: " + std::to_string(reactionIndex) + ".");
        rateMultipliers[reactionIndex] = value;
    }

    void scaleRateMultiplier(size_t reactionIndex, double factor)
    {
        if (reactionIndex >= numReactions)
            console::error("Rate multiplier index out of bounds: " + std::to_string(reactionIndex) + ".");
        rateMultipliers[reactionIndex] *= factor;
    }

    void resetRateMultipliers()
    {
        std::fill(rateMultipliers.begin(), rateMultipliers.end(), 1.0);
    }

private:
    size_t numReactions = 0;
    std::vector<Reaction *> reactions;
    std::vector<RateConstant> rateConstants;

    double totalReactionRate = 0;
    std::vector<double> reactionRates;
    std::vector<double> reactionProbabilities;
    std::vector<double> reactionCumulativeProbabilities;
    std::vector<double> rateMultipliers;

    double NAV;

    /**
     * @brief Calculate and update reaction rates for all reactions.
     * Also updates total reaction rate.
     */
    void updateReactionRates()
    {
        totalReactionRate = 0;
        for (size_t i = 0; i < numReactions; ++i)
        {
            reactionRates[i] = reactions[i]->calculateRate(NAV) * rateMultipliers[i];
            totalReactionRate += reactionRates[i];
        }
    }
};
