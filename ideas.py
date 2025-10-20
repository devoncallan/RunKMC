from dataclasses import dataclass


"""
Outputting chain length distributions from simulations?

- For sequences, we are doing the total monomer count, total sequence count, total sequence length2 for each monomer type across all polymers.
--- Also, discretized into buckets along the polymer chain to get a "spatially resolved" sequence length distribution.
--- Could output



At each analysis timepoint, mcPolymer specifices:

numMonomers
MW of monomer 1: X
Total Mn, Mw, and Mw/Mn (Molecular weight distribution)
numChains, numSCB, numLCB
Total Pn, Pw, and Pw/Pn (Chain length distribution)

chainLength, numChains, numSCB, numLCB, total count of monomer 1, ...
-, -, -, -, -, ...
-, -, -, -, -, ...


This is a histogram for each chain length.


We could for each polymer:
Output monomer count A, monomer count B, numSCB, numLCB

Or could store as histogram mode:
- 


Polymer

time, ChainLength, 


"""


@dataclass
class ChainLengthDistribution:
    """
    I suppose would be expanded for branching as well.
    -
    """

    pass


@dataclass
class DiscreteChainLengthDistribution(ChainLengthDistribution):
    """
    Obtained from stochastic simulations.
    - RunKMC
    - Sparks Hybrid
    """

    pass


@dataclass
class ChainLengthDistributionStats:
    """
    Obtained from averaging of discrete distributions
    or from method of moments equations.

    """

    pass


@dataclass
class SequenceLengthDistribution:
    """
    Distribution of monomer sequence lengths in a polymer.
    """

    pass


@dataclass
class DiscreteSequenceLengthDistribution(SequenceLengthDistribution):
    """
    Or discrete "average" SLD?
    For each monomer
    - Total monomer count across all polymers
    - Total sequence count across all polymers
    - Total sequence length2 across all polymers

    """

    pass


@dataclass
class ContinuousSequenceLengthDistribution(SequenceLengthDistribution):
    """"""

    pass
