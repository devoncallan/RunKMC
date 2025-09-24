from typing import Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np

from runkmc.results import SimulationResult


def plot_monomer_conversion(result: SimulationResult, ax: Optional[Axes] = None):
    """Plot monomer conversion over time from simulation result."""

    if ax is None:
        fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)


    state = result.results
    monomer_names = result.metadata.get_monomer_names()

    for monomer_name in monomer_names:
        monomer_conv = state.unit_convs[monomer_name]
        ax.plot(state.kmc_time, monomer_conv, label=monomer_name)

    ax.set_xlabel("KMC Time")
    ax.set_ylabel("Monomer Conversion")
    ax.legend(frameon=False)

    plt.tight_layout()

    return ax


def plot_monomer_counts(result: SimulationResult, ax: Optional[Axes] = None):
    """Plot monomer counts over time from simulation result."""

    if ax is None:
        fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)

    state = result.results
    monomer_names = result.metadata.get_monomer_names()

    for monomer_name in monomer_names:
        monomer_count = state.unit_counts[monomer_name]
        ax.plot(state.kmc_time, monomer_count, label=monomer_name)

    ax.set_xlabel("KMC Time")
    ax.set_ylabel("Monomer Count")
    ax.legend(frameon=False)

    plt.tight_layout()
    return ax
