from typing import Dict, List, Optional, Union

import numpy as np
from numpy.typing import NDArray

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import ListedColormap, BoundaryNorm

from runkmc.results.polymers import SpeciesID, PolymerMatrix


DEFAULT_COLORS = {
    1: "tab:blue",
    2: "tab:orange",
    3: "tab:green",
    4: "tab:red",
    5: "tab:purple",
}


def plot_polymers_as_image(
    polymers: PolymerMatrix,
    colors: Optional[Dict[SpeciesID, str]] = None,
    ax: Optional[Axes] = None,
) -> Axes:

    if ax is None:
        fig, ax = plt.subplots(figsize=(3.5, 3.5), dpi=300)

    colors = DEFAULT_COLORS if colors is None else colors

    if colors.get(0) is None:
        colors[0] = "#FFFFFF"

    ids = sorted(colors.keys())

    cmap = ListedColormap([colors[i] for i in ids])
    norm = BoundaryNorm(
        boundaries=[i - 0.5 for i in ids] + [ids[-1] + 0.5], ncolors=len(cmap.colors)
    )

    ax.imshow(polymers, cmap=cmap, norm=norm)
    ax.set_ylabel("Polymer Chains")
    ax.set_xlabel("Chain Position")

    plt.tight_layout()

    return ax
