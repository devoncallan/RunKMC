from pathlib import Path
from typing import List, Optional, TypeAlias

import numpy as np
from numpy.typing import NDArray


SpeciesID: TypeAlias = np.uint8
PolymerSequence: TypeAlias = NDArray[SpeciesID]
PolymerMatrix: TypeAlias = NDArray[SpeciesID]


def read_polymer_file(filepath: Path | str) -> List[PolymerSequence]:
    """
    Reads a polymer .dat file and returns a list of NumPy arrays.

    Args:
        filepath: Path to the polymer chain file

    Returns:
        List of NumPy arrays, each representing one polymer chain
    """

    with open(filepath, "r") as file:
        data = file.read().splitlines()

    polymers = [
        np.array([num for num in line.split()], dtype=SpeciesID)
        for line in data
        if line.strip()
    ]

    return polymers


def create_polymer_matrix(
    polymers: List[PolymerSequence],
    max_length: Optional[int] = None,
    max_polymers: Optional[int] = None,
) -> PolymerMatrix:
    """
    Creates a 2D NumPy array from polymer chains with padding for variable lengths.

    Args:
        polymers: List of polymer arrays
        max_length: Maximum length to pad to (uses longest polymer if None)
        max_polymers: Maximum number of polymers to include (None for all)

    Returns:
        2D array where each row is a polymer chain
    """

    if not polymers:
        return np.array([])

    if max_polymers is not None:
        polymers = polymers[:max_polymers]

    num_polymers = len(polymers)
    polymer_lengths = [len(p) for p in polymers]

    full_length = max(polymer_lengths)

    if max_length is None:
        max_length = int(np.ceil(np.mean(polymer_lengths)))
    elif max_length == -1:
        max_length = full_length
    elif max_length > full_length:
        max_length = full_length
    elif max_length < 1:
        raise ValueError(
            "max_length must be positive or -1 for full length. Default is average polymer length."
        )

    # Create padded array
    result = np.full(shape=(num_polymers, max_length), fill_value=0, dtype=SpeciesID)

    for i, polymer in enumerate(polymers):
        if len(polymer) > max_length:
            result[i, :] = polymer[:max_length]
        else:
            result[i, : len(polymer)] = polymer[:max_length]

    return result
