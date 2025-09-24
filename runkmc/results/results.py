from __future__ import annotations
from pathlib import Path
from typing import List, Optional
from dataclasses import dataclass

from .paths import SimulationPaths
from .state import StateData, Metadata, SequenceData
from .polymers import read_polymer_file, create_polymer_matrix, PolymerSequence


@dataclass
class SimulationResult:

    paths: SimulationPaths
    metadata: Metadata
    results: StateData
    sequence_data: Optional[SequenceData]
    polymer_data: Optional[List[PolymerSequence]] = None

    @staticmethod
    def load(output_dir: Path | str) -> SimulationResult:
        """Load a simulation result from the specified output directory."""

        output_dir = Path(output_dir)
        paths = SimulationPaths(output_dir)

        # Load metadata
        if not paths.metadata_filepath.exists():
            raise FileNotFoundError(
                f"Metadata file {paths.metadata_filepath} not found."
            )
        metadata = Metadata.load(paths.metadata_filepath)

        # Load results
        if not paths.results_filepath.exists():
            raise FileNotFoundError(f"Results file {paths.results_filepath} not found.")

        results = StateData.from_csv(paths.results_filepath, metadata)

        # Load sequence data if it exists
        sequence_data = None
        if paths.sequence_filepath.exists():
            sequence_data = SequenceData.from_csv(paths.sequence_filepath, metadata)

        # Load polymer data if it exists
        polymer_data = None
        if paths.polymers_filepath.exists():
            polymer_data = read_polymer_file(paths.polymers_filepath)

        return SimulationResult(paths, metadata, results, sequence_data, polymer_data)
