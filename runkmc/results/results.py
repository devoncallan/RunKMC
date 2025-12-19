from __future__ import annotations
from pathlib import Path
from typing import List, Optional, TYPE_CHECKING
from dataclasses import dataclass

from .paths import SimulationPaths
from .state import StateData, SequenceData, ChainHistogramData, SegmentHistogramData
from ..core.species import SpeciesRegistry
from .polymers import read_polymer_file, PolymerSequence

if TYPE_CHECKING:
    from .registry import SimulationRecord


@dataclass
class SimulationResult:

    paths: SimulationPaths
    species: SpeciesRegistry
    results: StateData
    chain_data: Optional[ChainHistogramData] = None
    segment_hist: Optional[SegmentHistogramData] = None
    sequence_data: Optional[SequenceData] = None
    polymer_data: Optional[List[PolymerSequence]] = None

    @staticmethod
    def load(output_dir: Path | str) -> SimulationResult:
        """Load a simulation result from the specified output directory."""

        output_dir = Path(output_dir)
        paths = SimulationPaths(output_dir)

        # Load species registry
        if not paths.species_filepath.exists():
            raise FileNotFoundError(f"Species file {paths.species_filepath} not found.")
        species = SpeciesRegistry.from_yaml(paths.species_filepath)

        # Load results
        if not paths.results_filepath.exists():
            raise FileNotFoundError(f"Results file {paths.results_filepath} not found.")

        results = StateData.from_csv(paths.results_filepath, species)

        # Load sequence data if it exists
        sequence_data = None
        if paths.sequence_filepath.exists():
            sequence_data = SequenceData.from_csv(paths.sequence_filepath, species)
            
        chain_data = None
        if paths.chain_stats_filepath.exists():
            chain_data = ChainHistogramData.load(paths.chain_stats_filepath, species)

        segment_hist = None
        if paths.segment_hist_filepath.exists():
            segment_hist = SegmentHistogramData.load(paths.segment_hist_filepath, species)

        # Load polymer data if it exists
        polymer_data = None
        if paths.polymers_filepath.exists():
            polymer_data = read_polymer_file(paths.polymers_filepath)

        return SimulationResult(
            paths, species, results, chain_data, segment_hist, sequence_data, polymer_data
        )

    @staticmethod
    def from_record(record: SimulationRecord) -> SimulationResult:
        """Load a simulation result from a registry record.

        Args:
            record: SimulationRecord containing the output directory path

        Returns:
            SimulationResult loaded from the record's output directory
        """
        return SimulationResult.load(record.dir)
