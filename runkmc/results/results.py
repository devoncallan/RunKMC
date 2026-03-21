from __future__ import annotations
from pathlib import Path
from typing import Dict, List, Optional, TYPE_CHECKING
from dataclasses import dataclass

from .paths import SimulationPaths
from .state import StateData, PosChainData, ChainRecordData, SegmentHistogramData
from ..core.species import SpeciesRegistry
from .polymers import read_polymer_file, PolymerSequence

if TYPE_CHECKING:
    from .registry import SimulationRecord


@dataclass
class SimulationResult:

    paths: SimulationPaths
    species: SpeciesRegistry
    results: StateData
    chain_data: Optional[Dict[str, ChainRecordData]] = None
    segment_hist: Optional[Dict[str, SegmentHistogramData]] = None
    sequence_data: Optional[Dict[str, PosChainData]] = None
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

        # Load per-container files
        chain_map: Dict[str, ChainRecordData] = {}
        segment_hist_map: Dict[str, SegmentHistogramData] = {}
        sequence_data_map: Dict[str, PosChainData] = {}

        for name in species.get_polymer_names():
            p = paths.chains_filepath(name)
            if p.exists():
                r = ChainRecordData.load(p, species)
                if r is not None:
                    chain_map[name] = r

            p = paths.segment_hist_filepath(name)
            if p.exists():
                r = SegmentHistogramData.load(p, species)
                if r is not None:
                    segment_hist_map[name] = r

            p = paths.pos_chain_filepath(name)
            if p.exists():
                r = PosChainData.load(p, species)
                if r is not None:
                    sequence_data_map[name] = r

        # Load polymer data if it exists
        polymer_data = None
        if paths.polymers_filepath.exists():
            polymer_data = read_polymer_file(paths.polymers_filepath)

        return SimulationResult(
            paths, species, results,
            chain_data=chain_map or None,
            segment_hist=segment_hist_map or None,
            sequence_data=sequence_data_map or None,
            polymer_data=polymer_data,
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
