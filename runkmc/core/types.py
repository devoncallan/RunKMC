from dataclasses import dataclass
from typing import List


@dataclass
class CommandLineConfig:
    input_filepath: str
    output_dir: str
    report_polymers: bool = False
    report_sequences: bool = False
    parse_only: bool = False
    debug: bool = False


@dataclass
class SimulationConfig:
    num_particles: int
    termination_time: float
    analysis_time: float


@dataclass
class SpeciesRead:
    name: str
    type: str


@dataclass
class UnitRead(SpeciesRead):
    C0: float = 0.0
    FW: float = 0.0
    efficiency: float = 1.0


@dataclass
class PolymerTypeRead(SpeciesRead):
    end_group_unit_names: List[str]


@dataclass
class PolymerLabelsRead(SpeciesRead):
    polymer_names: List[str]


@dataclass
class SpeciesSetRead:
    units: List[UnitRead]
    polymer_types: List[PolymerTypeRead]
    polymer_labels: List[PolymerLabelsRead]


@dataclass
class RateConstantRead:
    name: str
    k: float


@dataclass
class ReactionRead:
    type: str
    rate_constant_name: str
    reactant_names: List[str]
    product_names: List[str]


@dataclass
class KMCInputRead:
    config: SimulationConfig
    species: SpeciesSetRead
    rate_constants: List[RateConstantRead]
    reactions: List[ReactionRead]
