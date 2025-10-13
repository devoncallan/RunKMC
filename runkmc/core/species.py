from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, List

import yaml

from runkmc.core import C


@dataclass
class RegisteredSpecies:

    name: str
    type: str
    ID: int


@dataclass
class SpeciesRegistry:

    species: List[RegisteredSpecies]
    units: Dict[str, RegisteredSpecies]
    monomers: Dict[str, RegisteredSpecies]
    polymers: Dict[str, RegisteredSpecies]

    @staticmethod
    def from_yaml(filepath: Path | str) -> SpeciesRegistry:

        filepath = Path(filepath)
        with open(filepath, "r") as f:
            data = yaml.safe_load(f)

        assert data is not None and isinstance(
            data, dict
        ), f"Species file {filepath} is empty or invalid YAML."

        required_keys = [
            C.io.SPECIES_KEY,
            C.io.UNITS_KEY,
            C.io.MONOMERS_KEY,
            C.io.POLYMERS_KEY,
        ]
        missing_keys = [k for k in required_keys if k not in data]
        if missing_keys:
            raise ValueError(
                f"Species file {filepath} is missing required keys: {', '.join(missing_keys)}"
            )

        species = [
            RegisteredSpecies(
                name=s[C.io.NAME_KEY], type=s[C.io.TYPE_KEY], ID=s[C.io.ID_KEY]
            )
            for s in data[C.io.SPECIES_KEY]
        ]

        units = {s.name: s for s in species if s.name in data[C.io.UNITS_KEY]}
        monomers = {s.name: s for s in species if s.name in data[C.io.MONOMERS_KEY]}
        polymers = {s.name: s for s in species if s.name in data[C.io.POLYMERS_KEY]}

        return SpeciesRegistry(
            species=species, units=units, monomers=monomers, polymers=polymers
        )

    def get_monomer_names(self) -> List[str]:
        return list(self.monomers.keys())

    def get_unit_names(self) -> List[str]:
        return list(self.units.keys())

    def get_polymer_names(self) -> List[str]:
        return list(self.polymers.keys())
