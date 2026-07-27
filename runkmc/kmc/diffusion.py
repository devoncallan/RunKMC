from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path

import yaml


@dataclass
class MonomerDiffusionData:
    name: str
    Tg_m: float        # Glass transition temp of monomer (K)
    Tg_p: float        # Glass transition temp of polymer (K)
    Vf_m: float        # Free volume of pure monomer
    Vf_p: float        # Free volume of pure polymer
    alpha_m: float     # Thermal expansion coefficient monomer (1/K)
    alpha_p: float     # Thermal expansion coefficient polymer (1/K)
    Vf_c: float        # Critical free volume for kp/kfm diffusion control
    B: float           # Translational diffusion decay for kp
    Vf_i: float        # Critical free volume for f diffusion control
    C: float           # Initiation diffusion decay
    delta: float       # Segmental diffusion reaction radius (L/g)
    m: float           # Gel-effect model exponent (Mw)
    n: float           # Post-gel translational diffusion exponent
    A: float           # Gel-effect model parameter
    K3: float          # Gel-point onset criterion
    n_s: float         # Avg monomer units per chain segment
    l0: float          # Monomer unit length (cm)
    rho_m: float       # Density of monomer (g/cm^3)
    rho_p: float       # Density of polymer (g/cm^3)


def _to_plain(obj):
    """Recursively convert numpy scalars/arrays to plain Python types for YAML serialization."""
    if isinstance(obj, dict):
        return {k: _to_plain(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_to_plain(v) for v in obj]
    try:
        return float(obj)
    except (TypeError, ValueError):
        return obj


def write_diffusion_data_file(monomers: list[MonomerDiffusionData], path: Path) -> None:
    data = {"monomers": [_to_plain(asdict(m)) for m in monomers]}
    with open(path, "w") as f:
        yaml.dump(data, f, default_flow_style=False)
