from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("runkmc")
except PackageNotFoundError:
    __version__ = "0.1.0-dev"


class PATHS:
    PACKAGE_ROOT = Path(__file__).parent
    PROJECT_ROOT = PACKAGE_ROOT.parent
    CPP_DIR = PROJECT_ROOT / "cpp"
    BUILD_DIR = CPP_DIR / "build"
    EXECUTABLE_PATH = PACKAGE_ROOT / "build" / "RunKMC"
    if not EXECUTABLE_PATH.exists():
        EXECUTABLE_PATH = BUILD_DIR / "RunKMC"

    TEMPLATE_DIR = PACKAGE_ROOT / "models" / "templates"


from .kmc import RunKMC, SimulationConfig, KMCConfig
from .results import (
    SimulationResult,
    StateData,
    Metadata,
    SequenceData,
    SimulationPaths,
)

__all__ = [
    "__version__",
    "PATHS",
    "RunKMC",
    "SimulationConfig",
    "KMCConfig",
    "SimulationResult",
    "StateData",
    "Metadata",
    "SequenceData",
    "SimulationPaths",
]
