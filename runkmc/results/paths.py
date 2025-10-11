from pathlib import Path
from runkmc.core import C


# Should mirror c++ implementation
class SimulationPaths:

    def __init__(self, base_dir: Path | str):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

    @property
    def parsed_input_filepath(self) -> Path:
        return self.base_dir / C.paths.PARSED_INPUT_FILE

    @property
    def species_filepath(self) -> Path:
        return self.base_dir / C.paths.SPECIES_FILE

    @property
    def results_filepath(self) -> Path:
        return self.base_dir / C.paths.RESULTS_FILE

    @property
    def polymers_filepath(self) -> Path:
        return self.base_dir / C.paths.POLYMERS_FILE

    @property
    def sequence_filepath(self) -> Path:
        return self.base_dir / C.paths.SEQUENCES_FILE
