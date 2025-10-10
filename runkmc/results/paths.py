from pathlib import Path
from runkmc.core import C


# Should mirror c++ implementation
class SimulationPaths:

    def __init__(self, data_dir: Path | str):
        self.data_dir = Path(data_dir)

        self.data_dir.mkdir(parents=True, exist_ok=True)

    @property
    def input_filepath(self) -> Path:
        """Original input file (copied to output dir)"""
        return self.data_dir / "input.txt"

    @property
    def parsed_input_filepath(self) -> Path:
        """Parsed and normalized input file (YAML format) - used for hashing"""
        # Find any .kmc.yaml file in the directory
        kmc_yaml_files = list(self.data_dir.glob("*.kmc.yaml"))
        if kmc_yaml_files:
            return kmc_yaml_files[0]
        # If not found, construct expected name based on input.txt
        return self.data_dir / "input.kmc.yaml"

    @property
    def metadata_filepath(self) -> Path:
        return self.data_dir / "metadata.yaml"

    @property
    def results_filepath(self) -> Path:
        return self.data_dir / "results.csv"

    @property
    def sequence_filepath(self) -> Path:
        return self.data_dir / "sequences.csv"

    @property
    def polymers_filepath(self) -> Path:
        return self.data_dir / "polymers.dat"
