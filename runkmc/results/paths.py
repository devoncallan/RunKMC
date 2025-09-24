from pathlib import Path


# Should mirror c++ implementation
class SimulationPaths:

    def __init__(self, data_dir: Path | str):
        self.data_dir = Path(data_dir)

        self.data_dir.mkdir(parents=True, exist_ok=True)

    @property
    def input_filepath(self) -> Path:
        return self.data_dir / "input.txt"

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
