from typing import Optional
from pathlib import Path

from uuid import uuid4

from .config import SimulationConfig
from .execution import execute_simulation, compile_run_kmc
from runkmc.results import SimulationResult
from runkmc.models import create_input_file


class RunKMC:

    def __init__(self, base_dir: Path | str, compile: bool = False):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        compile_run_kmc(compile)

    def run_from_config(
        self, config: SimulationConfig, sim_id: Optional[str] = None
    ) -> SimulationResult:

        if sim_id is None:
            sim_id = f"sim_{uuid4()}"

        output_dir = self.base_dir / sim_id
        output_dir.mkdir(parents=True, exist_ok=True)

        self.input_filepath = output_dir / "input.txt"
        create_input_file(config.model_name, config.kmc_inputs, self.input_filepath)

        return self.run_from_file(
            self.input_filepath,
            config.report_polymers,
            config.report_sequences,
            sim_id=sim_id,
        )

    def run_from_file(
        self,
        input_filepath: Path | str,
        report_polymers: bool = False,
        report_sequences: bool = False,
        sim_id: Optional[str] = None,
    ) -> SimulationResult:

        if sim_id is None:
            sim_id = f"sim_{uuid4()}"
        output_dir = self.base_dir / sim_id

        self.input_filepath = Path(input_filepath)
        execute_simulation(
            self.input_filepath,
            output_dir,
            report_polymers,
            report_sequences,
        )

        results = SimulationResult.load(output_dir)

        return results
