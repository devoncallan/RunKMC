from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Any
import tempfile

from .execution import execute_simulation, parse_only
from runkmc.results import SimulationResult, SimulationPaths
from runkmc.models import create_input_file
from runkmc.results.registry import SimulationRegistry


@dataclass
class SimulationConfig:

    model_name: str
    kmc_inputs: Dict[str, Any]
    report_polymers: bool = False
    report_sequences: bool = False


@dataclass
class KMCConfig:

    num_units: int
    termination_time: float
    analysis_time: float


class RunKMC:

    def __init__(self, base_dir: Path | str, compile: bool = False):
        from .build import ensure_binary_exists

        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)
        ensure_binary_exists(force_rebuild=compile)

        self.registry = SimulationRegistry(self.base_dir)

    # def run_from_config(
    #     self, config: SimulationConfig, sim_id: Optional[str] = None
    # ) -> SimulationResult:

    #     if sim_id is None:
    #         sim_id = f"sim_{uuid4()}"

    #     output_dir = self.base_dir / sim_id
    #     output_dir.mkdir(parents=True, exist_ok=True)

    #     self.input_filepath = output_dir / "input.txt"
    #     create_input_file(config.model_name, config.kmc_inputs, self.input_filepath)

    #     return self.run_from_file(
    #         self.input_filepath,
    #         config.report_polymers,
    #         config.report_sequences,
    #         sim_id=sim_id,
    #     )

    def run_simulation(
        self,
        input_filepath: Path | str,
        overwrite: bool = False,
        **kwargs,
    ) -> SimulationResult:
        """Run a simulation with automatic caching based on input hash.

        Args:
            input_filepath: Path to the input file
            overwrite: If True, ignore cache and re-run simulation
            **kwargs: Additional arguments passed to execute_simulation
                     (report_polymers, report_sequences, etc.)

        Returns:
            SimulationResult object containing the results
        """
        input_filepath = Path(input_filepath)

        # Pre-parse simulation inputs
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)

            # Generate and load parsed KMC input
            parse_only(input_filepath, temp_dir_path)

            paths = SimulationPaths(temp_dir_path)
            parsed_input = paths.parsed_input_filepath

            if not parsed_input.exists():
                raise FileNotFoundError(
                    f"Parsed input file not found: {parsed_input}. "
                    "The C++ parser may not have generated the expected output."
                )

            # Hash the parsed KMC input file
            input_hash = SimulationRegistry.hash_file(parsed_input)

        # Check registry for existing simulation
        record = self.registry.latest_completed(input_hash)

        if record and not overwrite:
            print(f"✓ Found cached simulation results (hash: {input_hash[:8]}...)")
            print(f"  Loading from: {record.dir}")
            return SimulationResult.from_record(record)
        elif record and overwrite:
            print(f"↻ Overwriting existing simulation (hash: {input_hash[:8]}...)")
        else:  # Simulation does not exist
            print(f"➤ Running new simulation (hash: {input_hash[:8]}...)")
            record = self.registry.insert(input_hash)

        try:
            # Run the simulation
            execute_simulation(input_filepath, record.dir, **kwargs)

            # Mark as completed
            self.registry.update_completion(input_hash, completed=True)

            # Load and return results
            return SimulationResult.load(record.dir)

        except Exception as e:
            # If simulation fails, keep it marked as incomplete
            print(f"✗ Simulation failed: {e}")
            self.registry.update_completion(input_hash, completed=False)
            raise
