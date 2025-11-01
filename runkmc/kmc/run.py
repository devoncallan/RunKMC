from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Any
import tempfile

from .execution import execute_simulation, parse_only
from runkmc.results import SimulationResult, SimulationPaths
from runkmc.templates import create_input_file
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

    def run_from_template(
        self,
        template_name: str,
        template_params: Dict[str, Any],
        use_existing: bool = True,
        overwrite: bool = False,
        **kwargs,
    ) -> SimulationResult:

        temp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
        input_path = Path(temp_file.name)
        temp_file.close()
        create_input_file(template_name, template_params, input_path)

        try:
            result = self.run_simulation(
                input_filepath=input_path,
                use_existing=use_existing,
                overwrite=overwrite,
                **kwargs,
            )
            return result
        finally:
            if input_path.exists():
                input_path.unlink()

    def run_simulation(
        self,
        input_filepath: Path | str,
        use_existing: bool = True,
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
        if not input_filepath.exists():
            raise FileNotFoundError(f"Input file not found: {input_filepath}")

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
            input_hash = SimulationRegistry.hash_input(parsed_input)

        record = self.registry.get_latest(input_hash, completed=True)

        if record:
            print(f"✓ Found existing simulation (hash: {record.input_hash[:8]}...)")

            if not overwrite and use_existing:
                print(f"  Loading from: {record.dir}")
                return SimulationResult.load(record.dir)
            elif overwrite:
                print(f"  Overwriting existing simulation.")
            elif not use_existing:
                print(f"  Running new simulation despite existing record.")
                record = None

        if not record:
            record = self.registry.new_record(input_hash)

        output_dir = record.dir
        output_dir.mkdir(parents=True, exist_ok=True)

        try:
            # Run the simulation
            print(f"➤ Running new simulation (hash: {input_hash[:8]}...)")
            execute_simulation(input_filepath, record.dir, **kwargs)
            self.registry.update_completion(record, completed=True)
            print(f"✓ Simulation completed: {output_dir}")

        except Exception as e:
            print(f"✗ Simulation failed: {e}")
            self.registry.update_completion(record, completed=False)
            raise
        
        return SimulationResult.load(record.dir)
