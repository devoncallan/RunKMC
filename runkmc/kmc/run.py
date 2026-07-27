from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Any, List, Optional
import shutil
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
        extra_files: Optional[Dict[str, Path]] = None,
        **kwargs,
    ) -> SimulationResult:

        input_path = create_input_file(template_name, template_params)
        side_files: List[Path] = []

        try:
            # Write any extra files (e.g. diffusion data) alongside the input file
            if extra_files:
                for filename, src_path in extra_files.items():
                    dest = input_path.parent / filename
                    shutil.copy2(src_path, dest)
                    side_files.append(dest)

            # Include the full rendered input (which contains plugin config) in the
            # hash, since parsed_input.yaml strips the plugins section.
            # Also include any extra files (e.g. diffusion data) whose content
            # affects simulation behaviour but is not captured in parsed_input.yaml.
            extra_content = input_path.read_text()
            for f in side_files:
                if f.exists():
                    extra_content += f.read_text()

            result = self.run_simulation(
                input_filepath=input_path,
                use_existing=use_existing,
                overwrite=overwrite,
                input_hash_extra=extra_content,
                **kwargs,
            )
            return result
        finally:
            if input_path.exists():
                input_path.unlink()
            for f in side_files:
                if f.exists():
                    f.unlink()

    def run_ensemble_from_template(
        self,
        n_ensembles: int,
        template_name: str,
        template_params: Dict[str, Any],
        use_existing: bool = True,
        overwrite: bool = False,
        extra_files: Optional[Dict[str, Path]] = None,
        **kwargs,
    ) -> List[SimulationResult]:

        input_path = create_input_file(template_name, template_params)
        side_files: List[Path] = []

        try:
            if extra_files:
                for filename, src_path in extra_files.items():
                    dest = input_path.parent / filename
                    shutil.copy2(src_path, dest)
                    side_files.append(dest)

            extra_content = input_path.read_text()
            for f in side_files:
                if f.exists():
                    extra_content += f.read_text()
            input_hash = self._compute_input_hash(input_path, extra=extra_content)

            # Find existing completed runs
            existing = [r for r in self.registry.find_all(input_hash) if r.completed]
            results = []

            # Load existing if requested
            if use_existing and existing:
                n_load = min(len(existing), n_ensembles)
                print(f"Loading {n_load} existing simulation(s)...")
                for record in existing[:n_load]:
                    results.append(SimulationResult.load(record.dir))

            # Run additional if needed
            n_needed = n_ensembles - len(results)
            if n_needed > 0:
                print(f"Running {n_needed} new simulation(s)...")
                for _ in range(n_needed):
                    result = self.run_simulation(
                        input_filepath=input_path,
                        use_existing=False,
                        overwrite=overwrite,
                        input_hash=input_hash,
                        **kwargs,
                    )
                    results.append(result)

            return results
        finally:
            if input_path.exists():
                input_path.unlink()
            for f in side_files:
                if f.exists():
                    f.unlink()

    def _compute_input_hash(self, input_filepath: Path | str, extra: str = "") -> str:

        input_filepath = Path(input_filepath)
        if not input_filepath.exists():
            raise FileNotFoundError(f"Input file not found: {input_filepath}")

        # Pre-parse simulation inputs
        try:
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

                # Hash the parsed KMC input file plus any extra content (e.g. plugin
                # config) that the C++ parser strips from parsed_input.yaml.
                return SimulationRegistry.hash_input(parsed_input, params=extra or None)
        except Exception as e:
            raise RuntimeError(f"Failed to compute input hash: {e}") from e

    def run_simulation(
        self,
        input_filepath: Path | str,
        use_existing: bool = True,
        overwrite: bool = False,
        input_hash: str | None = None,
        input_hash_extra: str = "",
        **kwargs,
    ) -> SimulationResult:
        """Run a simulation with automatic caching based on input hash.

        Args:
            input_filepath: Path to the input file
            overwrite: If True, ignore cache and re-run simulation
            input_hash_extra: Additional content to mix into the hash (e.g. plugin config)
            **kwargs: Additional arguments passed to execute_simulation
                     (report_polymers, report_sequences, etc.)

        Returns:
            SimulationResult object containing the results
        """
        input_filepath = Path(input_filepath)
        if input_hash is None:
            input_hash = self._compute_input_hash(input_filepath, extra=input_hash_extra)
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
