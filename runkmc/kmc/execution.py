from dataclasses import dataclass
import subprocess
from pathlib import Path
from typing import List

from runkmc import PATHS


@dataclass
class CommandLineConfig:

    input_filepath: Path
    output_dir: Path
    report_polymers: bool = False
    report_sequences: bool = False

    def to_command(self) -> List[str]:
        cmd = [
            str(PATHS.EXECUTABLE_PATH.absolute()),
            str(self.input_filepath.absolute()),
            str(self.output_dir.absolute()),
        ]
        if self.report_polymers:
            cmd.append("--report-polymers")
        if self.report_sequences:
            cmd.append("--report-sequences")
        return cmd


def _execute_simulation(config: CommandLineConfig) -> None:

    if not config.input_filepath.exists():
        raise FileNotFoundError(f"Input file '{config.input_filepath}' does not exist")
    if not config.input_filepath.is_file():
        raise ValueError(f"'{config.input_filepath}' is not a file")

    config.output_dir.mkdir(parents=True, exist_ok=True)

    cmd = config.to_command()

    try:
        process = subprocess.Popen(cmd, cwd=PATHS.PROJECT_ROOT, text=True)
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            error_msg = f"Simulation failed with return code {process.returncode}\n"
            if stderr:
                error_msg += f"Error output:\n{stderr}"
            if stdout:
                error_msg += f"Standard output:\n{stdout}"
            raise RuntimeError(error_msg)

    except KeyboardInterrupt:
        process.terminate()
        try:
            process.wait(timeout=3)
        except subprocess.TimeoutExpired:
            process.kill()
        raise KeyboardInterrupt("Simulation interrupted by user.")

    print(f"RunKMC executed successfully.")


def execute_simulation(
    input_filepath: Path | str,
    output_dir: Path | str,
    report_polymers: bool = False,
    report_sequences: bool = False,
) -> None:

    config = CommandLineConfig(
        input_filepath=Path(input_filepath),
        output_dir=Path(output_dir),
        report_polymers=report_polymers,
        report_sequences=report_sequences,
    )
    _execute_simulation(config)
