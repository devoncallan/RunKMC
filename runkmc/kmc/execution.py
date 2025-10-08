import subprocess
from pathlib import Path

from runkmc import PATHS

from .build import ensure_binary_exists


def execute_simulation(
    input_filepath: Path | str,
    output_dir: Path | str,
    report_polymers: bool = False,
    report_sequences: bool = False,
) -> None:

    input_filepath = Path(input_filepath)
    output_dir = Path(output_dir)

    cmd = [
        str(PATHS.EXECUTABLE_PATH.absolute()),
        str(input_filepath.absolute()),
        str(output_dir.absolute()),
    ]

    if report_polymers:
        cmd.append("--report-polymers")
    if report_sequences:
        cmd.append("--report-sequences")

    print(f"Executing command: {' '.join(cmd)}")

    try:
        process = subprocess.Popen(
            cmd,
            cwd=PATHS.PROJECT_ROOT,
            text=True,
        )

        print("Running KMC simulation...")
        print(f"Results: {str(output_dir.absolute())}")

        stdout, stderr = process.communicate()
        if process.returncode != 0:
            error_msg = f"Simulation failed with return code {process.returncode}\n"
            if stderr:
                error_msg += f"Error output:\n{stderr}"
            if stdout:
                error_msg += f"Standard output:\n{stdout}"
            raise RuntimeError(error_msg)

        print(f"RunKMC executed successfully.")

    except KeyboardInterrupt:
        process.terminate()
        try:
            process.wait(timeout=3)
        except subprocess.TimeoutExpired:
            process.kill()
        raise KeyboardInterrupt("Simulation interrupted by user.")


# def compile_run_kmc(force: bool = False) -> None:
#     """Ensure RunKMC binary exists, building if necessary.

#     This function is maintained for backward compatibility.
#     New code should use runkmc.kmc.build.ensure_binary_exists() instead.

#     Args:
#         force: If True, rebuild from source even if binary exists.
#     """
#     from .build import ensure_binary_exists

#     ensure_binary_exists(force_rebuild=force, verbose=True)
