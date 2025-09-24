import os
import shutil
import subprocess
from pathlib import Path

from runkmc import PATHS, __version__


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


def compile_run_kmc(force: bool = False) -> None:

    # Check for precompiled binary
    precompiled_path = (
        PATHS.PACKAGE_ROOT / "bin" / ("RunKMC.exe" if os.name == "nt" else "RunKMC")
    )
    if not force and precompiled_path.exists():

        PATHS.BUILD_DIR.mkdir(exist_ok=True)
        shutil.copy2(precompiled_path, PATHS.EXECUTABLE_PATH)
        return

    # No precompiled binary, compile from source
    _compile_from_source(force)


def _compile_from_source(force: bool = False) -> None:

    if not force and PATHS.EXECUTABLE_PATH.exists():
        return

    # Ensure build directory exists
    PATHS.BUILD_DIR.mkdir(parents=True, exist_ok=True)

    # Use CMake to build
    configure_cmd = [
        "cmake",
        "-B",
        str(PATHS.BUILD_DIR),
        "-S",
        str(PATHS.CPP_DIR),
        f"-DRUNKMC_VERSION={__version__}",
    ]
    build_cmd = ["cmake", "--build", str(PATHS.BUILD_DIR)]

    try:
        subprocess.run(configure_cmd, check=True, capture_output=True, text=True)
        subprocess.run(build_cmd, check=True, capture_output=True, text=True)
        print("RunKMC compiled successfully.")

    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Failed to compile RunKMC executable: {e.stderr}")
