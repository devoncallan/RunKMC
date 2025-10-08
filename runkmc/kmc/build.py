"""Build and compilation logic for the RunKMC binary."""

import glob
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from runkmc import PATHS, __version__


def find_precompiled_binary() -> Optional[Path]:
    """Look for precompiled binary in the package bin directory.

    Returns:
        Path to precompiled binary if found, None otherwise.
    """
    binary_name = "RunKMC.exe" if os.name == "nt" else "RunKMC"
    precompiled_path = PATHS.PACKAGE_ROOT / "bin" / binary_name

    if precompiled_path.exists():
        return precompiled_path

    return None


def build_binary(force: bool = False, verbose: bool = True) -> Path:
    """Build the RunKMC binary from source using CMake.

    This function:
    1. Configures the CMake build
    2. Builds the binary
    3. Finds the binary (handles Windows Debug/ subfolder)
    4. Copies it to runkmc/build/

    Args:
        force: If True, rebuild even if binary already exists.
        verbose: If True, print progress messages.

    Returns:
        Path to the built binary.

    Raises:
        RuntimeError: If CMake configuration or build fails.
        FileNotFoundError: If built binary cannot be found.
    """
    package_build_dir = PATHS.PACKAGE_ROOT / "build"
    binary_name = "RunKMC.exe" if os.name == "nt" else "RunKMC"
    dest_path = package_build_dir / binary_name

    # Check if already built
    if not force and dest_path.exists():
        if verbose:
            print(f"Binary already exists at {dest_path}")
        return dest_path

    # Ensure build directories exist
    PATHS.BUILD_DIR.mkdir(parents=True, exist_ok=True)
    package_build_dir.mkdir(parents=True, exist_ok=True)

    # Configure CMake
    configure_cmd = [
        "cmake",
        "-B",
        str(PATHS.BUILD_DIR),
        "-S",
        str(PATHS.CPP_DIR),
        f"-DRUNKMC_VERSION={__version__}",
        f"-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
    ]

    # Build command
    build_cmd = ["cmake", "--build", str(PATHS.BUILD_DIR)]

    try:
        if verbose:
            print("Configuring CMake...")
        result = subprocess.run(
            configure_cmd, check=True, capture_output=True, text=True
        )
        if verbose and result.stdout:
            print(result.stdout)

        if verbose:
            print("Building RunKMC...")
        result = subprocess.run(build_cmd, check=True, capture_output=True, text=True)
        if verbose and result.stdout:
            print(result.stdout)

        # Find the binary (Windows puts it in Debug/ or Release/ subfolder)
        if os.name == "nt":
            pattern = str(PATHS.BUILD_DIR / "**" / "RunKMC.exe")
        else:
            pattern = str(PATHS.BUILD_DIR / "**" / "RunKMC")

        matches = glob.glob(pattern, recursive=True)
        if not matches:
            raise FileNotFoundError(
                f"Could not find binary matching {pattern} after build"
            )

        binary_path = Path(matches[0])
        if verbose:
            print(f"Found binary: {binary_path}")

        # Copy to package build directory
        shutil.copy2(binary_path, dest_path)
        if verbose:
            print(f"RunKMC compiled successfully and copied to {dest_path}")

        return dest_path

    except subprocess.CalledProcessError as e:
        error_msg = f"Failed to build RunKMC: {e.stderr if e.stderr else str(e)}"
        raise RuntimeError(error_msg)


def ensure_binary_exists(force_rebuild: bool = False, verbose: bool = True) -> Path:
    """Ensure the RunKMC binary exists, building if necessary.

    This function checks for the binary in the following order:
    1. Package build directory (runkmc/build/)
    2. Precompiled binary (runkmc/bin/)
    3. CPP build directory (cpp/build/)
    4. Build from source

    Args:
        force_rebuild: If True, rebuild from source even if binary exists.
        verbose: If True, print progress messages.

    Returns:
        Path to the binary.
    """
    # If forcing rebuild, go straight to building
    if force_rebuild:
        if verbose:
            print("Force rebuild requested, compiling from source...")
        return build_binary(force=True, verbose=verbose)

    # Check if binary already exists in package build dir
    if PATHS.EXECUTABLE_PATH.exists():
        if verbose:
            print(f"Using existing binary at {PATHS.EXECUTABLE_PATH}")
        return PATHS.EXECUTABLE_PATH

    # Check for precompiled binary
    precompiled = find_precompiled_binary()
    if precompiled is not None:
        if verbose:
            print(f"Found precompiled binary at {precompiled}")

        # Copy to package build directory for consistency
        package_build_dir = PATHS.PACKAGE_ROOT / "build"
        package_build_dir.mkdir(parents=True, exist_ok=True)
        dest_path = package_build_dir / precompiled.name
        shutil.copy2(precompiled, dest_path)

        if verbose:
            print(f"Copied to {dest_path}")
        return dest_path

    # No precompiled binary found, build from source
    if verbose:
        print("No precompiled binary found, compiling from source...")
    return build_binary(force=False, verbose=verbose)
