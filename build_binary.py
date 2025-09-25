#!/usr/bin/env python3
import subprocess
import sys
import os
import glob
import shutil


def main():
    # Build with cmake
    subprocess.run(
        [
            "cmake",
            "-B",
            "build",
            "-S",
            "cpp",
            "-DRUNKMC_VERSION=0.1.1",
            "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
        ],
        check=True,
    )

    subprocess.run(["cmake", "--build", "build"], check=True)

    # Create destination directory
    os.makedirs("runkmc/build", exist_ok=True)

    # Find the binary (Windows puts it in Debug/ subfolder)
    if sys.platform == "win32":
        pattern = "build/**/RunKMC.exe"
    else:
        pattern = "build/**/RunKMC"

    matches = glob.glob(pattern, recursive=True)
    if not matches:
        raise FileNotFoundError(f"Could not find binary matching {pattern}")

    src = matches[0]
    print(f"Found binary: {src}")

    # Copy to destination
    shutil.copy2(src, "runkmc/build/")
    print(f"Copied to runkmc/build/")


if __name__ == "__main__":
    main()
