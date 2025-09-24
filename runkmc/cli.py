#!/usr/bin/env python3
"""
Command-line interface for RunKMC.

This provides a direct command-line tool that users can run as:
    runkmc input.txt output_dir [options]

Instead of having to use the Python API.
"""

import sys
import argparse
from pathlib import Path
from typing import Optional

from .kmc.execution import execute_simulation


def create_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="runkmc",
        description="Run Kinetic Monte Carlo simulations for polymerization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  runkmc input.txt output/
  runkmc input.txt output/ --report-polymers
  runkmc input.txt output/ --report-polymers --report-sequences
        """,
    )

    parser.add_argument(
        "input_file", type=Path, help="Path to the input configuration file"
    )

    parser.add_argument(
        "output_dir", type=Path, help="Directory where simulation results will be saved"
    )

    parser.add_argument(
        "--report-polymers",
        action="store_true",
        help="Generate detailed polymer structure reports",
    )

    parser.add_argument(
        "--report-sequences",
        action="store_true",
        help="Generate sequence analysis reports",
    )

    parser.add_argument(
        "--version", action="version", version=f"runkmc {get_version()}"
    )

    return parser


def get_version() -> str:
    """Get the package version."""
    try:
        from . import __version__

        return __version__
    except ImportError:
        return "unknown"


def main() -> None:
    """Main entry point for the runkmc command-line tool."""
    parser = create_parser()
    args = parser.parse_args()

    # Validate input file exists
    if not args.input_file.exists():
        print(f"Error: Input file '{args.input_file}' does not exist", file=sys.stderr)
        sys.exit(1)

    if not args.input_file.is_file():
        print(f"Error: '{args.input_file}' is not a file", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)

    try:
        print(f"Running KMC simulation...")
        print(f"Input:  {args.input_file.absolute()}")
        print(f"Output: {args.output_dir.absolute()}")

        execute_simulation(
            input_filepath=args.input_file,
            output_dir=args.output_dir,
            report_polymers=args.report_polymers,
            report_sequences=args.report_sequences,
        )

        print("Simulation completed successfully!")

    except Exception as e:
        print(f"Error: Simulation failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
