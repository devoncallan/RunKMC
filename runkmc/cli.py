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

from .kmc.execution import execute_simulation, CommandLineConfig


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
            runkmc input.txt output/ --parse-only
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
        "--parse-only",
        action="store_true",
        help="Only parse the input file and generate input files without running the simulation",
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

    execute_simulation(
        input_filepath=args.input_file,
        output_dir=args.output_dir,
        report_polymers=args.report_polymers,
        report_sequences=args.report_sequences,
        parse_only=args.parse_only,
    )


if __name__ == "__main__":
    main()
