#!/usr/bin/env python3
"""Standalone script to build the RunKMC binary.

This is a thin wrapper around runkmc.kmc.build.build_binary() for
convenience when building from the command line or in CI/CD.
"""

from runkmc.kmc.build import build_binary


def main():
    """Build the RunKMC binary from source."""
    build_binary(force=True, verbose=True)


if __name__ == "__main__":
    main()
