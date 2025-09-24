# RunKMC - A High-Performance Kinetic Monte Carlo Simulator for Polymerization

**RunKMC** is a kinetic Monte Carlo simulation engine written in C++ with a thin Python wrapper for accessibility. The core engine is also accessible through the command line.

[![PyPI](https://badge.fury.io/py/runkmc.svg)](https://badge.fury.io/py/runkmc) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17172027.svg)](https://doi.org/10.5281/zenodo.17172027) [![Python](https://img.shields.io/badge/Python-%3E%3D3.10-blue?logo=python&logoColor=yellow)](https://python.org)

## Quick install
The easiest way to install **RunKMC** is via pip:

```shell
pip3 install runkmc
```

For local development, we recommend using cmake:

```shell
pip3 install cmake
make setup
```

We recommend setting a python virtual environment before installing the package.
```shell
python3 -m venv .venv
source .venv/bin/activate
```

## Examples

Documentation for RunKMC concepts can be found [here](docs/). Example input files can be found in [examples](docs/examples/README.md). More relevant examples and integrations with [**SPaRKS**🔗](https://github.com/devoncallan/sparks) can be found at the supporting data for the manuscript below. This can be found at [this repository](https://github.com/devoncallan/ReversibleCopolymerizations): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17172075.svg)](https://doi.org/10.5281/zenodo.17172075)

## How to Cite

**Citeable DOI for this version of RunKMC:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17172027.svg)](https://doi.org/10.5281/zenodo.17172027)

If you use RunKMC in your project, please cite
* Callan, D. H., and Bates, C. M. Efficient Deterministic Modeling of Reversile Copolymerizations, *Macromolecules*, **2025**, [doi:10.1021/acs.macromol.5c01421](https://doi.org/10.1021/acs.macromol.5c01421)
