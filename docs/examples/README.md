# Examples

### Running on the command line
After pip installing **RunKMC**, you can run the core C++ engine on the commandline. There are two required arguments:
* `input.txt`: The path to the input file
* `output/`: The output directory for the results

There are also two optional command line arguments:

* `--report-polymers`: outputs as `polymers.dat`. Contains the full polymer sequence.
* `--report-sequences`: outputs as `sequences.csv`. Contains positional sequence information along the polymer chain. 

**Example:**
* `runkmc CRP3_Example.txt output/ --report-sequences`

### Running through python

The file `run.py` in this directory has example code for how to run the above code in python.

### Running with SPaRKS

**SPaRKS** imports **RunKMC** into a more unified interface with deterministic models. Examples can be found in the **SPaRKS** [repository](https://github.com/devoncallan/sparks) and the supporting manuscript [repository](https://github.com/devoncallan/ReversibleCopolymerizations).