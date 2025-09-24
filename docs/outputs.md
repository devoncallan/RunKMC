
## Simulation Outputs

Every RunKMC simulation will output these files:
- results.csv
- input.txt
- metadata.yaml
- sequence.csv (optional)
- polymers.dat (optional)

`input.txt` is a copy of the input file used for the simulation.

`results.csv` contains information about the KMC state:
* Time
* Counts of all species
* Conversion of all unit species
* Chain/molecular weight distribution averages
* Sequence length distribution averages

`metadata.yaml` contains information about species and reactions and the information that RunKMC assigns to them. This helps with the processing of the results.

`sequences.csv` contains detailed sequence statistics across all polymer chains over the course of the simulation. The sequence statistics are discretized along the polymer chain into `Buckets`.

`polymers.dat` contains the full sequence information at the end of simulation. Each monomer is represented by its ID which can be found in the metadata.