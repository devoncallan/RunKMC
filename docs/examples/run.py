from pathlib import Path

import matplotlib.pyplot as plt

from runkmc.kmc import RunKMC
from runkmc.results import SimulationResult
from runkmc.vis.conversion import plot_monomer_conversion

EXAMPLES_DIR = Path(__file__).parent.resolve()

if __name__ == "__main__":

    from runkmc import RunKMC

    input_filepath = EXAMPLES_DIR / "CRP1_Example.txt"
    output_dir = EXAMPLES_DIR / "output_dir/"

    kmc = RunKMC(base_dir=output_dir)
    result = kmc.run_from_file(
        input_filepath, report_polymers=True, report_sequences=True
    )

    ax = plot_monomer_conversion(result)
    fig = plt.gcf()

    IMG_PATH = output_dir / "conversion.png"
    fig.savefig(IMG_PATH, dpi=300)
