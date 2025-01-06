# %%
import click

from pathlib import Path
from tartlet.utils import plotting
from tartlet.utils.utils import print


@click.command()
@click.option(
    "-p",
    "--peak-log",
    required=True,
    help="Path to the peak_log output file.",
)
@click.option(
    "-c",
    "--cluster-stats",
    required=True,
    help="Path to the cluster_stats output file.",
)
@click.option(
    "--name",
    required=True,
    help="Name used to label the plot.",
)
@click.option(
    "-o",
    "--output-path",
    required=True,
    help="Path and filename of the generated plot. Must have the '.png' suffix.",
)
def main(peak_log, cluster_stats, name, output_path):
    plotting.PeakPlot(plotting.__file__, peak_log, cluster_stats, output_path, name)
    print(f"Saved peak plot for {name} to {output_path}")


# %%
