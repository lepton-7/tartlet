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
@click.option(
    "--low-lim",
    default=-1.2,
    show_default=True,
    help="Lower bound for the peak plot y-axes.",
)
@click.option(
    "--up-lim",
    default=0.5,
    show_default=True,
    help="Upper bound for the peak plot y-axes.",
)
def exec_main(peak_log, cluster_stats, name, output_path, low_lim, up_lim):
    main(peak_log, cluster_stats, name, output_path, low_lim, up_lim)


def main(peak_log, cluster_stats, name, output_path, low_lim, up_lim):
    plotting.PeakPlot(
        abs_path=plotting.__file__,
        plog_path=peak_log,
        cstats_path=cluster_stats,
        out_path=output_path,
        name=name,
        plot_lim_low=low_lim,
        plot_lim_high=up_lim,
    )
    print(f"Saved peak plot for {name} to {output_path}")


# %%
