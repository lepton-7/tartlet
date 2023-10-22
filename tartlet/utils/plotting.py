import numpy as np
import matplotlib.pyplot as plt

from math import ceil
from matplotlib.axes import Axes
from tart.utils.read_parsing import AlignDat


class CoveragePlot:
    """Handle plotting of alignment data."""

    def __init__(self, alignDat: AlignDat, end_buffers: list[int]) -> None:
        self._dat = alignDat

        self.lbuff = end_buffers[0]
        self.rbuff = end_buffers[1]

        # x-axis array used in plots from the object
        self.x = [i for i in range(len(self._dat.readcov))]
        # Handle case where bin size is 1nt
        self.bin_x = self._dat.bin_ax if self._dat.bin_size > 1 else self.x

        # Start of the plot window
        self.buffstart = self._dat.switch_start - self.lbuff
        self.buffstart = 0 if self.buffstart < 0 else self.buffstart

        # End of the plot window
        self.buffend = self._dat.switch_end + self.rbuff
        self.buffend = self.buffend if self.buffend < len(self.x) else len(self.x) - 1

        # Plot window bounds used for binned arrays
        self.buffstart_bin = self.buffstart // self._dat.bin_size
        self.buffend_bin = ceil(self.buffend / self._dat.bin_size)

        # Use reasonable x ticks
        self.xticks = (
            self._dat.bin_ax
            if self._dat.bin_size >= 10
            else [i for i in range(0, len(self._dat.readcov), 10)]
        )

        # Set ticks that are in the plot frame
        self.xticks = (
            self.xticks[self.buffstart_bin : self.buffend_bin]
            if self._dat.bin_size >= 10
            else self.xticks[self.buffstart // 10 : ceil(self.buffend / 10)]
        )

    def _binned_ends_panel(self, ax: Axes):
        """Add the binned raw ends panel to the figure.

        Args:
            ax (Axes): Panel Axes.
        """
        ax.bar(
            self.bin_x[self.buffstart_bin : self.buffend_bin],
            self._dat.rawends[self.buffstart_bin : self.buffend_bin],
            color="slateblue",
            width=float(self._dat.bin_size),
            align="edge",
        )
        ax.set_title(f"Inferred fragment ends ({self._dat.bin_size}nt bins)")
        ax.set_xticks(self.xticks)
        ax.set_xlabel("Nucleotide position (bp)")

        bott, top = ax.get_ylim()
        ax.set_ylabel("Count")

        a_height = (top - bott) * 0.05

        ax.annotate(
            "",
            xy=(self._dat.switch_start, 0),
            xytext=(self._dat.switch_start, a_height),
            arrowprops=dict(facecolor="black"),
            annotation_clip=False,
        )
        ax.annotate(
            "",
            xy=(self._dat.switch_end, 0),
            xytext=(self._dat.switch_end, a_height),
            arrowprops=dict(facecolor="black"),
            annotation_clip=False,
        )

    def _coverage_panel(self, ax: Axes):
        """Add the coverage panel to the figure.

        Args:
            ax (Axes): Panel Axes.
        """
        coverage_counts = {
            "Read": self._dat.readcov,
            "Inferred": self._dat.infercov,
            "Clipped": self._dat.clipcov,
        }
        coverage_colours = {
            "Read": "slateblue",
            "Inferred": "crimson",
            "Clipped": "mediumseagreen",
        }
        bottom = np.zeros(len(self.x))

        for type, count in coverage_counts.items():
            ax.bar(
                self.x[self.buffstart : self.buffend],
                count[self.buffstart : self.buffend],
                label=type,
                bottom=bottom[self.buffstart : self.buffend],
                color=coverage_colours[type],
                width=1,
                align="edge",
            )

            bottom += count

        ax.set_title("Fragment coverage")
        ax.legend(loc="upper right")
        ax.set_ylabel("Count")

    def default(self, save_path: str):
        """Generate the default reference alignment plot and save to file.

        Args:
            save_path (str): Save path. Parent directories must exist.
        """
        fig, ax = plt.subplots(
            2, 1, sharex=True, figsize=(20, 10), dpi=100, constrained_layout=True
        )
        fig.suptitle(f"{self._dat.ref}")

        self._coverage_panel(ax[0])
        self._binned_ends_panel(ax[1])

        fig.savefig(f"{save_path}")

        plt.close()
