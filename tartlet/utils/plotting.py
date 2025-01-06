from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import rpy2.robjects as robjects

from math import ceil
from os import getcwd
from matplotlib.axes import Axes
from rpy2.robjects import pandas2ri
from tartlet.utils.utils import print
from matplotlib.patches import Rectangle
from tartlet.utils.read_parsing import AlignDat, SegregatedAlignDat
from tartlet.utils.activity_inference import _gen_kernel


class PeakPlot:
    """Handle plotting tool outputs"""

    def __init__(
        self,
        abs_path: Path | str,
        plog_path: Path | str,
        cstats_path: Path | str,
        out_path: Path | str,
        name: str,
    ):
        r = robjects.r
        r["source"](f"{Path(abs_path).parent}/peak_plotting.r")

        self.peakplotmaker_r = robjects.globalenv["peakplotmaker"]

        self.plog_path = Path(plog_path)
        self.cstats_path = Path(cstats_path)
        self.out_path = Path(out_path)
        self.name = name

        # FIXME: add input validation for paths and filetypes
        self.out_path.parent.mkdir(parents=True, exist_ok=True)

        self._make_peak_plot()

    def _make_peak_plot(self):
        self.peakplotmaker_r(
            str(self.plog_path), str(self.cstats_path), str(self.out_path), self.name
        )


class CoveragePlot:
    """Handle plotting of alignment data."""

    def __init__(
        self,
        alignDat: SegregatedAlignDat,
        end_buffers: list[int],
        palette: str = "heavypastel",
    ) -> None:
        self.palette_picker = {
            "apricots": {
                "Read": "#68904D",
                "Overlapped": "#14471E",
                "Inferred": "#DA6A00",
                "Clipped": "#C8D2D1",
                "axback": "#FFFFFF",
                "figback": "#FFFFFF",
            },
            "pastelbluepeach": {
                "Read": "#006D77",
                "Overlapped": "#FFDDD2",
                "Inferred": "#E29578",
                "Clipped": "#FFDDD2",
                "axback": "#FFFFFF",
                "figback": "#FFFFFF",
            },
            "heavypastel": {
                "Read": "#3D5A80",
                "Overlapped": "#293241",
                "Inferred": "#EE6C4D",
                "Clipped": "#98C1D9",
                # "Termini": "#698F3F",
                "Termini": "#F00699",
                "axback": "#FFFFFF",
                "figback": "#FFFFFF",
            },
            "figs": {
                "Read": "#483948",
                "Overlapped": "#180C0C",
                "Inferred": "#ED413E",
                "Clipped": "#B2AB2E",
                "axback": "#FFFFFF",
                "figback": "#FFFFFF",
            },
            "asterpink": {
                "Read": "#ab3777",
                "Overlapped": "#721836",
                "Inferred": "#efb420",
                "Clipped": "#373933",
                "axback": "#FFFFFF",
                "figback": "#FFFFFF",
            },
            "warmpastel": {
                "Read": "#C5998C",
                "Overlapped": "#4E4035",
                "Inferred": "#D67F54",
                "Clipped": "#DDE2DC",
                "axback": "#FFFFFF",
                "figback": "#DDE2DC",
            },
        }

        # Named colours:
        # https://matplotlib.org/stable/_images/sphx_glr_named_colors_003_2_00x.png
        self.palette = self.palette_picker[palette]

        self._dat = alignDat.total
        self._segregated = alignDat

        self.lbuff = end_buffers[0]
        self.rbuff = end_buffers[1]

        # Gaps between the tick labels if bins are small
        self.tick_sep = 20

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
            if self._dat.bin_size >= self.tick_sep
            else [i for i in range(0, len(self._dat.readcov), self.tick_sep)]
        )

        # Set ticks that are in the plot frame
        self.xticks = (
            self.xticks[self.buffstart_bin : self.buffend_bin]
            if self._dat.bin_size >= self.tick_sep
            else self.xticks[
                self.buffstart // self.tick_sep : ceil(self.buffend / self.tick_sep)
            ]
        )

        self.axis_label_args = {"fontsize": "xx-large"}
        self.axis_label_coverage = "Mapped bases"

        self.panel_title_size = 22
        self.tick_lab_size = 16

    def _binned_ends_panel(self, ax: Axes):
        """Add the binned raw ends panel to the figure.

        Args:
            ax (Axes): Panel Axes.
        """
        ax.bar(
            self.bin_x[self.buffstart_bin : self.buffend_bin],
            self._dat.binned_ends[self.buffstart_bin : self.buffend_bin],
            color=self.palette["Read"],
            width=float(self._dat.bin_size),
            align="edge",
        )
        ax.set_facecolor(self.palette["axback"])
        ax.set_title(
            f"Inferred fragment ends ({self._dat.bin_size}nt bins)",
            fontsize=self.panel_title_size,
        )
        ax.set_xticks(self.xticks)
        ax.xaxis.set_tick_params(labelsize=self.tick_lab_size)
        ax.yaxis.set_tick_params(labelsize=self.tick_lab_size)
        ax.set_xlabel("Nucleotide position (bp)", self.axis_label_args)

        ax.set_ylabel(self.axis_label_coverage, self.axis_label_args)

        self._add_switch_arrows(ax)
        self._add_boxes(ax)

    def _conv_ends_panel(self, ax: Axes):
        """Add the raw end convolution panel to the figure.

        Args:
            ax (Axes): Panel Axes.
        """
        try:
            conv_ends = self._dat.convends
        except AttributeError:
            k = _gen_kernel(51, 1.5)
            self._dat.convolve_rawends(k)

        ax.bar(
            self.x[self.buffstart : self.buffend],
            self._dat.convends[self.buffstart : self.buffend],
            color=self.palette["Read"],
            width=1,
            align="edge",
        )
        # ax.set_facecolor(self.palette["axback"])
        ax.set_title(
            f"Convolved inferred fragment ends", fontsize=self.panel_title_size
        )
        ax.set_xticks(self.xticks)
        ax.xaxis.set_tick_params(labelsize=self.tick_lab_size)
        ax.yaxis.set_tick_params(labelsize=self.tick_lab_size)
        ax.set_xlabel("Nucleotide position (bp)", self.axis_label_args)

        # The y axis doesn't really mean anything?
        # ax.set_ylabel("Count")

        self._add_switch_arrows(ax)
        self._add_boxes(ax)

    def _add_switch_arrows(self, ax: Axes):
        """Add little arrows to mark the bounds of the riboswitch.

        Args:
            ax (Axes): Panel Axes.
        """
        bott, top = ax.get_ylim()
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

    def _add_boxes(self, ax: Axes):
        """Add boxes to indicate functional unit bounds.

        Args:
            ax (Axes): Panel Axes.
        """
        bott, top = ax.get_ylim()
        switch_y = top * 0.8
        orf_y = top * 0.65
        box_height = top * 0.1

        # Riboswitch box
        ax.add_patch(
            Rectangle(
                xy=(self._dat.switch_start, switch_y),
                width=self._dat.switch_end - self._dat.switch_start,
                height=box_height,
                facecolor="xkcd:crimson",
                alpha=0.8,
                zorder=0.0,
            )
        )

        # ORF box
        ax.add_patch(
            Rectangle(
                xy=(self._dat.orf_start, orf_y),
                width=self._dat.orf_end - self._dat.orf_start,
                height=box_height,
                facecolor="green",
                alpha=0.8,
                zorder=0.0,
            )
        )

    def _coverage_panel(self, ax: Axes):
        """Add the coverage panel to the figure.

        Args:
            ax (Axes): Panel Axes.
        """
        coverage_counts = {
            "Inferred": self._dat.infercov,
            "Overlapped": self._dat.overlapcov,
            "Read": self._dat.readcov,
            "Termini": self._dat.terminalcov,
            "Clipped": self._dat.clipcov,
        }

        bottom = np.zeros(len(self.x))

        for type, count in coverage_counts.items():
            ax.bar(
                self.x[self.buffstart : self.buffend],
                count[self.buffstart : self.buffend],
                label=type,
                bottom=bottom[self.buffstart : self.buffend],
                color=self.palette[type],
                width=1,
                align="edge",
            )

            bottom += count

        ax.set_facecolor(self.palette["axback"])
        ax.set_title("Fragment coverage", fontsize=self.panel_title_size)
        ax.legend(loc="upper left")
        ax.set_ylabel(self.axis_label_coverage, self.axis_label_args)
        ax.yaxis.set_tick_params(labelsize=self.tick_lab_size)

    def default_depr(self, save_path: str | Path):
        """Generate the default reference alignment plot and save to file.

        Args:
            save_path (str): Save path. Parent directories must exist.
        """
        fig, ax = plt.subplots(
            2,
            1,
            sharex=True,
            figsize=(20, 10),
            dpi=100,
            constrained_layout=True,
            facecolor=self.palette["figback"],
        )
        fig.suptitle(f"{self._dat.ref}")

        self._coverage_panel(ax[0])
        self._binned_ends_panel(ax[1])

        fig.savefig(f"{save_path}")

        plt.close()

    def default(self, save_path: str | Path, with_conv=False):
        # TODO: docstring; Implement this properly
        fig, ax = plt.subplots(
            4,
            1,
            sharex=True,
            figsize=(20, 10),
            dpi=100,
            constrained_layout=True,
            facecolor=self.palette["figback"],
        )
        fig.suptitle(f"{self._dat.ref}", fontsize=self.panel_title_size)

        self._coverage_panel(ax[2])
        self._split_coverage_panels(ax[0], ax[1])
        self._conv_ends_panel(ax[-1]) if with_conv else self._binned_ends_panel(ax[-1])

        fig.savefig(f"{save_path}")

        plt.close()

    def distribution_plots(self, save_path: str | Path):
        """Statistic plot for a riboswitch alignment data object.

        Args:
            save_path (str | Path): Save path. Parent directories must exist..
        """
        fig, ax = plt.subplots(figsize=(20, 10), dpi=100)

        fig.suptitle(f"{self._dat.ref} fragments in switch")

        self._frags_in_switch_hist(ax, 20)

        fig.savefig(f"{save_path}")

        plt.close()

    def _frags_in_switch_hist(self, ax: Axes, bw: int):
        """Generate the binned distribution of fragments aligning to the riboswitch sequence.

        Args:
            ax (Axes): Panel Axes.
            bw (int): bin width.
        """
        ax.hist(
            self._dat.tlen_in_switch,
            bins=list(
                np.arange(
                    min(self._dat.tlen_in_switch),
                    max(self._dat.tlen_in_switch) + bw,
                    bw,
                )
            ),
            density=True,
        )

        ax.set_xlim((0, 1000))

    def _split_coverage_panels(self, ax3p: Axes, axnot3p: Axes):
        counts3p = {
            "Inferred": self._segregated.of_3prime.infercov,
            "Overlapped": self._segregated.of_3prime.overlapcov,
            "Read": self._segregated.of_3prime.readcov,
            "Termini": self._segregated.of_3prime.terminalcov,
            "Clipped": self._segregated.of_3prime.clipcov,
        }

        countsno3p = {
            "Inferred": self._segregated.except_3prime.infercov,
            "Overlapped": self._segregated.except_3prime.overlapcov,
            "Read": self._segregated.except_3prime.readcov,
            "Termini": self._segregated.except_3prime.terminalcov,
            "Clipped": self._segregated.except_3prime.clipcov,
        }

        bottom3p = np.zeros(len(self.x))

        for type, count in counts3p.items():
            ax3p.bar(
                self.x[self.buffstart : self.buffend],
                count[self.buffstart : self.buffend],
                label=type,
                bottom=bottom3p[self.buffstart : self.buffend],
                color=self.palette[type],
                width=1,
                align="edge",
            )

            bottom3p += count

        ax3p.set_facecolor(self.palette["axback"])
        ax3p.set_title("Fragment coverage", fontsize=self.panel_title_size)
        ax3p.legend(loc="upper left")
        ax3p.set_ylabel(self.axis_label_coverage, self.axis_label_args)
        ax3p.yaxis.set_tick_params(labelsize=self.tick_lab_size)

        bottomno3p = np.zeros(len(self.x))

        for type, count in countsno3p.items():
            axnot3p.bar(
                self.x[self.buffstart : self.buffend],
                count[self.buffstart : self.buffend],
                label=type,
                bottom=bottomno3p[self.buffstart : self.buffend],
                color=self.palette[type],
                width=1,
                align="edge",
            )

            bottomno3p += count

        axnot3p.set_facecolor(self.palette["axback"])
        axnot3p.set_title("Fragment coverage", fontsize=self.panel_title_size)
        axnot3p.legend(loc="upper left")
        axnot3p.set_ylabel(self.axis_label_coverage, self.axis_label_args)
        axnot3p.yaxis.set_tick_params(labelsize=self.tick_lab_size)

    def set_palette(self, palette):
        self.palette = self.palette_picker[palette]
