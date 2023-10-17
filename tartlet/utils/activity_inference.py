import numpy as np
import operator

import matplotlib.pyplot as plt

from scipy.stats import norm, ks_2samp
from math import ceil, isclose


class Peak:
    """Handles peaks in convolved data."""

    def __init__(
        self,
        from_switch_end,
        center=None,
        height=None,
        half_width=None,
        switch_end=None,
    ):
        self.switch_end = switch_end

        self.from_switch_end = None
        self.from_switch_end = (
            center - switch_end if switch_end is not None else self.from_switch_end
        )
        self.from_switch_end = (
            from_switch_end if from_switch_end is not None else self.from_switch_end
        )

        self.center = center
        self.height = height
        self.half_width = half_width

    def find_half_width(self, source_arr):
        """Find the half-max width of the peak from the source array.

        Args:
            source_arr (list or array): The source array containing the peak.

        Raises:
            AttributeError: If the Peak object is not instantiated with the center.
            AttributeError: If the Peak object is not instantiated with the height.

        Returns:
            bool: Whether the half-max width was successfully determined.
        """
        if self.center is None:
            raise AttributeError("Peak center needs to be defined")

        if self.height is None:
            raise AttributeError("Peak height needs to be defined")

        halfmax = abs(self.height / 2)

        i = 1
        minReached = abs(source_arr[self.center])
        checkLeft = True
        checkRight = True
        while checkLeft or checkRight:
            # Check if left descent valid
            if checkLeft:
                try:
                    cur = abs(source_arr[self.center - i])
                    prev = abs(source_arr[self.center - i + 1])

                    if cur - prev <= 0:
                        minReached = cur if cur < minReached else minReached

                    # Not descending anymore
                    else:
                        checkLeft = False

                except IndexError:
                    checkLeft = False

            # Check if right descent valid
            if checkRight:
                try:
                    cur = abs(source_arr[self.center + i])
                    prev = abs(source_arr[self.center + i - 1])

                    if cur - prev <= 0:
                        minReached = cur if cur < minReached else minReached

                    # Not descending anymore
                    else:
                        checkLeft = False

                except IndexError:
                    checkRight = False

            if minReached <= halfmax:
                self.half_width = i * 2
                return True

            i += 1

        return False

    def compare(self, peak: "Peak", heightMarg: float = 0.1, widthMarg: int = 4):
        """DEPRECATED

        Compares two peaks. Returns true if peaks match within the defined margins.

        Args:
            peak (Peak): Comparison subject Peak object.
            heightMarg (float, optional): Proportion of this Peak's margin within which to match peak. Defaults to 0.1.
            widthMarg (int, optional): Nucleotide margin within which to match peak's width. Defaults to 4.

        Returns:
            bool: Whether the height and widths of the peak being compared fall wlthin the margins of this peak.
        """
        hprop = peak.height / self.height
        wdiff = abs(self.half_width - peak.half_width)

        return (
            hprop > -1 - heightMarg and hprop < -1 + heightMarg and wdiff <= widthMarg
        )

    def compare_ks(self, peak: "Peak", source_arr):
        """Uses the 2-sample KS test to compare fragment end distributions between peaks.

        Args:
            peak (Peak): Comparison subject peak.
            source_arr (np.array): Array from which to extract fragment end distributions for self and subject peak.

        Returns:
            tuple: (KS_stats, self ends, peak ends, peak center)
        """
        this_ends = self.find_absolute_ends(source_arr)
        peak_ends = peak.find_absolute_ends(source_arr)

        return (ks_2samp(this_ends, peak_ends), this_ends, peak_ends, peak.center)

    def find_absolute_ends(self, source_arr):
        """Finds the distribution of raw, unconvolved ends across the range of the peak.

        Args:
            source_arr (list or np.array): The unconvolved ends array.

        Returns:
            list: A slice of the source array covered by the peak center ± half-width.
        """
        try:
            return self.abs_ends

        except AttributeError:
            left = self.center - int(self.half_width)
            right = self.center + int(self.half_width)

            # Validate
            left = 0 if left < 0 else left
            right = len(source_arr) - 1 if right >= len(source_arr) else right

            self.abs_ends = np.absolute(source_arr[left:right])
            return self.abs_ends

    def find_absolute_coverage_delta(self, sumcov: list) -> int:
        """Computes the absolute change in the summed coverage between the start and end of the peak (center ± half-width).

        Returns:
            int: Absolute change in summed coverage across the peak.
        """
        try:
            return self.abs_cov_delta

        except AttributeError:
            left = self.center - int(self.half_width)
            right = self.center + int(self.half_width)

            # Validate
            left = 0 if left < 0 else left
            right = len(sumcov) - 1 if right >= len(sumcov) else right

            self.abs_cov_delta = sumcov[right] - sumcov[left]
            return self.abs_cov_delta

    def find_relative_coverage_delta(self, sumcov: list) -> float:
        """Computes the change in summed coverage across the peak (center ± half-width) relative to the summed coverage at the lower bound (center - half-width) of the peak.

        Args:
            sumcov (list): Coverage summed across actual, inferred, and clipped reads per base.

        Returns:
            float: Change in summed coverage across the peak relative to coverage at the start of the peak.
        """
        try:
            return self.rel_cov_delta
        except AttributeError:
            l = self.center - int(self.half_width)
            l = 0 if l < 0 else l

            self.rel_cov_delta = self.find_absolute_coverage_delta(sumcov) / sumcov[l]
            return self.rel_cov_delta


class Candidate(Peak):
    """Wraps the relevant information to return as a candidate for transcription termination activity.

    Args:
        Peak (Peak): Parent class that handles most of the attribute setting.
    """

    def __init__(
        self,
        cand: Peak,
        switch_size: int,
        nocand_cov_delta: list,
        sumcov: list,
    ) -> None:
        super().__init__(
            from_switch_end=cand.from_switch_end,
            center=cand.center,
            height=cand.height,
            half_width=cand.half_width,
        )

        self.abs_cov_delta = cand.find_absolute_coverage_delta(sumcov)
        self.rel_cov_delta = cand.find_relative_coverage_delta(sumcov)

        self.switch_size = switch_size
        self.from_switch_end_relative = self.from_switch_end / self.switch_size

        self.coverage_delta_noise = nocand_cov_delta


def bin_counts(alignTup: tuple, bin_size: int = 1):
    """Bins then sums fragment ends within bins.

    Args:
        alignTup (tuple): Alignment tuple: ([*coverage], frag ends, (switch start, switch end)).
            The coverage list should contain the read coverage, inferred fragment coverage, and
            clipped coverage arrays in that order
        bin_size (int, optional): Binning size. Defaults to 1.

    Returns:
        _type_: _description_
    """
    cov, ends, (switch_start, switch_end) = alignTup

    readcov = cov[0]

    bin_pos = [i for i in range(0, len(readcov), bin_size)]
    numbins = len(bin_pos)

    binned_ends = np.ones(numbins)
    for i in range(numbins - 1):
        bstart = bin_pos[i]
        bend = bin_pos[i + 1]

        binned_ends[i] = sum(ends[bstart:bend])

    binned_ends[numbins - 1] = sum(ends[bin_pos[numbins - 1] :])

    return (cov, binned_ends, (switch_start, switch_end)), bin_pos


def plot_gen(
    ref: str,
    alignTup: tuple,
    save_path: str,
    buff: int = 40,
    bin_size: int = 1,
    bin_ax=None,
):
    """Generate a switch alignment plot and save to file.

    Args:
        ref (str): Name of the alignment reference for the plot.
        alignTup (tuple): Alignment tuple with the read coverage, inferred frag coverage,
            clipped coverage, fragment ends, and switch start and end
        save_path (str): Full save path for the plot
        buff (int, optional): Buffer around the riboswitch to show in the plot. Defaults to 40.
    """
    # The bin start (inclusive) values
    if bin_size > 1 and bin_ax is None:
        raise ValueError("Need to pass x values for binned plot")

    cov, frag_ends, (start, end) = alignTup
    readcov, infercov, clipcov = cov

    x = [i for i in range(len(readcov))]

    bin_x = bin_ax if bin_size > 1 else x

    # Use reasonable x ticks
    xticks = bin_ax if bin_size >= 10 else [i for i in range(0, len(readcov), 10)]

    fig, ax = plt.subplots(
        2, 1, sharex=True, figsize=(20, 10), dpi=100, constrained_layout=True
    )

    fig.suptitle(f"{ref}")

    buffstart = start - buff
    buffstart = 0 if buffstart < 0 else buffstart

    buffend = end + buff

    buffstart_bin = buffstart // bin_size
    buffend_bin = ceil(buffend / bin_size)

    # Select ticks that are in the plot frame
    xticks = (
        xticks[buffstart_bin:buffend_bin]
        if bin_size >= 10
        else xticks[buffstart // 10 : ceil(buffend / 10)]
    )

    # Add the coverage panel ----------------------------------------
    coverage_counts = {"Read": readcov, "Inferred": infercov, "Clipped": clipcov}
    coverage_colours = {
        "Read": "slateblue",
        "Inferred": "crimson",
        "Clipped": "mediumseagreen",
    }
    bottom = np.zeros(len(x))

    for type, count in coverage_counts.items():
        ax[0].bar(
            x[buffstart:buffend],
            count[buffstart:buffend],
            label=type,
            bottom=bottom[buffstart:buffend],
            color=coverage_colours[type],
            width=1,
            align="edge",
        )

        bottom += count

    ax[0].set_title("Fragment coverage")
    ax[0].legend(loc="upper right")
    ax[0].set_ylabel("Count")

    # Add the ends panel -------------------------------------------
    ax[1].bar(
        bin_x[buffstart_bin:buffend_bin],
        frag_ends[buffstart_bin:buffend_bin],
        color="slateblue",
        width=float(bin_size),
        align="edge",
    )
    ax[1].set_title(f"Inferred fragment ends ({bin_size}nt bins)")
    ax[1].set_xticks(xticks)
    ax[1].set_xlabel("Nucleotide position (bp)")

    bott, top = ax[1].get_ylim()
    ax[1].set_ylabel("Count")

    a_height = (top - bott) * 0.05

    ax[1].annotate(
        "",
        xy=(start, 0),
        xytext=(start, a_height),
        arrowprops=dict(facecolor="black"),
        annotation_clip=False,
    )
    ax[1].annotate(
        "",
        xy=(end, 0),
        xytext=(end, a_height),
        arrowprops=dict(facecolor="black"),
        annotation_clip=False,
    )

    fig.savefig(f"{save_path}")

    plt.close()


def find_peaks(arr, switch_end, l: int = 0, u: int = -1):
    """Finds peaks in the given convolved array.

    Args:
        arr (list or np.array): Data array/list convolved with a normal-distr weighted kernel.
        switch_end (int): Position of the riboswitch end.
        l (int, optional): Left bound of the region to find peaks in. Defaults to 0.
        u (int, optional): Right bound of the region to find peaks in. Setting to -1 indicates no bound.
                            Defaults to -1.

    Returns:
        list: List of valid Peak objects found within the specified bounds.
    """
    ret = []

    u = len(arr) if u < 0 else u

    for i in range(l + 1, u - 1):
        if abs(arr[i]) > abs(arr[i - 1]) and abs(arr[i]) > abs(arr[i + 1]):
            peak = Peak(center=i, height=arr[i], switch_end=switch_end)
            if peak.find_half_width(arr):
                ret.append(peak)

    return ret


def gen_kernel(kernel_size: int = 21, std_dev: float = 3.0):
    """Generates a convolution kernel for a normal pdf

    Args:
        kernel_size (int, optional): Size of the kernel. Only odd numbers. Defaults to 21.
        std_dev (float, optional): Standard deviation of the distribution with mean = 0. Defaults to 3.0.

    Returns:
        list: Normal PDF convolution kernel
    """
    num = kernel_size // 2
    k_in = [x for x in range(-num, num + 1)]

    kernel = norm.pdf(k_in, loc=0, scale=std_dev)

    return kernel


def check_cand_drop(
    cov: list, switchreg: tuple, cand: Peak, stdev: int, minDrop: float
):
    """DEPRECATED

    Checks whether coverage drop across (bounds determined by the standard
    deviation of the convolution kernel used to process the raw ends array)
    a candidate peak relative to the max coverage value in the riboswitch
    region is over the given threshold.

    Args:
        cov (list): Coverage list.
        switchreg (tuple): Riboswith bounds: (start, end)
        cand (Peak): Candidate Peak.
        stdev (int): Ends convolution kernel std dev.
        minDrop (float): Minimum relative drop threshold.

    Returns:
        bool: True if relative drop across the Peak exceeds the drop threshold.
    """
    # np.add can ONLY add two arrays.
    # The third param is the output array obj
    arr = np.add(cov[0], cov[1])
    arr = np.add(arr, cov[2])

    sstart, send = switchreg

    interv = 2 * stdev
    istart, iend = cand.center - interv, cand.center + interv

    istart = 0 if istart < 0 else istart
    iend = len(arr) - 1 if iend >= len(arr) else iend

    diff = arr[istart] - arr[iend]
    maxreads = max(arr[sstart:send])

    return diff / maxreads > minDrop


def is_interesting(
    alignTup: tuple,
    windowfrac: float = 0.15,
    threshtol: float = 0.15,
):
    """
    DEPRECATED

    Determines whether the read alignment is indicative of transcriptionally active riboswitches.

    Args:
        alignTup (tuple): Alignment tuple for this reference:
            ([read coverage, inferred frag cov, clipped cov], ends, (switch start, switch end))
        windowfrac (float, optional): Window size on either side of riboswitch 3' as a fraction
            of the total riboswitch size to check for fragment end peaks. Defaults to 15%.
        threshtol (float, optional): Allowable fractional percentage margin when checking
            if a peak height is close to the maximum peak height in the full
            riboswitch + window region. Defaults to 0.15.

    Returns:
        bool: True if alignTuple is determined to be interesting. False if not.
    """
    cov, rawends, (switch_left, switch_right) = alignTup
    readcov, infercov, clipcov = cov

    # OPTIONS -----------------------------------------------------------------
    kernel_size = 51
    kernel_stdev = 5

    minReadDrop = 0.3

    peakHeightTol = 0.15
    peakWidthMarg = 4
    # -------------------------------------------------------------------------

    kernel = gen_kernel(kernel_size=kernel_size, std_dev=kernel_stdev)
    ends = np.convolve(rawends, kernel, "same")

    # - strand riboswitches are reverse complemented during the reference generation,
    # so the right end in the reference is still the 3' end
    switch_end = switch_right

    peaks = find_peaks(ends)

    window = (int)((switch_right - switch_left) * windowfrac)
    left, right = switch_end - window, switch_end + window

    # Max peak from riboswitch 5' to the 3' window end
    maxpeak = max(ends[switch_left:right])

    for i in range(len(peaks)):
        cand: Peak = peaks[i]
        keepcand = True
        cand_stats = None

        # The candidate is within the window and within the margin of the tallest peak
        if (
            cand.center >= left
            and cand.center <= right
            and isclose(cand.height, maxpeak, rel_tol=threshtol)
            and check_cand_drop(
                cov,
                (switch_left, switch_right),
                cand,
                kernel_stdev,
                minDrop=minReadDrop,
            )
        ):
            cand_stats = []
            # Check if there is a similar mirrored peak
            # anywhere in the region of interest.
            # Only the frag start peaks prior to the candidate
            # frag end peak needs to be checked.
            for peak in peaks[:i]:
                # If the preceding peak is far out of the read range
                if peak.center < cand.center - 200:
                    continue
                # If there is a mirrored peak within margins,
                # immediately disqualify
                cand_stats.append(cand.compare_ks(peak, ends))
                if cand.compare(
                    peak, heightMarg=peakHeightTol, widthMarg=peakWidthMarg
                ):
                    keepcand = False
                    break

            # This candidate does not have a mirrored peak
            if keepcand:
                return (True, cand_stats)

    return (False, cand_stats)


def convolve_array(rawends, kernel_size, std_dev):
    """Helper function that generates and applies the normally-weighted kernal to the raw ends data.

    Args:
        rawends (list or np.array): The raw, unprocessed ends data array.
        kernel_size (int): Size of the 1-D convolution kernel. Must be odd?
        std_dev (float): Standard deviation used to generate the normally-weighted kernel.

    Returns:
        list: Convolved list of ends.
    """
    kernel = gen_kernel(kernel_size=kernel_size, std_dev=std_dev)
    ends = np.convolve(rawends, kernel, "same")

    return ends


def coverage_delta_per_peak(peaks: list, sumcov: list):
    """Helper function to find coverage deltas across a list of Peaks.

    Args:
        peaks (list): List of Peak objects.
        sumcov (list): Summed coverage list.

    Returns:
        list: Absolute coverage drop for each peak in the same order as peaks.
    """
    raw_cov_drop = []

    for peak in peaks:
        peak: Peak
        raw_cov_drop.append(peak.find_absolute_coverage_delta(sumcov))

    return raw_cov_drop


def peak_out_of_cov_delta(sorteddelta: list, i: int) -> (bool, list):
    """Checks whether the coverage delta of a given element is outside the
    range of deltas constituted by the rest of the list without the element being tested.

    Args:
        sorteddelta (list): List of deltas, ideally sorted increasingly by
                            absolute distance pf the source Peak from the
                            riboswitch end.
        i (int): Subject delta index in sorteddelta.

    Returns:
        tuple: (bool, list) -> boolean is True if the subject coverage delta
               is negative and the minumum delta in sorteddelta.
    """
    # Find coverage drops across peaks without cand
    cov_nocand = []
    cov_nocand.extend(sorteddelta[0:i])
    cov_nocand.extend(sorteddelta[i + 1 :])

    return (sorteddelta[i] <= min(cov_nocand) and sorteddelta[i] < 0, cov_nocand)


def has_candidate_peak(
    alignTup: tuple,
    kernel_size: int = 51,
    kernel_stdev: int = 5,
    left_margin=1.0,
    right_margin=1.0,
):
    cov, rawends, (switch_left, switch_right) = alignTup
    readcov, infercov, clipcov = cov

    # np.add can ONLY add two arrays.
    # The third param is the output array obj
    sumcov = np.add(readcov, infercov)
    sumcov = np.add(sumcov, clipcov)

    ends = convolve_array(rawends, kernel_size, kernel_stdev)

    # - strand riboswitches are reverse complemented during the reference generation,
    # so the right end in the reference is still the 3' end
    switch_end = switch_right
    switch_size = switch_right - switch_left

    # Set relevant regions to compare candidates within
    relevant_l = switch_left - switch_size * left_margin
    relevant_r = switch_right + switch_size * right_margin

    peaks = find_peaks(ends, switch_end, relevant_l, relevant_r)

    # Sort peaks in order of increasing absolute distance from riboswitch end
    close_peaks = sorted(peaks, key=abs(operator.attrgetter("from_switch_end")))

    # Record how coverage is changed by identified peaks in the region of interest
    cov_delta = coverage_delta_per_peak(close_peaks, sumcov)

    for i, cand in enumerate(close_peaks):
        is_sig, nocand_cov_delta = peak_out_of_cov_delta(cov_delta, i)
        if is_sig:
            cand_obj = Candidate(cand, switch_size, nocand_cov_delta, sumcov)
            return (True, cand_obj)

        else:
            continue

    return (False, None)
