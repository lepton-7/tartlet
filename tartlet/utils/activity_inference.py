import numpy as np
import operator

from math import ceil, isclose
from scipy.stats import norm, ks_2samp
from tart.utils.read_parsing import AlignDat


class Peak:
    """Handles peaks in convolved data."""

    def __init__(
        self,
        from_switch_end=None,
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

        # Only really used when sorting in one downstream step
        self.abs_from_switch_end = abs(self.from_switch_end)

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

    def set_frags_ending_at_peak(self, fragments: list[list[tuple]]):
        """Set attributes for the distribution of start positions and end positions of fragments that end in this Peak.

        Args:
            fragments (list[list[tuple]]): Fragment tracking list from an AlignDat object.
                Index of the outer list represents positions. Inner list collects fragments
                that have an end location at that index.
        """
        # self.ending_frags: list[tuple] = []
        starts = []
        ends = []
        left = self.center - int(self.half_width)
        right = self.center + int(self.half_width)

        # Keep track of all the start and end positions
        for col_list in fragments[left:right]:
            for tup in col_list:
                if tup[0] is not None:  # None-types break finding index offset
                    starts.append(tup[0])
                if tup[1] is not None:
                    ends.append(tup[1])

        # Find what the index offset is to form the base data of the start
        # and end histograms independent of their index in reference
        if len(starts) > 0 and len(ends) > 0:
            s_off = min(starts)
            e_off = min(ends)

            # Histogram array size
            s_size = max(starts) - s_off + 1
            e_size = max(ends) - e_off + 1

            self.fragment_starts = np.zeros(s_size)
            self.fragment_ends = np.zeros(e_size)

            for i in starts:
                self.fragment_starts[i - s_off] += 1

            for i in ends:
                self.fragment_ends[i - e_off] += 1

        else:
            self.fragment_starts = [0]
            self.fragment_ends = [0]


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
        alignDat: AlignDat,
    ) -> None:
        super().__init__(
            from_switch_end=cand.from_switch_end,
            center=cand.center,
            height=cand.height,
            half_width=cand.half_width,
        )

        self.abs_cov_delta = cand.find_absolute_coverage_delta(alignDat.summedcov)
        self.rel_cov_delta = cand.find_relative_coverage_delta(alignDat.summedcov)

        self.switch_size = switch_size
        self.from_switch_end_relative = self.from_switch_end / self.switch_size

        self.coverage_delta_noise = nocand_cov_delta

        # Find start position and end position distributions for fragments
        # that ended in this candidate.
        self.set_frags_ending_at_peak(alignDat.fragments)


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


def _gen_kernel(kernel_size: int = 21, std_dev: float = 3.0):
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
    summedcov: list, switchreg: tuple, cand: Peak, stdev: int, minDrop: float
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

    sstart, send = switchreg

    interv = 2 * stdev
    istart, iend = cand.center - interv, cand.center + interv

    istart = 0 if istart < 0 else istart
    iend = len(summedcov) - 1 if iend >= len(summedcov) else iend

    diff = summedcov[istart] - summedcov[iend]
    maxreads = max(summedcov[sstart:send])

    return diff / maxreads > minDrop


def is_interesting(
    alignDat: AlignDat,
    windowfrac: float = 0.15,
    threshtol: float = 0.15,
):
    """
    DEPRECATED

    Determines whether the read alignment is indicative of transcriptionally active riboswitches.

    Args:
        alignDat (tuple): Alignment tuple for this reference:
            ([read coverage, inferred frag cov, clipped cov], ends, (switch start, switch end))
        windowfrac (float, optional): Window size on either side of riboswitch 3' as a fraction
            of the total riboswitch size to check for fragment end peaks. Defaults to 15%.
        threshtol (float, optional): Allowable fractional percentage margin when checking
            if a peak height is close to the maximum peak height in the full
            riboswitch + window region. Defaults to 0.15.

    Returns:
        bool: True if alignTuple is determined to be interesting. False if not.
    """

    # OPTIONS -----------------------------------------------------------------
    kernel_size = 51
    kernel_stdev = 5

    minReadDrop = 0.3

    peakHeightTol = 0.15
    peakWidthMarg = 4
    # -------------------------------------------------------------------------

    kernel = _gen_kernel(kernel_size=kernel_size, std_dev=kernel_stdev)
    alignDat.convolve_rawends(kernel)

    # - strand riboswitches are reverse complemented during the reference generation,
    # so the right end in the reference is still the 3' end
    peaks = find_peaks(alignDat.convends)

    switch_left = alignDat.switch_start
    switch_right = alignDat.switch_end
    switch_end = alignDat.switch_end

    window = (int)((switch_right - switch_left) * windowfrac)
    left, right = switch_end - window, switch_end + window

    # Max peak from riboswitch 5' to the 3' window end
    maxpeak = max(alignDat.convends[switch_left:right])

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
                alignDat.summedcov,
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
                cand_stats.append(cand.compare_ks(peak, alignDat.convends))
                if cand.compare(
                    peak, heightMarg=peakHeightTol, widthMarg=peakWidthMarg
                ):
                    keepcand = False
                    break

            # This candidate does not have a mirrored peak
            if keepcand:
                return (True, cand_stats)

    return (False, cand_stats)


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


def peak_out_of_cov_delta(sorteddelta: list, i: int):
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
    alignDat: AlignDat,
    kernel_size: int = 51,
    kernel_stdev: int = 5,
    left_margin=1.0,
    right_margin=1.0,
):
    kernel = _gen_kernel(kernel_size, kernel_stdev)
    alignDat.convolve_rawends(kernel)

    # - strand riboswitches are reverse complemented during the reference
    # generation, so the right end in the reference is still the 3' end
    switch_end = alignDat.switch_end
    switch_size = alignDat.switch_end - alignDat.switch_start

    # Set relevant regions to compare candidates within ± switch size
    # proportion outside the switch
    relevant_l = alignDat.switch_start - int(switch_size * left_margin)
    relevant_r = alignDat.switch_end + int(switch_size * right_margin)

    peaks = find_peaks(alignDat.convends, switch_end, relevant_l, relevant_r)

    # Sort peaks in order of increasing absolute distance from riboswitch end
    close_peaks = sorted(peaks, key=operator.attrgetter("abs_from_switch_end"))

    # Record how coverage is changed by identified peaks in the region of interest
    cov_delta = coverage_delta_per_peak(close_peaks, alignDat.summedcov)

    for i, cand in enumerate(close_peaks):
        is_sig, nocand_cov_delta = peak_out_of_cov_delta(cov_delta, i)
        if is_sig:
            cand_obj = Candidate(cand, switch_size, nocand_cov_delta, alignDat)
            return cand_obj

        else:
            continue

    return None
