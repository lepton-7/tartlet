import operator
import numpy as np
import numpy.typing as npt

from typing import Optional
from tart.utils.read_parsing import AlignDat
from scipy.stats import norm, ks_2samp, multivariate_normal


class Peak:
    """Handles peaks in convolved data."""

    def __init__(
        self,
        summit: int,
        height: float,
        from_switch_end: Optional[int] = None,
        left: Optional[int] = None,
        right: Optional[int] = None,
        switch_end: Optional[int] = None,
    ):
        self.summit = summit
        self.height = height

        if from_switch_end is not None:
            self.from_switch_end = from_switch_end
        elif switch_end is not None:
            self.switch_end = switch_end
            self.from_switch_end = summit - switch_end

        if from_switch_end is not None or switch_end is not None:
            # Only really used when sorting in one downstream step
            self.abs_from_switch_end = abs(self.from_switch_end)

        if left is not None and right is not None:
            self.left = left
            self.right = right
            self.width = right - left

    def find_bounds(self, source_arr):
        """Find the region spanned by the peak in the source array.

        Args:
            source_arr (list or array): The source array containing the peak.

        Raises:
            AttributeError: If the Peak object is not instantiated with its summit.
            AttributeError: If the Peak object is not instantiated with its height.

        Returns:
            bool: Whether the bounds were successfully determined.
        """
        if self.summit is None:
            raise AttributeError("Peak center needs to be defined")

        if self.height is None:
            raise AttributeError("Peak height needs to be defined")

        il = self.summit
        ir = self.summit
        checkLeft = True
        checkRight = True

        # Descend down the sides of the peak. The bounds include the
        # monotonic region without crossing zero. This way, shoulders are
        # split across neighbouring Peaks
        while checkLeft:
            il -= 1
            # Outside array bounds
            if il < 0:
                break
            # Monotonic
            checkLeft = abs(source_arr[il]) < abs(source_arr[il + 1])
            # Zero-crossing
            checkLeft = checkLeft and source_arr[il] * source_arr[il + 1] > 0
            # Value essentially 0
            checkLeft = checkLeft and abs(source_arr[il]) > 0.1

        # Since left is the inclusive bound,
        # account for the descent going too far by 1
        il += 1

        while checkRight:
            ir += 1
            # Outside array bounds
            if il >= len(source_arr):
                break
            # Monotonic
            checkRight = abs(source_arr[ir]) < abs(source_arr[ir - 1])
            # Zero-crossing
            checkRight = checkRight and source_arr[ir] * source_arr[ir - 1] > 0
            # Value essentially 0
            checkRight = checkRight and abs(source_arr[ir]) > 0.1

        # Only continue if the peak has valid bounds
        if il < self.summit and ir > self.summit:
            self.left = il
            self.right = ir
            self.width = self.right - self.left

            return True

        else:
            return False

    def __compare(self, peak: "Peak", heightMarg: float = 0.1, widthMarg: int = 4):
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
        wdiff = abs(self.width - peak.width)

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

        return (ks_2samp(this_ends, peak_ends), this_ends, peak_ends, peak.summit)

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
            self.abs_ends = np.absolute(source_arr[self.left : self.right])
            return self.abs_ends

    def find_absolute_coverage_delta(self, sumcov: npt.NDArray[np.float64]) -> int:
        """Computes the absolute change in the summed coverage between the start and end of the peak (center ± half-width).

        Returns:
            int: Absolute change in summed coverage across the peak.
        """

        # Bounds are already validated when determined
        return sumcov[self.right] - sumcov[self.left]

    def find_relative_coverage_delta(self, sumcov: npt.NDArray[np.float64]) -> float:
        """Computes the change in summed coverage across the peak (center ± half-width) relative to the summed coverage at the lower bound (center - half-width) of the peak.

        Args:
            sumcov (list): Coverage summed across actual, inferred, and clipped reads per base.

        Returns:
            float: Change in summed coverage across the peak relative to coverage at the start of the peak.
        """

        return self.find_absolute_coverage_delta(sumcov) / sumcov[self.left]

    def set_frags_ending_at_peak(self, fragments: list[list[tuple]]):
        """Set attributes for the distribution of start positions and end positions of fragments that end in this Peak.

        Args:
            fragments (list[list[tuple]]): Fragment tracking list from an AlignDat object.
                Index of the outer list represents positions. Inner list collects fragments
                that have an end location at that index.
        """
        # self.ending_frags: list[tuple] = []
        starts: list[int] = []
        ends: list[int] = []

        # Keep track of all the start and end positions
        for col_list in fragments[self.left : self.right]:
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

    def end_symmetry_stat(self):
        """Do a 2-sample ks test between the fragment start and end positions
        to check for symmetry."""

        ksRes = ks_2samp(self.fragment_starts, self.fragment_ends)
        self.symks_stat = ksRes.statistic  # type: ignore
        self.symks_pval = ksRes.pvalue  # type: ignore


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
        cov_drop_pval: float,
        alignDat: AlignDat,
    ) -> None:
        super().__init__(
            from_switch_end=cand.from_switch_end,
            summit=cand.summit,
            height=cand.height,
            left=cand.left,
            right=cand.right,
        )

        self.abs_cov_delta = self.find_absolute_coverage_delta(alignDat.summedcov)
        self.rel_cov_delta = self.find_relative_coverage_delta(alignDat.summedcov)

        self.switch_size = switch_size
        self.from_switch_end_relative = self.from_switch_end / self.switch_size

        self.coverage_delta_noise = nocand_cov_delta
        self.coverage_drop_pvalue = cov_drop_pval

        self.note: str = ""

        # Find start position and end position distributions for fragments
        # that ended in this candidate.
        self.set_frags_ending_at_peak(alignDat.fragments)

        self.end_symmetry_stat()


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
    ret: list[Peak] = []

    u = len(arr) if u < 0 else u

    for i in range(l + 1, u - 1):
        if abs(arr[i]) > abs(arr[i - 1]) and abs(arr[i]) > abs(arr[i + 1]):
            peak = Peak(summit=i, height=arr[i], switch_end=switch_end)
            if peak.find_bounds(arr):
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


def coverage_delta_per_peak(peaks: list, sumcov: npt.NDArray[np.float64]):
    """Helper function to find coverage deltas across a list of Peaks.

    Args:
        peaks (list): List of Peak objects.
        sumcov (list): Summed coverage list.

    Returns:
        list: Absolute coverage drop for each peak in the same order as peaks.
    """
    raw_cov_drop: list[list[int]] = []

    for peak in peaks:
        peak: Peak
        raw_cov_drop.append([peak.find_absolute_coverage_delta(sumcov), peak.width])

    return raw_cov_drop


def peak_out_of_cov_delta(sorteddelta: list[list[int]], i: int):
    """DEPRECATED

    Checks whether the coverage delta of a given element is outside the
    range of deltas constituted by the rest of the list without the element being tested.

    Args:
        sorteddelta (list): List of deltas, ideally sorted increasingly by
                            absolute distance of the source Peak from the
                            riboswitch end.
        i (int): Subject delta index in sorteddelta.

    Returns:
        tuple: (bool, list) -> boolean is True if the subject coverage delta
               is negative and the minumum delta in sorteddelta.
    """
    # Find coverage drops across peaks without cand
    cov_nocand = []
    just_deltas = []
    cov_nocand.extend(sorteddelta[0:i])
    cov_nocand.extend(sorteddelta[i + 1 :])

    return (
        sorteddelta[i][0] <= min([x[0] for x in cov_nocand]) and sorteddelta[i][0] < 0,
        cov_nocand,
    )


def peak_significance(sorteddelta: list[list[int]], i: int) -> tuple[float, list]:
    """Sets up a multivariate Gaussian using the coverage drop and peak width,
    and calculate the p-value of sampling the peak from the given distribution.

    Args:
        sorteddelta (list[list[int]]): List of deltas, ideally sorted increasingly by
                            absolute distance of the source Peak from the
                            riboswitch end.
        i (int): Subject delta index in sorteddelta.

    Returns:
        _type_: _description_
    """
    cov_nocand = []
    cov_nocand.extend(sorteddelta[0:i])
    cov_nocand.extend(sorteddelta[i + 1 :])

    m = np.mean(cov_nocand, axis=0)
    sd = np.cov(cov_nocand, rowvar=0)  # type: ignore

    try:
        pval = multivariate_normal(mean=m, cov=sd, allow_singular=True).cdf(
            sorteddelta[i]
        )
    except ValueError:
        # some weird scipy error that occasionaly happens
        # check 7.all.filter_BAM_plots.out.2024-01-08_00-17-30.25487369 for details
        pval = 1.0

    return (
        pval,
        cov_nocand,
    )


def has_candidate_peak(
    alignDat: AlignDat,
    kernel_size: int = 51,
    kernel_stdev: float = 1.5,
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

    toRet: list[Candidate] = []
    for i, cand in enumerate(close_peaks):
        try:
            pval, nocand_cov_delta = peak_significance(cov_delta, i)
        except ValueError:  # Haven't yet figured out why multivariate fails
            pval, nocand_cov_delta = 1, [None]
            print()
        if pval <= 0.05:
            cand_obj = Candidate(cand, switch_size, nocand_cov_delta, pval, alignDat)
            toRet.append(cand_obj)

        else:
            continue

    return toRet
