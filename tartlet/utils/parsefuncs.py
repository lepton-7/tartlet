import pysam
import numpy as np

import matplotlib.pyplot as plt

from scipy.stats import norm
from collections import defaultdict
from math import ceil, isclose


class Peak:
    def __init__(self, center=None, height=None, half_width=None):
        self.center = center
        self.height = height
        self.half_width = half_width

    def find_half_width(self, source_arr):
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
        hprop = peak.height / self.height
        wdiff = abs(self.half_width - peak.half_width)

        return (
            hprop > -1 - heightMarg and hprop < -1 + heightMarg and wdiff <= widthMarg
        )


def tally_cigar(cigar: str):
    """Returns a dict with a tally for CIGAR markers from a string

    Args:
        cigar (str): CIGAR string for a read

    Returns:
        dict: Dictionary keyed by CIGAR opcodes storing tallies for the read
    """

    tally_dict = {
        "M": 0,
        "I": 0,
        "D": 0,
        "N": 0,
        "S": 0,
        "H": 0,
        "P": 0,
        "=": 0,
        "X": 0,
    }

    top = 0
    for i in range(len(cigar)):
        val = cigar[top:i]

        try:
            tally_dict[cigar[i]] = tally_dict[cigar[i]] + int(val)
            top = i + 1

        except KeyError:
            pass  # Haven't hit a letter yet

    return tally_dict


def split_cigar(cigar: str):
    """Turns a CIGAR string into a list by separating between numbers and CIGAR opcodes.

    Args:
        cigar (str): CIGAR string for a read

    Returns:
        list: list representation of the CIGAR string
    """
    cig_list = []

    options = {
        "M": "M",
        "I": "I",
        "D": "D",
        "N": "N",
        "S": "S",
        "H": "H",
        "P": "P",
        "=": "=",
        "X": "X",
    }

    top = 0
    for i in range(len(cigar)):
        val = cigar[top:i]

        try:
            letter = options[cigar[i]]
            cig_list.append(int(val))
            cig_list.append(letter)
            top = i + 1

        except KeyError:
            pass  # Haven't hit a letter yet

    return cig_list


def make_pair_dict(bam: pysam.AlignmentFile, refname: str):
    """Fetches reads for a contig (refname) and groups them into pairs based read name

    Args:
        bam (pysam.AlignmentFile): Sorted BAM file
        refname (str): Reference sequence name to look for in the sorted BAM

    Returns:
        collections.defaultdict(list): A dictionary keyed by read name for read pair lists
    """
    pair_dict = defaultdict(list)
    for read in bam.fetch(refname):
        pair_dict[read.query_name].append(read)

    return pair_dict


def paired_reads_orientation(pair: list):
    """Computes the orientation of a given read pair based on the ALIGNED segments
    without considering soft clipped ends.

    In the case of single reads, determines whether it is oriented in the forward (F)
    or reverse (R) direction relative to the reference sequence.

    In the case of a pair of reads, "FR" denotes that the left-most read is forward
    aligned and "RF" denotes that the left-most read is reverse aligned.

    "EQ" is returned if the reads align in opposite directions but the start and ends
    overlap with each other.

    "TANDEM" denotes a pair that is aligned in the same direction.

    Args:
        pair (list): A list of reads in a pair. List should have a length of 2.

    Returns:
        tuple: (no. reads in pair, orientation, f_read, r_read). f_read or r_read is None if that
        orientation is not present in the pair. Both are None if orientation is "TANDEM"
    """
    num = len(pair)
    orientation = ""
    fread = None
    revread = None

    if len(pair) > 1:  # Is the full pair mapped?
        # Check if the reads are tandem
        if pair[0].is_reverse == pair[1].is_reverse:
            return (num, "TANDEM", None, None)

        # Do the reads in this pair fully overlap?
        if (
            pair[0].reference_start == pair[1].reference_start
            and pair[0].reference_end == pair[1].reference_end
        ):
            orientation = "EQ"
            fread = pair[0]
            revread = pair[1]

        else:
            fread = pair[0] if not pair[0].is_reverse else pair[1]
            revread = pair[0] if pair[0].is_reverse else pair[1]

            orientation = (
                "FR" if fread.reference_start <= revread.reference_start else "RF"
            )

    else:  # just one read in this pair
        if pair[0].is_reverse:
            return (num, "R", None, pair[0])

        if not pair[0].is_reverse:
            return (num, "F", pair[0], None)

    return (num, orientation, fread, revread)


def process_read_ends(
    read: pysam.libcalignedsegment.AlignedSegment,
    orientation: str,
    ends_array: list,
    allowSoftClips=True,
):
    """Adds the tail of a read (the start of read synthesis) to an array
    that counts the number of read termini at each position in the reference.
    If the terminus is a likely fragment start, the count at its position is decremented.
    If the terminys is a liekly fragment end, the cout at its position is incremented.

    If the read is in the "F" orientation, the left end of the read is
    counted as a likely fragment end.

    In the case of an "R" read, the right end of the read is counted.

    Args:
        read (pysam.libcalignedsegment.AlignedSegment): An aligned read.
        orientation (str): Orientation of the read relative to the reference. Either "F" or "R".
        ends_array (list): The list in which likely fragment ends are tallied.
        allowSoftClips (bool, optional): If true, considers the start/end of the soft clipped region when
            determining where the likely fragment end is for this read. If set to False, only considers the
            aligned portion of the read regardless of whether the soft clipped region extends beyond.
            Defaults to True.

    Returns:
        bool: True if the read end was accounted for in the fragment ends list. False otherwise.
    """
    cig_arr = split_cigar(read.cigarstring)

    # Do nothing with this read if the soft clip is at the critical
    # end of the read and soft clipped reads are not considered
    if orientation == "F" and not allowSoftClips and cig_arr[1] == "S":
        return False

    if orientation == "R" and not allowSoftClips and cig_arr[-1] == "S":
        return False

    bounds_list = read.get_blocks()

    if orientation == "F":
        if allowSoftClips and cig_arr[1] == "S":
            # First calculate how far before the
            # aligned portion the clip extends
            try:
                clip_len = cig_arr[0]

                # start of the aligned section in the ref.
                # Clips are preceeded/succeeded by M regions so this
                # should not do arithmetic on None
                ref_clp_start = bounds_list[0][0] - clip_len

                ends_array[ref_clp_start] -= 1

            except IndexError:
                # The calculated clip start would be before the ref start
                ends_array[0] -= 1

        else:
            refloc = bounds_list[0][0]
            ends_array[refloc] -= 1

    elif orientation == "R":
        if allowSoftClips and cig_arr[-1] == "S":
            # First calculate how far before the
            # aligned portion the clip extends
            try:
                clip_len = cig_arr[-2]

                # end of the aligned section in the ref.
                # Clips are preceeded/succeeded by M regions so this
                # should not do arithmetic on None

                # The -1 is because the upper bound is not inclusive
                ref_clp_start = bounds_list[-1][1] + clip_len - 1

                ends_array[ref_clp_start] += 1

            except IndexError:
                # The calculated clip end would fall beyond the ref
                ends_array[-1] += 1

        else:
            refloc = bounds_list[-1][1] - 1
            ends_array[refloc] += 1

    return True


def add_to_frags(
    readtup: tuple,
    orientation: str,
    frag_array: list,
    infer_array: list,
    clip_array: list,
    allowSoftClips=True,
):
    """Increments elememts in a fragment coverage array that correspond to regions
        that are likely spanned by a given pair of reads. Also handles the presence of
        only one mapped read in a pair.

    Args:
        readtup (tuple): A tuple of reads: (forward read, reverse read).
            In the case of only one read being mapped, based on the mapped read's orientation,
            the other element in the tuple is None.
        orientation (str): "FR" or "EQ" for mapped pairs; "F" or "R" for one mapped read.
            Refers to the orientation of the read relative to the reference.
        frag_array (list): The likely fragment coverage array whose elements are incremented
            by 1 across the likely fragment span indicated by the read/read pair.
        allowSoftClips (bool, optional): Determines whether soft clipped regions at
            relevant ends of reads should be counted as part of the fragment. Defaults to True.
    """

    # The two types of coverage
    read_regions = []
    inferred_regions = []
    clipped_regions = []

    # Handle unpaired F or R reads:
    if orientation == "F" or orientation == "R":
        read = readtup[0] if readtup[1] is None else readtup[1]

        # get_blocks() does not include soft clipped regions
        bounds_list = read.get_blocks()

        fstart = bounds_list[0][0]
        rend = bounds_list[-1][1]

        read_regions.append((fstart, rend))

        # If soft clipped regions should be considered, compute the length of
        # the clip and increment the fragment coverage array over that region as well
        cig_list = split_cigar(read.cigarstring)

        if allowSoftClips:
            if cig_list[1] == "S":  # clipped start
                clip_len = cig_list[0]
                fclip_start = bounds_list[0][0] - clip_len

                # Validate bounds
                fclip_start = 0 if fclip_start < 0 else fclip_start
                # clipped region to left
                clipped_regions.append((fclip_start, fstart))

            if cig_list[-1] == "S":  # clipped end
                clip_len = cig_list[-2]
                rclip_end = bounds_list[-1][1] + clip_len

                # Validate bounds
                rclip_end = (
                    len(frag_array) if rclip_end > len(
                        frag_array) else rclip_end
                )
                # clipped region to right
                clipped_regions.append((rend, rclip_end))

    elif orientation == "FR" or orientation == "EQ":
        fread, revread = readtup

        fbounds = fread.get_blocks()
        rbounds = revread.get_blocks()

        fstart = fbounds[0][0]
        fend = fbounds[-1][1]

        rstart = rbounds[0][0]
        rend = rbounds[-1][1]

        # The regions actually covered by the reads, but we dont wan't to double count
        # bases that have overlapping reads of the same pair
        if fend < rstart:
            read_regions.extend([(fstart, fend), (rstart, rend)])
        else:
            read_regions.extend([(fstart, fend), (fend, rend)])

        # The regions between the heads of the reads that is
        # inferred to be part of the sequenced fragment. If
        # fend is greater than rstart, the tuple will ultimately not
        # do anything, which is very convenient
        inferred_regions.append((fend, rstart))

        if allowSoftClips:
            fcigs = split_cigar(fread.cigarstring)
            rcigs = split_cigar(revread.cigarstring)

            if fcigs[1] == "S":
                clip_len = fcigs[0]
                fclip_start = fbounds[0][0] - clip_len

                # Validate bounds
                fclip_start = 0 if fclip_start < 0 else fclip_start
                # clipped region to left
                clipped_regions.append((fclip_start, fstart))

            if rcigs[-1] == "S":
                clip_len = rcigs[-2]
                rclip_end = rbounds[-1][1] + clip_len

                # Validate bounds
                rclip_end = (
                    len(frag_array) if rclip_end > len(
                        frag_array) else rclip_end
                )
                # clipped region to right
                clipped_regions.append((rend, rclip_end))

    # Add actual read coverage
    for start, end in read_regions:
        for i in range(start, end):
            frag_array[i] += 1

    # Add inferred fragment coverage
    for start, end in inferred_regions:
        for i in range(start, end):
            infer_array[i] += 1

    # Add inferred clipped bases coverage
    for start, end in clipped_regions:
        for i in range(start, end):
            clip_array[i] += 1


def generate_plot_data(
    bam: pysam.AlignmentFile, refBounds: dict, allowSoftClips: bool = True
):
    """Generates the inferred fragment coverage and fragment ends arrays along with
    the switch start and stop indices within the ref sequence for each reference
    in the input sorted bam.

    Args:
        bam (pysam.AlignmentFile): Sorted and indexed BAM to process.
        refBounds (dict): A dictionary mapping the bounds of the
            reference sequence within its source contig.
        allowSoftClips (bool): Whether to allow soft clipped reads/regions to be considered.

    Returns:
        dict: Keyed by the reference sequence name.
            Each entry is a tuple: ([read coverage, inferred coverage, clipped coverage],
            fragment ends, (switch start idx, switch stop idx)).
    """
    outDict = {}
    # Iterate over each reference in the bam
    for idxstats in bam.get_index_statistics():
        # Only continue if there are aligned reads to this reference
        if idxstats.total == 0:
            continue

        ref = idxstats.contig
        ref_length = bam.get_reference_length(ref) + 1

        bounds = refBounds[ref]
        splits = ref.split("#")

        # Determine riboswitch bounds in the reference
        if splits[-1] == "+":
            switch_start = int(splits[-3]) - bounds[0]
            switch_end = int(splits[-2]) - bounds[0]

        else:
            switch_start = int(splits[-2]) - bounds[0]
            switch_end = int(splits[-3]) - bounds[0]

        frag_readcoverage = np.zeros(ref_length)
        frag_infercoverage = np.zeros(ref_length)
        frag_clipcoverage = np.zeros(ref_length)
        frag_ends = np.zeros(ref_length)

        # A dictionary with read pairs in lists keyed to read name
        pair_dict = make_pair_dict(bam, ref)

        # Iterate over each paired read for the current reference seq
        for name, pair in pair_dict.items():
            readIsProcessed = False

            num, orientation, fread, revread = paired_reads_orientation(pair)

            # Handle different read orientations
            if orientation == "TANDEM":  # Uh oh
                continue

            elif orientation == "RF":  # Uh oh
                continue

            elif orientation == "F":
                readIsProcessed = process_read_ends(
                    fread, "F", frag_ends, allowSoftClips
                )
                if readIsProcessed:
                    add_to_frags(
                        (fread, None),
                        "F",
                        frag_readcoverage,
                        frag_infercoverage,
                        frag_clipcoverage,
                        allowSoftClips,
                    )

            elif orientation == "R":
                readIsProcessed = process_read_ends(
                    revread, "R", frag_ends, allowSoftClips
                )
                if readIsProcessed:
                    add_to_frags(
                        (None, revread),
                        "F",
                        frag_readcoverage,
                        frag_infercoverage,
                        frag_clipcoverage,
                        allowSoftClips,
                    )

            # There's no functional difference between "FR" and "EQ" orientations
            elif orientation == "FR" or orientation == "EQ":
                fProcessed = process_read_ends(
                    fread, "F", frag_ends, allowSoftClips)
                rProcessed = process_read_ends(
                    revread, "R", frag_ends, allowSoftClips)

                if fProcessed and rProcessed:
                    add_to_frags(
                        (fread, revread),
                        orientation,
                        frag_readcoverage,
                        frag_infercoverage,
                        frag_clipcoverage,
                        allowSoftClips,
                    )

        outDict[ref] = (
            [frag_readcoverage, frag_infercoverage, frag_clipcoverage],
            frag_ends,
            (switch_start, switch_end),
        )
    return outDict


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

    binned_ends[numbins - 1] = sum(ends[bin_pos[numbins - 1]:])

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
    xticks = bin_ax if bin_size >= 10 else [
        i for i in range(0, len(readcov), 10)]

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
        else xticks[buffstart // 10: ceil(buffend / 10)]
    )

    # Add the coverage panel ----------------------------------------
    coverage_counts = {"Read": readcov,
                       "Inferred": infercov, "Clipped": clipcov}
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


def find_peaks(arr):
    ret = []

    for i in range(1, len(arr) - 1):
        if abs(arr[i]) > abs(arr[i - 1]) and abs(arr[i]) > abs(arr[i + 1]):
            peak = Peak(center=i, height=arr[i])

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
    """Determines whether the read alignment is indicative of transcriptionally active riboswitches.

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
            # Check if there is a similar mirrored peak
            # anywhere in the region of interest.
            # Only the frag start peaks prior to the candidate
            # frag end peak needs to be checked.
            for peak in peaks[:i]:
                # If there is a mirrored peak within margins,
                # immediately disqualify
                if cand.compare(
                    peak, heightMarg=peakHeightTol, widthMarg=peakWidthMarg
                ):
                    keepcand = False
                    break

            # This candidate does not have a mirrored peak
            if keepcand:
                return True

    return False
