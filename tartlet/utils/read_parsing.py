import pysam
import numpy as np

np.random.randint

from collections import defaultdict


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
                    len(frag_array) if rclip_end > len(frag_array) else rclip_end
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
                    len(frag_array) if rclip_end > len(frag_array) else rclip_end
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
                fProcessed = process_read_ends(fread, "F", frag_ends, allowSoftClips)
                rProcessed = process_read_ends(revread, "R", frag_ends, allowSoftClips)

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
