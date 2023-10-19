import pysam
import numpy as np

from collections import defaultdict


class SortedBAM:
    """Wrapper around a sorted BAM file to abstract away feature extraction"""

    def __init__(self, bam_path: str) -> None:
        self.bam: pysam.AlignmentFile = pysam.AlignmentFile(bam_path, "rb")
        self._set_reference_list()

    def _set_reference_list(self):
        """Make a list of references in this sorted BAM; only considering references with aligned reads"""
        self.ref_list = [
            idxstats.contig
            for idxstats in self.bam.get_index_statistics()
            if idxstats.total > 0
        ]


class Read:
    def __init__(self, read: pysam.AlignedSegment) -> None:
        self._read = read

        self.cigarstring = read.cigarstring
        self.cigarlist = self._split_cigar()
        self.cigartally = self._tally_cigar()

        self.orientation = "R" if read.is_reverse else "F"

        self.ref_start: int = read.reference_start
        if read.reference_end is not None:
            self.ref_end: int = read.reference_end

    def _split_cigar(self):
        """Turns the read CIGAR string into a list by separating between numbers and CIGAR opcodes.

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
        for i in range(len(self.cigarstring)):
            val = self.cigarstring[top:i]

            try:
                letter = options[self.cigarstring[i]]
                cig_list.append(int(val))
                cig_list.append(letter)
                top = i + 1

            except KeyError:
                pass  # Haven't hit a letter yet

        return cig_list

    def _tally_cigar(self):
        """Returns a dict with a tally of CIGAR markers.

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
        for i in range(len(self.cigarstring)):
            val = self.cigarstring[top:i]

            try:
                tally_dict[self.cigarstring[i]] = tally_dict[self.cigarstring[i]] + int(
                    val
                )
                top = i + 1

            except KeyError:
                pass  # Haven't hit a letter yet

        return tally_dict


class ReadPair:
    def __init__(self, *reads: Read) -> None:
        self.num_in_pair = len(reads)
        self.reads = reads

        self._compute_pair_orientation()

    def _compute_pair_orientation(self):
        """Computes read pair orientation based on the ALIGNED segments
        without considering soft clipped ends. Assumes sequencing context is paired-end.

        In the case of a fully mapped pair, "TANDEM" denotes the errant orientation where paired reads aligned
        in the same direction.

        "EQ" denotes reads aligned in opposite directions but overlapped across
        terminal base alignments.

        "FR" denotes that the left-most aligned
        base belongs to the forward aligned read and "RF" denotes that the
        left-most aligned base belongs to the reverse aligned read.

        "FFR" is an errant orientation denoting that the reverse read is fully
        aligned within the aligned termini of the forward read. Similarly,
        "RRF" is the errant orientation denoting that the forward read is fully
        aligned within the aligned termini of the reverse read.

        In the case of single reads with unmapped mates, determines whether it is oriented
        in the forward (F) or reverse (R) direction relative to the reference sequence.
        """

        if self.num_in_pair == 2:  # Is the full pair mapped?
            # Check if the reads are tandem
            if self.reads[0].orientation == self.reads[1].orientation:
                self.f, self.r = None, None
                self.orientation = "TANDEM"
                return

            # Do the reads in this pair fully overlap at the terminal aligned?
            if (
                self.reads[0].ref_start == self.reads[1].ref_start
                and self.reads[0].ref_end == self.reads[1].ref_end
            ):
                self.f, self.r = self.reads
                self.orientation = "EQ"
                return

            else:
                self.f = (
                    self.reads[0] if self.reads[0].orientation is "F" else self.reads[1]
                )
                self.r = (
                    self.reads[0] if self.reads[0].orientation is "R" else self.reads[1]
                )

                self.orientation = (
                    "FR" if self.f.ref_start <= self.r.ref_start else "RF"
                )

                # Check the case where r is exclusively within f alignment
                if self.orientation is "FR" and self.r.ref_end < self.f.ref_end:
                    self.orientation = "FFR"

                elif self.orientation is "RF" and self.f.ref_end <= self.r.ref_end:
                    self.orientation = "RRF"

        elif self.num_in_pair == 1:
            self.orientation = self.reads[0].orientation

            self.f = self.reads[0] if self.orientation is "F" else None
            self.r = self.reads[0] if self.orientation is "R" else None

        else:
            raise ValueError("More than two reads cannot be processed as a pair.")


class AlignDat:
    """
    Organises all the necessary alignment data required to infer riboswitch
    transcriptional activity for one reference in one sorted BAM.
    """

    def __init__(
        self,
        readcov: np.ndarray,
        infercov: np.ndarray,
        clipcov: np.ndarray,
        rawends: np.ndarray,
    ) -> None:
        self.readcov = readcov
        self.infercov = infercov
        self.clipcov = clipcov
        self.rawends = rawends

    @classmethod
    def base_init(cls, ref_length: int) -> "AlignDat":
        """Initialise an empty AlignDat object from reference length.

        Args:
            ref_length (int): Size of the reference in base pairs.

        Returns:
            AlignDat: Object with empty coverage and ends arrays initialised.
        """
        return cls(
            np.zeros(ref_length),
            np.zeros(ref_length),
            np.zeros(ref_length),
            np.zeros(ref_length),
        )


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


def process_read_ends(
    read: pysam.libcalignedsegment.AlignedSegment,
    orientation: str,
    ends_array: list,
    allowSoftClips=True,
):
    """Adds the tail of a read (the start of read synthesis) to an array
    that counts the number of read termini at each position in the reference.
    If the terminus is likely a fragment start, the count at its position is decremented.
    If the temrinus is likely a fragment end, the count at its position is incremented.

    If the read is in the "F" orientation, the left end of the read is
    counted as a likely fragment start.

    In the case of an "R" read, the right end of the read is counted as a likely fragment
    start.

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
    """Increments elements in a fragment coverage array that correspond to regions
        that are likely spanned by a given pair of reads. Also handles the presence of
        singular mapped reads.

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
        # Only proceed if there are aligned reads to this reference
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

        outDat = AlignDat.base_init(ref_length)
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
