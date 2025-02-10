import json
import pysam
import numpy as np
import numpy.typing as npt

# import tart.utils.filter_functions as fil

from click import echo
from typing import Optional
from collections import defaultdict


class Read:
    """Wraps pysam.AlignedSegment to simplify access to relevant features."""

    def __init__(self, read: pysam.AlignedSegment) -> None:
        self._read = read

        # Reads in a pair will have the same query_name.
        # Not sure if this is specific to HISAT
        self.name = read.query_name

        # Template length of the aligned read pair from the BAM
        # This can be positive or negative
        self.template_length = read.template_length

        self.cigarstring = read.cigarstring
        self.cigarlist = self._split_cigar()
        self.cigartally = self._tally_cigar()

        self.orientation = "R" if read.is_reverse else "F"

        self.ref_start: int = read.reference_start
        if read.reference_end is not None:
            self.ref_end: int = read.reference_end

        # get_blocks() does not include soft clipped regions
        self.block_loci = read.get_blocks()

        # The head is the leading end of the read
        self.head = (
            self.block_loci[-1][-1]
            if self.orientation == "F"
            else self.block_loci[0][0]
        )

        # The tail is the start end of the read
        self.tail = (
            self.block_loci[-1][-1]
            if self.orientation == "R"
            else self.block_loci[0][0]
        )

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
        if self.cigarstring is None:
            return []

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
        if self.cigarstring is None:
            return tally_dict

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
        """Aligned paired reads.

        Object attributes vary between fully mapped
        pairs and pairs with only one mapped mate. ReadPair.f and ReadPair.r
        are only defined if pair is fully mapped (ReadPair.in_pair == 2).
        ReadPair.read is only defined when one mate in pair is mapped
        (ReadPair.in_num == 1).
        """

        self.in_pair = len(reads)
        self._reads = reads

        # Make sure reads have the same name
        self._check_read_names(reads)
        self.name = reads[0].name

        # Absolute value for template length
        self.template_length = abs(reads[0].template_length)

        self._compute_pair_orientation()
        self._assign_frag_termini()
        self._compute_frag_length()

    def _compute_pair_orientation(self):
        """Computes read pair orientation based on the ALIGNED segments
        without considering soft clipped ends. Assumes sequencing context is paired-end.

        In the case of a fully mapped pair, "TANDEM" denotes the errant orientation where paired reads aligned
        in the same direction.

        "EQ" denotes reads aligned in opposite directions but overlapped across
        terminal base alignments.

        "FR" denotes that the left-most aligned base belongs to the forward
        aligned read and "RF" denotes that the left-most aligned base belongs
        to the reverse aligned read.

        "FFR" is an errant orientation denoting that the reverse read is fully
        aligned within the aligned termini of the forward read. Similarly,
        "RRF" is the errant orientation denoting that the forward read is fully
        aligned within the aligned termini of the reverse read.

        In the case of single reads with unmapped mates, determines whether it is oriented
        in the forward (F) or reverse (R) direction relative to the reference sequence.
        """

        if self.in_pair == 2:  # Is the full pair mapped?
            # Check if the reads are tandem
            if self._reads[0].orientation == self._reads[1].orientation:
                self.orientation = "TANDEM"
                return

            # Non-TANDEM case can be assigned at once
            self.f = (
                self._reads[0] if self._reads[0].orientation == "F" else self._reads[1]
            )
            self.r = (
                self._reads[0] if self._reads[0].orientation == "R" else self._reads[1]
            )
            # Do the reads in this pair fully overlap at aligned termini?
            if (
                self._reads[0].ref_start == self._reads[1].ref_start
                and self._reads[0].ref_end == self._reads[1].ref_end
            ):
                self.orientation = "EQ"
                return

            else:
                self.orientation = (
                    "FR" if self.f.ref_start <= self.r.ref_start else "RF"
                )

                # Check the case where r is exclusively within f alignment
                if self.orientation == "FR" and self.r.ref_end < self.f.ref_end:
                    self.orientation = "FFR"

                # Check the case where f is exclusively within r alignment
                elif self.orientation == "RF" and self.f.ref_end <= self.r.ref_end:
                    self.orientation = "RRF"

        elif self.in_pair == 1:
            self.orientation = self._reads[0].orientation

            # Don't need to specify f or r since there's just one mapped read in pair
            self.read = self._reads[0]

        else:
            raise ValueError("More than two reads cannot be processed as a pair.")

    def _check_read_names(self, readlist: tuple[Read, ...]):
        """Check if the reads passed to the instance have the same read name.

        Args:
            readlist (list[Read]): List of reads passed as a pair.

        Raises:
            ValueError: Read names do not match.
        """
        name = readlist[0].name
        if name is None:
            raise ValueError("Read does not have a name.")

        for read in readlist:
            if read.name != name:
                raise ValueError("Reads supplied do not seem to have the same name.")

    def _assign_frag_termini(self):
        """Determines the bounds [inclusive, exclusive) of the fragment spanned.

        Only records likely valid bounds. One or both bounds may be None if orientation
        is "F" or "R", or errant."""
        if self.in_pair == 1:
            self.fragment_termini = (
                (self.read.tail, -1)
                if self.orientation == "F"
                else (-1, self.read.tail)
            )

        elif self.in_pair == 2:
            if self.orientation == "FR" or self.orientation == "EQ":
                self.fragment_termini = (self.f.tail, self.r.tail)

            elif self.orientation in ["RF", "TANDEM", "FFR", "RRF"]:
                self.fragment_termini = (-1, -1)

    def _compute_frag_length(self):
        """Determines fragment size spanned by this fragment"""
        if self.in_pair == 1:
            self.fragment_size = None

        elif self.in_pair == 2:
            if self.orientation == "FR" or self.orientation == "EQ":
                self.fragment_size = self.fragment_termini[1] - self.fragment_termini[0]

            elif self.orientation in ["RF", "TANDEM", "FFR", "RRF"]:
                self.fragment_size = None


class AlignDat:
    """
    Organises all the necessary alignment data required to infer riboswitch
    transcriptional activity for one reference from one sorted BAM.
    """

    def __init__(self, ref: str, reflength: int, boundsDict: dict[str, list]) -> None:
        self.ref = ref
        self.ref_length = reflength

        self.__coverage_types = [
            "readcov",
            "infercov",
            "overlapcov",
            "clipcov",
            "terminalcov",
        ]

        # Appease checkers
        self.readcov: npt.NDArray[np.float64]
        self.infercov: npt.NDArray[np.float64]
        self.overlapcov: npt.NDArray[np.float64]
        self.clipcov: npt.NDArray[np.float64]
        self.terminalcov: npt.NDArray[np.float64]

        self.mapped_in_switch: int = 0
        self.unmapped_in_switch: int = 0

        # This is for easy future splitting of coverage types
        for covtype in self.__coverage_types:
            self.__setattr__(covtype, np.zeros(reflength))

        # Define now; populate later
        self.summedcov = np.zeros(reflength)

        # Stores terminal positions for the likely fragment spanned by each ReadPair.
        # Tuples with the same right terminal position are collected into lists.
        self.fragments: list[list[tuple]] = [[] for _ in range(reflength)]

        self.rawends = np.zeros(reflength)

        self.ref_bounds = boundsDict["ref"]
        self.orf_bounds = boundsDict["orf"]

        self._set_refrelative_switch_bounds()

    def _set_refrelative_switch_bounds(self) -> None:
        """Starting from the riboswitch locus relative to source genome
        indexing encoded in the ref string, calculates the riboswitch
        locus indices relative to the start of the alignment reference.
        """
        splits = self.ref.split("#")

        # Validate orf bounds
        try:
            _ = int(self.orf_bounds[0])
            _ = int(self.orf_bounds[1])
        except ValueError:
            self.orf_bounds = [0, 0]

        # Switch bounds are inclusive
        if splits[-1] == "+":
            self.switch_start = int(splits[-3]) - self.ref_bounds[0]
            self.switch_end = int(splits[-2]) - self.ref_bounds[0]

            self.orf_start = int(self.orf_bounds[0] - self.ref_bounds[0])
            self.orf_end = int(self.orf_bounds[1] - self.ref_bounds[0])

        else:
            # - strands need to be processed separately because the seq
            # is revcomp'd during ref generation itself
            self.switch_start = self.ref_bounds[1] - int(splits[-3])
            self.switch_end = self.ref_bounds[1] - int(splits[-2])

            self.orf_start = int(self.ref_bounds[1] - self.orf_bounds[0])
            self.orf_end = int(self.ref_bounds[1] - self.orf_bounds[1])

    def _intersect_switch(self, readstart: int, readend: int) -> bool:
        """Check whether the mapped portion of a read intersects with the riboswitch region.

        Args:
            readstart (int): Start index for the read.
            readend (int): End index for the read.

        Returns:
            bool: True if any part of the read intersects with the riboswitch.
        """
        return (
            (readstart >= self.switch_start and readstart <= self.switch_end)
            or (readend >= self.switch_start and readend <= self.switch_end)
            or (readstart <= self.switch_start and readend >= self.switch_end)
        )

    def _coalesce_into_cov(self, pair: ReadPair, allowSoftClips):
        """Extracts fragment coverage information from a ReadPair object.

        Args:
            pair (ReadPair): Pair of reads to process coverage.
            allowSoftClips (bool, optional): Consider soft-clipped regions for coverage processing. Defaults to True.
        """
        read_regions = []
        inferred_regions = []
        overlapped_regions = []
        clipped_regions = []
        terminal_regions = []

        # Handle unpaired F or R reads first
        if pair.in_pair == 1:
            read = pair.read

            fstart = read.block_loci[0][0]
            rend = read.block_loci[-1][1]

            if self._intersect_switch(fstart, rend):
                self.unmapped_in_switch += 1

            if pair.orientation == "F":
                read_regions.append((fstart + 1, rend))
                terminal_regions.append((fstart, fstart + 1))

            elif pair.orientation == "R":
                read_regions.append((fstart, rend - 1))
                terminal_regions.append((rend - 1, rend))

            # If soft clipped regions should be considered, compute the length of
            # the clip and increment the fragment coverage array over that region as well
            cig_list = read.cigarlist

            if allowSoftClips:
                if cig_list[1] == "S":  # clipped start
                    clip_len = cig_list[0]
                    fclip_start = read.block_loci[0][0] - clip_len

                    # Validate bounds
                    fclip_start = 0 if fclip_start < 0 else fclip_start
                    # clipped region to left
                    clipped_regions.append((fclip_start, fstart))

                if cig_list[-1] == "S":  # clipped end
                    clip_len = cig_list[-2]
                    rclip_end = read.block_loci[-1][1] + clip_len

                    # Validate bounds
                    rclip_end = (
                        self.ref_length if rclip_end > self.ref_length else rclip_end
                    )
                    # clipped region to right
                    clipped_regions.append((rend, rclip_end))

        elif pair.orientation == "FR" or pair.orientation == "EQ":
            fread = pair.f
            revread = pair.r

            fstart = fread.block_loci[0][0]
            fend = fread.block_loci[-1][1]

            rstart = revread.block_loci[0][0]
            rend = revread.block_loci[-1][1]

            if self._intersect_switch(fstart, fend) or self._intersect_switch(
                rstart, rend
            ):
                self.mapped_in_switch += 1

            # The regions actually covered by the reads, but we dont want to double count
            # bases that have overlapping reads of the same pair. These are processed separately.
            if fend < rstart:
                read_regions.extend([(fstart + 1, fend), (rstart, rend - 1)])
            else:  # reads in the pair overlap
                read_regions.extend([(fstart + 1, rstart), (fend, rend - 1)])

            # The regions between the heads of the reads that is
            # inferred to be part of the sequenced fragment. If
            # fend is greater than rstart, the tuple will ultimately not
            # do anything, which is very convenient
            inferred_regions.append((fend, rstart))

            # Overlap region between reads in a pair
            adj_rstart = rstart + 1 if rstart == fstart else rstart
            adj_fend = fend - 1 if fend == rend else fend

            overlapped_regions.append((adj_rstart, adj_fend))

            # terminal bases
            terminal_regions.extend([(fstart, fstart + 1), (rend - 1, rend)])

            if allowSoftClips:
                fcigs = fread.cigarlist
                rcigs = revread.cigarlist

                if fcigs[1] == "S":
                    clip_len = fcigs[0]
                    fclip_start = fread.block_loci[0][0] - clip_len

                    # Validate bounds
                    fclip_start = 0 if fclip_start < 0 else fclip_start
                    # clipped region to left
                    clipped_regions.append((fclip_start, fstart))

                if rcigs[-1] == "S":
                    clip_len = rcigs[-2]
                    rclip_end = revread.block_loci[-1][1] + clip_len

                    # Validate bounds
                    rclip_end = (
                        self.ref_length if rclip_end > self.ref_length else rclip_end
                    )

                    # clipped region to right
                    clipped_regions.append((rend, rclip_end))

        # TODO: Do this using np slicing to speed it up
        # Increment coverage arrays
        for start, end in read_regions:
            for i in range(start, end):
                self.readcov[i] += 1

        for start, end in inferred_regions:
            for i in range(start, end):
                self.infercov[i] += 1

        for start, end in overlapped_regions:
            for i in range(start, end):
                self.overlapcov[i] += 1

        for start, end in clipped_regions:
            for i in range(start, end):
                self.clipcov[i] += 1

        for start, end in terminal_regions:
            for i in range(start, end):
                self.terminalcov[i] += 1

    def _process_ends(self, read: Read, allowSoftClips):
        """Adds the tail of a read (the start of read synthesis) to an array
        that counts the number of read termini at each position in the reference.
        If the terminus is likely a fragment start, the count at its position is decremented.
        If the temrinus is likely a fragment end, the count at its position is incremented.

        If the read is in the "F" orientation, the left end of the read is
        counted as a likely fragment start.

        In the case of an "R" read, the right end of the read is counted as a likely fragment
        start.

        Args:
            read (Read): Aligned read.
            allowSoftClips (bool): If true, reads that have any soft-clipping
            are processed. Soft-clipped reads are thrown out if set to False.

        Returns:
            bool: True if the read end was accounted for in the fragment ends list. False otherwise.
        """

        # Throw out read if soft clips are not enabled
        if not allowSoftClips and read.cigartally["S"] > 0:
            return False

        # Do nothing with this read if the soft clip is at the critical
        # end of the read and soft clipped reads are not considered
        if read.orientation == "F" and not allowSoftClips and read.cigarlist[1] == "S":
            return False

        if read.orientation == "R" and not allowSoftClips and read.cigarlist[-1] == "S":
            return False

        if read.orientation == "F":
            if allowSoftClips and read.cigarlist[1] == "S":
                # First calculate how far the clip extends
                # before the aligned region
                try:
                    clip_len = read.cigarlist[0]

                    # start of the aligned section in the ref.
                    # Clips are preceeded/succeeded by M regions so this
                    # should not do arithmetic on None
                    ref_clp_start = read.block_loci[0][0] - clip_len

                    self.rawends[ref_clp_start] -= 1

                except IndexError:
                    # The calculated clip start would be before the ref start
                    self.rawends[0] -= 1

            else:
                refloc = read.block_loci[0][0]
                self.rawends[refloc] -= 1

        elif read.orientation == "R":
            if allowSoftClips and read.cigarlist[-1] == "S":
                # First calculate how far the clip extends
                # before the aligned region
                try:
                    clip_len = read.cigarlist[-2]

                    # end of the aligned section in the ref.
                    # Clips are preceeded/succeeded by M regions so this
                    # should not do arithmetic on None

                    # The -1 is because the upper bound is not inclusive
                    ref_clp_start = read.block_loci[-1][1] + clip_len - 1

                    self.rawends[ref_clp_start] += 1

                except IndexError:
                    # The calculated clip end would fall beyond the ref
                    self.rawends[-1] += 1

            else:
                refloc = read.block_loci[-1][1] - 1
                self.rawends[refloc] += 1

        return True

    def _process_frag_locus(self, pair: ReadPair):
        """Store the terminal positions in self.fragment indexed by the right end of pair fragment.

        Args:
            pair (ReadPair): Pair spanning the likely fragment.
        """
        l, r = pair.fragment_termini

        if r is not None:
            self.fragments[r].append(pair.fragment_termini)

    def process_pairs(
        self, pairs: list[ReadPair], allowSoftClips: bool, allowSingleReads: bool = True
    ) -> "AlignDat":
        """Process a list of ReadPair objects to extract alignment coverage and ends data.

        Args:
            pairs (list[ReadPair]): List of paired reads.
            allowSoftClips (bool): Consider soft clipped regions.
            allowSingleReads (bool, default True): Consider unmapped reads.

        Returns:
            AlignDat: Same object after incorporating all ReadPair objects in pairs.
        """
        self.allowSoftClips = allowSoftClips
        self.allowSingleReads = allowSingleReads

        # Filtered read pairs that span across any section of the riboswitch sequence
        self.tlen_in_switch: list[int] = []

        for pair in pairs:
            readIsProcessed = False

            # Skip pair if in RF or an errant orientation
            if pair.orientation in ["RF", "TANDEM", "RRF", "FFR"]:  # Uh oh
                continue

            elif pair.in_pair == 1:
                if not self.allowSingleReads:
                    continue

                readIsProcessed = self._process_ends(pair.read, self.allowSoftClips)
                if not readIsProcessed:
                    continue

            # There's no functional difference between handling "FR" and "EQ"
            # orientations.
            elif pair.orientation == "FR" or pair.orientation == "EQ":
                fProcessed = self._process_ends(pair.f, self.allowSoftClips)
                rProcessed = self._process_ends(pair.r, self.allowSoftClips)

                if not fProcessed or not rProcessed:
                    continue

            # Check if fragment intersects with the riboswitch
            if (
                self.switch_start in range(*pair.fragment_termini)
                or self.switch_end in range(*pair.fragment_termini)
                or pair.fragment_termini[0] in range(self.switch_start, self.switch_end)
            ):
                self.tlen_in_switch.append(pair.template_length)

            self._coalesce_into_cov(pair, self.allowSoftClips)
            self._process_frag_locus(pair)

        self.sum_cov()
        return self

    def sum_cov(self) -> None:
        """Sum the coverage arrays into an attribute. Also computes average fragment coverage across the riboswitch."""
        for covtype in self.__coverage_types:
            self.summedcov = np.add(self.summedcov, self.__getattribute__(covtype))

        # TODO: Check whether self.switch_start and self.switch_end are typical [) bounds.
        self.avg_switch_frag_cov: float = sum(
            self.summedcov[self.switch_start : self.switch_end]
        ) / abs(self.switch_end - self.switch_start)

    def is_coverage_threshold(
        self,
        covtype: str,
        thresh: int,
        l: Optional[int] = None,
        r: Optional[int] = None,
    ) -> bool:
        covpicker = {
            "read": self.readcov,
            "inferred": self.infercov,
            "clipped": self.clipcov,
            "summed": self.summedcov,
        }

        # Set default values
        l = self.switch_start if l is None else l
        r = self.switch_end if r is None else r

        return max(covpicker[covtype][l:r]) >= thresh

    def bin_rawends(self, bin_size: int = 1) -> None:
        """Bins then sums raw fragment ends within bins.

        Set attributes for binned ends and binned axis

        Args:
            bin_size (int, optional): Binning size. Defaults to 1.
        """
        self.bin_size = bin_size

        self.bin_ax = [i for i in range(0, len(self.readcov), bin_size)]
        numbins = len(self.bin_ax)

        self.binned_ends = np.ones(numbins)
        for i in range(numbins - 1):
            bstart = self.bin_ax[i]
            bend = self.bin_ax[i + 1]

            self.binned_ends[i] = sum(self.rawends[bstart:bend])

        self.binned_ends[numbins - 1] = sum(self.rawends[self.bin_ax[numbins - 1] :])

    def convolve_rawends(self, kernel) -> None:
        """Convolve the raw ends array with the normally-weighted kernal.

        Args:
            kernel (np.ndarray): The 1-D convolution kernel with normally distributed weights.
        """
        self.convends = np.convolve(self.rawends, kernel, "same")


class SegregatedAlignDat(AlignDat):
    def __init__(
        self, ref: str, reflength: int, boundsDict: dict[str, list], roi: float
    ) -> None:
        self.ref = ref
        self.ref_length = reflength
        self.bounds_dict = boundsDict
        self.roi = abs(roi)

        self.of_3prime: AlignDat
        self.except_3prime: AlignDat
        self.total: AlignDat

        self.ref_bounds = boundsDict["ref"]
        self.orf_bounds = boundsDict["orf"]

        self.switch_class = ref.split("#")[0]

        self._set_refrelative_switch_bounds()

    def process_pairs(
        self, pairs: list[ReadPair], allowSoftClips: bool, allowSingleReads: bool
    ) -> "SegregatedAlignDat":

        #  FIX LATER HOLY SHIT
        # marg = int(
        #     abs(self.switch_end - self.switch_start)
        #     * fil.DefaultThresholds.relative_size_bounds
        # )

        marg = int(abs(self.switch_end - self.switch_start) * self.roi)
        endmargins = (self.switch_end - marg, self.switch_end + marg)

        _pairs_of3prime: list[ReadPair] = []
        _pairs_except3prime: list[ReadPair] = []

        for pair in pairs:
            end_pos = pair.fragment_termini[1]

            if end_pos >= endmargins[0] and end_pos <= endmargins[1]:
                _pairs_of3prime.append(pair)
            else:
                _pairs_except3prime.append(pair)

        self.of_3prime = AlignDat(
            self.ref, self.ref_length, self.bounds_dict
        ).process_pairs(_pairs_of3prime, allowSoftClips, allowSingleReads)

        self.except_3prime = AlignDat(
            self.ref, self.ref_length, self.bounds_dict
        ).process_pairs(_pairs_except3prime, allowSoftClips, allowSingleReads)

        self.total = AlignDat(
            self.ref, self.ref_length, self.bounds_dict
        ).process_pairs(pairs, allowSoftClips, allowSingleReads)

        return self

    def sum_cov(self) -> None:
        self.total.sum_cov()
        self.of_3prime.sum_cov()
        self.except_3prime.sum_cov()

    def is_coverage_threshold(
        self,
        covtype: str,
        thresh: int,
        l: Optional[int] = None,
        r: Optional[int] = None,
    ) -> bool:
        return self.total.is_coverage_threshold(covtype, thresh, l, r)

    def bin_rawends(self, bin_size: int = 1) -> None:
        self.total.bin_rawends(bin_size)
        self.of_3prime.bin_rawends(bin_size)
        self.except_3prime.bin_rawends(bin_size)

    def convolve_rawends(self, kernel) -> None:
        self.total.convolve_rawends(kernel)
        self.of_3prime.convolve_rawends(kernel)
        self.except_3prime.convolve_rawends(kernel)


class SortedBAM:
    """Wrapper around a sorted BAM file to abstract away feature extraction"""

    def __init__(self, bam_path: str, ref_locus_dict_path: str) -> None:
        self._bam: pysam.AlignmentFile = pysam.AlignmentFile(bam_path, "rb")

        # Map storing the locus bounds of the full reference sequence
        # in the genome/MAG it was sliced from. These bounds are used
        # to offset the reference start on position axes in plots down to 0.
        with open(ref_locus_dict_path, "r") as f:
            self.ref_loc_dict: dict[str, dict[str, list]] = json.load(f)

        self._set_reference_list()

    def _set_reference_list(self):
        """Make a list of references in this sorted BAM; only considering references with aligned reads"""
        self.ref_list = [
            idxstats.contig
            for idxstats in self._bam.get_index_statistics()
            if idxstats.total > 0
        ]

    def reference_length(self, ref: str) -> int:
        """Wrapper to return the size of a reference in BAM.

        Args:
            ref (str): Reference name.

        Raises:
            ValueError: No reference in BAM with more than 0 aligned reads that matches ref.

        Returns:
            int: Size of the reference in base pairs.
        """
        if ref not in self.ref_list:
            raise ValueError(f"Reference '{ref}' with aligned reads not found in BAM.")

        return self._bam.get_reference_length(ref)

    def fetch_pairs(self, ref: str) -> list[ReadPair]:
        """Fetch read pairs for a reference in the BAM and return as a list.

        Args:
            ref (str): Reference name.

        Returns:
            list[ReadPair]: List of ReadPair objects of alignments to ref.
        """
        temp = defaultdict(list[Read])

        # Arrange pairs by matchig reads with same names
        for pys_read in self._bam.fetch(contig=ref, multiple_iterators=True):
            read = Read(pys_read)
            temp[read.name].append(read)

        toRet = []
        for _, pairlist in temp.items():
            try:
                toRet.append(ReadPair(*pairlist))
            except ValueError:
                echo(f"More than 2 reads assigned to pair: {pairlist}")

        return toRet

    def generate_ref_alignment_data(
        self, allowSoftClips: bool, allowSingleReads: bool, roi: float
    ) -> list[SegregatedAlignDat]:
        """Generate alignment data for each reference found in the sorted BAM using all aligned reads/read pairs and return as a list.

        Args:
            allowSoftClips (bool): If true, reads with soft-clipping are processed. Soft-clipped reads are thrown out otherwise.
            allowSingleReads (bool): If true, reads with unmapped mates are processed. Single reads with unmapped mates are thrown out otherwise.

        Returns:
            list[AlignDat]: Alignment data for each reference.
        """
        data: list[SegregatedAlignDat] = []

        for ref in self.ref_list:
            # The + 1 is to offset downstream plot x-axis to seem 1-indexed
            reflength = self.reference_length(ref) + 1
            readpairs = self.fetch_pairs(ref)
            boundDict = self.ref_loc_dict[ref]

            data.append(
                SegregatedAlignDat(ref, reflength, boundDict, roi).process_pairs(
                    readpairs, allowSoftClips, allowSingleReads
                )
            )

        return data
