import pysam
from dataclasses import dataclass
from patch import GenomicInterval, PatchPaths


@dataclass
class ReadPairStatus:
    """
    Tracks the primary alignment coordinates for mate 1 and mate 2.

    For deletion editing, a read pair is removed only if the entire
    paired-end fragment (or full singleton fragment) is fully contained inside the
    deletion interval.

    In pysam:
        reference_start is 0-based inclusive
        reference_end is 0-based exclusive

    A fragment is fully internal if:
        pair_start >= interval.start0
        pair_end   <= interval.end0

    where:
        pair_start = min(read1_start, read2_start)
        pair_end   = max(read1_end, read2_end)
    """

    read1_chrom: str | None = None
    read1_start: int | None = None
    read1_end: int | None = None

    read2_chrom: str | None = None
    read2_start: int | None = None
    read2_end: int | None = None

    def add_read(self, read: pysam.AlignedSegment) -> None:
        """
        Add coordinates from one alignment record.

        Only primary mapped read1/read2 records are used to define the
        fragment span.

        Secondary and supplementary records are ignored for span definition,
        but if the qname is later selected for removal, all records with that
        qname will be removed from the edited patch BAM.
        """
        if read.is_unmapped:
            return

        if read.is_secondary or read.is_supplementary:
            return

        if read.reference_name is None:
            return

        if read.reference_start is None or read.reference_end is None:
            return

        if read.is_read1:
            self.read1_chrom = read.reference_name
            self.read1_start = read.reference_start
            self.read1_end = read.reference_end

        elif read.is_read2:
            self.read2_chrom = read.reference_name
            self.read2_start = read.reference_start
            self.read2_end = read.reference_end

    @property
    def has_both_mates(self) -> bool:
        """
        Return True if both primary mate coordinates were found.
        """
        return (
            self.read1_chrom is not None
            and self.read1_start is not None
            and self.read1_end is not None
            and self.read2_chrom is not None
            and self.read2_start is not None
            and self.read2_end is not None
        )

    def both_mates_on_chrom(self, chrom: str) -> bool:
        """
        Return True if both primary mates map to the requested chromosome.
        """
        if not self.has_both_mates:
            return False

        return self.read1_chrom == chrom and self.read2_chrom == chrom

    @property
    def pair_start(self) -> int | None:
        """
        Leftmost coordinate of the paired-end fragment.
        """
        if not self.has_both_mates:
            return None

        return min(self.read1_start, self.read2_start)

    @property
    def pair_end(self) -> int | None:
        """
        Rightmost coordinate of the paired-end fragment.
        """
        if not self.has_both_mates:
            return None

        return max(self.read1_end, self.read2_end)

    def is_fully_inside_interval(self, interval: GenomicInterval) -> bool:
        """
        Return True if this qname should be removed for deletion editing.

        Removal rules:

        1. If both primary mates are present and both map to the CNV chromosome:
           remove only if the full paired-end fragment span is inside the interval.

        2. If one primary mate maps to the CNV chromosome and the other maps elsewhere:
           remove if the CNV-chromosome mate is fully inside the interval.

        3. If only one primary mate is present:
           remove if that single primary alignment is fully inside the interval.

        This keeps reads that cross the deletion breakpoints.
        """

        read1_on_target = (
                self.read1_chrom == interval.chrom
                and self.read1_start is not None
                and self.read1_end is not None
        )

        read2_on_target = (
                self.read2_chrom == interval.chrom
                and self.read2_start is not None
                and self.read2_end is not None
        )

        read1_inside = (
                read1_on_target
                and self.read1_start >= interval.start0
                and self.read1_end <= interval.end0
        )

        read2_inside = (
                read2_on_target
                and self.read2_start >= interval.start0
                and self.read2_end <= interval.end0
        )

        read1_crosses_left = (
                read1_on_target
                and self.read1_start < interval.start0
                and self.read1_end > interval.start0
        )

        read1_crosses_right = (
                read1_on_target
                and self.read1_start < interval.end0
                and self.read1_end > interval.end0
        )

        read2_crosses_left = (
                read2_on_target
                and self.read2_start < interval.start0
                and self.read2_end > interval.start0
        )

        read2_crosses_right = (
                read2_on_target
                and self.read2_start < interval.end0
                and self.read2_end > interval.end0
        )

        crosses_breakpoint = (
                read1_crosses_left
                or read1_crosses_right
                or read2_crosses_left
                or read2_crosses_right
        )

        if crosses_breakpoint:
            return False

        # Case 1: both mates are present and both are on the CNV chromosome.
        if self.has_both_mates and self.both_mates_on_chrom(interval.chrom):
            return (
                    self.pair_start >= interval.start0
                    and self.pair_end <= interval.end0
            )

        # Case 2: at least one observed primary alignment is fully inside
        # the deletion interval, and none crosses a breakpoint.
        return read1_inside or read2_inside