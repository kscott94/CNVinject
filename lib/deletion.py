#!/usr/bin/env python3

from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict

import pysam

from patch import GenomicInterval, PatchPaths
from helpers import sort_and_index_bam

@dataclass
class ReadPairStatus:
    """
    Tracks the primary alignment coordinates for mate 1 and mate 2.

    For CN0 deletion editing, a read pair is removed only if the entire
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
        Return True if this qname should be removed for CN0 deletion editing.

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


@dataclass
class CN0DeletionResult:
    """
    Summary of CN0 deletion patch editing.
    """

    input_patch_bam: Path
    edited_patch_bam: Path
    internal_qnames: Path
    n_internal_qnames: int
    n_input_records: int
    n_removed_records: int
    n_kept_records: int


class CN0DeletionEditor:
    """
    Edit a raw patch BAM for a CN0 deletion.

    This editor removes read pairs whose full paired-end fragment span is
    contained inside the deletion interval.

    It keeps:
    - breakpoint-crossing fragments
    - one-mate-inside / one-mate-outside fragments when the pair span extends outside the interval
    - singleton reads that cross a breakpoint or extend outside the interval
    - flanking reads
    - pairs mapping across chromosomes

    If a qname is selected for removal, all records with that qname are removed
    from the edited patch BAM, including secondary/supplementary records.
    """

    def __init__(
        self,
        input_patch_bam: str | Path,
        output_prefix: str | Path,
        interval: GenomicInterval,
        threads: int = 1,
    ):
        self.input_patch_bam = Path(input_patch_bam)
        self.paths = PatchPaths(output_prefix)
        self.interval = interval
        self.threads = threads

        self.edited_patch_bam = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.edited.patch.bam"
        )

        self.internal_qnames = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.internal.qnames.txt"
        )

        self.status_by_qname: dict[str, ReadPairStatus] = defaultdict(ReadPairStatus)
        self.qnames_to_remove: set[str] = set()

        self.n_input_records = 0
        self.n_removed_records = 0
        self.n_kept_records = 0

    def classify_patch_reads(self) -> None:
        """
        Classify read pairs by their full fragment span.

        A qname is marked for removal only if:
            - both primary mates are present
            - both primary mates map to the CNV chromosome
            - the full pair span is inside the deletion interval
        """
        with pysam.AlignmentFile(self.input_patch_bam, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.query_name is None:
                    continue

                status = self.status_by_qname[read.query_name]
                status.add_read(read)

        self.qnames_to_remove = {
            qname
            for qname, status in self.status_by_qname.items()
            if status.is_fully_inside_interval(self.interval)
        }

    def write_internal_qnames(self) -> None:
        """
        Write qnames that will be removed from the patch.
        """
        with open(self.internal_qnames, "w") as handle:
            for qname in sorted(self.qnames_to_remove):
                handle.write(qname + "\n")

    def write_edited_patch_bam(self) -> None:
        """
        Write edited patch BAM, excluding fully internal read pairs.
        """
        with pysam.AlignmentFile(self.input_patch_bam, "rb") as bam_in:
            with pysam.AlignmentFile(self.edited_patch_bam, "wb", template=bam_in) as bam_out:
                for read in bam_in.fetch(until_eof=True):
                    self.n_input_records += 1

                    if read.query_name in self.qnames_to_remove:
                        self.n_removed_records += 1
                        continue

                    bam_out.write(read)
                    self.n_kept_records += 1

    def sort_and_index_edited_patch(self) -> None:
        """
        Sort and index the edited patch BAM.
        """
        unsorted_bam = self.edited_patch_bam.with_name(
            f"{self.paths.prefix.name}.edited.patch.unsorted.tmp.bam"
        )

        # Rename the just-written BAM to a temporary unsorted name.
        self.edited_patch_bam.rename(unsorted_bam)

        sort_and_index_bam(
            input_bam=unsorted_bam,
            output_bam=self.edited_patch_bam,
            threads=self.threads,
            remove_input=True,
        )

    def run(self) -> CN0DeletionResult:
        """
        Run full CN0 deletion patch editing.
        """
        self.classify_patch_reads()
        self.write_internal_qnames()
        self.write_edited_patch_bam()
        self.sort_and_index_edited_patch()

        result = CN0DeletionResult(
            input_patch_bam=self.input_patch_bam,
            edited_patch_bam=self.edited_patch_bam,
            internal_qnames=self.internal_qnames,
            n_internal_qnames=len(self.qnames_to_remove),
            n_input_records=self.n_input_records,
            n_removed_records=self.n_removed_records,
            n_kept_records=self.n_kept_records,
        )

        self.print_summary(result)
        return result

    def print_summary(self, result: CN0DeletionResult) -> None:
        print("CN0 deletion patch editing complete")
        print(f"Deletion interval: {self.interval}")
        print(f"Input patch BAM: {result.input_patch_bam}")
        print(f"Edited patch BAM: {result.edited_patch_bam}")
        print(f"Internal qnames removed: {result.internal_qnames}")
        print(f"Fully internal fragment qnames: {result.n_internal_qnames:,}")
        print(f"Input patch records: {result.n_input_records:,}")
        print(f"Removed records: {result.n_removed_records:,}")
        print(f"Kept records: {result.n_kept_records:,}")