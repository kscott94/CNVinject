#!/usr/bin/env python3

from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
import csv

import pysam

from patch import GenomicInterval, PatchPaths


@dataclass
class PrimaryAlignment:
    """
    Stores one primary mapped alignment for a read.
    """

    chrom: str
    start: int
    end: int
    is_read1: bool
    is_read2: bool
    is_reverse: bool
    query_name: str

    @classmethod
    def from_read(cls, read: pysam.AlignedSegment) -> "PrimaryAlignment | None":
        """
        Create a PrimaryAlignment from a pysam read.

        Only primary mapped read1/read2 records are used.
        """
        if read.query_name is None:
            return None

        if read.is_unmapped:
            return None

        if read.is_secondary or read.is_supplementary:
            return None

        if read.reference_name is None:
            return None

        if read.reference_start is None or read.reference_end is None:
            return None

        if not (read.is_read1 or read.is_read2):
            return None

        return cls(
            chrom=read.reference_name,
            start=read.reference_start,
            end=read.reference_end,
            is_read1=read.is_read1,
            is_read2=read.is_read2,
            is_reverse=read.is_reverse,
            query_name=read.query_name,
        )

    def overlaps_interval(self, interval: GenomicInterval) -> bool:
        if self.chrom != interval.chrom:
            return False

        return self.start < interval.end0 and self.end > interval.start0

    def fully_inside_interval(self, interval: GenomicInterval) -> bool:
        if self.chrom != interval.chrom:
            return False

        return self.start >= interval.start0 and self.end <= interval.end0

    def left_of_interval(self, interval: GenomicInterval) -> bool:
        if self.chrom != interval.chrom:
            return False

        return self.end <= interval.start0

    def right_of_interval(self, interval: GenomicInterval) -> bool:
        if self.chrom != interval.chrom:
            return False

        return self.start >= interval.end0

    def crosses_left_breakpoint(self, interval: GenomicInterval) -> bool:
        if self.chrom != interval.chrom:
            return False

        return self.start < interval.start0 and self.end > interval.start0

    def crosses_right_breakpoint(self, interval: GenomicInterval) -> bool:
        if self.chrom != interval.chrom:
            return False

        return self.start < interval.end0 and self.end > interval.end0

    def position_class(self, interval: GenomicInterval) -> str:
        """
        Classify this single alignment relative to the deletion interval.
        """
        if self.chrom != interval.chrom:
            return "off_chrom"

        if self.fully_inside_interval(interval):
            return "inside"

        if self.crosses_left_breakpoint(interval):
            return "crosses_left_breakpoint"

        if self.crosses_right_breakpoint(interval):
            return "crosses_right_breakpoint"

        if self.left_of_interval(interval):
            return "left_flank"

        if self.right_of_interval(interval):
            return "right_flank"

        if self.overlaps_interval(interval):
            return "partial_overlap"

        return "outside"


@dataclass
class ReadPairBreakpointStatus:
    """
    Stores primary mate alignments for one qname and classifies its
    breakpoint relevance.
    """

    qname: str
    read1: PrimaryAlignment | None = None
    read2: PrimaryAlignment | None = None

    def add_alignment(self, aln: PrimaryAlignment) -> None:
        if aln.is_read1:
            self.read1 = aln
        elif aln.is_read2:
            self.read2 = aln

    @property
    def has_read1(self) -> bool:
        return self.read1 is not None

    @property
    def has_read2(self) -> bool:
        return self.read2 is not None

    @property
    def has_both_mates(self) -> bool:
        return self.read1 is not None and self.read2 is not None

    def any_on_target_overlap(self, interval: GenomicInterval) -> bool:
        """
        Return True if either primary alignment overlaps the CNV interval.
        """
        return (
            (self.read1 is not None and self.read1.overlaps_interval(interval))
            or (self.read2 is not None and self.read2.overlaps_interval(interval))
        )

    def any_crosses_breakpoint(self, interval: GenomicInterval) -> bool:
        """
        Return True if either primary alignment crosses either deletion boundary.
        """
        for aln in [self.read1, self.read2]:
            if aln is None:
                continue

            if aln.crosses_left_breakpoint(interval) or aln.crosses_right_breakpoint(interval):
                return True

        return False

    def get_read_position_classes(self, interval: GenomicInterval) -> tuple[str, str]:
        read1_class = "missing"
        read2_class = "missing"

        if self.read1 is not None:
            read1_class = self.read1.position_class(interval)

        if self.read2 is not None:
            read2_class = self.read2.position_class(interval)

        return read1_class, read2_class

    def classify(self, interval: GenomicInterval) -> str:
        """
        Classify qname-level breakpoint pattern.

        This assumes fully internal qnames have already been removed from
        the edited patch BAM.
        """
        r1_class, r2_class = self.get_read_position_classes(interval)

        # Singleton cases.
        if self.read1 is not None and self.read2 is None:
            if r1_class in ["crosses_left_breakpoint", "crosses_right_breakpoint", "partial_overlap"]:
                return f"singleton_read1_{r1_class}"
            if r1_class == "inside":
                return "singleton_read1_inside_remaining"
            return "not_breakpoint_candidate"

        if self.read2 is not None and self.read1 is None:
            if r2_class in ["crosses_left_breakpoint", "crosses_right_breakpoint", "partial_overlap"]:
                return f"singleton_read2_{r2_class}"
            if r2_class == "inside":
                return "singleton_read2_inside_remaining"
            return "not_breakpoint_candidate"

        # Complete pair cases.
        if self.read1 is not None and self.read2 is not None:
            classes = {r1_class, r2_class}

            if "crosses_left_breakpoint" in classes or "crosses_right_breakpoint" in classes:
                return f"pair_with_breakpoint_crossing_read.r1_{r1_class}.r2_{r2_class}"

            if r1_class == "inside" and r2_class == "left_flank":
                return "pair_read1_inside_read2_left_flank"

            if r1_class == "inside" and r2_class == "right_flank":
                return "pair_read1_inside_read2_right_flank"

            if r2_class == "inside" and r1_class == "left_flank":
                return "pair_read2_inside_read1_left_flank"

            if r2_class == "inside" and r1_class == "right_flank":
                return "pair_read2_inside_read1_right_flank"

            if r1_class == "inside" and r2_class == "off_chrom":
                return "pair_read1_inside_read2_off_chrom_remaining"

            if r2_class == "inside" and r1_class == "off_chrom":
                return "pair_read2_inside_read1_off_chrom_remaining"

            if r1_class == "left_flank" and r2_class == "right_flank":
                return "pair_spans_deletion_left_to_right"

            if r2_class == "left_flank" and r1_class == "right_flank":
                return "pair_spans_deletion_right_to_left"

            if self.any_on_target_overlap(interval):
                return f"pair_other_interval_overlap.r1_{r1_class}.r2_{r2_class}"

        return "not_breakpoint_candidate"

    def is_breakpoint_candidate(self, interval: GenomicInterval) -> bool:
        classification = self.classify(interval)

        if classification == "not_breakpoint_candidate":
            return False

        # These should ideally have been removed during internal deletion editing.
        # If they appear, report them in the TSV but do not use them as breakpoint
        # synthesis candidates yet.
        if classification.endswith("_inside_remaining"):
            return False

        if "off_chrom_remaining" in classification:
            return False

        return True

    def to_tsv_row(self, interval: GenomicInterval, classification: str) -> dict[str, str | int]:
        r1_class, r2_class = self.get_read_position_classes(interval)

        return {
            "qname": self.qname,
            "classification": classification,
            "is_breakpoint_candidate": str(self.is_breakpoint_candidate(interval)),
            "read1_chrom": self.read1.chrom if self.read1 else "NA",
            "read1_start0": self.read1.start if self.read1 else "NA",
            "read1_end0": self.read1.end if self.read1 else "NA",
            "read1_class": r1_class,
            "read1_reverse": str(self.read1.is_reverse) if self.read1 else "NA",
            "read2_chrom": self.read2.chrom if self.read2 else "NA",
            "read2_start0": self.read2.start if self.read2 else "NA",
            "read2_end0": self.read2.end if self.read2 else "NA",
            "read2_class": r2_class,
            "read2_reverse": str(self.read2.is_reverse) if self.read2 else "NA",
        }


@dataclass
class BreakpointCandidateResult:
    """
    Summary of breakpoint candidate classification.
    """

    candidates_tsv: Path
    n_qnames_seen: int
    n_candidates: int
    n_tsv_rows: int


class BreakpointCandidateClassifier:
    """
    Classify breakpoint-candidate qnames from an edited patch BAM.

    This should be run after CN0DeletionEditor has removed fully internal
    reads from the patch.
    """

    def __init__(
        self,
        edited_patch_bam: str | Path,
        output_prefix: str | Path,
        interval: GenomicInterval,
    ):
        self.edited_patch_bam = Path(edited_patch_bam)
        self.paths = PatchPaths(output_prefix)
        self.interval = interval

        self.candidates_tsv = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.breakpoint.candidates.tsv"
        )

        self.status_by_qname: dict[str, ReadPairBreakpointStatus] = {}

    def collect_primary_alignments(self) -> None:
        """
        Read the edited patch BAM and collect primary read1/read2 alignments
        by qname.
        """
        with pysam.AlignmentFile(self.edited_patch_bam, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                aln = PrimaryAlignment.from_read(read)

                if aln is None:
                    continue

                if aln.query_name not in self.status_by_qname:
                    self.status_by_qname[aln.query_name] = ReadPairBreakpointStatus(
                        qname=aln.query_name
                    )

                self.status_by_qname[aln.query_name].add_alignment(aln)

    def write_outputs(self) -> tuple[int, int]:
        """
        Write:
            PREFIX.breakpoint.candidates.tsv
            PREFIX.breakpoint.qnames.txt

        Returns:
            n_candidates, n_tsv_rows
        """
        fieldnames = [
            "qname",
            "classification",
            "is_breakpoint_candidate",
            "read1_chrom",
            "read1_start0",
            "read1_end0",
            "read1_class",
            "read1_reverse",
            "read2_chrom",
            "read2_start0",
            "read2_end0",
            "read2_class",
            "read2_reverse",
        ]

        candidate_qnames: set[str] = set()
        n_tsv_rows = 0

        with open(self.candidates_tsv, "w", newline="") as tsv_handle:
            writer = csv.DictWriter(
                tsv_handle,
                fieldnames=fieldnames,
                delimiter="\t",
            )

            writer.writeheader()

            for qname in sorted(self.status_by_qname):
                status = self.status_by_qname[qname]
                classification = status.classify(self.interval)

                # Skip qnames with no relevance to the deletion interval.
                if classification == "not_breakpoint_candidate":
                    continue

                row = status.to_tsv_row(
                    interval=self.interval,
                    classification=classification,
                )
                writer.writerow(row)
                n_tsv_rows += 1

                if status.is_breakpoint_candidate(self.interval):
                    candidate_qnames.add(qname)

        return len(candidate_qnames), n_tsv_rows

    def run(self) -> BreakpointCandidateResult:
        self.collect_primary_alignments()

        n_candidates, n_tsv_rows = self.write_outputs()

        result = BreakpointCandidateResult(
            candidates_tsv=self.candidates_tsv,
            n_qnames_seen=len(self.status_by_qname),
            n_candidates=n_candidates,
            n_tsv_rows=n_tsv_rows,
        )

        self.print_summary(result)
        return result

    def print_summary(self, result: BreakpointCandidateResult) -> None:
        print("Breakpoint read classification complete")
        print(f"Edited patch BAM: {self.edited_patch_bam}")
        print(f"Deletion interval: {self.interval}")
        print(f"Candidate TSV: {result.candidates_tsv}")
        print(f"Primary alignments: {result.n_qnames_seen:,}")
        print(f"TSV rows written: {result.n_tsv_rows:,}")
        print(f"Breakpoint candidate qnames: {result.n_candidates:,}")