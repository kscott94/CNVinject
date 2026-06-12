#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
import random
import subprocess

import pysam

from patch import GenomicInterval, PatchPaths
from helpers import (
    sort_and_index_bam,
    align_fastq_with_bwa,
    fastq_has_records,
    reverse_complement,
)


@dataclass
class BoundaryReadTemplate:
    donor_bam: Path
    qname: str
    side: str  # "left" or "right"


@dataclass
class DuplicationBreakpointResult:
    added_outer_breakpoint_records_bam: Path
    added_outer_breakpoint_jittered_bam: Path
    tandem_junction_fastq: Path
    tandem_junction_bam: Path
    candidate_tsv: Path
    selected_outer_qnames_tsv: Path
    selected_junction_qnames_tsv: Path

    n_recipient_left_boundary_qnames: int
    n_recipient_right_boundary_qnames: int

    n_requested_left_outer_qnames: int
    n_requested_right_outer_qnames: int
    n_selected_left_outer_qnames: int
    n_selected_right_outer_qnames: int
    n_added_outer_records: int

    n_requested_left_junction_templates: int
    n_requested_right_junction_templates: int
    n_selected_left_junction_templates: int
    n_selected_right_junction_templates: int
    n_tandem_junction_fastq_records: int


def added_fraction_from_copy_number(copy_number: float) -> float:
    copy_number = float(copy_number)
    if copy_number <= 2:
        raise ValueError("Duplication breakpoint workflow requires copy_number > 2.")
    return (copy_number - 2.0) / 2.0


def read_is_primary_usable(read: pysam.AlignedSegment, mapq: int) -> bool:
    if read.is_unmapped:
        return False
    if read.is_secondary or read.is_supplementary:
        return False
    if read.query_name is None:
        return False
    if read.reference_start is None or read.reference_end is None:
        return False
    if read.mapping_quality < mapq:
        return False
    return True


def read_overlaps_pos(read: pysam.AlignedSegment, pos0: int) -> bool:
    if read.reference_start is None or read.reference_end is None:
        return False
    return read.reference_start <= pos0 < read.reference_end


def read_fully_inside_interval(
    read: pysam.AlignedSegment,
    interval: GenomicInterval,
) -> bool:
    if read.reference_start is None or read.reference_end is None:
        return False
    return (
        read.reference_name == interval.chrom
        and read.reference_start >= interval.start0
        and read.reference_end <= interval.end0
    )

def scaled_count(n: int, fraction: float) -> int:
    if n == 0:
        return 0

    value = n * fraction

    if value > 0 and value < 1:
        return 1

    return round(value)

def qualities_to_fastq_string(read: pysam.AlignedSegment, length: int) -> str:
    if read.query_qualities is None:
        return "I" * length

    quals = "".join(chr(min(max(q, 0), 93) + 33) for q in read.query_qualities)

    if len(quals) != length:
        return "I" * length

    return quals


def reference_sequence_for_shifted_read(
    read: pysam.AlignedSegment,
    reference: pysam.FastaFile,
    shifted_start: int,
) -> str:
    """
    Build a shifted read sequence from the reference while preserving CIGAR shape.

    This keeps soft-clipped and inserted sequence from the original read,
    advances over deletions/skips, and fills aligned match/mismatch blocks
    from the shifted reference position.

    This is intentionally simpler than the deletion synthetic-read mismatch replay:
    it produces a reference-consistent shifted molecule while preserving CIGAR artifacts.
    """
    if read.query_sequence is None:
        raise ValueError(f"Read has no sequence: {read.query_name}")

    if read.cigartuples is None:
        raise ValueError(f"Read has no CIGAR: {read.query_name}")

    chrom = read.reference_name
    ref_pos = shifted_start
    query_pos = 0
    out = []

    original_seq = read.query_sequence

    for op, length in read.cigartuples:
        # M, =, X: consume query and reference
        if op in (0, 7, 8):
            seq = reference.fetch(chrom, ref_pos, ref_pos + length).upper()
            out.append(seq)
            ref_pos += length
            query_pos += length

        # I: consume query only; preserve inserted bases
        elif op == 1:
            out.append(original_seq[query_pos : query_pos + length])
            query_pos += length

        # D or N: consume reference only
        elif op in (2, 3):
            ref_pos += length

        # S: consume query only; preserve soft-clipped bases
        elif op == 4:
            out.append(original_seq[query_pos : query_pos + length])
            query_pos += length

        # H or P: consume neither query nor reference
        elif op in (5, 6):
            continue

        else:
            raise ValueError(f"Unsupported CIGAR op {op} in {read.query_name}")

    seq = "".join(out).upper()

    # Preserve N positions from the original read.
    if len(seq) == len(original_seq):
        seq_list = list(seq)
        for i, base in enumerate(original_seq.upper()):
            if base == "N":
                seq_list[i] = "N"
        seq = "".join(seq_list)

    if read.is_reverse:
        seq = reverse_complement(seq)

    return seq


class DuplicationBreakpointEditor:
    """
    Add duplication breakpoint evidence.

    Outputs:
      - added.outer.breakpoint.records.bam:
          copied donor left/right outer-boundary reads, renamed only
      - tandem.junction.synthetic.reads.fastq:
          synthetic reads representing interval_end -> interval_start
      - tandem.junction.synthetic.reads.bam:
          BWA alignment of the synthetic FASTQ
    """

    def __init__(
        self,
        recipient_patch_bam: str | Path,
        donor_bam_dir: str | Path,
        recipient_input_bam: str | Path,
        output_prefix: str | Path,
        interval: GenomicInterval,
        copy_number: float,
        seed: int,
        reference_fasta: str | Path,
        mapq: int = 0,
        threads: int = 1,
        bwa_args: str | None = None,
        allow_replacement: bool = False,
        jitter_bp: int = 10,
    ):
        self.recipient_patch_bam = Path(recipient_patch_bam)
        self.donor_bam_dir = Path(donor_bam_dir)
        self.recipient_input_bam = Path(recipient_input_bam).resolve()
        self.paths = PatchPaths(output_prefix)
        self.interval = interval
        self.copy_number = float(copy_number)
        self.seed = int(seed)
        self.reference_fasta = Path(reference_fasta)
        self.mapq = int(mapq)
        self.threads = int(threads)
        self.bwa_args = bwa_args
        self.allow_replacement = bool(allow_replacement)
        self.jitter_bp = int(jitter_bp)

        self.rng = random.Random(self.seed + 99173)

        self.candidate_tsv = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.dup.breakpoint.candidates.tsv"
        )

        self.selected_outer_qnames_tsv = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.donor.selected.outer.breakpoint.qnames.tsv"
        )

        self.selected_junction_qnames_tsv = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.donor.selected.tandem.junction.qnames.tsv"
        )

        self.added_outer_breakpoint_records_bam = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.added.outer.breakpoint.records.bam"
        )

        self.tandem_junction_fastq = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.tandem.junction.synthetic.reads.fastq"
        )

        self.tandem_junction_bam = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.tandem.junction.synthetic.reads.bam"
        )

        self.recipient_left_boundary_qnames: set[str] = set()
        self.recipient_right_boundary_qnames: set[str] = set()

        self.n_added_outer_records = 0
        self.n_tandem_junction_fastq_records = 0

        self.added_outer_breakpoint_jittered_fastq = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.added.outer.breakpoint.jittered.fastq"
        )

        self.added_outer_breakpoint_jittered_bam = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.added.outer.breakpoint.jittered.bam"
        )

        self.n_jittered_outer_fastq_records = 0

    def added_fraction(self) -> float:
        return added_fraction_from_copy_number(self.copy_number)

    def donor_bams(self) -> list[Path]:
        bams = sorted(self.donor_bam_dir.glob("*.bam"))

        usable: list[Path] = []
        recipient_name = self.recipient_input_bam.name

        for bam in bams:
            if bam.name.endswith(".bai"):
                continue

            if bam.name == recipient_name:
                continue

            try:
                if bam.resolve() == self.recipient_input_bam:
                    continue
            except FileNotFoundError:
                pass

            usable.append(bam)

        if not usable:
            raise ValueError(f"No usable donor BAMs found in {self.donor_bam_dir}")

        self.rng.shuffle(usable)
        return usable

    def classify_recipient_boundary_qnames(self) -> None:
        """
        Count recipient qnames overlapping the left or right outer boundary.

        These counts determine how many donor outer-boundary reads and tandem-junction
        templates should be added.
        """
        left_pos0 = self.interval.start0
        right_pos0 = self.interval.end0 - 1

        with pysam.AlignmentFile(self.recipient_patch_bam, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if not read_is_primary_usable(read, self.mapq):
                    continue

                if read.reference_name != self.interval.chrom:
                    continue

                # Fully internal reads are handled by added.internal.records.bam.
                if read_fully_inside_interval(read, self.interval):
                    continue

                if read_overlaps_pos(read, left_pos0):
                    self.recipient_left_boundary_qnames.add(read.query_name)

                if read_overlaps_pos(read, right_pos0):
                    self.recipient_right_boundary_qnames.add(read.query_name)

        print(
            f"Recipient left-boundary qnames: "
            f"{len(self.recipient_left_boundary_qnames):,}"
        )
        print(
            f"Recipient right-boundary qnames: "
            f"{len(self.recipient_right_boundary_qnames):,}"
        )

    def requested_left_count(self) -> int:
        return scaled_count(
            len(self.recipient_left_boundary_qnames),
            self.added_fraction(),
        )

    def requested_right_count(self) -> int:
        return scaled_count(
            len(self.recipient_right_boundary_qnames),
            self.added_fraction(),
        )

    def collect_donor_boundary_templates(self) -> dict[str, list[BoundaryReadTemplate]]:
        """
        Collect donor qnames that overlap the left or right outer boundary.

        Fully internal reads are excluded because internal reads are handled separately.
        """
        left_pos0 = self.interval.start0
        right_pos0 = self.interval.end0 - 1

        templates: dict[str, list[BoundaryReadTemplate]] = {
            "left": [],
            "right": [],
        }

        seen: set[tuple[Path, str, str]] = set()

        with open(self.candidate_tsv, "w") as out:
            out.write("donor_bam\tqname\tside\n")

            for donor_bam in self.donor_bams():
                donor_left = 0
                donor_right = 0

                with pysam.AlignmentFile(donor_bam, "rb") as bam:
                    for read in bam.fetch(
                        self.interval.chrom,
                        max(0, self.interval.start0 - 1000),
                        self.interval.end0 + 1000,
                    ):
                        if not read_is_primary_usable(read, self.mapq):
                            continue

                        if read.reference_name != self.interval.chrom:
                            continue

                        if read_fully_inside_interval(read, self.interval):
                            continue

                        side = None

                        if read_overlaps_pos(read, left_pos0):
                            side = "left"
                        elif read_overlaps_pos(read, right_pos0):
                            side = "right"

                        if side is None:
                            continue

                        key = (donor_bam, read.query_name, side)
                        if key in seen:
                            continue

                        seen.add(key)
                        templates[side].append(
                            BoundaryReadTemplate(
                                donor_bam=donor_bam,
                                qname=read.query_name,
                                side=side,
                            )
                        )
                        out.write(f"{donor_bam}\t{read.query_name}\t{side}\n")

                        if side == "left":
                            donor_left += 1
                        else:
                            donor_right += 1

                print(
                    f"Donor boundary candidates: {donor_bam.name}: "
                    f"left={donor_left:,}, right={donor_right:,}"
                )

        self.rng.shuffle(templates["left"])
        self.rng.shuffle(templates["right"])

        print(f"Total donor left-boundary candidates: {len(templates['left']):,}")
        print(f"Total donor right-boundary candidates: {len(templates['right']):,}")

        return templates

    def sample_templates(
        self,
        templates: list[BoundaryReadTemplate],
        n_requested: int,
        label: str,
    ) -> list[BoundaryReadTemplate]:
        if n_requested == 0:
            return []

        if len(templates) >= n_requested:
            return templates[:n_requested]

        if not self.allow_replacement:
            raise ValueError(
                f"Not enough donor {label} breakpoint templates.\n"
                f"Needed: {n_requested:,}\n"
                f"Available: {len(templates):,}\n"
                "Rerun with --allow-replacement to permit resampling donor qnames."
            )

        if not templates:
            raise ValueError(
                f"No donor {label} breakpoint templates available, "
                "so replacement sampling is impossible."
            )

        selected = list(templates)
        deficit = n_requested - len(selected)

        for _ in range(deficit):
            selected.append(self.rng.choice(templates))

        return selected

    def donor_label(self, donor_bam: Path) -> str:
        """
        Stable donor label for qname suffixes.

        Uses the BAM stem and sanitizes characters that could be annoying
        inside read names.
        """
        label = donor_bam.stem
        label = label.replace(" ", "_")
        label = label.replace("/", "_")
        label = label.replace("\\", "_")
        label = label.replace(":", "_")
        return label

    def renamed_outer_qname(
            self,
            qname: str,
            side: str,
            donor_bam: Path,
            occurrence_index: int,
    ) -> str:
        donor_label = self.donor_label(donor_bam)

        return (
            f"{qname}:CNVinject_outer_{side}_"
            f"{donor_label}_copy{occurrence_index + 1}"
        )

    def renamed_junction_qname(
        self,
        qname: str,
        side: str,
        occurrence_index: int,
    ) -> str:
        return (
            f"{qname}:CNVinject_tandem_junction_"
            f"{side}_copy{occurrence_index + 1}"
        )

    def write_added_outer_breakpoint_records_bam(
        self,
        selected_left: list[BoundaryReadTemplate],
        selected_right: list[BoundaryReadTemplate],
    ) -> None:
        """
        Copy donor outer-boundary records unchanged except qname.

        If a qname is sampled multiple times because --allow-replacement is used,
        extra replacement copies are written separately as FASTQ by
        write_jittered_outer_breakpoint_fastq() and realigned with BWA.
        """
        selected = selected_left + selected_right

        if not selected:
            # Create an empty BAM using first donor as template.
            first_donor = self.donor_bams()[0]
        else:
            first_donor = selected[0].donor_bam

        unsorted_bam = self.added_outer_breakpoint_records_bam.with_name(
            f"{self.paths.prefix.name}.added.outer.breakpoint.records.unsorted.tmp.bam"
        )

        selected_by_donor_side: dict[tuple[Path, str], list[BoundaryReadTemplate]] = defaultdict(list)
        for template in selected:
            selected_by_donor_side[(template.donor_bam, template.side)].append(template)

        with pysam.AlignmentFile(first_donor, "rb") as template_bam:
            with pysam.AlignmentFile(unsorted_bam, "wb", template=template_bam) as bam_out:
                with open(self.selected_outer_qnames_tsv, "w") as qout:
                    qout.write(
                        "donor_bam\toriginal_qname\trenamed_qname\tside\toccurrence_index\n"
                    )

                    for donor_index, ((donor_bam, side), templates) in enumerate(
                        selected_by_donor_side.items()
                    ):
                        counts_by_qname: dict[str, int] = defaultdict(int)
                        selected_qnames = {template.qname for template in templates}

                        rename_by_qname_occurrence: dict[tuple[str, int], str] = {}

                        for template in templates:
                            occurrence = counts_by_qname[template.qname]
                            counts_by_qname[template.qname] += 1

                            renamed = self.renamed_outer_qname(
                                qname=template.qname,
                                side=side,
                                donor_bam=donor_bam,
                                occurrence_index=occurrence,
                            )

                            rename_by_qname_occurrence[(template.qname, occurrence)] = renamed
                            qout.write(
                                f"{donor_bam}\t{template.qname}\t{renamed}\t"
                                f"{side}\t{occurrence + 1}\n"
                            )

                        # For writing, if a qname was selected multiple times because of replacement,
                        # write that donor qname multiple times with different renamed qnames.
                        occurrence_total_by_qname = dict(counts_by_qname)

                        with pysam.AlignmentFile(donor_bam, "rb") as bam_in:
                            for read in bam_in.fetch(
                                self.interval.chrom,
                                max(0, self.interval.start0 - 1000),
                                self.interval.end0 + 1000,
                            ):
                                if read.query_name not in selected_qnames:
                                    continue

                                if read.mapping_quality < self.mapq:
                                    continue

                                n_occurrences = occurrence_total_by_qname[read.query_name]

                                # Only the first occurrence is copied unchanged.
                                # Additional replacement occurrences are written by
                                # write_jittered_outer_breakpoint_fastq() and realigned with BWA.
                                occurrence = 0

                                copied = pysam.AlignedSegment.fromstring(
                                    read.to_string(),
                                    bam_in.header,
                                )
                                copied.query_name = rename_by_qname_occurrence[
                                    (read.query_name, occurrence)
                                ]

                                bam_out.write(copied)
                                self.n_added_outer_records += 1

        sort_and_index_bam(
            input_bam=unsorted_bam,
            output_bam=self.added_outer_breakpoint_records_bam,
            threads=self.threads,
            remove_input=True,
        )

    def write_jittered_outer_breakpoint_fastq(
        self,
        selected_left: list[BoundaryReadTemplate],
        selected_right: list[BoundaryReadTemplate],
    ) -> None:
        """
        Write jittered replacement copies for outer-boundary reads.

        This only writes FASTQ records for donor qnames that were sampled more
        than once because --allow-replacement was used. The first occurrence is
        copied unchanged into added.outer.breakpoint.records.bam. Additional
        occurrences are regenerated from a nearby shifted reference position and
        later realigned with BWA.
        """
        selected = selected_left + selected_right

        selected_by_donor_side: dict[tuple[Path, str], list[BoundaryReadTemplate]] = defaultdict(list)
        for template in selected:
            selected_by_donor_side[(template.donor_bam, template.side)].append(template)

        with pysam.FastaFile(str(self.reference_fasta)) as reference:
            with open(self.added_outer_breakpoint_jittered_fastq, "w") as fastq_out:

                for (donor_bam, side), templates in selected_by_donor_side.items():
                    occurrence_total_by_qname: dict[str, int] = defaultdict(int)

                    for template in templates:
                        occurrence_total_by_qname[template.qname] += 1

                    # Only qnames sampled more than once need jittered replacement copies.
                    repeated_qnames = {
                        qname
                        for qname, n_occurrences in occurrence_total_by_qname.items()
                        if n_occurrences > 1
                    }

                    if not repeated_qnames:
                        continue

                    with pysam.AlignmentFile(donor_bam, "rb") as bam_in:
                        for read in bam_in.fetch(
                            self.interval.chrom,
                            max(0, self.interval.start0 - 1000),
                            self.interval.end0 + 1000,
                        ):
                            if read.query_name not in repeated_qnames:
                                continue

                            if not read_is_primary_usable(read, self.mapq):
                                continue

                            # Make sure we use reads from the relevant boundary side.
                            if side == "left" and not read_overlaps_pos(read, self.interval.start0):
                                continue

                            if side == "right" and not read_overlaps_pos(read, self.interval.end0 - 1):
                                continue

                            n_occurrences = occurrence_total_by_qname[read.query_name]

                            # occurrence 0 is already copied unchanged into
                            # added.outer.breakpoint.records.bam.
                            # occurrences 1+ are jittered and realigned.
                            for occurrence in range(1, n_occurrences):
                                renamed = self.renamed_outer_qname(
                                    qname=read.query_name,
                                    side=side,
                                    donor_bam=donor_bam,
                                    occurrence_index=occurrence,
                                )

                                shift = self.rng.randint(1, self.jitter_bp)
                                if self.rng.random() < 0.5:
                                    shift *= -1

                                shifted_start = max(0, read.reference_start + shift)

                                seq = reference_sequence_for_shifted_read(
                                    read=read,
                                    reference=reference,
                                    shifted_start=shifted_start,
                                )

                                qual = qualities_to_fastq_string(read, len(seq))

                                fastq_out.write(
                                    f"@{renamed}\n{seq}\n+\n{qual}\n"
                                )

                                self.n_jittered_outer_fastq_records += 1

    def bases_before_boundary_on_read(
        self,
        read: pysam.AlignedSegment,
        boundary_pos0: int,
    ) -> int:
        """
        Estimate how many query bases lie before the boundary.

        This uses aligned pairs and therefore respects CIGAR-level artifacts better than
        simply subtracting reference positions.
        """
        before = 0

        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False):
            if query_pos is None:
                continue

            if ref_pos is None:
                # insertion/soft-clipped query base: assign it to the side based
                # on current count; this is approximate but preserves artifacts.
                continue

            if ref_pos < boundary_pos0:
                before += 1

        return before

    def synthetic_junction_sequence_for_read(
        self,
        read: pysam.AlignedSegment,
        side: str,
        reference: pysam.FastaFile,
    ) -> str:
        """
        Build a synthetic read sequence spanning interval_end -> interval_start.

        The template read's CIGAR-derived breakpoint position determines where
        the synthetic sequence crosses the tandem junction.
        """
        query_len = read.query_length
        if query_len is None or query_len <= 0:
            query_len = len(read.query_sequence)

        if query_len is None or query_len <= 0:
            raise ValueError(f"Read has no query length: {read.query_name}")

        flank_len = max(query_len + 50, 250)

        right_context = reference.fetch(
            self.interval.chrom,
            max(0, self.interval.end0 - flank_len),
            self.interval.end0,
        ).upper()

        left_context = reference.fetch(
            self.interval.chrom,
            self.interval.start0,
            self.interval.start0 + flank_len,
        ).upper()

        junction_ref = right_context + left_context
        junction_offset = len(right_context)

        if side == "left":
            # original read crosses left flank -> interval start.
            # Bases before interval.start become sequence before junction,
            # now interval.end-side sequence.
            bases_before = self.bases_before_boundary_on_read(
                read,
                self.interval.start0,
            )
        elif side == "right":
            # original read crosses interval end -> right flank.
            # Bases before interval.end remain sequence before junction.
            bases_before = self.bases_before_boundary_on_read(
                read,
                self.interval.end0,
            )
        else:
            raise ValueError(f"Unexpected side: {side}")

        bases_before = max(1, min(query_len - 1, bases_before))
        start = junction_offset - bases_before
        end = start + query_len

        if start < 0 or end > len(junction_ref):
            raise ValueError(
                f"Could not build junction sequence for {read.query_name}: "
                f"start={start}, end={end}, junction_len={len(junction_ref)}"
            )

        seq = junction_ref[start:end].upper()

        # Preserve N positions from the original template read.
        original_seq = read.query_sequence or ""
        if len(original_seq) == len(seq):
            seq_list = list(seq)
            for i, base in enumerate(original_seq.upper()):
                if base == "N":
                    seq_list[i] = "N"
            seq = "".join(seq_list)

        if read.is_reverse:
            seq = reverse_complement(seq)

        return seq

    def write_tandem_junction_fastq(
        self,
        selected_left: list[BoundaryReadTemplate],
        selected_right: list[BoundaryReadTemplate],
    ) -> None:
        """
        Write synthetic reads representing the tandem duplication junction:
            interval_end -> interval_start
        """
        selected = selected_left + selected_right

        selected_by_donor_side: dict[tuple[Path, str], list[BoundaryReadTemplate]] = defaultdict(list)
        for template in selected:
            selected_by_donor_side[(template.donor_bam, template.side)].append(template)

        with pysam.FastaFile(str(self.reference_fasta)) as reference:
            with open(self.tandem_junction_fastq, "w") as fastq_out:
                with open(self.selected_junction_qnames_tsv, "w") as qout:
                    qout.write(
                        "donor_bam\toriginal_qname\trenamed_qname\tside\toccurrence_index\n"
                    )

                    for (donor_bam, side), templates in selected_by_donor_side.items():
                        counts_by_qname: dict[str, int] = defaultdict(int)
                        selected_qnames = {template.qname for template in templates}

                        occurrence_total_by_qname: dict[str, int] = defaultdict(int)
                        for template in templates:
                            occurrence_total_by_qname[template.qname] += 1

                        with pysam.AlignmentFile(donor_bam, "rb") as bam_in:
                            for read in bam_in.fetch(
                                self.interval.chrom,
                                max(0, self.interval.start0 - 1000),
                                self.interval.end0 + 1000,
                            ):
                                if read.query_name not in selected_qnames:
                                    continue

                                if not read_is_primary_usable(read, self.mapq):
                                    continue

                                # Make sure we use a read that actually overlaps the relevant edge.
                                if side == "left" and not read_overlaps_pos(read, self.interval.start0):
                                    continue
                                if side == "right" and not read_overlaps_pos(read, self.interval.end0 - 1):
                                    continue

                                n_occurrences = occurrence_total_by_qname[read.query_name]

                                for occurrence in range(n_occurrences):
                                    renamed = self.renamed_junction_qname(
                                        qname=read.query_name,
                                        side=side,
                                        occurrence_index=occurrence,
                                    )

                                    seq = self.synthetic_junction_sequence_for_read(
                                        read=read,
                                        side=side,
                                        reference=reference,
                                    )

                                    qual = qualities_to_fastq_string(read, len(seq))

                                    fastq_out.write(
                                        f"@{renamed}\n{seq}\n+\n{qual}\n"
                                    )

                                    qout.write(
                                        f"{donor_bam}\t{read.query_name}\t{renamed}\t"
                                        f"{side}\t{occurrence + 1}\n"
                                    )

                                    self.n_tandem_junction_fastq_records += 1


    def align_tandem_junction_fastq(self) -> None:
        if not fastq_has_records(self.tandem_junction_fastq):
            print("No tandem junction FASTQ records generated; skipping BWA alignment.")
            return

        align_fastq_with_bwa(
            fastq=self.tandem_junction_fastq,
            output_bam=self.tandem_junction_bam,
            reference=self.reference_fasta,
            threads=self.threads,
            interleaved=False,
            bwa_args=self.bwa_args,
        )

    def align_jittered_outer_breakpoint_fastq(self) -> None:
        if not fastq_has_records(self.added_outer_breakpoint_jittered_fastq):
            print("No jittered outer breakpoint FASTQ records generated; skipping BWA alignment.")
            return

        align_fastq_with_bwa(
            fastq=self.added_outer_breakpoint_jittered_fastq,
            output_bam=self.added_outer_breakpoint_jittered_bam,
            reference=self.reference_fasta,
            threads=self.threads,
            interleaved=False,
            bwa_args=self.bwa_args,
        )

    def run(self) -> DuplicationBreakpointResult:
        self.classify_recipient_boundary_qnames()

        n_left = self.requested_left_count()
        n_right = self.requested_right_count()

        print(f"Added fraction for breakpoint evidence: {self.added_fraction():.4f}")
        print(f"Requested left-boundary donor qnames: {n_left:,}")
        print(f"Requested right-boundary donor qnames: {n_right:,}")

        donor_templates = self.collect_donor_boundary_templates()

        selected_left_outer = self.sample_templates(
            templates=donor_templates["left"],
            n_requested=n_left,
            label="left outer-boundary",
        )
        selected_right_outer = self.sample_templates(
            templates=donor_templates["right"],
            n_requested=n_right,
            label="right outer-boundary",
        )

        # Independently sample junction templates. If no replacement is used and
        # donor pool is large enough, try to avoid reusing the same exact first block.
        self.rng.shuffle(donor_templates["left"])
        self.rng.shuffle(donor_templates["right"])

        selected_left_junction = self.sample_templates(
            templates=donor_templates["left"],
            n_requested=n_left,
            label="left tandem-junction",
        )
        selected_right_junction = self.sample_templates(
            templates=donor_templates["right"],
            n_requested=n_right,
            label="right tandem-junction",
        )

        print(f"Selected left outer-boundary qnames: {len(selected_left_outer):,}")
        print(f"Selected right outer-boundary qnames: {len(selected_right_outer):,}")
        print(f"Selected left junction templates: {len(selected_left_junction):,}")
        print(f"Selected right junction templates: {len(selected_right_junction):,}")

        self.write_added_outer_breakpoint_records_bam(
            selected_left=selected_left_outer,
            selected_right=selected_right_outer,
        )

        self.write_jittered_outer_breakpoint_fastq(
            selected_left=selected_left_outer,
            selected_right=selected_right_outer,
        )

        self.align_jittered_outer_breakpoint_fastq()

        self.write_tandem_junction_fastq(
            selected_left=selected_left_junction,
            selected_right=selected_right_junction,
        )

        self.align_tandem_junction_fastq()

        result = DuplicationBreakpointResult(
            added_outer_breakpoint_records_bam=self.added_outer_breakpoint_records_bam,
            added_outer_breakpoint_jittered_bam=self.added_outer_breakpoint_jittered_bam,
            tandem_junction_fastq=self.tandem_junction_fastq,
            tandem_junction_bam=self.tandem_junction_bam,
            candidate_tsv=self.candidate_tsv,
            selected_outer_qnames_tsv=self.selected_outer_qnames_tsv,
            selected_junction_qnames_tsv=self.selected_junction_qnames_tsv,
            n_recipient_left_boundary_qnames=len(self.recipient_left_boundary_qnames),
            n_recipient_right_boundary_qnames=len(self.recipient_right_boundary_qnames),
            n_requested_left_outer_qnames=n_left,
            n_requested_right_outer_qnames=n_right,
            n_selected_left_outer_qnames=len(selected_left_outer),
            n_selected_right_outer_qnames=len(selected_right_outer),
            n_added_outer_records=self.n_added_outer_records,
            n_requested_left_junction_templates=n_left,
            n_requested_right_junction_templates=n_right,
            n_selected_left_junction_templates=len(selected_left_junction),
            n_selected_right_junction_templates=len(selected_right_junction),
            n_tandem_junction_fastq_records=self.n_tandem_junction_fastq_records,
        )

        self.print_summary(result)
        return result

    def print_summary(self, result: DuplicationBreakpointResult) -> None:
        print("Duplication breakpoint editing complete")
        print(f"Added outer breakpoint records BAM: {result.added_outer_breakpoint_records_bam}")
        print(f"Tandem junction FASTQ: {result.tandem_junction_fastq}")
        print(f"Tandem junction BAM: {result.tandem_junction_bam}")
        print(f"Recipient left-boundary qnames: {result.n_recipient_left_boundary_qnames:,}")
        print(f"Recipient right-boundary qnames: {result.n_recipient_right_boundary_qnames:,}")
        print(f"Selected left outer-boundary qnames: {result.n_selected_left_outer_qnames:,}")
        print(f"Selected right outer-boundary qnames: {result.n_selected_right_outer_qnames:,}")
        print(f"Added outer breakpoint records: {result.n_added_outer_records:,}")
        print(f"Selected left junction templates: {result.n_selected_left_junction_templates:,}")
        print(f"Selected right junction templates: {result.n_selected_right_junction_templates:,}")
        print(f"Tandem junction FASTQ records: {result.n_tandem_junction_fastq_records:,}")
        print(f"Jittered outer breakpoint BAM: {result.added_outer_breakpoint_jittered_bam}")
        print(f"Jittered outer breakpoint FASTQ records: {self.n_jittered_outer_fastq_records:,}")


def merge_duplication_patch_with_breakpoints(
    edited_patch_bam: str | Path,
    added_outer_breakpoint_records_bam: str | Path,
    tandem_junction_bam: str | Path,
    output_prefix: str | Path,
    threads: int = 1,
    added_outer_breakpoint_jittered_bam: str | Path | None = None,
) -> Path:
    """
    Merge duplication edited patch with outer-boundary and tandem-junction BAMs.
    """
    prefix = Path(output_prefix)

    final_patch_bam = prefix.with_name(f"{prefix.name}.final.patch.bam")
    unsorted_bam = prefix.with_name(f"{prefix.name}.final.patch.unsorted.tmp.bam")

    merge_inputs = [
        str(edited_patch_bam),
        str(added_outer_breakpoint_records_bam),
    ]

    if added_outer_breakpoint_jittered_bam is not None and Path(added_outer_breakpoint_jittered_bam).exists():
        merge_inputs.append(str(added_outer_breakpoint_jittered_bam))

    if Path(tandem_junction_bam).exists():
        merge_inputs.append(str(tandem_junction_bam))

    cmd = [
        "samtools",
        "merge",
        "-@",
        str(threads),
        "-f",
        str(unsorted_bam),
        *merge_inputs,
    ]

    print("Merging duplication edited patch with breakpoint evidence:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)

    sort_and_index_bam(
        input_bam=unsorted_bam,
        output_bam=final_patch_bam,
        threads=threads,
        remove_input=True,
    )

    return final_patch_bam