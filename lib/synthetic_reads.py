#!/usr/bin/env python3

from dataclasses import dataclass
from pathlib import Path
import csv

import pysam

from patch import GenomicInterval, PatchPaths
from helpers import (
    parse_cigar_tuples,
    reverse_complement,
    force_mismatch,
    apply_n_mask,
    make_fastq_record,
)


@dataclass
class SyntheticFastqResult:
    paired_fastq: Path
    singleton_fastq: Path
    n_candidate_rows: int
    n_paired_records: int
    n_singleton_records: int
    n_qnames_written: int


class SyntheticBreakpointReadGenerator:
    """
    Generate synthetic breakpoint FASTQs from breakpoint.candidates.tsv.

    This module does not align reads. It only writes FASTQ files.

    Outputs:
        PREFIX.breakpoint.paired.interleaved.fastq
        PREFIX.breakpoint.singletons.fastq
    """

    def __init__(
        self,
        edited_patch_bam: str | Path,
        candidates_tsv: str | Path,
        output_prefix: str | Path,
        interval: GenomicInterval,
        reference_fasta: str | Path,
    ):
        self.edited_patch_bam = Path(edited_patch_bam)
        self.candidates_tsv = Path(candidates_tsv)
        self.paths = PatchPaths(output_prefix)
        self.interval = interval
        self.reference_fasta = Path(reference_fasta)

        self.paired_fastq = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.breakpoint.paired.interleaved.fastq"
        )

        self.singleton_fastq = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.breakpoint.singletons.fastq"
        )

    def load_candidate_rows(self) -> list[dict[str, str]]:
        """
        Load candidate rows from breakpoint.candidates.tsv.

        Only rows with is_breakpoint_candidate == True are used.
        """
        rows = []

        with open(self.candidates_tsv) as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            for row in reader:
                if row.get("is_breakpoint_candidate") == "True":
                    rows.append(row)

        return rows

    def load_reads_by_qname(self, qnames: set[str]) -> dict[str, dict[str, pysam.AlignedSegment]]:
        """
        Pull primary read1/read2 alignments from the edited patch BAM.

        Returns:
            {
                qname: {
                    "read1": pysam.AlignedSegment,
                    "read2": pysam.AlignedSegment,
                }
            }
        """
        reads_by_qname: dict[str, dict[str, pysam.AlignedSegment]] = {}

        with pysam.AlignmentFile(self.edited_patch_bam, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.query_name not in qnames:
                    continue

                if read.is_unmapped:
                    continue

                if read.is_secondary or read.is_supplementary:
                    continue

                if read.query_name not in reads_by_qname:
                    reads_by_qname[read.query_name] = {}

                if read.is_read1:
                    reads_by_qname[read.query_name]["read1"] = read

                elif read.is_read2:
                    reads_by_qname[read.query_name]["read2"] = read

        return reads_by_qname

    def original_read_fastq_sequence(self, read: pysam.AlignedSegment) -> str:
        """
        Return the original read sequence in FASTQ orientation.

        For reverse-strand BAM records, pysam stores the aligned sequence.
        To recreate FASTQ-like orientation, use get_forward_sequence().
        """
        seq = read.get_forward_sequence()

        if seq is None:
            raise ValueError(f"Read {read.query_name} has no query sequence.")

        return seq.upper()

    def reference_span_sequence(self, fasta: pysam.FastaFile, read: pysam.AlignedSegment) -> str:
        """
        Fetch the original reference sequence over the read's alignment span.
        """
        if read.reference_name is None:
            raise ValueError(f"Read {read.query_name} has no reference name.")

        if read.reference_start is None or read.reference_end is None:
            raise ValueError(f"Read {read.query_name} has no reference span.")

        return fasta.fetch(
            read.reference_name,
            read.reference_start,
            read.reference_end,
        ).upper()

    def fetch_new_reference_sequence(
        self,
        fasta: pysam.FastaFile,
        read: pysam.AlignedSegment,
        action: str,
    ) -> str:
        """
        Fetch the new reference-derived sequence for this read.

        The returned sequence is in reference-forward orientation and should have
        length equal to the original reference-consuming length of the read.

        Actions:
            keep_original:
                fetch the original reference span

            shift_to_right_flank:
                move the read from inside the deleted interval to the right flank

            shift_to_left_flank:
                move the read from inside the deleted interval to the left flank

            junction:
                create sequence across the deletion junction
        """
        ref_len = read.reference_length

        if ref_len is None:
            raise ValueError(f"Read {read.query_name} has no reference_length.")

        if action == "keep_original":
            return self.reference_span_sequence(fasta, read)

        if action == "shift_to_right_flank":
            # Preserve offset from the left breakpoint.
            offset = read.reference_start - self.interval.start0
            new_start = self.interval.end0 + offset
            new_end = new_start + ref_len

            return fasta.fetch(self.interval.chrom, new_start, new_end).upper()

        if action == "shift_to_left_flank":
            chrom_len = fasta.get_reference_length(self.interval.chrom)

            # If this is a left-terminal deletion, there is no left flank.
            # Fall back to placing the read in the right flank.
            if self.interval.start0 == 0:
                offset = read.reference_start - self.interval.start0
                new_start = self.interval.end0 + max(0, offset)
                new_end = new_start + ref_len

                if new_end > chrom_len:
                    raise ValueError(
                        f"Read {read.query_name} shift_to_left_flank fallback "
                        f"extends beyond chromosome end: {self.interval.chrom}:{new_start}-{new_end}"
                    )

                return fasta.fetch(self.interval.chrom, new_start, new_end).upper()

            offset = self.interval.end0 - read.reference_end
            new_end = self.interval.start0 - offset
            new_start = new_end - ref_len

            if new_start < 0:
                raise ValueError(
                    f"Read {read.query_name} shift_to_left_flank produced "
                    f"negative coordinate: {new_start}"
                )

            return fasta.fetch(self.interval.chrom, new_start, new_end).upper()

        if action == "junction":
            return self.fetch_junction_reference_sequence(fasta, read)

        raise ValueError(f"Unknown synthetic read action: {action}")

    def fetch_junction_reference_sequence(
            self,
            fasta: pysam.FastaFile,
            read: pysam.AlignedSegment,
    ) -> str:
        """
        Generate a reference-forward sequence for reads that cross a deletion boundary.

        Internal deletion:
            left flank + right flank junction

        Left-terminal deletion, e.g. chr1:1-1003800:
            no left flank exists
            breakpoint-crossing reads are converted to right-flank-only sequence

        Right-terminal deletion:
            no right flank exists
            breakpoint-crossing reads are converted to left-flank-only sequence
        """
        ref_len = read.reference_length

        if ref_len is None:
            raise ValueError(f"Read {read.query_name} has no reference_length.")

        chrom = self.interval.chrom
        chrom_len = fasta.get_reference_length(chrom)

        is_left_terminal_deletion = self.interval.start0 == 0
        is_right_terminal_deletion = self.interval.end0 >= chrom_len

        crosses_left = (
                read.reference_name == chrom
                and read.reference_start < self.interval.start0
                and read.reference_end > self.interval.start0
        )

        crosses_right = (
                read.reference_name == chrom
                and read.reference_start < self.interval.end0
                and read.reference_end > self.interval.end0
        )

        # ------------------------------------------------------------
        # Left-terminal deletion:
        # deleted interval starts at chr coordinate 1.
        #
        # There is no left flank, so a read crossing the right boundary
        # should become a right-flank-only read.
        # ------------------------------------------------------------
        if is_left_terminal_deletion:
            new_start = self.interval.end0
            new_end = new_start + ref_len

            if new_end > chrom_len:
                raise ValueError(
                    f"Read {read.query_name} left-terminal deletion replacement "
                    f"extends beyond chromosome end: {chrom}:{new_start}-{new_end}, "
                    f"chrom length={chrom_len}"
                )

            return fasta.fetch(chrom, new_start, new_end).upper()

        # ------------------------------------------------------------
        # Right-terminal deletion:
        # deleted interval extends to the end of the chromosome.
        #
        # There is no right flank, so a read crossing the left boundary
        # should become a left-flank-only read.
        # ------------------------------------------------------------
        if is_right_terminal_deletion:
            new_end = self.interval.start0
            new_start = new_end - ref_len

            if new_start < 0:
                raise ValueError(
                    f"Read {read.query_name} right-terminal deletion replacement "
                    f"extends before chromosome start: {chrom}:{new_start}-{new_end}"
                )

            return fasta.fetch(chrom, new_start, new_end).upper()

        # ------------------------------------------------------------
        # Standard internal deletion.
        # ------------------------------------------------------------
        if crosses_left:
            left_len = self.interval.start0 - read.reference_start
            right_len = ref_len - left_len

            left_seq = fasta.fetch(
                chrom,
                read.reference_start,
                self.interval.start0,
            ).upper()

            right_seq = fasta.fetch(
                chrom,
                self.interval.end0,
                self.interval.end0 + right_len,
            ).upper()

            return left_seq + right_seq

        if crosses_right:
            right_len = read.reference_end - self.interval.end0
            left_len = ref_len - right_len

            left_start = self.interval.start0 - left_len

            if left_start < 0:
                raise ValueError(
                    f"Read {read.query_name} junction sequence produced "
                    f"negative coordinate: {left_start}"
                )

            left_seq = fasta.fetch(
                chrom,
                left_start,
                self.interval.start0,
            ).upper()

            right_seq = fasta.fetch(
                chrom,
                self.interval.end0,
                read.reference_end,
            ).upper()

            return left_seq + right_seq

        # If the read does not itself cross a breakpoint, use original span.
        return self.reference_span_sequence(fasta, read)

    def build_synthetic_sequence(
        self,
        read: pysam.AlignedSegment,
        original_ref_seq: str,
        new_ref_seq: str,
    ) -> str:
        """
        Build a synthetic read sequence by replaying the original CIGAR onto
        a new reference sequence.

        CIGAR policy:
            M:
                if original read matched original reference:
                    use new reference base
                else:
                    use force_mismatch(new reference base)

            =:
                use new reference base

            X:
                use force_mismatch(new reference base)

            I:
                preserve original inserted query bases

            D, N:
                advance reference pointer, add no bases

            S:
                preserve original soft-clipped query bases

            H, P:
                ignore

        Finally:
            preserve original N bases at the same query positions
            convert to FASTQ/original-read orientation if read is reverse
        """
        if read.cigarstring is None:
            raise ValueError(f"Read {read.query_name} has no CIGAR string.")

        query_seq = read.query_sequence

        if query_seq is None:
            raise ValueError(f"Read {read.query_name} has no query sequence.")

        query_seq = query_seq.upper()
        original_ref_seq = original_ref_seq.upper()
        new_ref_seq = new_ref_seq.upper()

        qpos = 0
        rpos = 0
        synthetic_parts = []

        for length, op in parse_cigar_tuples(read.cigarstring):

            if op == "M":
                for _ in range(length):
                    original_read_base = query_seq[qpos]
                    original_ref_base = original_ref_seq[rpos]
                    new_ref_base = new_ref_seq[rpos]

                    if original_read_base == "N":
                        synthetic_base = "N"
                    elif original_read_base == original_ref_base:
                        synthetic_base = new_ref_base
                    else:
                        synthetic_base = force_mismatch(new_ref_base)

                    synthetic_parts.append(synthetic_base)

                    qpos += 1
                    rpos += 1

            elif op == "=":
                synthetic_parts.append(new_ref_seq[rpos : rpos + length])
                qpos += length
                rpos += length

            elif op == "X":
                for i in range(length):
                    synthetic_parts.append(force_mismatch(new_ref_seq[rpos + i]))

                qpos += length
                rpos += length

            elif op == "I":
                # Preserve original inserted query bases.
                synthetic_parts.append(query_seq[qpos : qpos + length])
                qpos += length

            elif op == "D":
                # Preserve deletion pattern by skipping reference bases.
                rpos += length

            elif op == "N":
                # Treat skipped reference like deletion for sequence generation.
                rpos += length

            elif op == "S":
                # Preserve original soft-clipped query bases.
                synthetic_parts.append(query_seq[qpos : qpos + length])
                qpos += length

            elif op in {"H", "P"}:
                # No query or reference bases to add.
                continue

            else:
                raise ValueError(
                    f"Unsupported CIGAR operation {op} in read {read.query_name}"
                )

        synthetic_seq_aligned_orientation = "".join(synthetic_parts).upper()

        # Preserve original read-base Ns at the same query positions.
        synthetic_seq_aligned_orientation = apply_n_mask(
            synthetic_seq_aligned_orientation,
            query_seq,
        )

        # Convert from BAM/aligned orientation back to FASTQ/original-read orientation.
        if read.is_reverse:
            return reverse_complement(synthetic_seq_aligned_orientation)

        return synthetic_seq_aligned_orientation

    def choose_actions_for_pair(self, row: dict[str, str]) -> tuple[str | None, str | None]:
        """
        Decide read1/read2 actions from the breakpoint candidate row.

        Returns:
            (read1_action, read2_action)

        Actions:
            keep_original
            shift_to_right_flank
            shift_to_left_flank
            junction
            None
        """
        r1_class = row["read1_class"]
        r2_class = row["read2_class"]

        r1_action = None
        r2_action = None

        # Breakpoint-crossing reads get junction sequence.
        if r1_class in {"crosses_left_breakpoint", "crosses_right_breakpoint", "partial_overlap"}:
            r1_action = "junction"

        if r2_class in {"crosses_left_breakpoint", "crosses_right_breakpoint", "partial_overlap"}:
            r2_action = "junction"

        # Internal mate paired with left flank should move to right flank.
        if r1_class == "inside" and r2_class == "left_flank":
            r1_action = "shift_to_right_flank"
            r2_action = "keep_original"

        if r2_class == "inside" and r1_class == "left_flank":
            r2_action = "shift_to_right_flank"
            r1_action = "keep_original"

        # Internal mate paired with right flank should move to left flank.
        if r1_class == "inside" and r2_class == "right_flank":
            r1_action = "shift_to_left_flank"
            r2_action = "keep_original"

        if r2_class == "inside" and r1_class == "right_flank":
            r2_action = "shift_to_left_flank"
            r1_action = "keep_original"

        # Pair spanning the deletion: both reads can be realigned unchanged.
        if r1_class == "left_flank" and r2_class == "right_flank":
            r1_action = "keep_original"
            r2_action = "keep_original"

        if r2_class == "left_flank" and r1_class == "right_flank":
            r1_action = "keep_original"
            r2_action = "keep_original"

        # If one read crosses a breakpoint and the other is flanking, keep the flanking mate.
        if r1_action == "junction" and r2_class in {"left_flank", "right_flank"}:
            r2_action = "keep_original"

        if r2_action == "junction" and r1_class in {"left_flank", "right_flank"}:
            r1_action = "keep_original"

        return r1_action, r2_action

    def make_record_for_read(
        self,
        read: pysam.AlignedSegment,
        action: str,
        fasta: pysam.FastaFile,
        mate_suffix: str,
    ) -> str:
        """
        Make one FASTQ record for a read using the selected action.
        """
        if action == "keep_original":
            seq = self.original_read_fastq_sequence(read)
        else:
            original_ref_seq = self.reference_span_sequence(fasta, read)
            new_ref_seq = self.fetch_new_reference_sequence(
                fasta=fasta,
                read=read,
                action=action,
            )

            seq = self.build_synthetic_sequence(
                read=read,
                original_ref_seq=original_ref_seq,
                new_ref_seq=new_ref_seq,
            )

        fastq_name = f"{read.query_name}/{mate_suffix}"

        return make_fastq_record(fastq_name, seq, phred_char="I")

    def run(self) -> SyntheticFastqResult:
        rows = self.load_candidate_rows()
        qnames = {row["qname"] for row in rows}
        reads_by_qname = self.load_reads_by_qname(qnames)

        n_paired_records = 0
        n_singleton_records = 0
        qnames_written: set[str] = set()

        with pysam.FastaFile(self.reference_fasta) as fasta:
            with open(self.paired_fastq, "w") as paired_handle, open(
                self.singleton_fastq, "w"
            ) as singleton_handle:

                for row in rows:
                    qname = row["qname"]

                    if qname not in reads_by_qname:
                        continue

                    read1 = reads_by_qname[qname].get("read1")
                    read2 = reads_by_qname[qname].get("read2")

                    read1_action, read2_action = self.choose_actions_for_pair(row)

                    has_read1_record = read1 is not None and read1_action is not None
                    has_read2_record = read2 is not None and read2_action is not None

                    # Complete pair: write interleaved read1/read2.
                    if has_read1_record and has_read2_record:
                        paired_handle.write(
                            self.make_record_for_read(
                                read=read1,
                                action=read1_action,
                                fasta=fasta,
                                mate_suffix="1",
                            )
                        )
                        paired_handle.write(
                            self.make_record_for_read(
                                read=read2,
                                action=read2_action,
                                fasta=fasta,
                                mate_suffix="2",
                            )
                        )

                        n_paired_records += 2
                        qnames_written.add(qname)
                        continue

                    # Singleton read1.
                    if has_read1_record:
                        singleton_handle.write(
                            self.make_record_for_read(
                                read=read1,
                                action=read1_action,
                                fasta=fasta,
                                mate_suffix="1",
                            )
                        )

                        n_singleton_records += 1
                        qnames_written.add(qname)
                        continue

                    # Singleton read2.
                    if has_read2_record:
                        singleton_handle.write(
                            self.make_record_for_read(
                                read=read2,
                                action=read2_action,
                                fasta=fasta,
                                mate_suffix="2",
                            )
                        )

                        n_singleton_records += 1
                        qnames_written.add(qname)
                        continue

        result = SyntheticFastqResult(
            paired_fastq=self.paired_fastq,
            singleton_fastq=self.singleton_fastq,
            n_candidate_rows=len(rows),
            n_paired_records=n_paired_records,
            n_singleton_records=n_singleton_records,
            n_qnames_written=len(qnames_written),
        )

        self.print_summary(result)
        return result

    def print_summary(self, result: SyntheticFastqResult) -> None:
        print("Synthetic breakpoint FASTQ generation complete")
        print(f"Edited patch BAM: {self.edited_patch_bam}")
        print(f"Candidates TSV: {self.candidates_tsv}")
        print(f"Reference FASTA: {self.reference_fasta}")
        print(f"Paired FASTQ: {result.paired_fastq}")
        print(f"Singleton FASTQ: {result.singleton_fastq}")
        print(f"Candidate reads: {result.n_candidate_rows:,}")
        print(f"Paired FASTQ records: {result.n_paired_records:,}")
        print(f"Singleton FASTQ records: {result.n_singleton_records:,}")
        print(f"Qnames written: {result.n_qnames_written:,}")