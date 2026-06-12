#!/usr/bin/env python3

from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict

import pysam
from read_pairs import ReadPairStatus
from patch import GenomicInterval, PatchPaths
from helpers import (sort_and_index_bam,
                     deletion_fraction,
                     sample_qnames_file_for_deletion,
                     load_qnames)


@dataclass
class DeletionResult:
    """
    Summary of deletion patch editing.
    """

    input_patch_bam: Path
    edited_patch_bam: Path
    internal_qnames: Path
    deleted_internal_qnames: Path
    n_internal_qnames: int
    n_deleted_internal_qnames: int
    n_input_records: int
    n_removed_records: int
    n_kept_records: int


class DeletionEditor:
    """
    Edit a raw patch BAM for a deletion.

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
        copy_number: float = 0.0,
        seed: int = 1,
        all_alignments: bool = False,
        threads: int = 1,
    ):
        self.input_patch_bam = Path(input_patch_bam)
        self.paths = PatchPaths(output_prefix)
        self.interval = interval
        self.copy_number = copy_number
        self.seed = seed
        self.all_alignments = all_alignments
        self.threads = threads

        self.edited_patch_bam = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.edited.patch.bam"
        )

        self.internal_qnames = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.internal.qnames.txt"
        )

        self.deleted_internal_qnames = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.deleted.internal.qnames"
        )

        self.status_by_qname: dict[str, ReadPairStatus] = defaultdict(ReadPairStatus)
        self.internal_qnames_all: set[str] = set()
        self.qnames_to_remove: set[str] = set()


        self.n_input_records = 0
        self.n_removed_records = 0
        self.n_deleted_internal_qnames = 0
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

        self.internal_qnames_all = {
            qname
            for qname, status in self.status_by_qname.items()
            if status.is_fully_inside_interval(self.interval)
        }

    def write_internal_qnames(self) -> None:
        """
        Write all fully internal qnames before copy-number sampling.
        """
        with open(self.internal_qnames, "w") as handle:
            for qname in sorted(self.internal_qnames_all):
                handle.write(qname + "\n")


    def write_deleted_internal_qnames(self) -> None:
        """
        Randomly sample internal qnames according to requested copy number.
        """
        n_selected = sample_qnames_file_for_deletion(
            input_qnames=self.internal_qnames,
            output_qnames=self.deleted_internal_qnames,
            copy_number=float(self.copy_number),
            seed=int(self.seed),
        )

        self.qnames_to_remove = load_qnames(self.deleted_internal_qnames)
        self.n_deleted_internal_qnames = n_selected



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

    def run(self) -> DeletionResult:
        """
        Run deletion patch editing.
        """
        self.classify_patch_reads()
        self.write_internal_qnames()
        self.write_deleted_internal_qnames()
        self.write_edited_patch_bam()
        self.sort_and_index_edited_patch()

        result = DeletionResult(
            input_patch_bam=self.input_patch_bam,
            edited_patch_bam=self.edited_patch_bam,
            internal_qnames=self.internal_qnames,
            deleted_internal_qnames=self.deleted_internal_qnames,
            n_internal_qnames=len(self.internal_qnames_all),
            n_deleted_internal_qnames=len(self.qnames_to_remove),
            n_input_records=self.n_input_records,
            n_removed_records=self.n_removed_records,
            n_kept_records=self.n_kept_records,
        )

        self.print_summary(result)
        return result

    def print_summary(self, result: DeletionResult) -> None:
        print("Deletion patch editing complete")
        print(f"Copy number: {self.copy_number}")
        print(f"Deletion fraction: {deletion_fraction(self.copy_number):.4f}")
        print(f"Deletion interval: {self.interval}")
        print(f"Input patch BAM: {result.input_patch_bam}")
        print(f"Edited patch BAM: {result.edited_patch_bam}")
        print(f"Input patch records: {result.n_input_records:,}")
        print(f"All internal qnames: {result.n_internal_qnames:,}")
        print(f"Deleted internal qnames sampled: {result.n_deleted_internal_qnames:,}")
        print(f"Removed records: {result.n_removed_records:,}")
        print(f"Kept records: {result.n_kept_records:,}")