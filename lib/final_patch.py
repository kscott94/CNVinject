#!/usr/bin/env python3

from dataclasses import dataclass
from pathlib import Path
import csv
import tempfile
import subprocess

import pysam

from patch import PatchPaths
from helpers import sort_and_index_bam

@dataclass
class FinalPatchResult:
    final_patch_bam: Path
    n_breakpoint_qnames: int
    n_records_removed_from_patch: int
    n_records_kept_from_patch: int


class FinalPatchBuilder:
    """
    Build final edited patch BAM.

    This removes original breakpoint-candidate qnames from the edited patch BAM,
    then merges in newly aligned synthetic breakpoint reads.

    Permanent output:
        PREFIX.final.patch.bam
        PREFIX.final.patch.bam.bai

    Temporary BAMs are created in a temp directory and deleted automatically.
    """

    def __init__(
        self,
        edited_patch_bam: str | Path,
        breakpoint_candidates_tsv: str | Path,
        output_prefix: str | Path,
        synthetic_bams: list[str | Path],
        threads: int = 1,
    ):
        self.edited_patch_bam = Path(edited_patch_bam)
        self.breakpoint_candidates_tsv = Path(breakpoint_candidates_tsv)
        self.paths = PatchPaths(output_prefix)
        self.synthetic_bams = [Path(bam) for bam in synthetic_bams]
        self.threads = threads

        self.final_patch_bam = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.final.patch.bam"
        )

        self.n_records_removed_from_patch = 0
        self.n_records_kept_from_patch = 0

    def load_breakpoint_qnames(self) -> set[str]:
        """
        Load qnames from breakpoint.candidates.tsv where is_breakpoint_candidate is True.
        """
        qnames: set[str] = set()

        with open(self.breakpoint_candidates_tsv) as handle:
            reader = csv.DictReader(handle, delimiter="\t")

            for row in reader:
                if row.get("is_breakpoint_candidate") == "True":
                    qnames.add(row["qname"])

        return qnames

    def write_patch_without_breakpoint_reads(
        self,
        breakpoint_qnames: set[str],
        output_bam: Path,
    ) -> None:
        """
        Remove original breakpoint qnames from edited patch BAM.
        """
        with pysam.AlignmentFile(self.edited_patch_bam, "rb") as bam_in:
            with pysam.AlignmentFile(output_bam, "wb", template=bam_in) as bam_out:
                for read in bam_in.fetch(until_eof=True):
                    if read.query_name in breakpoint_qnames:
                        self.n_records_removed_from_patch += 1
                        continue

                    bam_out.write(read)
                    self.n_records_kept_from_patch += 1

    def existing_nonempty_bams(self) -> list[Path]:
        """
        Return synthetic BAMs that exist and contain records.
        """
        usable_bams = []

        for bam in self.synthetic_bams:
            if not bam.exists():
                continue

            if bam.stat().st_size == 0:
                continue

            try:
                with pysam.AlignmentFile(bam, "rb") as handle:
                    has_record = False
                    for _ in handle.fetch(until_eof=True):
                        has_record = True
                        break

                if has_record:
                    usable_bams.append(bam)

            except ValueError:
                # Not a readable BAM.
                continue

        return usable_bams

    def merge_and_sort(
        self,
        cleaned_patch_bam: Path,
        synthetic_bams: list[Path],
        unsorted_output_bam: Path,
    ) -> None:
        """
        Merge cleaned patch BAM with synthetic breakpoint BAMs.
        """
        if synthetic_bams:
            merge_cmd = [
                "samtools",
                "merge",
                "-@",
                str(self.threads),
                "-f",
                str(unsorted_output_bam),
                str(cleaned_patch_bam),
                *[str(bam) for bam in synthetic_bams],
            ]

            print("Merging cleaned patch with synthetic breakpoint BAMs:")
            print(" ".join(merge_cmd))

            subprocess.run(merge_cmd, check=True)

        else:
            # No synthetic BAMs. The cleaned patch is the unsorted output.
            # Use pysam to copy it so downstream sort/index behavior is the same.
            with pysam.AlignmentFile(cleaned_patch_bam, "rb") as bam_in:
                with pysam.AlignmentFile(unsorted_output_bam, "wb", template=bam_in) as bam_out:
                    for read in bam_in.fetch(until_eof=True):
                        bam_out.write(read)

    def sort_and_index_final(self, unsorted_bam: Path) -> None:
        """
        Sort and index final patch BAM.
        """
        sort_and_index_bam(
            input_bam=unsorted_bam,
            output_bam=self.final_patch_bam,
            threads=self.threads,
            remove_input=False,
        )

    def run(self) -> FinalPatchResult:
        breakpoint_qnames = self.load_breakpoint_qnames()

        with tempfile.TemporaryDirectory(prefix="cnvinject_final_patch_") as tmpdir:
            tmpdir = Path(tmpdir)

            cleaned_patch_bam = tmpdir / "cleaned_patch.bam"
            unsorted_final_bam = tmpdir / "final_patch.unsorted.bam"

            self.write_patch_without_breakpoint_reads(
                breakpoint_qnames=breakpoint_qnames,
                output_bam=cleaned_patch_bam,
            )

            usable_synthetic_bams = self.existing_nonempty_bams()

            self.merge_and_sort(
                cleaned_patch_bam=cleaned_patch_bam,
                synthetic_bams=usable_synthetic_bams,
                unsorted_output_bam=unsorted_final_bam,
            )

            self.sort_and_index_final(unsorted_final_bam)

        result = FinalPatchResult(
            final_patch_bam=self.final_patch_bam,
            n_breakpoint_qnames=len(breakpoint_qnames),
            n_records_removed_from_patch=self.n_records_removed_from_patch,
            n_records_kept_from_patch=self.n_records_kept_from_patch,
        )

        self.print_summary(result)
        return result

    def print_summary(self, result: FinalPatchResult) -> None:
        print("Final patch construction complete")
        print(f"Edited patch BAM: {self.edited_patch_bam}")
        print(f"Breakpoint candidates TSV: {self.breakpoint_candidates_tsv}")
        print(f"Final patch BAM: {result.final_patch_bam}")
        print(f"Breakpoint reads removed from patch: {result.n_breakpoint_qnames:,}")
        print(f"Records removed from patch: {result.n_records_removed_from_patch:,}")
        print(f"Records kept from patch: {result.n_records_kept_from_patch:,}")