#!/usr/bin/env python3

from dataclasses import dataclass
from pathlib import Path
import tempfile
import subprocess
import pysam
from patch import PatchPaths
from helpers import sort_and_index_bam


@dataclass
class MergePatchResult:
    full_bam: Path
    n_patch_qnames: int
    n_records_removed_from_full_bam: int
    n_records_kept_from_full_bam: int


class MergePatchBuilder:
    """
    Merge a final edited patch BAM back into the original full BAM.

    This removes all original patch qnames from the full BAM, then merges in
    the final edited patch BAM.

    Permanent output:
        PREFIX.full.bam
        PREFIX.full.bam.bai

    Temporary intermediate BAMs are deleted automatically.
    """

    def __init__(
        self,
        full_bam: str | Path,
        final_patch_bam: str | Path,
        patch_qnames: str | Path,
        output_bam: str | Path,
        threads: int = 1,
    ):
        self.full_bam = Path(full_bam)
        self.final_patch_bam = Path(final_patch_bam)
        self.patch_qnames = Path(patch_qnames)
        self.output_full_bam = Path(output_bam)
        self.threads = threads

        self.n_records_removed_from_full_bam = -1
        self.n_records_kept_from_full_bam = -1

    def load_patch_qnames(self) -> set[str]:
        """
        Load all original patch qnames.

        These are the qnames originally extracted from the full BAM.
        They must be removed from the full BAM before the final edited patch
        is merged back in.
        """
        qnames: set[str] = set()

        with open(self.patch_qnames) as handle:
            for line in handle:
                qname = line.strip()
                if qname:
                    qnames.add(qname)

        return qnames

    def write_full_bam_without_patch_reads(
            self,
            patch_qnames: set[str],
            output_bam: Path,
    ) -> None:
        """
        Remove all original patch qnames from the full BAM using samtools.

        This is much faster than scanning the full BAM in Python.

        samtools logic:
            -N patch_qnames.txt selects reads with qnames in the file
            -U output_bam writes reads NOT selected by -N
            -o /dev/null discards the selected patch reads
        """
        cmd = [
            "samtools",
            "view",
            "-@",
            str(self.threads),
            "-b",
            "-N",
            str(self.patch_qnames),
            "-U",
            str(output_bam),
            "-o",
            "/dev/null",
            str(self.full_bam),
        ]

        print("Removing patch reads from input BAM...")
        print(" ".join(cmd))

        subprocess.run(cmd, check=True)

        # We no longer count records removed/kept here because that would require
        # another full-BAM scan. Keep qname count as the main summary metric.
        self.n_records_removed_from_full_bam = -1
        self.n_records_kept_from_full_bam = -1

    def merge_clean_full_with_patch(
        self,
        clean_full_bam: Path,
        final_patch_bam: Path,
        unsorted_output_bam: Path,
    ) -> None:
        """
        Merge the clean full BAM and final patch BAM.
        """
        merge_cmd = [
            "samtools",
            "merge",
            "-@",
            str(self.threads),
            "-f",
            str(unsorted_output_bam),
            str(clean_full_bam),
            str(final_patch_bam),
        ]

        print("Generating final BAM...")
        print(" ".join(merge_cmd))

        subprocess.run(merge_cmd, check=True)

    def run(self) -> MergePatchResult:
        if not self.full_bam.exists():
            raise FileNotFoundError(f"Full BAM does not exist: {self.full_bam}")

        if not self.final_patch_bam.exists():
            raise FileNotFoundError(f"Final patch BAM does not exist: {self.final_patch_bam}")

        if not self.patch_qnames.exists():
            raise FileNotFoundError(f"Patch qnames file does not exist: {self.patch_qnames}")

        patch_qnames = self.load_patch_qnames()

        tmp_parent = self.output_full_bam.parent
        tmp_parent.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory(
                prefix="temp_",
                dir=tmp_parent,
        ) as tmpdir:

            tmpdir = Path(tmpdir)

            clean_full_bam = tmpdir / "full_without_original_patch_reads.bam"
            unsorted_full_bam = tmpdir / "merged_full.unsorted.bam"

            self.write_full_bam_without_patch_reads(
                patch_qnames=patch_qnames,
                output_bam=clean_full_bam,
            )

            self.merge_clean_full_with_patch(
                clean_full_bam=clean_full_bam,
                final_patch_bam=self.final_patch_bam,
                unsorted_output_bam=unsorted_full_bam,
            )

            sort_and_index_bam(
                input_bam=unsorted_full_bam,
                output_bam=self.output_full_bam,
                threads=self.threads,
                remove_input=False,
            )

        result = MergePatchResult(
            full_bam=self.output_full_bam,
            n_patch_qnames=len(patch_qnames),
            n_records_removed_from_full_bam=self.n_records_removed_from_full_bam,
            n_records_kept_from_full_bam=self.n_records_kept_from_full_bam,
        )

        self.print_summary(result)
        return result

    def print_summary(self, result: MergePatchResult) -> None:
        print("Full BAM reconstruction complete")
        print(f"Original full BAM: {self.full_bam}")
        print(f"Final patch BAM: {self.final_patch_bam}")
        print(f"Patch qnames: {self.patch_qnames}")
        print(f"Output full BAM: {result.full_bam}")
        print(f"Original patch qnames removed: {result.n_patch_qnames:,}")
        print(f"Records removed from full BAM: {result.n_records_removed_from_full_bam:,}")
        print(f"Records kept from full BAM: {result.n_records_kept_from_full_bam:,}")