#!/usr/bin/env python3

from dataclasses import dataclass
from pathlib import Path
from typing import Set
import pysam
import subprocess

from helpers import index_bam


@dataclass(frozen=True)
class GenomicInterval:
    """
    Genomic interval using 1-based coordinates.

    Example:
        chr17:30780079-31936302

    pysam uses 0-based half-open coordinates internally, so this class
    provides start0 and end0 properties for pysam fetch().
    """

    chrom: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 1:
            raise ValueError("Interval start must be >= 1.")
        if self.end < self.start:
            raise ValueError("Interval end must be >= start.")

    @property
    def start0(self) -> int:
        return self.start - 1

    @property
    def end0(self) -> int:
        return self.end

    @classmethod
    def from_string(cls, interval_string: str) -> "GenomicInterval":
        """
        Parse interval string like:
            chr17:30780079-31936302
        """
        if ":" not in interval_string or "-" not in interval_string:
            raise ValueError(
                f"Interval must be in chr:start-end format. Got: {interval_string}"
            )

        chrom, coordinates = interval_string.split(":", 1)
        start_str, end_str = coordinates.split("-", 1)

        start = int(start_str.replace(",", ""))
        end = int(end_str.replace(",", ""))

        return cls(chrom=chrom, start=start, end=end)

    def with_buffer(self, buffer_size: int) -> "GenomicInterval":
        """
        Return a new interval expanded by buffer_size.
        """
        return GenomicInterval(
            chrom=self.chrom,
            start=max(1, self.start - buffer_size),
            end=self.end + buffer_size,
        )

    def to_region_string(self) -> str:
        """
        Return samtools-style region string.
        """
        return f"{self.chrom}:{self.start}-{self.end}"

    def __str__(self) -> str:
        return self.to_region_string()


@dataclass(frozen=True)
class PatchPaths:
    """
    Standard patch output paths generated from an output prefix.

    Example prefix:
        Sample2.chr17_30780079_31936302.CN0

    Generates:
        Sample2.chr17_30780079_31936302.CN0.patch.bam
        Sample2.chr17_30780079_31936302.CN0.patch.bam.bai
        Sample2.chr17_30780079_31936302.CN0.patch.qnames.txt
    """

    output_prefix: str | Path

    @property
    def prefix(self) -> Path:
        return Path(self.output_prefix)

    @property
    def patch_bam(self) -> Path:
        return self.prefix.with_name(f"{self.prefix.name}.patch.bam")

    @property
    def patch_bai(self) -> Path:
        return self.prefix.with_name(f"{self.prefix.name}.patch.bam.bai")

    @property
    def patch_qnames(self) -> Path:
        return self.prefix.with_name(f"{self.prefix.name}.patch.qnames.txt")


@dataclass
class PatchDissectionResult:
    """
    Summary of patch dissection.
    """

    patch_bam: Path
    patch_qnames: Path
    n_patch_qnames: int
    n_patch_records: int


class PatchDissector:
    """
    Extract a buffered patch from a full BAM.

    This class prepares the shared patch used by deletion and duplication workflows.

    It does two passes:

    Pass 1:
        Fetch reads overlapping the buffered patch interval.
        Collect their qnames.

    Pass 2:
        Scan the full BAM and write every alignment record whose qname was found
        in pass 1.

    Why two passes?
        pysam.fetch(region) only gets records overlapping the region.
        It does not automatically retrieve mate records outside the region.

        By collecting qnames first and then scanning the full BAM, this class can
        include mate records for reads whose partner overlaps the patch.
    """

    def __init__(
        self,
        input_bam: str | Path,
        output_prefix: str | Path,
        interval: GenomicInterval,
        buffer_size: int = 10000,
        mapq: int = 0,
        threads: int = 1,
    ):
        self.input_bam = Path(input_bam)
        self.paths = PatchPaths(output_prefix)
        self.interval = interval
        self.buffer_size = buffer_size
        self.mapq = mapq
        self.threads = threads

        self.patch_interval = self.interval.with_buffer(self.buffer_size)

        # Make sure the output directory exists.
        self.paths.prefix.parent.mkdir(parents=True, exist_ok=True)

    def extract_patch_bam_with_samtools(self) -> None:
        """
        Extract buffered patch reads and their mates using samtools view -P.

        This replaces the slow Python two-pass method.
        """
        region = self.patch_interval.to_region_string()

        cmd = [
            "samtools",
            "view",
            "-@",
            str(self.threads),
            "-b",
            "-P",
            "-q",
            str(self.mapq),
            "-o",
            str(self.paths.patch_bam),
            str(self.input_bam),
            region,
        ]

        print("Extracting patch BAM...")
        print(" ".join(cmd))

        subprocess.run(cmd, check=True)

    def collect_qnames_from_patch_bam(self) -> set[str]:
        """
        Collect all qnames present in the extracted patch BAM.
        """
        qnames: set[str] = set()

        with pysam.AlignmentFile(self.paths.patch_bam, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.query_name is not None:
                    qnames.add(read.query_name)

        return qnames

    def write_patch_qnames(self, patch_qnames: set[str]) -> None:
        """
        Write all patch qnames to PREFIX.patch.qnames.txt.
        """
        with open(self.paths.patch_qnames, "w") as handle:
            for qname in sorted(patch_qnames):
                handle.write(qname + "\n")

    def count_patch_records(self) -> int:
        """
        Count alignment records in the patch BAM.
        """
        n_records = 0

        with pysam.AlignmentFile(self.paths.patch_bam, "rb") as bam:
            for _ in bam.fetch(until_eof=True):
                n_records += 1

        return n_records

    def index_patch_bam(self) -> None:
        """
        Index the patch BAM.
        """
        index_bam(self.paths.patch_bam, threads=self.threads)

    def run(self) -> PatchDissectionResult:
        """
        Run fast patch dissection using samtools.
        """
        self.extract_patch_bam_with_samtools()

        patch_qnames = self.collect_qnames_from_patch_bam()

        self.write_patch_qnames(patch_qnames)

        self.index_patch_bam()

        result = PatchDissectionResult(
            patch_bam=self.paths.patch_bam,
            patch_qnames=self.paths.patch_qnames,
            n_patch_qnames=len(patch_qnames),
            n_patch_records=self.count_patch_records(),
        )

        self.print_summary(result)
        return result

    def print_summary(self, result: PatchDissectionResult) -> None:
        print("Patch dissection complete")
        print(f"Input BAM: {self.input_bam}")
        print(f"CNV interval: {self.interval}")
        print(f"Patch interval: {self.patch_interval}")
        print(f"MAPQ threshold: {self.mapq}")
        print(f"Patch BAM: {result.patch_bam}")
        print(f"Patch qnames: {result.patch_qnames}")
        print(f"Unique patch qnames: {result.n_patch_qnames:,}")
        print(f"Patch records written including mates: {result.n_patch_records:,}")