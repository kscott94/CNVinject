#!/usr/bin/env python3

from dataclasses import dataclass
from pathlib import Path
from collections import defaultdict
import random
import subprocess

import pysam

from patch import GenomicInterval, PatchPaths
from read_pairs import ReadPairStatus
from helpers import sort_and_index_bam


@dataclass
class DuplicationResult:
    input_patch_bam: Path
    added_internal_records_bam: Path
    edited_patch_bam: Path
    recipient_internal_qnames: Path
    donor_selected_qnames: Path
    n_recipient_internal_qnames: int
    n_qnames_requested: int
    n_qnames_selected: int
    n_added_records: int
    n_original_patch_records: int


class DuplicationEditor:
    """
    Build a duplication patch by adding internal reads sampled from donor BAMs.

    MVP behavior:
      - keep original recipient patch unchanged
      - identify recipient fully-internal qnames
      - compute how many donor qnames are needed for target CN
      - equally sample fully-internal qnames from donor BAMs
      - write sampled donor records to PREFIX.added.internal.records.bam
      - merge original patch + added internal records into PREFIX.edited.patch.bam
    """

    def __init__(
        self,
        input_patch_bam: str | Path,
        donor_bam_dir: str | Path,
        recipient_input_bam: str | Path,
        output_prefix: str | Path,
        interval: GenomicInterval,
        copy_number: float,
        seed: int,
        mapq: int = 0,
        threads: int = 1,
    ):
        self.input_patch_bam = Path(input_patch_bam)
        self.donor_bam_dir = Path(donor_bam_dir)
        self.recipient_input_bam = Path(recipient_input_bam).resolve()
        self.paths = PatchPaths(output_prefix)
        self.interval = interval
        self.copy_number = float(copy_number)
        self.seed = int(seed)
        self.mapq = int(mapq)
        self.threads = int(threads)

        self.rng = random.Random(self.seed)

        self.recipient_internal_qnames = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.recipient.internal.qnames.txt"
        )

        self.donor_selected_qnames = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.donor.selected.internal.qnames.tsv"
        )

        self.added_internal_records_bam = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.added.internal.records.bam"
        )

        self.edited_patch_bam = self.paths.prefix.with_name(
            f"{self.paths.prefix.name}.edited.patch.bam"
        )

        self.recipient_status_by_qname: dict[str, ReadPairStatus] = defaultdict(ReadPairStatus)
        self.recipient_internal_qnames_set: set[str] = set()

        self.n_qnames_requested = 0
        self.n_qnames_selected = 0
        self.n_added_records = 0
        self.n_original_patch_records = 0

    def added_fraction(self) -> float:
        if self.copy_number <= 2:
            raise ValueError("Duplication copy number must be > 2.")

        return (self.copy_number - 2.0) / 2.0

    def collect_recipient_internal_qnames(self) -> None:
        with pysam.AlignmentFile(self.input_patch_bam, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.query_name is None:
                    continue

                self.recipient_status_by_qname[read.query_name].add_read(read)

        self.recipient_internal_qnames_set = {
            qname
            for qname, status in self.recipient_status_by_qname.items()
            if status.is_fully_inside_interval(self.interval)
        }

        with open(self.recipient_internal_qnames, "w") as handle:
            for qname in sorted(self.recipient_internal_qnames_set):
                handle.write(qname + "\n")

    def calculate_requested_qnames(self) -> int:
        self.n_qnames_requested = round(
            len(self.recipient_internal_qnames_set) * self.added_fraction()
        )
        print(f"Requested donor qnames: {self.n_qnames_requested:,}")
        return self.n_qnames_requested

    def donor_bams(self) -> list[Path]:
        bams = sorted(self.donor_bam_dir.glob("*.bam"))

        usable = []
        recipient_name = self.recipient_input_bam.name

        for bam in bams:
            if bam.name.endswith(".bai"):
                continue

            # Exclude same basename as recipient. Resolve if possible, but do not
            # require every donor path to already exist cleanly.
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

    def collect_internal_qnames_from_donor(self, donor_bam: Path) -> set[str]:
        status_by_qname: dict[str, ReadPairStatus] = defaultdict(ReadPairStatus)

        with pysam.AlignmentFile(donor_bam, "rb") as bam:
            for read in bam.fetch(
                self.interval.chrom,
                self.interval.start0,
                self.interval.end0,
            ):
                if read.query_name is None:
                    continue

                if read.mapping_quality < self.mapq:
                    continue

                status_by_qname[read.query_name].add_read(read)

        return {
            qname
            for qname, status in status_by_qname.items()
            if status.is_fully_inside_interval(self.interval)
        }

    def sample_donor_qnames_evenly(self) -> dict[Path, set[str]]:
        """
        Sample donor qnames as evenly as possible across donor BAMs.

        Behavior:
          1. Assign each donor an equal target.
          2. If a donor has fewer eligible qnames than its target, take all available.
          3. Redistribute the remaining deficit across donors with extra eligible qnames.
          4. Fail only if the total donor pool is smaller than the total number needed.
        """

        donors = self.donor_bams()
        n_needed = self.calculate_requested_qnames()

        if n_needed == 0:
            return {donor: set() for donor in donors}

        # Collect eligible qnames from all donors first.
        available_by_donor: dict[Path, list[str]] = {}

        for donor_bam in donors:
            donor_internal = list(self.collect_internal_qnames_from_donor(donor_bam))
            self.rng.shuffle(donor_internal)
            available_by_donor[donor_bam] = donor_internal
            print(
                f"Donor eligible qnames: {donor_bam.name}: "
                f"{len(donor_internal):,}"
            )

        total_available = sum(len(qnames) for qnames in available_by_donor.values())

        if total_available < n_needed:
            raise ValueError(
                "Not enough eligible internal donor qnames across all donor BAMs.\n"
                f"Needed: {n_needed:,}\n"
                f"Available: {total_available:,}\n"
                "Try using more donor BAMs, lowering --mapq, increasing donor pool size, "
                "or using a lower duplication copy number."
            )

        # Initial equal target per donor.
        base = n_needed // len(donors)
        remainder = n_needed % len(donors)

        selected_by_donor: dict[Path, set[str]] = {donor: set() for donor in donors}
        leftovers_by_donor: dict[Path, list[str]] = {}

        selected_count = 0

        for donor_index, donor_bam in enumerate(donors):
            target = base + (1 if donor_index < remainder else 0)
            available = available_by_donor[donor_bam]

            n_take = min(target, len(available))

            selected = set(available[:n_take])
            leftovers = available[n_take:]

            selected_by_donor[donor_bam] = selected
            leftovers_by_donor[donor_bam] = leftovers

            selected_count += n_take

        # Fill any deficit from donors that had extra qnames.
        deficit = n_needed - selected_count

        if deficit > 0:
            donor_cycle = donors.copy()
            self.rng.shuffle(donor_cycle)

            while deficit > 0:
                made_progress = False

                for donor_bam in donor_cycle:
                    if deficit == 0:
                        break

                    leftovers = leftovers_by_donor[donor_bam]

                    if not leftovers:
                        continue

                    qname = leftovers.pop()
                    selected_by_donor[donor_bam].add(qname)
                    deficit -= 1
                    made_progress = True

                if not made_progress:
                    raise RuntimeError(
                        "Internal error while redistributing donor qnames. "
                        "Total available qnames should have been sufficient."
                    )

        return selected_by_donor

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

    def renamed_qname(self, qname: str, donor_bam: Path) -> str:
        return f"{qname}:CNVinject_donor_{self.donor_label(donor_bam)}"

    def write_added_internal_records_bam(
        self,
        selected_by_donor: dict[Path, set[str]],
    ) -> None:
        """
        Write all records from donor BAMs for selected donor qnames.

        Eligibility is determined at the qname level using primary alignments,
        but once a donor qname is selected, all interval-overlapping records
        with that qname are copied into added.internal.records.bam.
        """
        if not selected_by_donor:
            raise ValueError("No donor qnames were selected.")

        first_donor = next(iter(selected_by_donor))

        unsorted_bam = self.added_internal_records_bam.with_name(
            f"{self.paths.prefix.name}.added.internal.records.unsorted.tmp.bam"
        )

        with pysam.AlignmentFile(first_donor, "rb") as template_bam:
            with pysam.AlignmentFile(
                unsorted_bam,
                "wb",
                template=template_bam,
            ) as bam_out:

                with open(self.donor_selected_qnames, "w") as qname_out:
                    qname_out.write("donor_bam\toriginal_qname\trenamed_qname\n")

                    for donor_bam, selected_qnames in selected_by_donor.items():
                        print(
                            f"Donor selected qnames: {donor_bam.name}: "
                            f"{len(selected_qnames):,}"
                        )

                        self.n_qnames_selected += len(selected_qnames)

                        rename_map = {
                            qname: self.renamed_qname(qname, donor_bam)
                            for qname in selected_qnames
                        }

                        for original_qname, renamed_qname in sorted(rename_map.items()):
                            qname_out.write(
                                f"{donor_bam}\t{original_qname}\t{renamed_qname}\n"
                            )

                        with pysam.AlignmentFile(donor_bam, "rb") as bam_in:
                            for read in bam_in.fetch(
                                self.interval.chrom,
                                self.interval.start0,
                                self.interval.end0,
                            ):
                                if read.query_name not in selected_qnames:
                                    continue

                                if read.mapping_quality < self.mapq:
                                    continue

                                copied = pysam.AlignedSegment.fromstring(
                                    read.to_string(),
                                    bam_in.header,
                                )
                                copied.query_name = rename_map[read.query_name]

                                bam_out.write(copied)
                                self.n_added_records += 1

        sort_and_index_bam(
            input_bam=unsorted_bam,
            output_bam=self.added_internal_records_bam,
            threads=self.threads,
            remove_input=True,
        )


    def write_edited_patch_bam(self) -> None:
        """
        Merge original recipient patch + added donor internal records.
        """
        unsorted = self.edited_patch_bam.with_name(
            f"{self.paths.prefix.name}.edited.patch.unsorted.tmp.bam"
        )

        cmd = [
            "samtools",
            "merge",
            "-@",
            str(self.threads),
            "-f",
            str(unsorted),
            str(self.input_patch_bam),
            str(self.added_internal_records_bam),
        ]

        print("Merging recipient patch with added donor internal records:")
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)

        sort_and_index_bam(
            input_bam=unsorted,
            output_bam=self.edited_patch_bam,
            threads=self.threads,
            remove_input=True,
        )

        with pysam.AlignmentFile(self.input_patch_bam, "rb") as bam:
            for _ in bam.fetch(until_eof=True):
                self.n_original_patch_records += 1

    def run(self) -> DuplicationResult:
        self.collect_recipient_internal_qnames()
        selected_by_donor = self.sample_donor_qnames_evenly()
        self.write_added_internal_records_bam(selected_by_donor)
        self.write_edited_patch_bam()

        result = DuplicationResult(
            input_patch_bam=self.input_patch_bam,
            added_internal_records_bam=self.added_internal_records_bam,
            edited_patch_bam=self.edited_patch_bam,
            recipient_internal_qnames=self.recipient_internal_qnames,
            donor_selected_qnames=self.donor_selected_qnames,
            n_recipient_internal_qnames=len(self.recipient_internal_qnames_set),
            n_qnames_requested=self.n_qnames_requested,
            n_qnames_selected=self.n_qnames_selected,
            n_added_records=self.n_added_records,
            n_original_patch_records=self.n_original_patch_records,
        )

        self.print_summary(result)
        return result

    def print_summary(self, result: DuplicationResult) -> None:
        print("Duplication patch editing complete")
        print(f"Copy number: {self.copy_number}")
        print(f"Added fraction: {self.added_fraction():.4f}")
        print(f"Duplication interval: {self.interval}")
        print(f"Input patch BAM: {result.input_patch_bam}")
        print(f"Added internal records BAM: {result.added_internal_records_bam}")
        print(f"Edited patch BAM: {result.edited_patch_bam}")
        print(f"Recipient internal qnames: {result.n_recipient_internal_qnames:,}")
        print(f"Requested donor qnames: {result.n_qnames_requested:,}")
        print(f"Selected donor qnames: {result.n_qnames_selected:,}")
        print(f"Original patch records: {result.n_original_patch_records:,}")
        print(f"Added donor records: {result.n_added_records:,}")