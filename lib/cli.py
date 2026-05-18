#!/usr/bin/env python3

import argparse
from pathlib import Path
from parser import build_parser
from patch import GenomicInterval, PatchDissector
from deletion import CN0DeletionEditor
from breakpoints import BreakpointCandidateClassifier
from synthetic_reads import SyntheticBreakpointReadGenerator
from final_patch import FinalPatchBuilder
from mergepatch import MergePatchBuilder
from helpers import (align_fastq_with_bwa,
                     fastq_has_records,
                     make_output_prefix,
                     cleanup_intermediate_files)

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "del":
        try:
            run_deletion(args)
        except ValueError as e:
            parser.error(str(e))

    elif args.command == "dup":
        try:
            run_duplication(args)
        except ValueError as e:
            parser.error(str(e))

    elif args.command == "mergepatch":
        try:
            run_mergepatch(args)
        except ValueError as e:
            parser.error(str(e))

    else:
        parser.error(f"Unknown command: {args.command}")


def run_deletion(args: argparse.Namespace) -> None:

    if args.copy_number >= 2:
        raise ValueError("cnvinject del requires --copy-number < 2")

    seed_label = "none" if args.seed is None else str(args.seed)

    print("Running deletion workflow")
    print(f"Seed: {seed_label}")

    interval = GenomicInterval.from_string(args.interval)
    output_prefix = make_output_prefix(args)

    dissector = PatchDissector(
        input_bam=args.input,
        output_prefix=output_prefix,
        interval=interval,
        buffer_size=args.buffer,
        mapq=args.mapq,
        threads=args.threads,
    )

    patch_result = dissector.run()

    if args.copy_number == 0:
        editor = CN0DeletionEditor(
            input_patch_bam=patch_result.patch_bam,
            output_prefix=output_prefix,
            interval=interval,
            threads=args.threads,
        )

        deletion_result = editor.run()

        breakpoint_classifier = BreakpointCandidateClassifier(
            edited_patch_bam=deletion_result.edited_patch_bam,
            output_prefix=output_prefix,
            interval=interval,
        )

        breakpoint_result = breakpoint_classifier.run()

        synthetic_read_generator = SyntheticBreakpointReadGenerator(
            edited_patch_bam=deletion_result.edited_patch_bam,
            candidates_tsv=breakpoint_result.candidates_tsv,
            output_prefix=output_prefix,
            interval=interval,
            reference_fasta=args.reference,
        )

        synthetic_result = synthetic_read_generator.run()

        paired_breakpoint_bam = Path(f"{output_prefix}.breakpoint.paired.bam")
        singleton_breakpoint_bam = Path(f"{output_prefix}.breakpoint.singletons.bam")

        synthetic_bams = []

        if fastq_has_records(synthetic_result.paired_fastq):
            align_fastq_with_bwa(
                fastq=synthetic_result.paired_fastq,
                reference=args.reference,
                output_bam=paired_breakpoint_bam,
                threads=args.threads,
                interleaved=True,
                bwa_args=args.bwa_args,
            )
            synthetic_bams.append(paired_breakpoint_bam)

        if fastq_has_records(synthetic_result.singleton_fastq):
            align_fastq_with_bwa(
                fastq=synthetic_result.singleton_fastq,
                reference=args.reference,
                output_bam=singleton_breakpoint_bam,
                threads=args.threads,
                interleaved=False,
                bwa_args=args.bwa_args,
            )
            synthetic_bams.append(singleton_breakpoint_bam)

        final_patch_builder = FinalPatchBuilder(
            edited_patch_bam=deletion_result.edited_patch_bam,
            breakpoint_candidates_tsv=breakpoint_result.candidates_tsv,
            output_prefix=output_prefix,
            synthetic_bams=synthetic_bams,
            threads=args.threads,
        )

        final_patch_result = final_patch_builder.run()

        if not args.getpatch:
            merge_builder = MergePatchBuilder(
                full_bam=args.input,
                final_patch_bam=final_patch_result.final_patch_bam,
                patch_qnames=patch_result.patch_qnames,
                output_bam=Path(f"{output_prefix}.final.bam"),
                threads=args.threads,
            )

            merge_builder.run()
        else:
            print("--getpatch, skipping full BAM reconstruction.")

        if args.disable_cleanup:
            print("--disable-cleanup, keeping intermediate files.")
        else:
            cleanup_intermediate_files(
                output_prefix=output_prefix,
                keep_full_final=not args.getpatch,
            )


    elif args.copy_number == 1:
        raise NotImplementedError("CN1 deletion editing is not implemented yet.")


def run_duplication(args: argparse.Namespace) -> None:
    if args.copy_number <= 2:
        raise ValueError("cnvinject dup requires --copy-number > 2")

    seed_label = "none" if args.seed is None else str(args.seed)

    print("Running duplication workflow")
    print(f"Seed: {seed_label}")

    interval = GenomicInterval.from_string(args.interval)
    output_prefix = make_output_prefix(args)

    dissector = PatchDissector(
        input_bam=args.input,
        output_prefix=output_prefix,
        interval=interval,
        buffer_size=args.buffer,
        mapq=args.mapq,
        threads=args.threads,
    )

    dissector.run()


def run_mergepatch(args: argparse.Namespace) -> None:
    print("Running mergepatch workflow")

    merge_builder = MergePatchBuilder(
        full_bam=args.full_bam,
        final_patch_bam=args.patch_bam,
        patch_qnames=args.patch_reads,
        output_bam=args.output,
        threads=args.threads,
    )

    merge_builder.run()


if __name__ == "__main__":
    main()