#!/usr/bin/env python3

import argparse

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cnvinject",
        description="Inject synthetic CNVs into BAM files."
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True
    )

    # cnvinject del ----------------------------------------------
    del_parser = subparsers.add_parser(
        "del",
        help="Inject a deletion into a BAM."
    )

    add_common_cnv_args(del_parser)

    del_parser.add_argument(
        "--copy-number",
        dest="copy_number",
        type=int,
        required=True,
        choices=[0, 1],
        help="Deletion copy number. Supported: 0 or 1."
    )

    # cnvinject dup ----------------------------------------------
    dup_parser = subparsers.add_parser(
        "dup",
        help="Inject a duplication/amplification into a BAM."
    )

    add_common_cnv_args(dup_parser)

    dup_parser.add_argument(
        "--copy-number",
        dest="copy_number",
        type=int,
        required=True,
        help="Duplication copy number. Must be greater than 2."
    )

    dup_parser.add_argument(
        "--donor-bam-dir",
        dest="donor_bam_dir",
        required=True,
        help="Directory containing BAMs for duplication read sampling. Will not use a bam file with the same name as input bam file."
    )

    # cnvinject mergepatch ---------------------------------------
    merge_parser = subparsers.add_parser(
        "mergepatch",
        help="Merge an edited patch BAM back into the full BAM."
    )

    merge_parser.add_argument(
        "--full-bam",
        dest="full_bam",
        required=True,
        help="Original full BAM."
    )

    merge_parser.add_argument(
        "--patch-bam",
        dest="patch_bam",
        required=True,
        help="Edited patch BAM."
    )

    merge_parser.add_argument(
        "--patch-reads",
        dest="patch_reads",
        required=True,
        help="Text file containing line-separated patch qnames (i.e. <file>.patch.qnames.txt)."
    )

    merge_parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="Output edited BAM."
    )

    merge_parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=1,
        help="Number of threads. Default 1"
    )

    return parser


def add_common_cnv_args(parser: argparse.ArgumentParser) -> None:
    """
    Arguments shared by all cnvinject commands
    """

    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        help="Input BAM."
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        required=True,
        help="Output prefix."
    )

    parser.add_argument(
        "-r",
        dest="reference",
        required=True,
        help="Reference genome FASTA used for synthetic read generation and bwa mem realignment.",
    )

    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory. Default: current working directory.",
    )

    parser.add_argument(
        "--bwa",
        dest="bwa_args",
        default=None,
        help=(
            "Optional extra arguments passed to bwa mem, as a quoted string. "
            "Example: --bwa \"-B 4 -O 6 -E 1 -L 5\""
        ),
    )

    parser.add_argument(
        "--getpatch",
        action="store_true",
        help="Output edited patch only instead of full edited BAM."
    )

    parser.add_argument(
        "--disable-cleanup",
        action="store_true",
        help=("Keep intermediate files. By default, cnvinject removes intermediate files."
        ),
    )

    parser.add_argument(
        "--interval",
        dest="interval",
        required=True,
        help=("CNV interval in chr:start-end format. "
              "Example: chr17:30780079-31936302"
              )
    )

    parser.add_argument(
        "--mapq",
        dest="mapq",
        type=int,
        default=10,
        help="Minimum MAPQ for reads eligible for mutation. Default: 10."
    )

    parser.add_argument(
        "--buffer",
        dest="buffer",
        type=int,
        default=10000,
        help="Buffer size around CNV interval for patch extraction. Default 10000"
    )

    parser.add_argument(
        "-s",
        "--seed",
        dest="seed",
        type=int,
        default=None,
        required=False,
        help="Random seed for read sampling when copy number > 0."
    )

    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=1,
        help="Number of threads."
    )

