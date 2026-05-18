import re
import subprocess
import argparse
import shlex
from pathlib import Path


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1].upper()


def force_mismatch(base: str) -> str:
    """
    Deterministically convert a reference base into a mismatch base.

    Purines:
        A <-> G

    Pyrimidines:
        C <-> T
    """
    base = base.upper()

    if base == "A":
        return "G"
    if base == "G":
        return "A"
    if base == "C":
        return "T"
    if base == "T":
        return "C"

    return "N"


def parse_cigar_tuples(cigar: str) -> list[tuple[int, str]]:
    """
    Parse CIGAR string into a list of (length, op) tuples.
    """
    if cigar is None or cigar == "*":
        raise ValueError("Missing CIGAR string.")

    parts = re.findall(r"(\d+)([MIDNSHP=X])", cigar)

    if not parts or "".join(length + op for length, op in parts) != cigar:
        raise ValueError(f"Malformed CIGAR string: {cigar}")

    return [(int(length), op) for length, op in parts]


def new_cigar(cigar: str, read_seq: str, ref_seq: str) -> str:
    """
    Convert M operations in a CIGAR string into explicit = and X by comparing
    read_seq to ref_seq.

    If the CIGAR has no M operations, return the original CIGAR unchanged.

    Keeps I, D, N, S, H, P, existing =, and existing X as-is.

    Assumptions:
      - read_seq is the BAM read sequence.
      - ref_seq is the reference sequence over the read's reference span:
            fasta.fetch(read.reference_name,
                        read.reference_start,
                        read.reference_end)
    """
    if cigar is None or cigar == "":
        raise ValueError(
            "Cannot create new CIGAR because the read has no CIGAR string. "
            "This usually means the read is unmapped or the BAM record is malformed."
        )

    if cigar == "*":
        raise ValueError(
            "Cannot create new CIGAR from '*'. "
            "This usually means the read is unmapped."
        )

    # If there are no M operations, there is nothing ambiguous to convert.
    # Existing = and X are already explicit.
    if "M" not in cigar:
        return cigar

    read_seq = read_seq.upper()
    ref_seq = ref_seq.upper()

    parts = parse_cigar_tuples(cigar)

    qpos = 0
    rpos = 0
    output = []

    def add_op(length: int, op: str) -> None:
        if length == 0:
            return

        if output and output[-1][1] == op:
            old_len, old_op = output[-1]
            output[-1] = (old_len + length, old_op)
        else:
            output.append((length, op))

    for length, op in parts:
        if op == "M":
            for _ in range(length):
                if qpos >= len(read_seq):
                    raise ValueError(
                        f"Cannot convert CIGAR {cigar}: read sequence is shorter "
                        f"than expected. Stopped at query index {qpos}, "
                        f"read length is {len(read_seq)}."
                    )

                if rpos >= len(ref_seq):
                    raise ValueError(
                        f"Cannot convert CIGAR {cigar}: reference sequence is shorter "
                        f"than expected. Stopped at reference index {rpos}, "
                        f"reference length is {len(ref_seq)}."
                    )

                read_base = read_seq[qpos]
                ref_base = ref_seq[rpos]

                if read_base == ref_base and read_base != "N":
                    add_op(1, "=")
                else:
                    add_op(1, "X")

                qpos += 1
                rpos += 1

        elif op in {"=", "X"}:
            # Already explicit, so keep as-is.
            add_op(length, op)
            qpos += length
            rpos += length

        elif op in {"I", "S"}:
            add_op(length, op)
            qpos += length

        elif op in {"D", "N"}:
            add_op(length, op)
            rpos += length

        elif op in {"H", "P"}:
            add_op(length, op)

        else:
            raise ValueError(f"Unsupported CIGAR operation '{op}' in CIGAR: {cigar}")

    return "".join(f"{length}{op}" for length, op in output)


def apply_n_mask(new_seq: str, original_seq: str) -> str:
    """
    Preserve N bases from the original read sequence.

    If original_seq has N at a query position, force new_seq to also have N
    at that same position.
    """
    if len(new_seq) != len(original_seq):
        raise ValueError(
            f"Cannot apply N mask because sequence lengths differ: "
            f"new_seq={len(new_seq)}, original_seq={len(original_seq)}"
        )

    new_seq_list = list(new_seq.upper())
    original_seq = original_seq.upper()

    for i, base in enumerate(original_seq):
        if base == "N":
            new_seq_list[i] = "N"

    return "".join(new_seq_list)


def make_fastq_record(name: str, seq: str, phred_char: str = "I") -> str:
    """
    Create a FASTQ record with a constant quality character.

    Default phred_char='I' corresponds to Q40.
    """
    seq = seq.upper()
    qual = phred_char * len(seq)

    return f"@{name}\n{seq}\n+\n{qual}\n"

def fetch_segments(fasta, segments: str) -> str:
    """
    Fetch and concatenate reference segments.

    Input format:
        chr17:100-200
        chr17:100-200|chr17:500-550

    Coordinates are 0-based half-open.
    """
    seq_parts = []

    for segment in segments.split("|"):
        chrom, coords = segment.split(":")
        start, end = coords.split("-")

        start = int(start)
        end = int(end)

        seq_parts.append(fasta.fetch(chrom, start, end).upper())

    return "".join(seq_parts)


def fastq_has_records(fastq: str | Path) -> bool:
    """
    Return True if a FASTQ file exists and contains at least one complete record.

    FASTQ records are 4 lines each.
    """
    fastq = Path(fastq)

    if not fastq.exists():
        return False

    if fastq.stat().st_size == 0:
        return False

    with open(fastq) as handle:
        first_line = handle.readline().strip()

    return first_line.startswith("@")

def sort_bam(
    input_bam: str | Path,
    output_bam: str | Path,
    threads: int = 1,
) -> Path:
    """
    Sort a BAM file with samtools sort.
    """
    input_bam = Path(input_bam)
    output_bam = Path(output_bam)

    if not input_bam.exists():
        raise FileNotFoundError(f"Input BAM does not exist: {input_bam}")

    cmd = [
        "samtools",
        "sort",
        "-@",
        str(threads),
        "-o",
        str(output_bam),
        str(input_bam),
    ]

    print("Sorting BAM...")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)

    return output_bam


def index_bam(
    bam: str | Path,
    threads: int = 1,
) -> Path:
    """
    Index a BAM file with samtools index.
    """
    bam = Path(bam)

    if not bam.exists():
        raise FileNotFoundError(f"BAM does not exist: {bam}")

    cmd = [
        "samtools",
        "index",
        "-@",
        str(threads),
        str(bam),
    ]

    print("Indexing BAM...")
    print(" ".join(cmd))

    subprocess.run(cmd, check=True)

    return Path(str(bam) + ".bai")


def sort_and_index_bam(
    input_bam: str | Path,
    output_bam: str | Path,
    threads: int = 1,
    remove_input: bool = False,
) -> Path:
    """
    Sort input_bam into output_bam, then index output_bam.

    If remove_input=True, delete input_bam after successful sorting/indexing.
    """
    input_bam = Path(input_bam)
    output_bam = Path(output_bam)

    sort_bam(
        input_bam=input_bam,
        output_bam=output_bam,
        threads=threads,
    )

    index_bam(
        bam=output_bam,
        threads=threads,
    )

    if remove_input and input_bam != output_bam:
        input_bam.unlink()

    return output_bam


def align_fastq_with_bwa(
    fastq: str | Path,
    reference: str | Path,
    output_bam: str | Path,
    threads: int = 1,
    interleaved: bool = False,
    bwa_args: str | None = None,
) -> Path:
    """
    Align a FASTQ file to a reference using bwa mem, then sort and index with samtools.

    Parameters
    ----------
    fastq:
        Input FASTQ file.

    reference:
        Reference FASTA indexed for bwa.

    bwa:
        additional bwa arguments

    output_bam:
        Sorted BAM output path.

    threads:
        Number of threads for bwa mem and samtools sort.

    interleaved:
        If True, run bwa mem with -p for interleaved paired-end FASTQ.

    Returns
    -------
    Path to the sorted/indexed BAM.
    """
    fastq = Path(fastq)
    reference = Path(reference)
    output_bam = Path(output_bam)

    if not fastq.exists():
        raise FileNotFoundError(f"FASTQ does not exist: {fastq}")

    if not reference.exists():
        raise FileNotFoundError(f"Reference FASTA does not exist: {reference}")

    output_bam.parent.mkdir(parents=True, exist_ok=True)

    bwa_cmd = [
        "bwa",
        "mem",
        "-t",
        str(threads),
    ]

    if interleaved:
        bwa_cmd.append("-p")

    if bwa_args:
        bwa_cmd.extend(shlex.split(bwa_args))

    bwa_cmd.extend([
        str(reference),
        str(fastq),
    ])

    sort_cmd = [
        "samtools",
        "sort",
        "-@",
        str(threads),
        "-o",
        str(output_bam),
        "-",
    ]

    print("Aligning synthetic reads...")
    print(" ".join(bwa_cmd) + " | " + " ".join(sort_cmd))

    bwa_process = subprocess.Popen(
        bwa_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    )

    sort_process = subprocess.Popen(
        sort_cmd,
        stdin=bwa_process.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    )

    # Allow bwa to receive SIGPIPE if samtools exits early.
    if bwa_process.stdout is not None:
        bwa_process.stdout.close()

    sort_stdout, sort_stderr = sort_process.communicate()
    bwa_stderr = bwa_process.stderr.read() if bwa_process.stderr is not None else b""

    bwa_returncode = bwa_process.wait()
    sort_returncode = sort_process.returncode

    if bwa_returncode != 0:
        raise RuntimeError(
            "bwa mem failed.\n"
            f"Command: {' '.join(bwa_cmd)}\n"
            f"stderr:\n{bwa_stderr.decode(errors='replace')}"
        )

    if sort_returncode != 0:
        raise RuntimeError(
            "samtools sort failed.\n"
            f"Command: {' '.join(sort_cmd)}\n"
            f"stderr:\n{sort_stderr.decode(errors='replace')}"
        )

    index_cmd = [
        "samtools",
        "index",
        "-@",
        str(threads),
        str(output_bam),
    ]

    print("Indexing BAM...")
    print(" ".join(index_cmd))

    subprocess.run(index_cmd, check=True)

    return output_bam

def make_output_prefix(args: argparse.Namespace) -> Path:
    """
    Combine --outdir and -o/--output into one output prefix path.

    If args.output already includes a directory, it is still placed under --outdir
    unless it is absolute.
    """
    outdir = Path(args.outdir)

    output = Path(args.output)

    if output.is_absolute():
        return output

    return outdir / output


def remove_file_if_exists(path: str | Path) -> None:
    """
    Delete a file if it exists.
    """
    path = Path(path)

    if path.exists():
        print(f"Removing intermediate file: {path}")
        path.unlink()


def cleanup_intermediate_files(
    output_prefix: str | Path,
    keep_full_final: bool = True,
) -> None:
    """
    Remove intermediate files created by the CN0 workflow.

    Default retained files:
        PREFIX.final.patch.bam
        PREFIX.final.patch.bam.bai
        PREFIX.final.bam
        PREFIX.final.bam.bai

    If keep_full_final=False, retain only:
        PREFIX.final.patch.bam
        PREFIX.final.patch.bam.bai
    """
    prefix = Path(output_prefix)

    keep = {
        Path(f"{prefix}.final.patch.bam"),
        Path(f"{prefix}.final.patch.bam.bai"),
    }

    if keep_full_final:
        keep.add(Path(f"{prefix}.final.bam"))
        keep.add(Path(f"{prefix}.final.bam.bai"))

    # Known intermediate suffixes generated by the current pipeline.
    candidates = [
        Path(f"{prefix}.patch.bam"),
        Path(f"{prefix}.patch.bam.bai"),
        Path(f"{prefix}.edited.patch.bam"),
        Path(f"{prefix}.edited.patch.bam.bai"),
        Path(f"{prefix}.breakpoint.paired.bam"),
        Path(f"{prefix}.breakpoint.paired.bam.bai"),
        Path(f"{prefix}.breakpoint.singletons.bam"),
        Path(f"{prefix}.breakpoint.singletons.bam.bai"),
        Path(f"{prefix}.breakpoint.paired.interleaved.fastq"),
        Path(f"{prefix}.breakpoint.singletons.fastq"),
        #Path(f"{prefix}.patch.qnames.txt"), keep patch qnames for mergepatch
        Path(f"{prefix}.internal.qnames.txt"),
        Path(f"{prefix}.breakpoint.candidates.tsv"),
    ]

    for path in candidates:
        if path not in keep:
            remove_file_if_exists(path)