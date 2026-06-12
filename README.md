# CNVinject
current release: version 0.0.2

CNVinject is a command-line tool for injecting artificial copy number variants (CNVs) into existing BAM files while preserving the noise profile of the original sample.

The primary goal of CNVinject is to generate realistic positive-control BAMs for benchmarking sequencing depth- or breakpoint-based CNV detection workflows, especially in difficult sequencing contexts such as low-input, whole-genome-amplified, low-coverage, sparse, noisy, or unevenly amplified libraries. CNVinject was designed for situations where the **noise is part of the data** and should be retained.

CNVinject modifies an existing BAM by extracting aligned reads from a defined local genomic region (patch), editing the reads in that patch to represent the desired copy number state, and then replacing those reads in the input BAM to simulate a CNV. This approach preserves many properties of the original data and is compatible with read-depth and breakpoint-aware CNV callers. Read eligibility for modification is defined by the user.

What CNVinject preserves:

- breakpoints where possible;
- native sequencing depth;
- uneven and low coverage;
- whole-genome amplification artifacts;
- library-specific noise;
- duplicate reads;
- local mapping artifacts recorded in the CIGAR string (mismatches, deletions, insertions);
- incomplete read pairs;
- read-length, insert-size, and alignment behavior

---

> **Current implementation status**
>
> - `cnvinject del` (`--copy-number` 0 to <2) is implemented for simulating deletions, including fractional/mosaic deletions.
> - `cnvinject dup` (`--copy-number` > 2) is implemented for simulating duplications, including fractional/mosaic duplications.
> - `cnvinject mergepatch` is implemented for merging an edited patch back into the original full BAM.

---

## Intended use cases

CNVinject was originally designed for ultra-low-input, whole-genome-amplified sequencing libraries, where amplification bias, allelic dropout, uneven coverage, duplicate reads, and sparse genomic representation can make some simulation approaches unrealistic. CNVinject generates realistic artificial CNV-positive BAM files for:

- benchmarking read-depth CNV callers;
- when cell lines with desired CNVs are not accessible, including mosaic genotypes;
- testing CNV detection limits in data with low coverage, sparse coverage, amplification bias, sequencing and mapping artifacts;
- creating positive controls from real BAMs;
- comparing breakpoint sensitivity across defined genomic intervals;
- validating pipelines where retaining the original sample noise is important.

Intermediate FASTQ files (synthetic breakpoint/junction reads only) can be retained with `--disable-cleanup` where raw FASTQ files are necessary for benchmarking. Alternatively, users can convert the edited BAM back to FASTQ. It is important to note that the edited reads are designed to maintain sequencing and mapping artifacts present in the original alignments, and therefore small variants should be treated as artifacts.

CNVinject was developed and tested with bwa-mem alignments of ~150 bp Illumina paired-end reads. CNVinject parses alignments by query of qnames and retains both singleton and paired reads as a feature. CNVinject is compatible with both paired-end and single-end read workflows.

The current implementation of CNVinject excludes secondary and supplementary alignments from modification. Future implementations may incorporate user control over these alignments.

CNVinject uses a patch-based workflow for speed and memory optimization:

1. **Extract a local patch** around the target CNV interval (interval ± `--buffer`).
2. **Identify reads overlapping the target interval.**
3. **Edit the patch** to introduce the desired copy number state (remove reads for deletions; add donor reads for duplications).
4. **Synthesize FASTQ reads** spanning the breakpoints (deletions) or the novel tandem junction (duplications).
5. **Align synthetic reads** to the reference genome with bwa-mem.
6. **Assemble the final patch** from the edited reads and aligned synthetic reads.
7. **Patch the input BAM** with the edited reads (unless `--getpatch` is used).

To simulate deletions, CNVinject removes read pairs (or singletons) whose fragment lies fully inside the target interval and edits reads that overlap breakpoints. Breakpoint reads are edited such that their sequences match the reference genome adjacent to the target interval while preserving insert length and pair orientation. Soft-clipped sequences, single-nucleotide mutations, and indels recorded in the CIGAR string of the original unmodified read are perpetuated in the edited reads in order to preserve pre-existing sequencing artifacts as much as possible. Substitutions are introduced as follows: C>T, T>C, A>G, and G>A. This mutation scheme follows a purine>purine and pyrimidine>pyrimidine mutation rule. **Small variants in reads edited by CNVinject should be treated as artifacts.**

To simulate duplications, CNVinject adds reads rather than removing them; recipient reads are never deleted. Extra coverage is built from **donor BAMs** — other samples that, ideally, were prepared with the same library strategy as the input and therefore share its noise and coverage profile. Fully internal read pairs (or singletons) are sampled from the donor BAMs and added to the interval to raise depth in proportion to the requested copy number, spread as evenly as possible across the available donors. To represent the structural change, CNVinject also adds donor reads spanning the interval's outer boundaries (scaled to copy number) and generates synthetic reads spanning the novel tandem-duplication junction (interval end joined to interval start), which are aligned to the reference with bwa-mem. Donor BAMs whose filename matches the input BAM are ignored so the recipient is never sampled from itself.

It is important to note that the defined target interval of the input BAM is presumed diploid, and therefore reads will be added or removed (and breakpoint reads modified) at a rate relative to the input BAM. For example, if the specified copy number is 1, then half of the reads in the interval will be randomly selected for removal/breakpoint editing to reduce the coverage in the target interval by ~50%. For a copy number of 3, roughly 50% additional coverage is added from donor BAMs.

---
### Read eligibility
Currently, only primary alignments are eligible for modification. Secondary and supplementary alignments will not be touched since their alignment to the target interval is not certain and therefore may be interpreted as noise or artifact. A future `--all-alignments` flag is reserved for this purpose but is **not yet implemented** (currently a no-op placeholder).

`--mapq` controls which reads are eligible for CNV editing. Users can opt to retain low-MAPQ reads in the target interval at their original concentration to simulate low-quality/uncertain alignments in the target interval. The appropriate MAPQ threshold depends on the benchmarking goal. If the goal is to preserve the full noise profile of the original BAM, `--mapq 10` is preferable and is the default.

In whole-genome-amplified or ultra-low-input sequencing libraries, duplicates may be part of the observed read-depth structure. Removing all duplicates before simulation can make the injected CNV less representative of the original data. CNVinject is designed to preserve duplicated reads. After a CNV-injected BAM is generated, it is recommended to remove duplicate markings and/or remark duplicated reads.

---

### Limitations

CNVinject is under active development.
- Haplotype- and allele-specific injections are not supported at this time.
- The copy number of the target interval in the input BAM is presumed diploid. A later implementation may support haploid samples.
- The duplication workflow is newer than the deletion workflow and has had less validation.

---


# Installation

### Dependencies and recommended installation

CNVinject was developed and tested with:

- Python 3.12
- pysam 0.24.0
- samtools 1.19.2 using htslib 1.19
- bwa 0.7.17-r1188

The reference FASTA passed to `-r` must be indexed for bwa (`bwa index reference.fa`) and faidx (`samtools faidx reference.fa`).

```bash
# 1. Clone the repository
git clone https://github.com/kscott94/CNVinject.git
cd CNVinject
chmod +x bin/cnvinject

# 2. Create the conda environment
conda env create -f environment.yml
conda activate cnvinject

# 3. Add CNVinject/bin to your `PATH`. Replace `/path/to/CNVinject` with the full path to your cloned repository.
echo 'export PATH="/path/to/CNVinject/bin:$PATH"' >> ~/.bash_profile
source ~/.bash_profile

# 4. Test the installation
cnvinject --help
```

You should see the main help menu with available subcommands:

```text
cnvinject del
cnvinject dup
cnvinject mergepatch
```
---

## Quick start: full-copy deletion patch

Example: create a complete deletion patch to simulate NF1 microdeletion syndrome. Note, the outdir will be automatically created if it does not already exist.

```bash
cnvinject del \
  -i Sample1.bam \
  -o Sample1.NF1.CN0 \
  -r reference.fa \
  --getpatch \
  --copy-number 0 \
  --interval chr17:30780079-31936302 \
  --outdir ~/project/NF1_CN0
```
Output:
```text
Sample1.NF1.CN0.final.patch.bam
Sample1.NF1.CN0.final.patch.bam.bai
Sample1.NF1.CN0.patch.qnames.txt
```

This command extracts a patch around the requested interval, removes or modifies eligible reads from the deletion interval, and writes a coordinate-sorted and indexed `final.patch.bam`. With `--getpatch`, the program will stop after patch generation. If the user wants a whole-genome BAM, do not include `--getpatch`, or run `cnvinject mergepatch`.

To patch the original input BAM:

```bash
cnvinject mergepatch \
  --full-bam Sample1.bam \
  --patch-bam Sample1.NF1.CN0.final.patch.bam \
  --patch-reads Sample1.NF1.CN0.patch.qnames.txt \
  -o Sample1.NF1.CN0.final.bam
```

Adjust the `--patch-bam` and `--patch-reads` filenames to match the exact files produced by your `cnvinject del --getpatch` command.

---

## Quick start: tandem duplication patch

Example: create a CN3 tandem-duplication patch. Duplications require `--copy-number > 2` and a directory of donor BAMs supplied with `--donor-bam-dir`.

```bash
cnvinject dup \
  -i Sample2.bam \
  -o Sample2.NF1.CN3 \
  -r reference.fa \
  --copy-number 3 \
  --interval chr17:30780079-31936302 \
  --donor-bam-dir ~/project/donor_bams \
  --getpatch \
  --outdir ~/project/NF1_CN3
```
Output:
```text
Sample2.NF1.CN3.final.patch.bam
Sample2.NF1.CN3.final.patch.bam.bai
Sample2.NF1.CN3.patch.qnames.txt
```

Donor BAMs should be prepared with the same library strategy as the input BAM so they share its noise and coverage profile. Any donor BAM whose filename matches the input BAM is skipped. If the donor pool is too small to supply the requested boundary reads, rerun with `--allow-replacement`.

---

# Manual

Available commands:

| Command | Status | Purpose |
|---|---:|---|
| `del` | Implemented | Inject a deletion into a BAM. |
| `dup` | Implemented (active development) | Inject a duplication/amplification into a BAM. |
| `mergepatch` | Implemented | Merge an edited patch BAM back into the original input BAM. |


## Main command


```bash
cnvinject --help
cnvinject <command> [options] -i <input.bam> -o <output_prefix> -r <reference.fa> --copy-number <CN> --interval <chr:start-end>
```

`--input, -i`
```text
-i ~/project/bams/input.bam
```

`--reference, -r`
```text
-r ~/project/references/reference.fa
```
Required for `del` and `dup`. Used to fetch reference sequence and to realign synthetic breakpoint/junction reads with bwa-mem. Must be bwa- and faidx-indexed.

`--copy-number`
```text
--copy-number 0      → retain 0% reads     → homozygous deletion
--copy-number 0.25   → retain 12.5% reads  → mosaic deletion
--copy-number 1      → retain 50% reads    → hemizygous deletion
--copy-number 1.5    → retain 75% reads    → mosaic deletion
--copy-number 2      → does nothing (exits on error)
--copy-number 2.5    → subsamples donor bams to increase coverage by 25%    → mosaic duplication
--copy-number 3      → subsamples donor bams to increase coverage by 50%    → tandem duplication
--copy-number 3.5    → subsamples donor bams to increase coverage by 75%    → mosaic tandem duplication
--copy-number 4      → subsamples donor bams to increase coverage by 100%   → tandem duplication
```

Note, if `--copy-number 2` then the program will exit. Intervals are assumed to have a copy number of 2.


`--interval`

Takes samtools-style syntax.

chromosome:1_base_start_genomic_coordinate-end_genomic_coordinate. Be sure to use the chromosome labels in your BAM file.

```text
chr17:30780079-31936302
```

---

### `cnvinject del`

```bash
cnvinject del --help
cnvinject del [options] -i <input.bam> -o <OUTPUT_PREFIX> -r <reference.fa> --copy-number <0 to <2> --interval <chr:start-end>
```

#### options


| Argument | Required | Default | Description |
|---|---:|---:|---|
| `-i`, `--input` | Yes | none | Input BAM file. |
| `-o`, `--output` | Yes | none | Output prefix or output BAM path, depending on workflow. For patch generation, this is typically used as the prefix for patch-related output files. |
| `-r`, `--reference` | Yes | none | Reference FASTA (bwa- and faidx-indexed). Used for synthetic breakpoint read generation and bwa-mem realignment. |
| `--copy-number` | Yes | none | Target deletion copy number. Accepts decimals from `0` to `<2` (e.g. `0` = homozygous deletion, `1` = hemizygous, `1.5` = mosaic). |
| `--interval` | Yes | none | CNV interval in `chr:start-end` format. Example: `chr17:30780079-31936302`. |
| `--outdir` | No | current working directory | Output files will be directed to this directory. |
| `--bwa` | No | none | Extra arguments passed to `bwa mem`, as a quoted string. Example: `--bwa "-B 4 -O 6 -E 1 -L 5"`. |
| `--getpatch` | No | false | Write the edited patch BAM only instead of immediately producing a full edited BAM. This is useful for a two-step workflow where `mergepatch` is run separately. |
| `--disable-cleanup` | No | false | Keep all intermediate files. By default, cnvinject removes intermediates. It is recommended to disable cleanup for debugging purposes or if the user simply wants to better understand what is going on under the hood. |
| `--mapq` | No | `10` | Minimum mapping quality for reads eligible for mutation. Use `0` to allow all mapped reads regardless of mapping quality. |
| `--buffer` | No | `10000` | Number of bases to include upstream and downstream of the CNV interval when extracting the patch. |
| `-s`, `--seed` | No | none | Random seed. If none, a runtime seed will be randomly generated and printed to the console. |
| `-t`, `--threads` | No | `1` | Number of threads. |



#### Expected outputs

Exact filenames may depend on the implementation and the `-o/--output` value, but the deletion patch workflow is expected to produce files in the following categories:

| Output | Description |
|---|---|
| Edited patch BAM | BAM containing the extracted patch after CNV editing. |
| Edited patch BAM index | BAM index for the edited patch BAM, if indexing is performed. |
| Patch read-name list | Text file containing the original patch read names that should be removed from the full BAM before merging the edited patch. |
| Logs or progress messages | Run information useful for troubleshooting and reproducibility. |

A typical output set may look like:

```text
Sample1.NF1.CN0.final.patch.bam
Sample1.NF1.CN0.final.patch.bam.bai
Sample1.NF1.CN0.patch.qnames.txt
```

`patch.qnames.txt` is used by `cnvinject mergepatch`.

---

### `cnvinject dup`

Inject a duplication/amplification by adding reads sampled from donor BAMs.

```bash
cnvinject dup --help
cnvinject dup [options] -i <input.bam> -o <OUTPUT_PREFIX> -r <reference.fa> --copy-number <CN > 2> --interval <chr:start-end> --donor-bam-dir <dir>
```

Duplications require `--copy-number > 2` and a directory of donor BAMs (`--donor-bam-dir`). Donor BAMs should be prepared with the same library strategy as the input BAM so they share its noise and coverage profile. Any donor BAM whose filename matches the input BAM is ignored.

#### options

| Argument | Required | Default | Description |
|---|---:|---:|---|
| `-i`, `--input` | Yes | none | Input (recipient) BAM file. |
| `-o`, `--output` | Yes | none | Output prefix. |
| `-r`, `--reference` | Yes | none | Reference FASTA (bwa- and faidx-indexed). Used to build and realign synthetic tandem-junction reads. |
| `--copy-number` | Yes | none | Target duplication copy number. Must be `> 2`. Decimals supported (e.g. `3` = tandem duplication, `2.5` = mosaic gain). |
| `--interval` | Yes | none | CNV interval in `chr:start-end` format. |
| `--donor-bam-dir` | Yes | none | Directory of donor BAMs used to sample added reads. A BAM with the same name as the input is ignored. |
| `--allow-replacement` | No | false | Permit donor breakpoint reads to be sampled with replacement when the donor boundary pool is too small. By default, an insufficient pool raises an error. Replacement copies are jittered and realigned with bwa-mem rather than coordinate-shifted. |
| `--outdir` | No | current working directory | Output files will be directed to this directory. |
| `--bwa` | No | none | Extra arguments passed to `bwa mem`, as a quoted string. Example: `--bwa "-B 4 -O 6 -E 1 -L 5"`. |
| `--getpatch` | No | false | Write the edited patch BAM only instead of immediately producing a full edited BAM. |
| `--disable-cleanup` | No | false | Keep all intermediate files. |
| `--mapq` | No | `10` | Minimum mapping quality for reads eligible for sampling. Use `0` to allow all mapped reads. |
| `--buffer` | No | `10000` | Number of bases to include upstream and downstream of the CNV interval when extracting the patch. |
| `-s`, `--seed` | No | none | Random seed. If none, a runtime seed will be randomly generated and printed to the console. |
| `-t`, `--threads` | No | `1` | Number of threads. |

#### Expected outputs

| Output | Description |
|---|---|
| Edited patch BAM | BAM containing the recipient patch with added donor internal reads, boundary reads, and synthetic tandem-junction reads. |
| Edited patch BAM index | BAM index for the edited patch BAM, if indexing is performed. |
| Patch read-name list | Text file containing the original patch read names, used by `mergepatch`. |
| Logs or progress messages | Run information useful for troubleshooting and reproducibility. |

A typical output set may look like:

```text
Sample2.NF1.CN3.final.patch.bam
Sample2.NF1.CN3.final.patch.bam.bai
Sample2.NF1.CN3.patch.qnames.txt
```

---

## `cnvinject mergepatch`

Merge an edited patch BAM back into the original full BAM.

```bash
cnvinject mergepatch \
  --full-bam Sample2.bam \
  --patch-bam Sample2.final.patch.bam \
  --patch-reads Sample2.patch.qnames.txt \
  -o Sample2.final.bam \
  -t 1
```

### Purpose
By default, cnvinject generates a whole-genome BAM file with a modified target interval. The `--getpatch` flag will disable full BAM generation and instead produce only the patch BAM. This operation is relatively fast, whereas patching the full BAM (which calls samtools merge, sort, and index) is slow for large BAM files, like those produced from human libraries. It is recommended to apply `--getpatch` to quickly verify the pipeline is functioning properly before committing to a time-intensive full BAM generation.

For a modular approach, `mergepatch` was implemented to generate a whole-genome BAM containing the artificial CNV from the `final.patch.bam` and the `patch.qnames.txt` generated from `cnvinject del` and `cnvinject dup` operations. `mergepatch` removes the patch reads from the `--full-bam` and the result is merged with the `--patch-bam`.

#### Arguments

| Argument | Required | Default | Description |
|---|---:|---:|---|
| `--full-bam` | Yes | none | Original full BAM. |
| `--patch-bam` | Yes | none | Edited patch BAM produced by a CNVinject patch-generation command. |
| `--patch-reads` | Yes | none | Text file containing line-separated patch read names. This file is an output ending with `patch.qnames.txt`. |
| `-o`, `--output` | Yes | none | Output BAM file name including path. Example: `-o ~/project/output/final.bam`. |
| `-t`, `--threads` | No | `1` | Number of threads. |


#### Expected outputs

| Output | Description |
|---|---|
| Output BAM | Full BAM with the artificial CNV injected. |
| Output BAM index | BAM index, if indexing is performed by the workflow. |
| Logs or progress messages | Run information useful for troubleshooting. |

A typical final output is coordinate-sorted and indexed.

```text
Sample1.NF1.CN0.final.bam
Sample1.NF1.CN0.final.bam.bai
```



## Recommended validation after CNV injection

After generating a CNV-injected BAM, validate that the expected copy number change is present. Users should inspect the output using a genome browser and coverage summaries. For deletion simulations, users should observe reduced coverage across the target interval relative to nearby background regions and/or relative to the original BAM. For duplication simulations, users should observe increased coverage across the interval proportional to the requested copy number, and reads spanning the novel tandem junction at the interval boundaries.

BAM file integrity checks include:

```bash
samtools quickcheck -v Sample1.NF1.CN0.final.patch.bam
```

Coverage validation with `samtools coverage`:

```bash
samtools coverage Sample1.NF1.CN0.final.patch.bam
```

---

## Troubleshooting

**`error: the following arguments are required: -r`** — `del` and `dup` both require a reference FASTA. Pass `-r reference.fa`.

**bwa fails or the reference cannot be opened** — ensure the reference passed to `-r` is indexed: `bwa index reference.fa` and `samtools faidx reference.fa`.

**`--copy-number 2` exits with an error** — intervals are presumed diploid. Use `< 2` for a deletion (`cnvinject del`) or `> 2` for a duplication (`cnvinject dup`).

**`No usable donor BAMs found ...`** (dup) — the `--donor-bam-dir` is empty or contains only a BAM matching the input name. Add donor BAMs from other samples.

**`Not enough eligible internal donor qnames ...`** (dup) — the donor pool is too small for the requested copy number. Add more donor BAMs, lower `--mapq`, or use a lower copy number.

**`Not enough donor ... breakpoint templates ...`** (dup) — donor boundary coverage is too sparse. Rerun with `--allow-replacement` to permit resampling donor boundary reads.



## Development notes

CNVinject is currently in early development. The interface and output filenames may change.


## Citation

No formal citation is available yet. If you use CNVinject in a publication or internal benchmark, please cite the GitHub repository.
