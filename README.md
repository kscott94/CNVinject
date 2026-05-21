# CNVinject

CNVinject is a command-line tool for injecting artificial copy number variants (CNVs) into existing BAM files while preserving the noise profile of the original sample.

The primary goal of CNVinject is to generate realistic positive-control BAMs for benchmarking sequencing depth- or breakpoint-based CNV detection workflows, especially in difficult sequencing contexts such as low-input, whole-genome-amplified, low-coverage, sparsley, noisy, or unevenly amplified libraries. CNVinject was designed for situations where the **noise is part of the data** and should be retained.

CNVinject modifies an existing BAM by extracting aligned reads from a defined local genomic region (patch), editing the reads in that patch to represent the desired copy number state, and then replacing edited reads in original BAM to simulate a CNV. This approach preserves many properties of the original data and is compatible with read-depth and breakpoint aware CNV callers. The program does not rely on contig generation and read eligibility is defined by the user. 

What CNVinject preserves:

- breakpoints at defined interval
- native sequencing depth;
- uneven coverage;
- whole-genome amplification artifacts;
- sparse or low-coverage regions;
- library-specific noise;
- duplicate reads;
- local mapping artifacts recorded in CIGAR string (mismatches, deletions, insertions);
- incomplete read pairs;
- read-length, insert-size, and alignment behavior.

---

> **Current implementation status**
>
> - `cnvinject del --copy-number [0 to <2] is implemented for simulating deletions.
> - `cnvinject dup` is a placeholder for duplication simulation.
> - `cnvinject mergepatch` is implemented for merging an edited patch or alignments back into the original full BAM.

---

## Intended use cases

CNVinject is intended for generating realistic artificial CNV-positive BAM files for:

- benchmarking read-depth CNV callers;
- when cell lines with desired CNVs are not accessible;
- testing CNV detection limits in data with low coverage, sparse coverage, amplification bias, sequencing and mapping artifacts;
- creating positive controls from real BAMs;
- comparing breakpoint sensitivity across defined genomic intervals;
- validating pipelines where retaining the original sample noise is important.


CNVinject was originally designed for ultra-low-input, whole-genome-amplified sequencing libraries, where amplification bias, allelic dropout, uneven coverage, duplicate reads, and sparse genomic representation can make some simulation approaches unrealistic.

CNVinject was developed and tested with bwa-mem alignments of ~150 bp Illumina paired-end reads. Since CNVinject retaines singleton reads as a feature, this pipeline is compatible with single-end read samples. 

CNVinject uses a patch-based workflow:

1. **Extract a local patch** around the target CNV interval.
2. **Identify reads overlapping the target interval.**
3. **Edit the patch** to create the desired copynumber state.
4. **Synthesize artifical fastq reads** that overlap with breakpoints.
5. **Align synthetic reads** to reference genome.
6. **Replace original patch reads with synthic aligned read in input BAM** using `cnvinject mergepatch`.

To simulate deletions, CNVinject removes read pairs (or singletons) from the target interval and edits reads that overlap breakpoints. Breakpoint reads are edited such that their sequecences match that of the reference genome adjacent to the target interval while preserving insert length and pari orientation. Softclipped sequences, single nucleotide mutations, single nucleotide deletions, and indels recorded in the CIGAR string in the original unmodified read are perpetuated in the edited reads in order to preserve pre-existing sequencing artifacts as much as possible. Substitutions are introduced as follows: C>T, T>C, A>G, and G>A. This mutation scheme follows a purine>purine and pyrimidnie>pyrimidine mutation rule. **Small variants in reads edited by CNVinject should be treated as artifacts**. 

To simulate duplications, reads that overlap the target interval are randomly sampled from doner bam files that presumably were prepared using the same libarary prep strategy as the input bam and therefore share the same noise and coverage profile as the input bam. Internal and breakpoint reads from each doner bam are extracted and combined into a single file before read pairs (or singletons) are randomly subsampled to increase the coverage at target region proportional to copy number. 

It is important to note that the defined target inerval of the input bam is presumed diploid, and therefore reads will be added or removed (and breakpoint reads modified) at a rate relative to the imput bam. For example, if the specified copy number is 1, then half of the reads in the interval will be randomly selected for removal to reduce the coverage in the target interal by 50%.

---

## Limitations

CNVinject is under active development. 
- Haplotype and allele specific injections are not supported at this time.
- The copy number of target loci in the input BAM is presumed diploid. A later implemenation may support a haploid-state. 

---


# Installation

## Dependencies and recommended installation

CNVinject was developed and tested with:

- Python 3.12
- pysam 0.24.0
- samtools 1.19.2 using htslib 1.19
- bwa 0.7.17-r1188

### 1. Clone the repository

```bash
git clone https://github.com/kscott94/CNVinject.git
cd CNVinject
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate cnvinject
```

### 3. Add CNVinject/bin to your `PATH`

Replace `</path/to/CNVinject>` with the full path to your cloned repository.

```bash
echo 'export PATH="</path/to/CNVinject/bin>:$PATH"'
```

### 4. Test the installation

```bash
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

Example: create a complete deletion patch to simulate NF1 microdeletion syndrome. Note, the ourdir will be automatically created if it does not already exist. 

```bash
cnvinject del \
  -i Sample2.bam \
  -o Sample2.NF1.CN0 \
  --getpatch \
  --copy-number 0 \
  --interval chr17:30780079-31936302 \
  --outdir ~/project/out \
  -t 1
```

This command extracts a patch around the requested interval, removes or modifies eligible reads from the deletion interval, and writes a coordinate sorted and indexed final.patch.bam. With --getpatch, the program will stop after patch generation. If the user wants a whole genome BAM, do not include --getpatch or run cnvinject mergepatch.

To patch the original input BAM:

```bash
cnvinject mergepatch \
  --full-bam Sample2.bam \
  --patch-bam Sample2.NF1.CN0.final.patch.bam \
  --patch-reads Sample2.NF1.CN0.patch.qnames.txt \
  -o Sample2.NF1.CN0.full.bam \
  -t 8
```


> Adjust the `--patch-bam` and `--patch-reads` filenames to match the exact files produced by your `cnvinject del --getpatch` command.


---

# Manual 


## Main command


```bash
cnvinject --help
cnvinject <command> [options] -i <input.bam> -o <output_prefix> --copy-number <CN> --interval <chr:start-end>
```


Available commands:

| Command | Status | Purpose |
|---|---:|---|
| `del` | Partially implemented | Inject a deletion into a BAM. |
| `mergepatch` | Implemented | Merge an edited patch BAM back into the original input BAM. |
| `dup` | Placeholder | Future command for duplication/amplification simulation. |

`--input, -i`
```text
-i ~/project/bams/input.bam
```

`--copy-number [0-2]`
```text
--copy-number 0      → retain 0% reads    → homozygous deletion
--copy-number 0.25   → retain 12.5% reads  → mosaic deletion
--copy-number 1      → retain 50% reads    → hemizygous deletion
--copy-number 1.5    → retain 75% reads    → mosaic deletion
--copy-number 2      → does nothing
```

`--copy-number [>2]`
```text
--copy-number 2.5    → subsamples doner bams to increase coverage by 25%    → mosaic duplication
--copy-number 3      → subsamples doner bams to increase coverage by 50%    → duplication
--copy-number 3.5    → subsamples doner bams to increase coverage by 75%    → mosaic duplication
--copy-number 4      → subsamples doner bams to increase coverage by 100%   → duplication
```

Note, if --copy-number 2 then the program will exit. Intervals are assumed to have a copy number of 2. If the sample is haploid, you can simulate a diploid state by invoking cnvinject dup --copy-number 4, which will double the coverage at the target interval".




`--interval`
```text
#samtools-style syntax. 
#chromosome:1_base_start_genomic_coordinate-end_genomic_coordinate 
chr17:30780079-31936302
```


---


## `cnvinject del`

```bash
cnvinject del --help
```


```bash
cnvinject del \
  -i ~/project/input.bam \
  -o OUTPUT_PREFIX \
  --outdir ~/project/output \
  --copy-number 0 \
  --interval chr:start-end \
  --getpatch \
  --mapq 10
```


### options


| Argument | Required | Default | Description |
|---|---:|---:|---|
| `-i`, `--input` | Yes | none | Input BAM file. |
| `-o`, `--output` | Yes | none | Output prefix or output BAM path, depending on workflow. For patch generation, this is typically used as the prefix for patch-related output files. |
| `--copy-number` | Yes | none | Target deletion copy number. The CLI currently accepts `0` or `1`; only `0` is considered completed/validated at this time. |
| `--interval` | Yes | none | CNV interval in `chr:start-end` format. Example: `chr17:30780079-31936302`. |
| `--getpatch` | No | false | Write the edited patch BAM only instead of immediately producing a full edited BAM. This is useful for a two-step workflow where `mergepatch` is run separately. |
| `--disable-cleanup` | No | false |  Keep all intermediate files. By default, cnvinject removes intermediates. It is recommended to disable cleanup for debugging purposes or if the user simply wants to better understand what is going on under the hood. |
| `--mapq` | No | `0` | Minimum mapping quality for reads eligible for mutation. Use `0` to allow all mapped reads regardless of mapping quality. |
| `--buffer` | No | `10000` | Number of bases to include upstream and downstream of the CNV interval when extracting the patch. |
| `-s`, `--seed` | No | none | Random seed. This is mainly relevant for workflows involving random read sampling. |
| `-t`, `--threads` | No | `1` | Number of threads. |



### Expected outputs

Exact filenames may depend on the implementation and the `-o/--output` value, but the deletion patch workflow is expected to produce files in the following categories:

| Output | Description |
|---|---|
| Edited patch BAM | BAM containing the extracted patch after CNV editing. |
| Edited patch BAM index | BAM index for the edited patch BAM, if indexing is performed. |
| Patch read-name list | Text file containing the original patch read names that should be removed from the full BAM before merging the edited patch. |
| Logs or progress messages | Run information useful for troubleshooting and reproducibility. |

A typical output set may look like:

```text
Sample2.NF1.CN0.final.patch.bam
Sample2.NF1.CN0.final.patch.bam.bai
Sample2.NF1.CN0.patch.qnames.txt
```

patch.qnames.txt is used by `cnvinject mergepatch`.

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
By default, cnvnject generates a whole genome BAM file with a modified target interval. The `--getpatch` flag will disable full BAM generation and instead produce only the patch BAM. This operation is relatively fast whereas  patching the full BAM (which calls samtools merge, sort, and index) is slow for large BAM files, like those produced from human libraries. It is recommended to apply `--getpatch` to quickly verify the pipeline is functioning properly before committing to a time intensive full BAM generation. 

For a modular approach, `mergepatch` was implemented to generate a whole genome BAM containing the artificial CNV from the final.patch.bam and the patch.qnames.txt generated from `cnvinject del` and `cnvinject dup` operations. `mergepatch` removes the patch reads from the `full-bam` and the standard out is merged with the `patch-bam`. 

### Arguments

| Argument | Required | Default | Description |
|---|---:|---:|---|
| `--full-bam` | Yes | none | Original full BAM. |
| `--patch-bam` | Yes | none | Edited patch BAM produced by a CNVinject patch-generation command. |
| `--patch-reads` | Yes | none | Text file containing line-separated patch read names. This file is an output with ending with patch.qnames.txt |
| `-o`, `--output` | Yes | none | Output bam file name including path. Example: `-o ~/project/output/final.bam` |
| `-t`, `--threads` | No | `1` | Number of threads. |


### Expected outputs

| Output | Description |
|---|---|
| Output BAM | Full BAM with the artificial CNV injected. |
| Output BAM index | BAM index, if indexing is performed by the workflow. |
| Logs or progress messages | Run information useful for troubleshooting. |

A typical final output are coordinate sorted and indexed.

```text
Sample2.NF1.CN0.final.bam
Sample2.NF1.CN0.final.bam.bai
```



## Recommended validation after cnv injection

After generating a CNV-injected BAM, validate that the expected copy-number change is present. Users should inspect the output using a genome browser and/or coverage summaries. For deletion simulations, users should observe reduced coverage across the target interval relative to nearby background regions and/or relative to the original BAM.

BAM file integrity checks include:

```bash
samtools quickcheck -v NF1.CN0.final.patch.bam
```

Coverage validation with `samtools coverage`:

```bash
samtools coverage NF1.CN0.final.patch.bam
```

---

## Read eligibility
Currently, only primary alignments are eligible for modification. Secondary and supplemental alignments will not be touched since their alignment to the target interval is not certain and therefore may be interpreted as noise or atifact. 

`--mapq` controls which reads are eligible for CNV editing. Users can opt to retain low mapq reads in the target interval at their original concentration to simualte low quality/uncertain alignments in the target interval. For noisy whole-genome-amplified or sparse libraries, the appropriate MAPQ threshold depends on the benchmarking goal. If the goal is to preserve the full noise profile of the original BAM, `--mapq 10` is preferable and is the default.


In whole-genome-amplified or ultra-low-input sequencing libraries, duplicates may be part of the observed read-depth structure. Removing all duplicates before simulation can make the injected CNV less representative of the original data. CNVinject is designed to preserve duplicated reads. After a CNV-injected bam is generated, it is recommended to remove and/or remark duplicated reads if desired. 

---

## Troubleshooting



## Development notes

CNVinject is currently in early development. The interface and output filenames may change.


## Citation

No formal citation is available yet. If you use CNVinject in a publication or internal benchmark, please cite the GitHub repository.

