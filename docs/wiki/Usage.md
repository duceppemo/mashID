# Usage

```
usage: mashID [-h] -i PATH -o DIR [-d FILE.msh|NAME] [--db-metadata FILE.tsv]
              [--identity 0.9] [--p-value 0.05] [-n 10]
              [-s {identity,multiplicity}] [--no-winner-take-all]
              [--max-reads N] [--ambiguity-margin 0.005] [--min-ref-length BP]
              [--skip-stats] [-t 64] [-p 2] [--debug] [-v]
```

## Options

### input/output

| Option | Description |
| --- | --- |
| `-i, --input PATH` | Input directory (searched recursively) with fastq/fasta files, or a single file, gzipped or not. Paired-end files (R1/R2) are screened together. |
| `--sample-sheet FILE.tsv` | Instead of `-i`: a TSV/CSV with columns `sample` and `file` (one row per file, or files separated by `;`). Relative paths are resolved from the sheet's directory. Other columns are copied into the summary. See below. |
| `-o, --output DIR` | Output directory, created if needed. |
| `-d, --database FILE.msh\|NAME` | Mash sketch database: a `.msh` file, or the name of a downloaded pre-built database (`mycobacteriaceae`, `mycobacteriaceae-2025`, `listeria`, `listeria-2025`, `salmonella`, `escherichia`, `progenomes4`, `progenomes3`, `refseq_bacteria`). Default: `mycobacteriaceae`. |
| `--db-metadata FILE.tsv` | Table mapping reference accessions to organism names and TaxIDs. Default: the `<database>.metadata.tsv` sidecar. |

### screening

| Option | Description |
| --- | --- |
| `--identity 0.9` | Minimum identity to report, between 0 and 1. |
| `--p-value 0.05` | Maximum p-value to report. |
| `-n, --n-hits 10` | Number of top hits to report per sample. |
| `-s, --sort-by {identity,multiplicity}` | How to rank hits; determines the top hit in the summary. |
| `--no-winner-take-all` | Disable Mash's winner-take-all strategy (`-w`). Reports more redundant hits and disables the mixture note. |
| `--max-reads N` | Screen only the first N reads of each fastq sample, split across R1/R2 and streamed to Mash. `Sequences`/`Bases` then describe the screened reads. |
| `--ambiguity-margin 0.005` | Flag a sample as ambiguous when a different organism scores within this identity margin of the top hit. |
| `--min-ref-length 100000` | Ignore hits to references shorter than this many bp. `0` keeps all hits. |
| `--skip-stats` | Do not count reads/bases of the input files. |
| `--fail-on {none,no-hit,note}` | Exit with code 2 when any sample has no hit, or has any note or no hit. For pipelines. |

### performance

| Option | Description |
| --- | --- |
| `-t, --threads N` | Total threads, shared between parallel samples. Default: all CPUs. |
| `-p, --parallel 2` | Number of samples processed concurrently. |

`--debug` prints every external command; `-v` prints the version.

## Examples

For a runnable worked example with its expected output, see [Example](Example).

```bash
# Default (downloaded) Mycobacteriaceae database, 4 samples at a time
mashID -i /data/run42/fastq -o /data/run42/mashID -t 16 -p 4

# A downloaded database by name
mashID -i assemblies/ -o results -d listeria

# Custom database, quick check on the first 200k reads of each sample
mashID -i /data/run42/fastq -o /data/run42/mashID_quick -d /db/my_db.msh --max-reads 200000

# Rank by multiplicity to pick the dominant organism in a mixed sample
mashID -i mixed.fastq.gz -o results -s multiplicity
```

## Input files and sample names

Accepted extensions: `.fastq`, `.fq`, `.fasta`, `.fa`, `.fna`, optionally followed by `.gz`.
Directories are searched recursively; hidden files and directories are ignored, and symlinked
directories are followed once.

The sample name is the file name without its extension and without Illumina read designations
(`_S1_L001_R1_001`, `_R1`, `_1`). For example:

| File | Sample |
| --- | --- |
| `MBWGS440_S10_L001_R1_001.fastq.gz` and `..._R2_001.fastq.gz` | `MBWGS440` (one sample, both files screened) |
| `isolate_A_R1.fq.gz` | `isolate_A` |
| `GCF_000195955.2.fna` | `GCF_000195955.2` |
| `barcode01_pass.fastq.gz` | `barcode01_pass` |

Files from several lanes of the same sample are grouped the same way. A sample cannot mix fasta and
fastq files. Nothing is concatenated or copied: Mash reads all files of a sample directly.

## Sample sheets

When file names do not encode sample names cleanly (Nanopore barcodes, files from several runs, LIMS
identifiers), give `--sample-sheet` instead of `-i`:

```
sample     file                              expected
ISO-2024-1 runA/barcode01_pass.fastq.gz      Mycobacterium bovis
ISO-2024-2 runA/barcode02_pass.fastq.gz      Mycobacterium bovis
ISO-2024-3 runB/S3_R1.fastq.gz;runB/S3_R2.fastq.gz  Listeria monocytogenes
```

Tab- or comma-separated, one row per file or per sample with files joined by `;`. Paths are relative to
the sheet. Any other column (`expected` above) is carried into `summary_mashID.tsv` and
`summary_mashID.json`, which makes a "matches expectation" check a one-liner downstream.

## Quick screening of large runs

Mash screen reads every base of the input, so a 5 GB Nanopore run takes minutes. `--max-reads` streams
only the first N reads to Mash and typically gives the same identification in seconds:

| Input | Command | Time |
| --- | --- | --- |
| 250 MB Illumina fastq.gz | full file | 13 s |
| same file | `--max-reads 200000` | 3 s |

For identification, 100k–500k reads are plenty. Coverage-related columns (`Median_Multiplicity`) then
reflect the subsample.
