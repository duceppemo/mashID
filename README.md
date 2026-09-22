# mashID

Identify organisms from genome assemblies (fasta) or raw sequencing reads (fastq) with
[Mash](https://github.com/marbl/Mash).

mashID screens each sample against a Mash sketch database (`mash screen`), keeps the best hits and
reports a clean organism name (genus, species, subspecies/variant) parsed from the reference
headers. It works with Illumina (paired-end or single-end), Ion Torrent, Nanopore and PacBio reads,
and with assemblies. Inputs may be gzipped.

Pre-built databases (Mycobacteriaceae, *Listeria*, proGenomes v3, RefSeq bacteria) are downloaded on
demand with `mashID_download_db`, and any Mash sketch database built with `make_mashID_db` can be used.

## Installation

Requirements: Python ≥ 3.10 and Mash 2.3. No Python dependencies.

```bash
# 1. Environment with Mash (bioconda's mash 2.3 needs GSL 2.7.0: newer GSL builds break it)
conda create -n mashID -c conda-forge -c bioconda python=3.12 mash=2.3 "gsl==2.7" python-isal
conda activate mashID

# 2. Install mashID
git clone https://github.com/duceppemo/mashID
cd mashID
pip install .

# 3. Download the default database (Mycobacteriaceae, 27 MB)
mashID_download_db mycobacteriaceae

mashID -h
```

`python-isal` is optional and only speeds up gzip decompression when counting reads. You can also
create the environment from the repository's `environment.yml`.

Running from the source tree without installing still works: `python mashID.py -h`.

## Usage

```
usage: mashID [-h] -i PATH -o DIR [-d FILE.msh|NAME] [--db-metadata FILE.tsv]
              [--identity 0.9] [--p-value 0.05] [-n 10]
              [-s {identity,multiplicity}] [--no-winner-take-all]
              [--max-reads N] [--ambiguity-margin 0.005] [--skip-stats]
              [-t 64] [-p 2] [--debug] [-v]

input/output:
  -i, --input PATH      Input directory (searched recursively) with fastq/fasta files, or a single
                        fastq/fasta file, gzipped or not. Paired-end files (R1/R2) are screened together.
  -o, --output DIR      Output directory (created if needed).
  -d, --database FILE.msh|NAME
                        Mash sketch database: a .msh file, or the name of a downloaded pre-built
                        database (mycobacteriaceae, listeria, progenomes3, refseq_bacteria; see
                        mashID_download_db). Default: mycobacteriaceae
  --db-metadata FILE.tsv
                        Table mapping reference accessions to organism names (and TaxIDs). Default:
                        the <database>.metadata.tsv sidecar written by make_mashID_db, if present.

screening:
  --identity 0.9        Minimum identity to report, between 0 and 1. Default: 0.9
  --p-value 0.05        Maximum p-value to report. Default: 0.05
  -n, --n-hits 10       Number of top hits to report per sample. Default: 10
  -s, --sort-by {identity,multiplicity}
                        How to rank hits; determines the "top hit" in the summary. Default: identity
  --no-winner-take-all  Disable Mash's winner-take-all strategy (-w). Reports more redundant hits.
  --max-reads N         Screen only the first N reads of each fastq sample (split across R1/R2),
                        streamed to Mash. Much faster on large runs; Sequences/Bases then describe
                        the screened reads.
  --ambiguity-margin 0.005
                        Flag a sample as ambiguous when a different organism scores within this
                        identity margin of the top hit. Default: 0.005
  --min-ref-length BP   Ignore hits to references shorter than this many bp (partial records in a
                        database otherwise produce false top hits). 0 keeps all hits. Default: 100000
  --skip-stats          Do not count reads/bases of the input files (faster on very large fastq).

performance:
  -t, --threads 64      Total number of threads, shared between parallel samples. Default: all
  -p, --parallel 2      Number of samples to process concurrently. Default: 2
```

Examples:

```bash
# Default (downloaded) Mycobacteriaceae database, 4 samples at a time
mashID -i /data/run42/fastq -o /data/run42/mashID -t 16 -p 4

# Custom database, quick check on the first 200k reads of each sample
mashID -i /data/run42/fastq -o /data/run42/mashID_quick -d /db/my_db.msh --max-reads 200000
```

### Databases

```bash
mashID_download_db                    # list databases and whether they are installed
mashID_download_db listeria           # download one (MD5-verified) into the database directory
mashID_download_db --all --dir /db    # download everything somewhere else
```

Databases live in `$MASHID_DB_DIR` if set, else `$XDG_DATA_HOME/mashID/db`
(`~/.local/share/mashID/db`). Pass a name (`-d listeria`) or a path (`-d /db/x.msh`) to `mashID`.

### Input files and sample names

Accepted extensions: `.fastq`, `.fq`, `.fasta`, `.fa`, `.fna`, optionally followed by `.gz`.
Directories are searched recursively; hidden files are ignored.

The sample name is the file name without its extension and without Illumina read designations, so
`MBWGS440_S10_L001_R1_001.fastq.gz` and `MBWGS440_S10_L001_R2_001.fastq.gz` both become sample
`MBWGS440` and are screened together (Mash reads both files directly; nothing is concatenated or
copied). Files from several lanes of the same sample are grouped the same way. A sample cannot mix
fasta and fastq files.

## Outputs

- `<sample>_mashID.tsv`: the top `-n` hits for each sample, with columns
  `Rank`, `Identity`, `Shared_Hashes`, `Median_Multiplicity`, `P_Value`, `Accession`, `TaxID`,
  `Identification` (organism name) and `Description` (the full reference header).
- `summary_mashID.tsv`: one line per sample with the number of sequences (reads or contigs), total
  bases, the best hit and a `Note`. Samples without a hit above the thresholds are reported as
  `No significant hit in database`. The same table is printed to the terminal.

The `Note` column flags results that deserve a second look:

| Note | Meaning |
| --- | --- |
| `Possible mixture with: X` | With winner-take-all (default), another organism still scores ≥ 0.99 identity. A complete genome with ≥ 95% ANI to the top hit cannot do that once shared k-mers are credited to the winner, so X is supported by its own k-mers: expect contamination or a mixed culture. A "may be partial" hint marks references shorter than half the top hit's, whose k-mers can escape winner-take-all. |
| `Ambiguous: X at identity ...` | The best hit of a different organism is within `--ambiguity-margin` of the top hit. The database cannot separate them for this sample. |
| `Low coverage: median multiplicity N` | Reads only: the top hit's k-mers were seen fewer than 5 times on average. Identification is plausible but weakly supported. |

Organism names, TaxIDs and reference lengths come from the database's `<db>.metadata.tsv` sidecar.
`make_mashID_db` writes it, `mashID_download_db` generates it for downloaded databases, and `mashID`
creates it on first use for any other database (names parsed from the reference headers).

Hits to references shorter than `--min-ref-length` (100 kb) are ignored. Public genome collections
contain partial records, and a 1 kb "genome" is fully contained in any related sample, which makes it
the top hit at identity 1.0. The 2025 Mycobacteriaceae database, for example, contains nine such
records; without the filter an *M. bovis* sample is reported as *M. tuberculosis* from a 947-k-mer
fragment.

Within very close groups such as the *Mycobacterium tuberculosis* complex, Mash identities of the
members differ in the fourth decimal and the top hit is not a reliable variant call. Use mashID for
species-level identification and a SNP-based method for variant or lineage assignment.

`Identity` is Mash's containment estimate (0-1). `Median_Multiplicity` is the median number of times
the shared hashes were seen in the sample and approximates coverage for reads (it is 1 for
assemblies). With `-s multiplicity` hits are ranked by multiplicity, which can help pick the dominant
organism in a mixed sample.

## Building a custom database

```bash
make_mashID_db -i /path/to/genomes -o /path/to/db -p my_database -s 10000 -k 21 -t 16
# -> /path/to/db/my_database.msh and /path/to/db/my_database.metadata.tsv
```

`-i` accepts a directory (searched recursively for `.fna/.fa/.fasta[.gz]`) or a text file listing one
genome path per line. A sketch size of 10000 is recommended for species-level identification.

Organism names for the sidecar come from, in order of preference:

1. `--metadata FILE`: a TSV/CSV with columns `Accession`, `Organism` and optionally `TaxID`.
2. `--assembly-report assembly_data_report.jsonl`: the report inside every
   `datasets download genome` archive from NCBI, which carries organism names and TaxIDs.
3. The first fasta header of each file (NCBI-style headers such as
   `NZ_CP060409.1 Mycolicibacterium fortuitum strain W4 chromosome` parse well).

References shorter than `--min-length` (100 kb) are excluded with a warning; use `--min-length 0`
for plasmid or viral collections. Accessions are recognised in file names
(`GCF_000195955.2_ASM19595v2_genomic.fna.gz`). To add or refresh the sidecar of an existing database,
including the downloaded ones, without re-sketching:

```bash
make_mashID_db --annotate /path/to/db.msh --assembly-report assembly_data_report.jsonl
```

Run `make_mashID_db -h` for all options.

`scripts/build_mycobacteriaceae_db.sh` shows a full workflow: download all NCBI genomes of a taxon
with `datasets`, dereplicate them per species, and sketch the database.

### Pre-built databases

1. [Mycobacteriaceae (2025-02-20)](https://figshare.com/articles/dataset/Mycobacteriaceae_database_for_mashID_-_2025-02-20_update/28489304?file=52609139)
2. [Listeria spp. (2025-02-18)](https://figshare.com/articles/dataset/Listeria_database_for_mashID_-_2025-02-18_update/28489262?file=52609082)
3. [proGenomes v3](https://figshare.com/articles/dataset/progenomes3_msh/22312282)
4. [RefSeq bacteria, dereplicated at 0.01 (2023-01-19)](https://figshare.com/articles/dataset/Untitled_Item/22312240?file=39690568)

## Development

```bash
pip install -e ".[test]"
ruff check .
pytest
```

The integration tests build a small synthetic database and are skipped when `mash` is not on `PATH`.
GitHub Actions runs the same checks on every push (`.github/workflows/ci.yml`). A bioconda recipe
template is in `recipe/`.

## License

MIT. See `LICENSE`.
