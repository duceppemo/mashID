# Outputs

The output directory holds the summary in three forms, one table per sample, and a provenance file.
The summary table is also printed.

| File | Content |
| --- | --- |
| `summary_mashID.tsv` | One line per sample with its best hit and notes. |
| `summary_mashID.json` | The same plus every reported hit per sample, for pipelines. |
| `summary_mashID_mqc.tsv` | MultiQC custom-content table; `multiqc <dir>` picks it up. |
| `<sample>_mashID.tsv` | The top `-n` hits of each sample. |
| `mashID_run.json` | Provenance: mashID and Mash versions, command line, parameters, database path, size, MD5 and reference count, input files, timestamps, exit code. |

## `summary_mashID.tsv`

One line per sample with its best hit.

| Column | Meaning |
| --- | --- |
| `Sample` | Sample name derived from the file name(s). |
| `Sequences` | Number of reads (fastq) or contigs (fasta) across all files of the sample. `NA` with `--skip-stats`. |
| `Bases` | Total bases. |
| `Identity` | Mash containment estimate of the reference in the sample, between 0 and 1. |
| `Shared_Hashes` | Hashes of the reference sketch found in the sample, e.g. `994/1000`. |
| `Median_Multiplicity` | Median number of times the shared hashes were seen. Approximates k-mer coverage for reads; 1 for assemblies. |
| `Est_Depth` | Reads only: total bases divided by the top reference's length, i.e. the sequencing depth if the sample is that organism. `NA` for assemblies or when the reference length is unknown. |
| `P_Value` | Probability of the observed sharing by chance. |
| `Accession` | Reference accession, parsed from the reference file name (`GCF_…`/`GCA_…` recognised). |
| `TaxID` | NCBI TaxID from the database metadata, or `NA`. |
| `Identification` | Organism name from the database metadata, or `No significant hit in database`. |
| `Note` | Warnings; see below. |

Columns from a [sample sheet](Usage#sample-sheets) are appended after `Note`.

## `<sample>_mashID.tsv`

The top `-n` hits of a sample, ranked by `--sort-by`, with columns
`Rank`, `Identity`, `Shared_Hashes`, `Median_Multiplicity`, `P_Value`, `Accession`, `TaxID`,
`Identification` and `Description` (the reference's full fasta header). A sample with no hit has a
header-only file.

## The `Note` column

| Note | Meaning |
| --- | --- |
| `Possible mixture with: X` | With winner-take-all (default), a different organism still scores ≥ 0.99 identity. A complete genome with ≥ 95% ANI to the top hit cannot do that once shared k-mers are credited to the winner, so X is supported by its own k-mers: expect contamination or a mixed culture. A "may be partial" hint marks references shorter than half the top hit's, whose k-mers can escape winner-take-all. |
| `Ambiguous: X at identity …` | The best hit of a different organism is within `--ambiguity-margin` of the top hit. The database cannot separate them for this sample. |
| `Low coverage: median multiplicity N` | Reads only: the top hit's k-mers were seen fewer than 5 times on average. Identification is plausible but weakly supported. |

The mixture note is disabled with `--no-winner-take-all`, because without winner-take-all every
related reference scores high.

## Ignored references

Hits to references shorter than `--min-ref-length` (100 kb) are dropped before ranking, and the log
says how many. Public genome collections contain partial records, and a 1 kb "genome" is fully
contained in any related sample, which makes it the top hit at identity 1.0. The 2025 Mycobacteriaceae
database, for example, contains nine such records; without the filter an *M. bovis* sample is reported
as *M. tuberculosis* from a 947-k-mer fragment.

## What Mash identification can and cannot do

Mash identity is a k-mer containment estimate. It separates species reliably, and often subspecies,
but within very close groups such as the *Mycobacterium tuberculosis* complex the members differ in
the fourth decimal and the top hit is not a trustworthy variant call. Use mashID for species-level
identification and a SNP-based method for variant or lineage assignment.

Ranking by `-s multiplicity` orders hits by how often their k-mers were seen, which helps pick the
dominant organism in a mixed sample rather than the best-matching reference.
