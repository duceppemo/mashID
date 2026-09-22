# Databases

Any Mash sketch database works with mashID. Databases are stored outside the package and outside the
repository.

## Pre-built databases

```bash
mashID_download_db                    # list databases and whether they are installed
mashID_download_db mycobacteriaceae   # download one (MD5-verified)
mashID_download_db --all --dir /db    # download everything somewhere else
```

| Name | Content | Size |
| --- | --- | --- |
| `mycobacteriaceae` (default) | Mycobacteriaceae (NCBI taxon 1762): 16,189 GenBank assemblies as of 2026-09-22, dereplicated per species at 99.9% identity to 3534 references, NCBI TaxIDs for all. k=21, s=10000. Ships with its metadata sidecar. [Figshare](https://doi.org/10.6084/m9.figshare.33965176) | 284 MB |
| `mycobacteriaceae-2025` | The previous Mycobacteriaceae database (2025-02-20): 3309 references, k=21, s=1000, no TaxIDs. Kept for reproducibility. [Figshare](https://doi.org/10.6084/m9.figshare.28489304.v1) | 27 MB |
| `listeria` | *Listeria* (NCBI taxon 1637): 7453 RefSeq assemblies as of 2026-09-22, dereplicated per species at 99.9% identity to 1441 references covering 42 species and subspecies, NCBI TaxIDs for all. k=21, s=10000. Ships with its metadata sidecar. [Figshare](https://doi.org/10.6084/m9.figshare.33968620) | 116 MB |
| `listeria-2025` | The previous *Listeria* database (2025-02-18), no TaxIDs. Kept for reproducibility. [Figshare](https://doi.org/10.6084/m9.figshare.28489262.v1) | 46 MB |
| `progenomes3` | proGenomes v3 representative genomes, all bacteria and archaea. [Figshare](https://doi.org/10.6084/m9.figshare.22312282.v1) | 331 MB |
| `refseq_bacteria` | RefSeq bacteria as of 2023-01-19, dereplicated at 99% identity. [Figshare](https://doi.org/10.6084/m9.figshare.22312240.v2) | 653 MB |

Databases live in `$MASHID_DB_DIR` if set, else `$XDG_DATA_HOME/mashID/db`
(`~/.local/share/mashID/db`). Pass a name (`-d listeria`) or a path (`-d /db/x.msh`) to `mashID`.
A metadata sidecar (see below) is generated right after download.

## Building your own

```bash
make_mashID_db -i /path/to/genomes -o /path/to/db -p my_database -s 10000 -k 21 -t 16
# -> /path/to/db/my_database.msh and /path/to/db/my_database.metadata.tsv
```

| Option | Description |
| --- | --- |
| `-i PATH` | Directory searched recursively for `.fna/.fa/.fasta[.gz]`, or a text file listing one genome path per line. |
| `-o DIR`, `-p NAME` | Output directory and database name (`.msh` is appended). |
| `-s 10000` | Sketch size. 10000 is recommended for species-level identification. |
| `-k 21` | K-mer size (1–32). |
| `--min-length 100000` | Exclude references shorter than this; use `0` for plasmid or viral collections. |
| `--metadata FILE` | TSV/CSV with `Accession`, `Organism` and optionally `TaxID` columns. |
| `--assembly-report FILE.jsonl` | `assembly_data_report.jsonl` from an NCBI `datasets download genome` archive. |
| `--annotate DB.msh` | Only (re)write the metadata sidecar of an existing database. |
| `--check DB.msh` | Report quality issues of a database: references shorter than `--min-length`, unparsed organism names, missing TaxIDs, species listed under several genera (synonyms), single-reference organisms. |

`scripts/build_mycobacteriaceae_db.sh` in the repository shows a full workflow: download all NCBI
genomes of a taxon with `datasets`, dereplicate them per species with
[Assembly-dereplicator](https://github.com/rrwick/Assembly-dereplicator), and sketch the database.

## Checking a database

```bash
make_mashID_db --check /path/to/db.msh
```

Run this on any database before trusting it. On the 2025 Mycobacteriaceae database (`mycobacteriaceae-2025`) it reports seven
references under 100 kb (partial records that would otherwise become false top hits), no TaxIDs, and
40 species that appear under two or three genera because of the *Mycobacterium* →
*Mycolicibacterium* / *Mycobacteroides* / *Mycolicibacter* reclassification. The last point matters:
mashID counts *Mycobacterium abscessus* and *Mycobacteroides abscessus* as different organisms, so a
sample may show a spurious "Possible mixture" between the two names. A sidecar built from NCBI's
assembly report (`--annotate --assembly-report`) or a curated `--metadata` table harmonises them.

## Provenance and known issues of the pre-built databases

| Database | Built | Source and filters | Known issues |
| --- | --- | --- | --- |
| `mycobacteriaceae` | 2026-09-22 | `scripts/build_mycobacteriaceae_db.sh`: all 16,189 GenBank assemblies of taxon 1762 (atypical excluded) via NCBI Datasets, binned by species from the assembly report, dereplicated per species at Mash distance 0.001 with Assembly-dereplicator 0.3.2 (*M. tuberculosis*: 211 of 8722 kept), references < 100 kb excluded, sketched k=21, s=10000. 3534 references; names and TaxIDs from the assembly report, strain text removed. | `--check` reports one species under two genus spellings (`[Mycobacterium] chelonae`, NCBI's own naming). MTBC members are not separable (see Interpreting results). |
| `mycobacteriaceae-2025` | 2025-02-20 | All NCBI assemblies of taxon 1762 from the datasets web table (26,101), renamed and binned by the species in the first header, dereplicated per species at Mash distance 0.001 with Assembly-dereplicator 0.3.2, sketched k=21, s=1000. 3309 references. | Seven references < 100 kb; names parsed from headers, so genus synonyms are not merged; no TaxIDs; s=1000 gives identities in steps of 0.001. Superseded by `mycobacteriaceae`. |
| `listeria` | 2026-09-22 | `scripts/build_mycobacteriaceae_db.sh` with `ASSEMBLY_SOURCE=RefSeq` and taxon 1637: all 7453 RefSeq assemblies (atypical excluded; GenBank's 79,000 are mostly redundant *L. monocytogenes*), binned by species from the assembly report, dereplicated per species at Mash distance 0.001 (*L. monocytogenes*: 892 kept), references < 100 kb excluded, sketched k=21, s=10000. 1441 references with NCBI names and TaxIDs. | `--check` reports no issues. Validated on Nanopore *L. ivanovii* subsp. *londoniensis* and *L. monocytogenes* reads and an assembly. |
| `listeria-2025` | 2025-02-18 | All NCBI *Listeria* assemblies, dereplicated as above. | Names parsed from headers; no TaxIDs. Superseded by `listeria`. |
| `progenomes3` | 2023-03 | proGenomes v3 representative genomes, sketched from the proGenomes fasta. | Header format differs from NCBI's; names come from the sidecar generated on first use, check them with `--check`. |
| `refseq_bacteria` | 2023-01-19 | RefSeq bacteria, dereplicated at 0.01. | Large (653 MB); three years old. |

## The metadata sidecar

`<database>.metadata.tsv` sits next to the `.msh` file and holds, for every reference:
`Accession`, `Organism`, `TaxID`, `Length`, `Hashes`, `Source_File`, `Description`.

It is the source of `Identification` and `TaxID` in the outputs, and of the reference lengths used by
`--min-ref-length`. Organism names come from, in order of preference:

1. `--metadata FILE`: your own table.
2. `--assembly-report`: NCBI's report, which carries organism names and TaxIDs for every accession.
3. The first fasta header of each genome, parsed for genus, species and infraspecific ranks
   (`NZ_CP060409.1 Mycolicibacterium fortuitum strain W4 chromosome` → *Mycolicibacterium fortuitum*).

Accessions are recognised in file names (`GCF_000195955.2_ASM19595v2_genomic.fna.gz` →
`GCF_000195955.2`). If a database has no sidecar, `mashID` builds one on first use from `mash info`,
and writes it next to the database when the directory is writable. To add or refresh a sidecar,
including for downloaded databases:

```bash
make_mashID_db --annotate /path/to/db.msh --assembly-report assembly_data_report.jsonl
```

## Choosing a sketch size

The pre-built Mycobacteriaceae database uses s=1000, which resolves species but yields coarse
identities (steps of 0.001). For a custom database aimed at close relatives, s=10000 is worth the
tenfold larger file: identities become finer and the shared-hash counts more informative.
