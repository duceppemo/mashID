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
| `salmonella` | *Salmonella* (NCBI taxon 590): 21,133 RefSeq assemblies as of 2026-09-22, binned by subspecies and dereplicated to 1442 references covering *S. bongori* and all six *S. enterica* subspecies, NCBI TaxIDs for all. k=21, s=10000. Ships with its metadata sidecar. Species and subspecies level: the serovar in a reference name is that of the nearest reference, not a serotype determination; use SeqSero2 or SISTR for serotyping. [Figshare](https://doi.org/10.6084/m9.figshare.33968917) | 116 MB |
| `brucella` | *Brucella* (NCBI taxon 234, including the former *Ochrobactrum* species): all 1995 GenBank assemblies as of 2026-09-22, dereplicated per species at Mash distance 0.0003 to 387 references covering 29 species, NCBI TaxIDs for all. k=21, s=10000. Ships with its metadata sidecar. Held-out concordance 93.4%; the classical species are host-adapted lineages 0.2 to 0.5% apart, so their calls are nearest-reference calls and biovars are not resolved. [Figshare](https://doi.org/10.6084/m9.figshare.33970168) | 31 MB |
| `escherichia` | *Escherichia* and *Shigella* (NCBI taxa 561 and 620): 57,146 RefSeq assemblies as of 2026-09-22, dereplicated to 4503 references (*E. coli* 3316, five other *Escherichia* species, four *Shigella* species), NCBI TaxIDs for all. k=21, s=10000. Ships with its metadata sidecar. Held-out concordance 96.4% at species level; *S. boydii*, *S. dysenteriae* and *S. flexneri* share lineages and are not reliably separated from each other, nor enteroinvasive *E. coli* from *Shigella*: use ShigaTyper or ShigEiFinder. [Figshare](https://doi.org/10.6084/m9.figshare.33969760) | 362 MB |
| `progenomes4` | proGenomes 4 representative genomes: 32,887 species-cluster representatives across 4384 genera of bacteria and archaea, current NCBI names and TaxIDs for all. k=21, s=2000. Ships with its metadata sidecar. For species-level screening of anything. [Figshare](https://doi.org/10.6084/m9.figshare.33968806) | 536 MB |
| `progenomes3` | proGenomes v3 representative genomes (2023), names parsed from headers. Superseded by `progenomes4`. [Figshare](https://doi.org/10.6084/m9.figshare.22312282.v1) | 331 MB |
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

`scripts/build_taxon_db.sh` in the repository shows a full workflow: download all NCBI
genomes of a taxon with `datasets`, dereplicate them per species with
[Assembly-dereplicator](https://github.com/rrwick/Assembly-dereplicator), and sketch the database.

## Validating a database on held-out genomes

`scripts/validate_db_heldout.py <build_dir> <db.msh>` screens genomes that the build downloaded but
did not keep (dropped by the bin cap or by dereplication) and reports, per species, how often mashID
returns the NCBI name, how often a note fired, and every mismatch. Run it after each build; the
numbers in the provenance table above come from it.

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
| `mycobacteriaceae` | 2026-09-22 | `scripts/build_taxon_db.sh`: all 16,189 GenBank assemblies of taxon 1762 (atypical excluded) via NCBI Datasets, binned by species from the assembly report, dereplicated per species at Mash distance 0.001 with Assembly-dereplicator 0.3.2 (*M. tuberculosis*: 211 of 8722 kept), references < 100 kb excluded, sketched k=21, s=10000. 3534 references; names and TaxIDs from the assembly report, strain text removed. | `--check` reports one species under two genus spellings (`[Mycobacterium] chelonae`, NCBI's own naming). MTBC members are not separable (see Interpreting results). |
| `mycobacteriaceae-2025` | 2025-02-20 | All NCBI assemblies of taxon 1762 from the datasets web table (26,101), renamed and binned by the species in the first header, dereplicated per species at Mash distance 0.001 with Assembly-dereplicator 0.3.2, sketched k=21, s=1000. 3309 references. | Seven references < 100 kb; names parsed from headers, so genus synonyms are not merged; no TaxIDs; s=1000 gives identities in steps of 0.001. Superseded by `mycobacteriaceae`. |
| `listeria` | 2026-09-22 | `scripts/build_taxon_db.sh` with `ASSEMBLY_SOURCE=RefSeq` and taxon 1637: all 7453 RefSeq assemblies (atypical excluded; GenBank's 79,000 are mostly redundant *L. monocytogenes*), binned by species from the assembly report, dereplicated per species at Mash distance 0.001 (*L. monocytogenes*: 892 kept), references < 100 kb excluded, sketched k=21, s=10000. 1441 references with NCBI names and TaxIDs. | `--check` reports no issues. Validated on Nanopore *L. ivanovii* subsp. *londoniensis* and *L. monocytogenes* reads and an assembly. |
| `listeria-2025` | 2025-02-18 | All NCBI *Listeria* assemblies, dereplicated as above. | Names parsed from headers; no TaxIDs. Superseded by `listeria`. |
| `salmonella` | 2026-09-22 | `scripts/build_taxon_db.sh` with `ASSEMBLY_SOURCE=RefSeq BIN_RANK=subspecies MAX_BIN=3000` and taxon 590: all 21,133 RefSeq assemblies (GenBank's 630,000 are mostly redundant), binned by subspecies, each bin capped at 3000 keeping complete genomes first, dereplicated per bin at Mash distance 0.001, assemblies named only "Salmonella sp." excluded, sketched k=21, s=10000. 1442 references with NCBI names and TaxIDs. | Nearest-reference serovars appear in names but Mash does not serotype: serovars are defined by three loci, several are polyphyletic, and rare ones have few references. `--check` reports no issues. Validated on SKESA assemblies and MiSeq reads of *S. enterica*. |
| `brucella` | 2026-09-22 | `scripts/build_taxon_db.sh` with `ASSEMBLY_SOURCE=GenBank EXCLUDE_SP=1 DEREP_DISTANCE=0.0003` and taxon 234: all 1995 GenBank assemblies, binned by species, dereplicated per species at Mash distance 0.0003 (0.001 left only 21 *B. melitensis* and 9 *B. suis* references and scored 90.7% held-out; the classical species are only 0.2 to 0.5% apart), "sp." bin excluded, sketched k=21, s=10000. 387 references. | Held-out validation on 226 genomes: 93.4%; *B. abortus* 29/30, *B. melitensis* 29/30, *B. suis* 28/30, *B. canis*, *B. ovis*, *B. ceti*, *B. pinnipedialis*, *B. microti*, *B. neotomae* 100%. Misses match another species at > 0.9995, including a submission batch deposited as *B. intermedia* identical to *B. ciceri*. For genera whose species are this close, set `DEREP_DISTANCE` below the inter-species distance. |
| `escherichia` | 2026-09-22 | `scripts/build_taxon_db.sh` with `ASSEMBLY_SOURCE=RefSeq MAX_BIN=5000 EXCLUDE_SP=1` and taxa 561,620: all RefSeq *Escherichia* (53,845) and *Shigella* (3,301) assemblies, binned by species, each bin capped at 5000 keeping complete genomes first (for *E. coli* that is the 5595 complete genomes), dereplicated per bin at Mash distance 0.001, "sp." bins excluded, sketched k=21, s=10000. | Validated with `scripts/validate_db_heldout.py` on 251 genomes left out of the database: 96.4% species concordance, 100% for every *Escherichia* species but *E. coli* (28/30, both misses matching *S. sonnei* at 0.99); the other misses are *Shigella* genomes matching another *Shigella* species at > 0.999 identity, a labelling matter Mash cannot resolve. Complete genomes over-represent clinical and outbreak lineages. |
| `progenomes4` | 2026-09-22 | `scripts/build_progenomes_db.sh`: the proGenomes 4 representatives multi-fasta (45.9 GB) split per assembly accession, names and TaxIDs fetched from NCBI Datasets for every accession, sketched k=21, s=2000, no dereplication (one genome per species cluster already). | Single representative per species, so identities of 0.95 to 0.99 are normal; use a group-specific database for finer resolution. `--check` reports 2230 epithets shared across genera and 91 references named "bacterium" or similar: expected at this scale, not synonyms. |
| `progenomes3` | 2023-03 | proGenomes v3 representative genomes, sketched from the proGenomes fasta. | Header format differs from NCBI's; names come from the sidecar generated on first use, check them with `--check`. Superseded by `progenomes4`. |
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
