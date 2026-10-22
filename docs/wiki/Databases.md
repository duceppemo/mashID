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
| `refseq_viral` | All 19,625 viral genomes and genome segments of NCBI RefSeq release 237 (6427 phages), current ICTV species names and TaxIDs from the GenBank flat files. k=21, s=2000. Ships with its metadata sidecar. Median reference 6.9 kb, so mashID applies no reference-length filter. For virus identification and for prophage or viral contamination detection in any sample. [Figshare](https://doi.org/10.6084/m9.figshare.34018548) | 304 MB |
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
| `mycobacteriaceae` | 2026-09-22 | `scripts/build_taxon_db.sh`: all 16,189 GenBank assemblies of taxon 1762 (atypical excluded) via NCBI Datasets, binned by species from the assembly report, dereplicated per species at Mash distance 0.001 with Assembly-dereplicator 0.3.2 (*M. tuberculosis*: 211 of 8722 kept), references < 100 kb excluded, sketched k=21, s=10000. 3534 references; names and TaxIDs from the assembly report, strain text removed. | `--check` reports one species under two genus spellings (`[Mycobacterium] chelonae`, NCBI's own naming). MTBC members are not separable (see Interpreting results). Two references are mislabelled at NCBI and are corrected in `curation/label_corrections.tsv`: `GCA_013697125.1` ("Mycobacterium canetti") is *M. tuberculosis* variant *bovis*, and `GCA_025822535.1` ("Mycolicibacterium goodii") is *M. wolinskyi*. Databases downloaded before this was found carry the deposited names; apply the corrections with `scripts/apply_label_corrections.py`. |
| `mycobacteriaceae-2025` | 2025-02-20 | All NCBI assemblies of taxon 1762 from the datasets web table (26,101), renamed and binned by the species in the first header, dereplicated per species at Mash distance 0.001 with Assembly-dereplicator 0.3.2, sketched k=21, s=1000. 3309 references. | Seven references < 100 kb; names parsed from headers, so genus synonyms are not merged; no TaxIDs; s=1000 gives identities in steps of 0.001. Superseded by `mycobacteriaceae`. |
| `listeria` | 2026-09-22 | `scripts/build_taxon_db.sh` with `ASSEMBLY_SOURCE=RefSeq` and taxon 1637: all 7453 RefSeq assemblies (atypical excluded; GenBank's 79,000 are mostly redundant *L. monocytogenes*), binned by species from the assembly report, dereplicated per species at Mash distance 0.001 (*L. monocytogenes*: 892 kept), references < 100 kb excluded, sketched k=21, s=10000. 1441 references with NCBI names and TaxIDs. | `--check` reports no issues. Validated on Nanopore *L. ivanovii* subsp. *londoniensis* and *L. monocytogenes* reads and an assembly. |
| `listeria-2025` | 2025-02-18 | All NCBI *Listeria* assemblies, dereplicated as above. | Names parsed from headers; no TaxIDs. Superseded by `listeria`. |
| `salmonella` | 2026-09-22 | `scripts/build_taxon_db.sh` with `ASSEMBLY_SOURCE=RefSeq BIN_RANK=subspecies MAX_BIN=3000` and taxon 590: all 21,133 RefSeq assemblies (GenBank's 630,000 are mostly redundant), binned by subspecies, each bin capped at 3000 keeping complete genomes first, dereplicated per bin at Mash distance 0.001, assemblies named only "Salmonella sp." excluded, sketched k=21, s=10000. 1442 references with NCBI names and TaxIDs. | Nearest-reference serovars appear in names but Mash does not serotype: serovars are defined by three loci, several are polyphyletic, and rare ones have few references. `--check` reports no issues. Validated on SKESA assemblies and MiSeq reads of *S. enterica*. |
| `brucella` | 2026-09-22 | `scripts/build_taxon_db.sh` with `ASSEMBLY_SOURCE=GenBank EXCLUDE_SP=1 DEREP_DISTANCE=0.0003` and taxon 234: all 1995 GenBank assemblies, binned by species, dereplicated per species at Mash distance 0.0003 (0.001 left only 21 *B. melitensis* and 9 *B. suis* references and scored 90.7% held-out; the classical species are only 0.2 to 0.5% apart), "sp." bin excluded, sketched k=21, s=10000. 387 references. | Held-out validation on 226 genomes: 93.4%; *B. abortus* 29/30, *B. melitensis* 29/30, *B. suis* 28/30, *B. canis*, *B. ovis*, *B. ceti*, *B. pinnipedialis*, *B. microti*, *B. neotomae* 100%. Misses match another species at > 0.9995, including a submission batch deposited as *B. intermedia* identical to *B. ciceri*. For genera whose species are this close, set `DEREP_DISTANCE` below the inter-species distance. |
| `escherichia` | 2026-09-22 | `scripts/build_taxon_db.sh` with `ASSEMBLY_SOURCE=RefSeq MAX_BIN=5000 EXCLUDE_SP=1` and taxa 561,620: all RefSeq *Escherichia* (53,845) and *Shigella* (3,301) assemblies, binned by species, each bin capped at 5000 keeping complete genomes first (for *E. coli* that is the 5595 complete genomes), dereplicated per bin at Mash distance 0.001, "sp." bins excluded, sketched k=21, s=10000. | Validated with `scripts/validate_db_heldout.py` on 251 genomes left out of the database: 96.4% species concordance, 100% for every *Escherichia* species but *E. coli* (28/30, both misses matching *S. sonnei* at 0.99); the other misses are *Shigella* genomes matching another *Shigella* species at > 0.999 identity, a labelling matter Mash cannot resolve. Complete genomes over-represent clinical and outbreak lineages. |
| `refseq_viral` | 2026-09-28 | `scripts/build_virus_db.sh`: the viral division of RefSeq release 237 from the NCBI FTP release files, split into one file per record, names and TaxIDs from the GenBank flat file (`scripts/refseq_gbff_metadata.py`), sketched k=21, s=2000, no length filter, no dereplication (one reference per species or segment). 19,625 references. | Names are RefSeq's ORGANISM field, i.e. ICTV species names ("limnipivirus B1"); the familiar name is in `Description`. Segmented viruses appear as separate segments. Divergent viruses without a close RefSeq relative are not detected at k=21. Tested: RefSeq phages self-match at identity 1; two *S.* Infantis assemblies show a complete Salmonella phage MET_P1_103_31 prophage at 0.9996; *M. bovis* assemblies and reads give no hit. |
| `progenomes4` | 2026-09-22 | `scripts/build_progenomes_db.sh`: the proGenomes 4 representatives multi-fasta (45.9 GB) split per assembly accession, names and TaxIDs fetched from NCBI Datasets for every accession, sketched k=21, s=2000, no dereplication (one genome per species cluster already). | Single representative per species, so identities of 0.95 to 0.99 are normal; use a group-specific database for finer resolution. `--check` reports 2230 epithets shared across genera and 91 references named "bacterium" or similar: expected at this scale, not synonyms. |
| `progenomes3` | 2023-03 | proGenomes v3 representative genomes, sketched from the proGenomes fasta. | Header format differs from NCBI's; names come from the sidecar generated on first use, check them with `--check`. Superseded by `progenomes4`. |
| `refseq_bacteria` | 2023-01-19 | RefSeq bacteria, dereplicated at 0.01. | Large (653 MB); three years old. |

`scripts/build_virus_db.sh` builds a database of every RefSeq viral genome from the NCBI RefSeq release
files, with names and TaxIDs from the GenBank flat files (`scripts/refseq_gbff_metadata.py`).

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

## Organisms the database cannot separate

Whether two organisms can be told apart is a property of the database, not of a sample: the same two
references are the same distance apart in every run. Computing it once makes the per-sample
annotation free, and correct whether or not winner-take-all is in use — unlike `--ambiguity-margin`,
which cannot see a close relative under winner-take-all at all (see Outputs).

```bash
make_mashID_db --confusable /path/to/db.msh -t 16          # writes db.confusable.tsv
make_mashID_db --confusable /path/to/db.msh --confusable-distance 0.002
```

For each reference the table records the closest reference of every *other* organism within
`--confusable-distance` (default 0.005, matching `--ambiguity-margin`), nearest pair first:

```
Accession        Organism                                  Confusable_Organism             Distance  Nearest_Accession
GCA_025821985.1  Mycobacterium timonense                   Mycobacterium bouchedurhonense  0.000041  GCA_025822125.1
GCA_013697125.1  Mycobacterium tuberculosis variant bovis  Mycobacterium tuberculosis      0.000972  GCA_021325765.1
```

mashID picks the sidecar up automatically and notes affected calls; `--db-confusable` points at a
different table. A pre-built database that publishes the table alongside it installs the table too,
so the O(n^2) computation only has to happen once, for whoever built the database. A database
downloaded before its table was published picks the table up on the next
`mashID_download_db <NAME>` for it, fetching only the sidecars and leaving the sketch alone — the
listing calls a database installed on the strength of its `.msh`, so it will say `yes` while the
sidecars are still missing. Runs that finished before the table existed do not have to be screened
again either:

```bash
scripts/reannotate_summary.py /path/to/results --dry-run   # then without --dry-run
```

rebuilds their summary files from the per-sample tables and the current sidecars, which is also how
to apply a corrected organism name to results that predate the correction. Names are compared as
organisms, not as strings, so `Genus sp.` placeholders and the
same epithet under two genera (*Mycobacterium*/*Mycobacteroides abscessus*) are not reported as
conflicts.

A name and a rank below it *are* reported, which is two thirds of the rows for `mycobacteriaceae`
(1306 of 1993). Some of those pairs are the most important thing the table has to say —
*M. tuberculosis* against *M. tuberculosis* variant *bovis* is the second row above — while others
are merely redundant, such as *Mycobacteroides abscessus* against *M. abscessus* subsp. *abscessus*.
The two are the same shape and cannot be told apart without a taxonomy, so both are kept: a
redundant line in a note costs little, and suppressing the shape would hide the bovis case.

It is an all-pairs comparison, so cost grows with the square of the reference count: 3,534
references of `mycobacteriaceae` take under two minutes and 850 MB and flag 945 references, but a
database with hundreds of thousands of references is a different proposition. This is why it is a
separate step rather than part of every build.

Read the table top-down and it doubles as the sharpest mislabel detector available: a reference
named for one species sitting a hair from another species' references is either a genuinely
inseparable pair or a bad label. Both mislabelled `mycobacteriaceae` references were in the first 20
of 2,102 rows, `GCA_025822535.1` at rank 4.

## Correcting a mislabelled reference

Some assemblies are deposited under the wrong organism name, and a reference named for a species it
does not belong to misidentifies every sample of the species it really is. Because Mash's
winner-take-all strategy credits each shared hash to a single reference, such a reference only has to
win by a few hashes to take them all, which also pushes the correct references far down the hit list:
a runner-up at a much lower identity is not evidence that it matched poorly.

NCBI's own ANI check is where to start looking, but it does not settle the question on its own. The
`average_nucleotide_identity` block of a NCBI Datasets assembly report carries `match_status`:

- `species_match` — the best-matching type assembly is the declared species.
- `below_threshold_match` — the best match *is* the declared species but below the species threshold,
  usually a divergent or undescribed lineage rather than a mislabel.
- `below_threshold_mismatch` — common and usually benign where species sit above 98% ANI or the
  declared species has no type assembly, as across the *M. tuberculosis* complex.
- `mismatch` — the best-matching type assembly is a *different* species. Worth investigating, though
  most hits are genus reclassifications (*Mycobacterium* → *Mycolicibacterium*) or reciprocal synonym
  pairs: of 3,534 references in `mycobacteriaceae`, 38 were flagged and 2 were real.

Confirm before relabelling, with evidence that suits the group. ANI against the type strain of both
candidate species decides it for most genera; within the *M. tuberculosis* complex, where ANI is
uninformative, the regions of difference do (`GCA_013697125.1` was settled by the bovis-exclusive RD4
deletion). Record the correction, the deposited name and the evidence in `curation/label_corrections.tsv`,
then apply it.

On an existing database, rewrite the sidecar and `bins.tsv` in place — no re-sketching, and
`--dry-run` first:

```bash
scripts/apply_label_corrections.py curation/label_corrections.tsv \
    ~/.local/share/mashID/db/<db>.metadata.tsv <build_dir>/bins.tsv
```

`Identification` and `TaxID` are then correct. `Description` comes from the reference's own fasta
header by way of the sketch, so it keeps the deposited name until the database is rebuilt.

On a rebuild, pass the same file to both steps. The binner needs it too: binning is what decides
which genomes a reference is dereplicated against, and a mislabelled genome put in the wrong bin is
compared only with the species it is named for, so it survives as a near-duplicate of the cluster it
actually belongs to — which is how both `mycobacteriaceae` references got in.

```bash
scripts/bin_by_species.py assembly_data_report.jsonl ncbi_dataset/data binned/ \
    --corrections curation/label_corrections.tsv
make_mashID_db -i binned/ -o /path/to/db -p my_database \
    --assembly-report assembly_data_report.jsonl \
    --metadata curation/label_corrections.tsv
```

## Choosing a sketch size

The pre-built Mycobacteriaceae database uses s=1000, which resolves species but yields coarse
identities (steps of 0.001). For a custom database aimed at close relatives, s=10000 is worth the
tenfold larger file: identities become finer and the shared-hash counts more informative.
