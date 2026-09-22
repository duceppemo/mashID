# FAQ and troubleshooting

## Installation

**`mash: error while loading shared libraries: libgsl.so.25`**
An older bioconda build of Mash 2.3 was installed with a newer GSL. Install the current build:
`conda install -c conda-forge -c bioconda "mash=2.3=*_11"`.

**`The 'mash' executable was not found on PATH`**
Activate the environment that has Mash, or `conda install -c bioconda mash`.

**Does it run on macOS or Windows?**
Linux and macOS (CI tests both; bioconda has Mash for Intel and Apple silicon). Not Windows, except
through WSL.

## Databases

**`No database given and the default 'mycobacteriaceae' database is not installed`**
`mashID_download_db mycobacteriaceae`, or pass `-d /path/to/db.msh`.

**Where are downloaded databases stored?**
`$MASHID_DB_DIR` if set, else `~/.local/share/mashID/db`. `mashID_download_db` with no arguments lists
them.

**Can I use a Mash sketch made elsewhere?**
Yes, any `.msh` works. Organism names are parsed from the sketch comments (the first fasta header of
each genome) into a metadata sidecar on first use. Run `make_mashID_db --check` to see how well that
went, and `--annotate` with a `--metadata` table to fix names.

**How big should the sketch size be?**
1000 resolves species; 10000 gives finer identities and better shared-hash counts, and is worth it for
close relatives. The database file is ten times larger.

## Results

**Every sample of a run is a "possible mixture" with the same organism.**
Either the run is contaminated, or the database lists that species under two genus names (check with
`make_mashID_db --check`), or a partial reference is involved (look at the per-sample table's
`Shared_Hashes` denominator and `Description`). See [Interpreting results](Interpreting-results).

**The identification alternates between *M. tuberculosis* and *M. bovis*.**
Expected. Mash cannot separate members of the *M. tuberculosis* complex; report the complex and use a
SNP-based method for the variant.

**`No significant hit in database` for a sample I know is bacterial.**
The organism is not in the database (try `refseq_bacteria` or `progenomes3`), the reads are mostly
something else (host, adapters), or the file is empty or truncated. `--identity 0.8` shows weaker hits.

**`Sequences` and `Bases` are `NA`.**
`--skip-stats` was given, or counting failed (see the log).

**`Est_Depth` is `NA` for reads.**
The database has no reference lengths, which happens with a `--db-metadata` table that lacks a
`Length` column. Let mashID use the sidecar, or add the column.

**Why is `Median_Multiplicity` lower than `Est_Depth`?**
Multiplicity counts k-mer observations and is reduced by sequencing errors and read ends (a 150 bp read
holds 130 21-mers), so it is typically 20 to 40% below the base-level depth. A much larger gap means
much of the sample is not the reported organism.

**Two files of the same sample were treated as two samples.**
Names must match after stripping the extension and Illumina read designations (`_S1_L001_R1_001`,
`_R1`, `_1`). For other conventions use a [sample sheet](Usage#sample-sheets).

**Two different samples were merged into one.**
Their names collapse to the same string after suffix stripping (for example `isolate_1.fastq` and
`isolate_2.fastq` both become `isolate`). Rename them or use a sample sheet.

## Performance

**It is slow on large Nanopore runs.**
Use `--max-reads 500000`. Screening 500k reads takes seconds and gives the same identification.

**How much memory does it need?**
Mash loads the database sketches (a few hundred MB for the largest pre-built databases) plus a
hash table of the sample's k-mers. A bacterial isolate needs well under 2 GB per concurrent sample.

**Reads and bases counting takes long.**
It reads every gzipped file once in Python. Install `python-isal` for a two- to three-fold speedup,
use `--max-reads` (counts only what is screened), or `--skip-stats`.
