# Changelog

## Unreleased

### Added
- `progenomes4` database: proGenomes 4 representatives (32,887 genomes, s=2000, NCBI names and
  TaxIDs, [Figshare](https://doi.org/10.6084/m9.figshare.33968806)) with its metadata sidecar;
  `scripts/build_progenomes_db.sh` and `scripts/split_multifasta_by_genome.py` build it.
- Build helpers: bin by subspecies or serovar, cap bin size preferring complete assemblies;
  assembly reports in `datasets summary` (snake_case) form are accepted.

## 0.2.5 (2026-09-22)

### Fixed
- Symbolic links to missing files in the input are reported up front with their targets, instead of a
  raw Mash error after read counting.

### Changed
- The `listeria` database is the 2026-09-22 rebuild from RefSeq (1441 references, s=10000, NCBI TaxIDs,
  [Figshare](https://doi.org/10.6084/m9.figshare.33968620)) with its metadata sidecar. The 2025 one
  remains available as `listeria-2025`.

### Added
- `ASSEMBLY_SOURCE` option (GenBank or RefSeq) in the database build script.

## 0.2.4 (2026-09-22)

### Fixed
- The MultiQC table's header no longer embeds the mashID version; the example's expected table went
  stale on every version bump, which made the 0.2.3 tarball's example test fail.

## 0.2.3 (2026-09-22)

### Changed
- The example ships and compares its MultiQC table; the Example wiki page lists every file produced.
- Workflow diagram redrawn vertically with larger text, in the README and on the wiki Home page.

## 0.2.2 (2026-09-22)

### Changed
- The default `mycobacteriaceae` database is the 2026-09-22 rebuild (3534 references, s=10000, NCBI
  TaxIDs, no partial records, [Figshare](https://doi.org/10.6084/m9.figshare.33965176)); it ships its
  metadata sidecar, which `mashID_download_db` now fetches alongside the sketch. The 2025 database
  remains available as `mycobacteriaceae-2025`.

### Added
- `Est_Depth` column: bases divided by the top reference's length (reads only).
- `summary_mashID.json`, `summary_mashID_mqc.tsv` (MultiQC custom content) and `mashID_run.json`
  (provenance: versions, command line, parameters, database MD5, inputs, timestamps, exit code).
- `--sample-sheet`: explicit sample-to-file mapping; extra columns are carried into the summary.
- `--fail-on {no-hit,note}` exit code 2 for pipeline gating.
- `make_mashID_db --check`: database quality report (short references, unparsed names, missing TaxIDs,
  genus synonyms).
- `scripts/build_mycobacteriaceae_db.sh` rewritten around NCBI `datasets` dehydrated downloads and
  `scripts/bin_by_species.py`; TaxIDs from the assembly report.
- Header-parser regression corpus of 3309 real NCBI headers; macOS job in CI.
- Wiki pages: Interpreting results, Pipelines, FAQ, Comparison; database provenance and known issues.

## 0.2.1 (2026-09-22)

### Added
- `example/`: a self-contained synthetic dataset (three fictional genomes, paired-end reads, a mixed
  culture, an assembly, an unrelated sequence) with expected results and `run_example.sh` to verify
  an installation offline. Checked in CI.
- Minimal README with logo and badges; documentation moved to the wiki, whose sources live in
  `docs/wiki/` and are published by a workflow. Release workflow, citation file, contributing guide,
  code of conduct, security policy, issue and pull request templates, dependabot, pre-commit.

## 0.2.0 (2026-09-21)

### Breaking changes
- The 2022 Mycobacteria database is no longer bundled. Databases are downloaded on demand with
  `mashID_download_db` (MD5-verified) into `$MASHID_DB_DIR` or `~/.local/share/mashID/db`; the
  default database is the 2025-02-20 Mycobacteriaceae one, which supersedes the bundled file.
  `-d` accepts a path or a database name.
- `summary_mashID.tsv` replaces `topID.tsv` (the name the README always documented).
- Column names no longer contain `%` or `-`: `Sequences`, `Bases`, `Identity`, `Shared_Hashes`,
  `Median_Multiplicity`, `P_Value`, `Accession`, `TaxID`, `Identification`, `Note`. The old
  `%-Identity` column was a fraction, not a percentage. Per-sample tables gain `TaxID`,
  `Identification` and `Description` columns.
- Sample names are derived by stripping the extension and Illumina read designations
  (`_S1_L001_R1_001`, `_R1`, `_1`), instead of everything after the first underscore.
  `isolate_A_R1.fq.gz` is now sample `isolate_A` and no longer collides with `isolate_B_R1.fq.gz`.
- `-m/--memory` is accepted but ignored (BBMap is no longer used).
- Installable package: `pip install .` provides the `mashID` and `make_mashID_db` commands.
  `python mashID.py` and `python make_mashID_db.py` keep working.

### Dependencies
- Dropped pandas, psutil, BBMap and `pkg_resources`. Only Python ≥ 3.10 and Mash 2.3 are required.
  Reads/bases are counted in pure Python (optionally accelerated by `python-isal`).

### Fixes
- Paired-end files are passed to `mash screen` directly instead of being concatenated into a copy
  under `<output>/tmp`, which was then deleted with `shutil.rmtree` (destroying any pre-existing
  `tmp` folder in the output directory).
- Argument validation was ineffective (`if 0 < x > 1` never triggers for negatives; thread/memory
  checks returned values that were discarded). Ranges are now validated by argparse.
- `mash` errors were silenced (`stderr` discarded, exit code ignored) and surfaced as
  "No significant hit". Failures now abort with Mash's error message; a missing `mash` executable or
  database is reported up front.
- `-p` larger than `-t` produced `-p 0` for Mash; threads per sample are now at least 1.
- Organism-name parsing no longer crashes on short or unusual headers and handles `Candidatus`,
  bracketed genera, `sp.`, serovars, multiple infraspecific ranks, and headers without an accession.
  The former check for `sub` matched unrelated words such as "substrain".
- Accessions are derived correctly from `.fna.gz` references.
- Read counts from BBMap `stats.sh` were slightly low on some fastq files; counts now match `awk`.
- `make_mashID_db`: k-mer size range check was inverted; silent failures now raise; a `.msh` suffix
  given in `-p` no longer yields `name.msh.msh`; the temporary file list is no longer left in the
  output directory; input may be a directory or a list file.

### Security
- The NCBI API key that was hard-coded in `mashID_Mycobactriaceae_DB.sh` has been removed; the
  rewritten `scripts/build_mycobacteriaceae_db.sh` reads `NCBI_API_KEY` from the environment.
  The old key remains in git history and must be regenerated at NCBI.

### Added
- `Note` column flagging possible mixtures, ambiguous calls and low coverage
  (`--ambiguity-margin`).
- Hits to references shorter than `--min-ref-length` (100 kb) are ignored, and `make_mashID_db`
  excludes such references (`--min-length`). Partial records otherwise become false top hits.
- `--max-reads N`: stream only the first N reads of each fastq sample to `mash screen` for quick
  checks on large runs.
- Metadata sidecar `<db>.metadata.tsv` (accession, organism, TaxID) written by `make_mashID_db` from
  `--metadata` tables or NCBI `--assembly-report` files, or parsed headers, and including reference
  lengths; `--annotate` adds one to an existing database and `mashID` generates it on first use.
  `mashID` uses it (or `--db-metadata`) instead of parsing headers.
- `mashID_download_db` and a registry of the Figshare databases with sizes and checksums.
- `--no-winner-take-all`, `--skip-stats`, `--debug` options.
- Test suite (`pytest`), `ruff` lint, GitHub Actions CI, bioconda recipe, `pyproject.toml`,
  `environment.yml`, `.gitignore`.

## 0.1.1
- Previous release.
