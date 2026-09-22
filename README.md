<p align="center">
  <img src="docs/images/mashID_logo.png" alt="mashID logo" width="200">
</p>

<h1 align="center">mashID</h1>

<p align="center">
  Species identification from genome assemblies or raw sequencing reads, powered by
  <a href="https://github.com/marbl/Mash">Mash</a>.
</p>

<p align="center">
  <a href="https://github.com/duceppemo/mashID/actions/workflows/ci.yml"><img src="https://github.com/duceppemo/mashID/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://github.com/duceppemo/mashID/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/mashID?label=release" alt="Latest release"></a>
  <!-- Enable once the bioconda recipe is merged (https://github.com/bioconda/bioconda-recipes/pull/69474):
  <a href="https://anaconda.org/bioconda/mashid"><img src="https://img.shields.io/conda/vn/bioconda/mashid?label=bioconda" alt="Bioconda"></a>
  -->
  <img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10+">
  <a href="LICENSE"><img src="https://img.shields.io/github/license/duceppemo/mashID" alt="License: MIT"></a>
  <!-- Enable once a Zenodo DOI exists for the release:
  <a href="https://doi.org/10.5281/zenodo.XXXXXXX"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg" alt="DOI"></a>
  -->
</p>

mashID screens each sample against a Mash sketch database, reports the best-matching organisms with
identity, coverage and p-value, and flags mixtures, ambiguous calls and low coverage. It takes fasta
or fastq files, gzipped or not, from Illumina, Ion Torrent, Nanopore or PacBio. Paired-end files are
screened together.

## Quick start

```bash
conda create -n mashID -c conda-forge -c bioconda python=3.12 mash=2.3 python-isal
conda activate mashID
pip install https://github.com/duceppemo/mashID/archive/refs/tags/v0.2.1.tar.gz

mashID_download_db mycobacteriaceae          # 27 MB, MD5-verified
mashID -i /path/to/fastq_or_fasta -o results
```

```
Sample    Sequences  Bases      Identity  Shared_Hashes  Median_Multiplicity  P_Value  Accession        TaxID  Identification                            Note
MBWGS440  2741728    625918949  0.999713  994/1000       122                  0        NAZK01000001.1   NA     Mycobacterium tuberculosis variant bovis
```

Results are written to `summary_mashID.tsv` (one line per sample) and `<sample>_mashID.tsv` (top hits).

To check an installation without downloading anything, run the bundled [example](example/):

```bash
bash example/run_example.sh
```

## Documentation

Everything else lives in the [wiki](https://github.com/duceppemo/mashID/wiki), whose sources are
maintained in [`docs/wiki`](docs/wiki):

- [Installation](https://github.com/duceppemo/mashID/wiki/Installation) — conda, pip, troubleshooting
- [Usage](https://github.com/duceppemo/mashID/wiki/Usage) — options, input files, sample naming, quick screening of large runs
- [Outputs](https://github.com/duceppemo/mashID/wiki/Outputs) — columns, the `Note` flags, limits of Mash identification
- [Databases](https://github.com/duceppemo/mashID/wiki/Databases) — pre-built databases, building your own, metadata sidecar
- [Development](https://github.com/duceppemo/mashID/wiki/Development) — tests, CI, releases, bioconda

## Citation

If mashID is useful in your work, please cite it (see [`CITATION.cff`](CITATION.cff)):

> Duceppe, M.-O. (2026). mashID: species identification from genome assemblies and raw reads using Mash (v0.2.1). https://github.com/duceppemo/mashID

and Mash: Ondov et al. (2016) *Mash: fast genome and metagenome distance estimation using MinHash*, Genome Biology 17:132. https://doi.org/10.1186/s13059-016-0997-x

## License

[MIT](LICENSE)
