<p align="center"><img src="https://raw.githubusercontent.com/duceppemo/mashID/master/docs/images/mashID_logo.png" alt="mashID logo" width="180"></p>

# mashID

Species identification from genome assemblies (fasta) or raw sequencing reads (fastq) with
[Mash](https://github.com/marbl/Mash).

mashID runs `mash screen` for every sample against a sketch database, keeps the best hits, names the
organism from the database metadata, and flags results that deserve a second look: possible
mixtures, ambiguous calls, low coverage. It works with Illumina (paired-end or single-end),
Ion Torrent, Nanopore and PacBio reads, and with assemblies; inputs may be gzipped.

## Pages

- [Installation](Installation)
- [Usage](Usage)
- [Example](Example)
- [Interpreting results](Interpreting-results)
- [Pipelines](Pipelines)
- [FAQ](FAQ)
- [Comparison](Comparison)
- [Outputs](Outputs)
- [Databases](Databases)
- [Development](Development)

## In one minute

```bash
conda create -n mashID -c conda-forge -c bioconda python=3.12 mash=2.3 python-isal
conda activate mashID
pip install https://github.com/duceppemo/mashID/archive/refs/tags/v0.2.3.tar.gz
mashID_download_db mycobacteriaceae
mashID -i /path/to/reads -o results
```

## How it works

<p align="center"><img src="https://raw.githubusercontent.com/duceppemo/mashID/master/docs/images/workflow_v3.svg" alt="mashID workflow" width="760"></p>

1. Input files are grouped into samples by name; R1/R2 pairs and multi-lane files belong to one sample.
2. Reads and bases are counted (pure Python, optionally accelerated by `python-isal`).
3. `mash screen -w` estimates, for every reference in the database, how much of its sketch is
   contained in the sample. Winner-take-all credits shared hashes to the best reference.
4. Hits to references shorter than 100 kb are ignored (partial records), the rest are ranked, and the
   organism name, TaxID and reference length are read from the database's metadata sidecar.
5. Notes are added from the hit list, and the summary (TSV, JSON, MultiQC), per-sample tables and a
   run-provenance file are written.

> These pages are generated from [`docs/wiki`](https://github.com/duceppemo/mashID/tree/master/docs/wiki)
> in the repository. Edit them there; a workflow publishes changes to this wiki.
