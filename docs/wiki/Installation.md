# Installation

Requirements: Python ≥ 3.10 and Mash 2.3. mashID has no Python dependencies.

## conda + pip (recommended)

```bash
conda create -n mashID -c conda-forge -c bioconda python=3.12 mash=2.3 python-isal
conda activate mashID
pip install https://github.com/duceppemo/mashID/archive/refs/tags/v0.2.7.tar.gz
mashID -h
```

`python-isal` is optional; it makes gzip decompression two to three times faster when counting reads.

Three commands are installed: `mashID`, `make_mashID_db` and `mashID_download_db`.

## From a clone

```bash
git clone https://github.com/duceppemo/mashID
cd mashID
conda env create -f environment.yml
conda activate mashID
pip install .
```

Running from the source tree without installing also works: `python mashID.py -h`.

## Bioconda

A recipe is [submitted to bioconda](https://github.com/bioconda/bioconda-recipes/pull/69474). Once
merged, `conda install -c conda-forge -c bioconda mashid` installs everything, Mash included.

## Verify the installation

The repository ships a small synthetic [example](https://github.com/duceppemo/mashID/tree/master/example)
with its own database and expected results, so nothing needs to be downloaded:

```bash
git clone https://github.com/duceppemo/mashID && cd mashID   # if not already cloned
bash example/run_example.sh
```

It screens four samples (paired-end reads, a mixed culture, an assembly, and an unrelated sequence)
and reports `OK` when the summary matches `example/expected/summary_mashID.tsv`. The [Example](Example)
page walks through the results.

## Databases

No database is bundled. Download the default one (284 MB) after installing:

```bash
mashID_download_db mycobacteriaceae
```

See [Databases](Databases) for the other pre-built databases and how to build your own.

## Troubleshooting

**`mash: error while loading shared libraries: libgsl.so.25`**
The solver picked an older bioconda build of Mash 2.3 linked against GSL 2.6/2.7. Request the current
build, which links against GSL 2.8:

```bash
conda install -c conda-forge -c bioconda "mash=2.3=*_11"
```

**`The 'mash' executable was not found on PATH`**
Activate the conda environment that contains Mash, or install it: `conda install -c bioconda mash`.

**`No database given and the default 'mycobacteriaceae' database is not installed`**
Run `mashID_download_db mycobacteriaceae`, or pass `-d /path/to/your.msh`.
