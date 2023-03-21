# mashID
Identify organisms using Mash.

## Description
MashID was developed to quickly identify bacterial species from WGS data. Users can provide their own mash sketches database to identify other organisms.

MashID works with Illumina (paired-end or single-end), Ion Torrent, Nanopore and PacBio data. Actually any fastq (sequencing data) or fasta (genome assemblies), gzipped or not, can be used as input.

MashID comes with a pre-compiled database to identify mycobacteria species.
## Installation
```commandline
# Create and activate virtual envrionment
conda create -n mashID -y -c bioconda python=3.10 mash=2.3 pandas=1.4.4 psutil=5.9.2 bbmap=38.99
conda activate mashID

# Clone repo and test mashID
git clone https://github.com/duceppemo/mashID
cd mashID
python mashID.py -h
```

## Usage
```
usage: python mashID.py [-h] -i /input/folder/ -o /output/folder/ [-d /path/to/mash_databse.msh] [--identity 0.9] [--p-value 0.05] [-n 10] [-t 16] [-p 2] [-m 57] [-v]

Species identification from NGS data using Mash.

options:
  -h, --help            show this help message and exit
  -i /input/folder/, --input /input/folder/
                        Input directory with fastq/fasta files or a single fastq/fasta file, gzipped or not.
  -o /output/folder/, --output /output/folder/
                        Output directory
  -d /path/to/mash_databse.msh, --database /path/to/mash_databse.msh
                        Mash sketch database. Will run a Mycobacteria mash database by default.
  --identity 0.9        Minimum identity to report. [0-1]. Default is 0.9
  --p-value 0.05        Maximum p-value to report
  -n 10, --n-hits 10    Number of top-hits to report (sorted by % identity). Default is 10.
  -t 16, --threads 16   Number of threads. Default is maximum available(16). Optional.
  -p 2, --parallel 2    Number of samples to process in parallel. Default is 2. Optional.
  -m 57, --memory 57    Memory in GB. Default is 85% of total memory (57)
  -v, --version         show program's version number and exit
```

## Outputs
- sample1_mashID.tsv: individual sample `mash screen` output table.
- summary_mashID.tsv: if more than one sample, this file will hold the top identification result for each sample.

## Building custom database
You can use the `make_mashID_db.py` script to build a suitable database for mashID.
```commandline
python make_mashID_db.py \
    -i /media/data/sample_folder \
    -o /media/database/mycoID' \
```
Run `python make_mashID_db.py -h` for detailed help.

- Pre-compiled databases for Refseq bacteria and proGenomes v3 can be found [here](https://figshare.com/account/home#/projects/162688)
