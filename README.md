# mashID
Identify organisms using Mash.

## Description
MashID was developed to quickly identify Mycobacteria species from WGS data. Users can provide their own mash sketches database to identify other organisms.

MashID works with Illumina (paired-end or single-end), Ion Torrent, Nanopore and PacBio data. Actually any fastq (sequencing data) or fasta (genome assemblies), gzipped or not, can be used as input.

## Installation
```commandline
git clone 
```

## Usage
```commandline
usage: python mashID.py [-h] -i /input/folder/ -o /output/folder/
                 [-m /path/to/file.msh] [-v]

Species identification from NGS data using Mash.

optional arguments:
  -h, --help            Show this help message and exit
  -i /input/folder/, --input /input/folder/
                        Input directory with fastq files or a single fastq
                        file, gzipped or not.
  -o /output/folder/, --output /output/folder/
                        Output directory
  -m /path/to/file.msh, --mash-db /path/to/file.msh
                        Mash sketch database. Will run a Mycobacteria mash
                        database by default.
  -v, --version         Show program's version number and exit

```