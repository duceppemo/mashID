# mashID
Identify organisms using Mash.

## Description
MashID was developed to quickly identify Mycobacteria species from WGS data. Users can provide their own mash sketches database to identify other organisms.

MashID works with Illumina (paired-end or single-end), Ion Torrent, Nanopore and PacBio data. Actually any fastq (sequencing data) or fasta (genome assemblies), gzipped or not, can be used as input.

## Installation
```commandline
# Create and activate virtual envrionment
conda create -n mashID -y -c bioconda python=3.10 mash=2.3 pandas=1.4.4 
conda activate mashID

# Clone repo and test mashID
git clone https://github.com/duceppemo/mashID
cd mashID
python mashID.py -h
```

## Usage
```
usage: python mashID.py [-h] -i /input/folder/ -o /output/folder/
                 [-m /path/to/file.msh] [-v]

Species identification from NGS data using Mash.

optional arguments:
  -h, --help            show this help message and exit
  -i /input/folder/, --input /input/folder/
                        Input directory with fastq/fasta files or a single
                        fastq/fasta file, gzipped or not.
  -o /output/folder/, --output /output/folder/
                        Output directory
  -m /path/to/file.msh, --mash-db /path/to/file.msh
                        Mash sketch database. Will run a Mycobacteria mash
                        database by default.
  -v, --version         show program's version number and exit
```

## Outputs
- sample1_mashID.tsv: individual sample `mash screen` output table.
- summary_mashID.tsv: if more than one sample, this file will hold the top identification result for each sample.

## Building custom databse

