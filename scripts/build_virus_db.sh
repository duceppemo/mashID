#!/usr/bin/env bash
# Build a mashID database of all RefSeq viral genomes from the RefSeq release files on the NCBI FTP site.
#
# Usage:  build_virus_db.sh <work_dir> [db_name]
# Env:    THREADS      default: all CPUs
#         SKETCH_SIZE  default 2000 (viral genomes are small; a sketch cannot exceed the genome's k-mers)
#         KMER_SIZE    default 21
# Needs on PATH: curl, zcat, mash, make_mashID_db, python3.
#
# Steps (each skipped when its output exists, so the script can be re-run):
#   1. download viral.*.genomic.fna.gz and viral.*.genomic.gbff.gz of the current RefSeq release;
#   2. split the fasta into one file per record (scripts/split_multifasta_by_genome.py), so every
#      RefSeq genome or genome segment becomes one reference;
#   3. extract organism names and TaxIDs from the GenBank flat files (scripts/refseq_gbff_metadata.py);
#   4. sketch with make_mashID_db, keeping references of every length (--min-length 0).
# RefSeq viral holds one reference per species (or segment), so there is no dereplication.
set -euo pipefail

work="${1:?Usage: $0 <work_dir> [db_name]}"
release="$(curl -sL --retry 3 https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER | tr -d '[:space:]')"
name="${2:-refseq_viral_r${release}_$(date +%F)}"
threads="${THREADS:-$(nproc)}"
sketch_size="${SKETCH_SIZE:-2000}"
kmer_size="${KMER_SIZE:-21}"
base="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
log() { echo "[$(date '+%F %T')] $*"; }

for tool in curl zcat mash make_mashID_db python3; do
    command -v "$tool" > /dev/null || { echo "$tool not found on PATH" >&2; exit 1; }
done

mkdir -p "$work"
cd "$work"
log "RefSeq release $release"

log "1. Downloading the viral release files"
files=$(curl -sL --retry 3 "$base/" | grep -oE 'href="viral\.[0-9]+\.(1\.genomic\.fna|genomic\.gbff)\.gz"' | sed 's/href=//; s/"//g' | sort -u)
[[ -n "$files" ]] || { echo "No viral release files found at $base" >&2; exit 1; }
for f in $files; do
    [[ -f "$f.done" ]] && continue
    curl -sL -C - --retry 5 --retry-delay 30 -o "$f" "$base/$f"
    touch "$f.done"
done
echo "$release" > release.txt

if [[ ! -f genomes/accessions.txt ]]; then
    log "2. Splitting into one fasta per record"
    zcat viral.*.1.genomic.fna.gz | python3 "$here/split_multifasta_by_genome.py" genomes
fi
log "   $(wc -l < genomes/accessions.txt) records"

if [[ ! -s metadata.tsv ]]; then
    log "3. Extracting names and TaxIDs from the GenBank flat files"
    python3 "$here/refseq_gbff_metadata.py" metadata.tsv viral.*.genomic.gbff.gz
fi

log "4. Sketching (k=$kmer_size, s=$sketch_size) and writing the metadata sidecar"
make_mashID_db -i genomes -o . -p "$name" -s "$sketch_size" -k "$kmer_size" -t "$threads" \
    --min-length 0 --metadata metadata.tsv
log "Done: $work/$name.msh  (+ $name.metadata.tsv)"
