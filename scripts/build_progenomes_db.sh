#!/usr/bin/env bash
# Build a mashID database from proGenomes representative genomes (https://progenomes.embl.de).
#
# Usage:  build_progenomes_db.sh <work_dir> [db_name] [url]
# Env:    NCBI_API_KEY   optional; used when fetching TaxIDs from NCBI
#         THREADS        default: all CPUs
#         SKETCH_SIZE    default 10000
# Needs on PATH: curl, zcat, datasets (ncbi-datasets-cli), mash, make_mashID_db, python3.
#
# Steps (each skipped when its output exists, so the script can be re-run):
#   1. download the representatives multi-fasta (resumable);
#   2. split it into one fasta per genome (scripts/split_multifasta_by_genome.py);
#   3. fetch organism names and TaxIDs for the assembly accessions from NCBI Datasets;
#   4. sketch with make_mashID_db (references < 100 kb excluded) and write the metadata sidecar.
# proGenomes representatives are already one genome per species cluster, so there is no dereplication.
set -euo pipefail

work="${1:?Usage: $0 <work_dir> [db_name] [url]}"
name="${2:-progenomes4_$(date +%F)}"
url="${3:-https://progenomes.embl.de/data/repGenomes/pg4_genomes_representatives.fna.gz}"
threads="${THREADS:-$(nproc)}"
sketch_size="${SKETCH_SIZE:-10000}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
api=()
[[ -n "${NCBI_API_KEY:-}" ]] && api=(--api-key "$NCBI_API_KEY")
log() { echo "[$(date '+%F %T')] $*"; }

for tool in curl zcat datasets mash make_mashID_db python3; do
    command -v "$tool" > /dev/null || { echo "$tool not found on PATH" >&2; exit 1; }
done

mkdir -p "$work"
cd "$work"
archive="$(basename "$url")"

if [[ ! -f "$archive.done" ]]; then
    log "1. Downloading $url"
    curl -L -C - --retry 5 --retry-delay 30 -o "$archive" "$url"
    touch "$archive.done"
fi

if [[ ! -f genomes/accessions.txt ]]; then
    log "2. Splitting into one fasta per genome"
    zcat "$archive" | python3 "$here/split_multifasta_by_genome.py" genomes
fi
log "   $(wc -l < genomes/accessions.txt) genomes"

if [[ ! -s assembly_data_report.jsonl ]]; then
    log "3. Fetching organism names and TaxIDs from NCBI Datasets"
    datasets summary genome accession --inputfile genomes/accessions.txt --as-json-lines "${api[@]}" \
        > assembly_data_report.jsonl
fi
log "   $(wc -l < assembly_data_report.jsonl) assembly records"

log "4. Sketching (k=21, s=$sketch_size) and writing the metadata sidecar"
make_mashID_db -i genomes -o . -p "$name" -s "$sketch_size" -k 21 -t "$threads" \
    --assembly-report assembly_data_report.jsonl
log "Done: $work/$name.msh  (+ $name.metadata.tsv)"
