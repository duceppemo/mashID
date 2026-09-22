#!/usr/bin/env bash
# Build a mashID database for an NCBI taxon (default: Mycobacteriaceae, taxid 1762).
#
# Usage:  build_mycobacteriaceae_db.sh <work_dir> [taxid] [db_name]
# Env:    NCBI_API_KEY     optional, raises NCBI rate limits (never hard-code it)
#         THREADS          default: all CPUs
#         DEREP_DISTANCE   Mash distance below which assemblies of a species are collapsed (default 0.001)
#         SKETCH_SIZE      default 10000
#         DEREPLICATOR     path to dereplicator.py (https://github.com/rrwick/Assembly-dereplicator)
# Needs on PATH: datasets (ncbi-datasets-cli), unzip, mash, make_mashID_db, python3.
#
# Steps (each is skipped when its output already exists, so the script can be re-run):
#   1. one dehydrated `datasets` archive for the taxon (GenBank, atypical assemblies excluded), then
#      parallel rehydration into gzipped fasta files, with assembly_data_report.jsonl for names/TaxIDs;
#   2. bin the genomes by species (scripts/bin_by_species.py);
#   3. dereplicate each species bin with Assembly-dereplicator;
#   4. sketch the survivors with make_mashID_db, excluding references < 100 kb, and write the
#      metadata sidecar from the assembly report.
set -euo pipefail

work="${1:?Usage: $0 <work_dir> [taxid] [db_name]}"
taxid="${2:-1762}"
name="${3:-mycobacteriaceae_$(date +%F)}"
threads="${THREADS:-$(nproc)}"
derep_distance="${DEREP_DISTANCE:-0.001}"
sketch_size="${SKETCH_SIZE:-10000}"
derep="${DEREPLICATOR:-$HOME/prog/Assembly-dereplicator/dereplicator.py}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
api=()
[[ -n "${NCBI_API_KEY:-}" ]] && api=(--api-key "$NCBI_API_KEY")
log() { echo "[$(date '+%F %T')] $*"; }

for tool in datasets unzip mash make_mashID_db python3; do
    command -v "$tool" > /dev/null || { echo "$tool not found on PATH" >&2; exit 1; }
done
[[ -f "$derep" ]] || { echo "dereplicator.py not found at $derep (set DEREPLICATOR)" >&2; exit 1; }

mkdir -p "$work"
cd "$work"
report=ncbi/ncbi_dataset/data/assembly_data_report.jsonl

if [[ ! -f "$report" ]]; then
    log "1. Downloading the dehydrated archive for taxon $taxid"
    datasets download genome taxon "$taxid" --assembly-source GenBank --exclude-atypical \
        --include genome --dehydrated --filename ncbi.zip "${api[@]}"
    unzip -q -o ncbi.zip -d ncbi
fi
log "1. Rehydrating genomes (gzipped, 10 workers)"
datasets rehydrate --directory ncbi --gzip --max-workers 10 "${api[@]}"
log "   $(find ncbi/ncbi_dataset/data -name '*.fna.gz' | wc -l) genome files"

if [[ ! -f binned/bins.tsv ]]; then
    log "2. Binning by species"
    python3 "$here/bin_by_species.py" "$report" ncbi/ncbi_dataset/data binned
fi

log "3. Dereplicating each species at Mash distance $derep_distance"
mkdir -p derep
derep_bin() {
    local bin="$1" out="derep/$(basename "$1")"
    [[ -d "$out" ]] && return 0
    python3 "$DEREP" --distance "$DEREP_DISTANCE" --sketch_size "$SKETCH_SIZE" --threads 4 "$bin" "$out.tmp" \
        > "derep/$(basename "$1").log" 2>&1 && mv "$out.tmp" "$out"
}
export -f derep_bin
export DEREP="$derep" DEREP_DISTANCE="$derep_distance" SKETCH_SIZE="$sketch_size"
# largest bins first so they do not end up last on a single worker
find binned -mindepth 1 -maxdepth 1 -type d -print0 | xargs -0 -I{} sh -c 'echo "$(ls {} | wc -l) {}"' \
    | sort -rn | cut -d' ' -f2- | tr '\n' '\0' \
    | xargs -0 -P "$(( threads / 4 > 0 ? threads / 4 : 1 ))" -I{} bash -c 'derep_bin "$@"' _ {}
log "   $(find derep -mindepth 2 -name '*.fna.gz' | wc -l) genomes kept from $(find binned -mindepth 2 -name '*.fna.gz' | wc -l)"

log "4. Sketching (k=21, s=$sketch_size) and writing the metadata sidecar"
make_mashID_db -i derep -o . -p "$name" -s "$sketch_size" -k 21 -t "$threads" --assembly-report "$report"
log "Done: $work/$name.msh  (+ $name.metadata.tsv)"
