#!/usr/bin/env bash
# Build a mashID database for any NCBI taxon from its GenBank or RefSeq assemblies
# (default taxon: Mycobacteriaceae, taxid 1762).
#
# Usage:  build_taxon_db.sh <work_dir> [taxid[,taxid...]] [db_name]
# Examples:
#   build_taxon_db.sh /db/myco 1762 mycobacteriaceae_2026-09-22
#   ASSEMBLY_SOURCE=RefSeq build_taxon_db.sh /db/listeria 1637 listeria_2026-09-22
#   ASSEMBLY_SOURCE=RefSeq BIN_RANK=subspecies MAX_BIN=3000 build_taxon_db.sh /db/salmonella 590 salmonella_2026-09-22
#   ASSEMBLY_SOURCE=RefSeq MAX_BIN=5000 EXCLUDE_SP=1 build_taxon_db.sh /db/escherichia 561,620 escherichia_shigella_2026-09-22
# Env:    NCBI_API_KEY     optional, raises NCBI rate limits (never hard-code it)
#         ASSEMBLY_SOURCE  GenBank (default) or RefSeq. Use RefSeq for heavily sequenced taxa
#                          (e.g. Listeria: 79,000 GenBank vs 7,500 RefSeq assemblies), since
#                          dereplication cost grows with the square of the largest species bin
#         BIN_RANK         species (default), subspecies or serovar: how genomes are grouped before
#                          dereplication (scripts/bin_by_species.py --rank)
#         MAX_BIN          cap per bin, best assembly levels kept first (default 0 = no cap); use for
#                          taxa with tens of thousands of assemblies per species
#         EXCLUDE_SP       1 to leave out bins of unnamed species ("Genus sp.", "uncultured ...") so the
#                          database never answers "Genus sp." (default 0)
#         THREADS          default: all CPUs
#         DEREP_DISTANCE   Mash distance below which assemblies of a species are collapsed (default 0.001)
#         SKETCH_SIZE      default 10000
#         DEREPLICATOR     path to dereplicator.py (https://github.com/rrwick/Assembly-dereplicator)
# Needs on PATH: datasets (ncbi-datasets-cli), unzip, mash, make_mashID_db, python3.
#
# Steps (each is skipped when its output already exists, so the script can be re-run):
#   1. one dehydrated `datasets` archive for the taxon (GenBank or RefSeq, atypical assemblies excluded), then
#      parallel rehydration into gzipped fasta files, with assembly_data_report.jsonl for names/TaxIDs;
#   2. bin the genomes by species (scripts/bin_by_species.py);
#   3. dereplicate each species bin with Assembly-dereplicator;
#   4. sketch the survivors with make_mashID_db, excluding references < 100 kb, and write the
#      metadata sidecar from the assembly report.
set -euo pipefail

work="${1:?Usage: $0 <work_dir> [taxid] [db_name]}"
taxa="${2:-1762}"
name="${3:-mycobacteriaceae_$(date +%F)}"
threads="${THREADS:-$(nproc)}"
source="${ASSEMBLY_SOURCE:-GenBank}"
bin_rank="${BIN_RANK:-species}"
max_bin="${MAX_BIN:-0}"
exclude_sp="${EXCLUDE_SP:-0}"
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

# One dehydrated archive per taxon (ncbi_<taxid>/), merged into ncbi/ through symbolic links.
for taxid in ${taxa//,/ }; do
    if [[ ! -f "ncbi_$taxid/ncbi_dataset/data/assembly_data_report.jsonl" ]]; then
        log "1. Downloading the dehydrated archive for taxon $taxid ($source assemblies)"
        datasets download genome taxon "$taxid" --assembly-source "$source" --exclude-atypical \
            --include genome --dehydrated --filename "ncbi_$taxid.zip" "${api[@]}"
        unzip -q -o "ncbi_$taxid.zip" -d "ncbi_$taxid"
    fi
    log "1. Rehydrating genomes of taxon $taxid (gzipped, 10 workers)"
    datasets rehydrate --directory "ncbi_$taxid" --gzip --max-workers 10 "${api[@]}"
done
mkdir -p ncbi/ncbi_dataset/data
: > "$report"
for taxid in ${taxa//,/ }; do
    cat "ncbi_$taxid/ncbi_dataset/data/assembly_data_report.jsonl" >> "$report"
    for d in "ncbi_$taxid"/ncbi_dataset/data/GC*; do
        [[ -d "$d" ]] && ln -sfn "$(readlink -f "$d")" "ncbi/ncbi_dataset/data/$(basename "$d")"
    done
done
log "   $(find -L ncbi/ncbi_dataset/data -name '*.fna.gz' | wc -l) genome files"

if [[ ! -f binned/bins.tsv ]]; then
    log "2. Binning by species"
    python3 "$here/bin_by_species.py" --rank "$bin_rank" --max-bin "$max_bin" "$report" ncbi/ncbi_dataset/data binned
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
if [[ "$exclude_sp" == "1" ]]; then
    find "$PWD/derep" -mindepth 2 -name '*.fna.gz' | grep -vE '/derep/([^/]+_sp|uncultured_[^/]+|unknown)/' > genomes_named.txt
    log "   excluding unnamed species: $(wc -l < genomes_named.txt) of $(find derep -mindepth 2 -name '*.fna.gz' | wc -l) genomes kept"
    input=genomes_named.txt
else
    input=derep
fi
make_mashID_db -i "$input" -o . -p "$name" -s "$sketch_size" -k 21 -t "$threads" --assembly-report "$report"
log "Done: $work/$name.msh  (+ $name.metadata.tsv)"
