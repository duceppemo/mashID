#!/usr/bin/env bash
# Example: build a Mycobacteriaceae mashID database from all NCBI genomes of taxon 1762.
#
# Requirements (e.g. in a "ncbi" conda env): ncbi-datasets-cli, unzip, GNU parallel or xargs,
# Assembly-dereplicator (https://github.com/rrwick/Assembly-dereplicator), and mashID.
#
# Usage:
#   export NCBI_API_KEY=...            # optional but strongly recommended (higher rate limits)
#   bash build_mycobacteriaceae_db.sh /path/to/work_dir /path/to/ncbi_dataset.tsv
#
# The TSV is the accession table exported from https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=1762
# (first column = assembly accession).

set -euo pipefail

base_dir="${1:?Usage: $0 <work_dir> <ncbi_dataset.tsv>}"
accession_table="${2:?Usage: $0 <work_dir> <ncbi_dataset.tsv>}"
db_prefix="mycobacteriaceae_$(date +%F)"
threads="$(nproc)"

# Never hard-code your NCBI API key in a script: pass it through the environment.
api_key_args=()
if [[ -n "${NCBI_API_KEY:-}" ]]; then
    api_key_args=(--api-key "$NCBI_API_KEY")
fi

mkdir -p "$base_dir"/{ncbi_dataset,fasta}

# 1. Accession list (skip header, keep first column)
tail -n +2 "$accession_table" | cut -f 1 > "$base_dir/acc_list.txt"

# 2. Download each genome as a zip archive
get_genome() {
    local acc="$1"
    datasets download genome accession "$acc" "${api_key_args[@]}" \
        --filename "$base_dir/ncbi_dataset/$acc.zip"
}
export -f get_genome
export base_dir
export api_key_args
xargs -a "$base_dir/acc_list.txt" -P "$threads" -I {} bash -c 'get_genome "$@"' _ {}

# 3. Extract the fasta files into a single folder
decompress() {
    local sample
    sample="$(basename "$1" .zip)"
    unzip -q "$1" -d "$base_dir/ncbi_dataset/$sample"
}
export -f decompress
find "$base_dir/ncbi_dataset" -type f -name "*.zip" -print0 \
    | xargs -0 -P "$threads" -I {} bash -c 'decompress "$@"' _ {}

find "$base_dir/ncbi_dataset" -type f -name "*.fna" -exec mv -t "$base_dir/fasta" {} +
rm -rf "$base_dir/ncbi_dataset"

# 4. (Optional) bin genomes by species and dereplicate each bin at 99.9% identity, so that the
#    database stays small while keeping one representative per closely related group.
#    rename_and_bin_fasta.py is a helper that names/bins files from the species in the first header.
if command -v rename_and_bin_fasta.py >/dev/null 2>&1; then
    rename_and_bin_fasta.py -i "$base_dir/fasta" -o "$base_dir/binned_by_species"
    derep_bin() {
        local species
        species="$(basename "$1")"
        dereplicator.py --distance 0.001 --sketch_size 10000 --threads 4 "$1" "$2/$species"
    }
    export -f derep_bin
    find "$base_dir/binned_by_species" -mindepth 1 -maxdepth 1 -type d -print0 \
        | xargs -0 -P "$((threads / 4 > 0 ? threads / 4 : 1))" -I {} bash -c 'derep_bin "$@"' _ {} "$base_dir/derep0.001"
    db_input="$base_dir/derep0.001"
else
    db_input="$base_dir/fasta"
fi

# 5. Sketch the database
make_mashID_db -i "$db_input" -o "$base_dir" -p "$db_prefix" -s 10000 -k 21 -t "$threads"
echo "Database ready: $base_dir/$db_prefix.msh"
