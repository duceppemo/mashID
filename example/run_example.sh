#!/usr/bin/env bash
# Run mashID on the example dataset and compare the summary with the expected result.
# Usage: bash example/run_example.sh [--update-expected]
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="${TMPDIR:-/tmp}/mashID_example_$$"

mashID -i "$here/reads" -o "$out" -d "$here/db/example.msh" -t 2 -p 2

if [[ "${1:-}" == "--update-expected" ]]; then
    mkdir -p "$here/expected"
    cp "$out"/summary_mashID.tsv "$out"/summary_mashID_mqc.tsv "$out"/*_mashID.tsv "$here/expected/"
    echo "Expected output updated in $here/expected"
elif diff "$here/expected/summary_mashID.tsv" "$out/summary_mashID.tsv" \
        && diff "$here/expected/summary_mashID_mqc.tsv" "$out/summary_mashID_mqc.tsv"; then
    echo "OK: results match $here/expected/summary_mashID.tsv"
else
    echo "MISMATCH: compare $out with $here/expected" >&2
    exit 1
fi
rm -rf "$out"
