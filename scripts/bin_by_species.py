#!/usr/bin/env python3
"""Bin NCBI ``datasets`` genomes into one directory per species, using assembly_data_report.jsonl.

Usage: bin_by_species.py assembly_data_report.jsonl ncbi_dataset/data binned/

Creates binned/<Genus_species>/<accession>_....fna[.gz] as symbolic links (or copies with --copy).
"Genus sp." assemblies of one genus go to a single "Genus_sp" bin; "Candidatus" is kept in the name.
A summary of bins and sizes is printed; the accession -> species table is written to binned/bins.tsv.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
from collections import Counter
from pathlib import Path


def species_of(organism: str) -> str:
    tokens = organism.replace("[", "").replace("]", "").split()
    if not tokens:
        return "unknown"
    if tokens[0] == "Candidatus" and len(tokens) >= 3:
        tokens = ["Candidatus_" + tokens[1]] + tokens[2:]
    genus = tokens[0]
    species = tokens[1] if len(tokens) > 1 else "sp"
    if species in ("sp.", "sp", "spp.", "bacterium"):
        species = "sp"
    name = f"{genus}_{species}"
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("report", type=Path, help="assembly_data_report.jsonl from the datasets archive")
    ap.add_argument("data_dir", type=Path, help="ncbi_dataset/data directory holding <accession>/ folders")
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--copy", action="store_true", help="Copy files instead of symlinking")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    missing = 0
    with open(args.report) as fh, open(args.out_dir / "bins.tsv", "w") as table:
        table.write("Accession\tOrganism\tTaxID\tBin\n")
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            acc = rec["accession"]
            organism = (rec.get("organism") or {}).get("organismName", "unknown")
            taxid = (rec.get("organism") or {}).get("taxId", "")
            files = sorted(p for p in (args.data_dir / acc).glob("*")
                           if re.search(r"\.(fna|fa|fasta)(\.gz)?$", p.name))
            if not files:
                missing += 1
                continue
            bin_name = species_of(organism)
            bin_dir = args.out_dir / bin_name
            bin_dir.mkdir(exist_ok=True)
            for src in files:
                dst = bin_dir / src.name
                if dst.exists() or dst.is_symlink():
                    continue
                if args.copy:
                    shutil.copy2(src, dst)
                else:
                    os.symlink(src.resolve(), dst)
            counts[bin_name] += 1
            table.write(f"{acc}\t{organism}\t{taxid}\t{bin_name}\n")

    total = sum(counts.values())
    print(f"{total} assemblies in {len(counts)} bins; {missing} accession(s) without a sequence file",
          file=sys.stderr)
    for name, n in counts.most_common(15):
        print(f"  {n:6d}  {name}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
