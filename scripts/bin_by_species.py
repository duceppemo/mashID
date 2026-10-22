#!/usr/bin/env python3
"""Bin NCBI ``datasets`` genomes into one directory per species, using assembly_data_report.jsonl.

Usage: bin_by_species.py assembly_data_report.jsonl ncbi_dataset/data binned/

Creates binned/<Genus_species>/<accession>_....fna[.gz] as symbolic links (or copies with --copy).
"Genus sp." assemblies of one genus go to a single "Genus_sp" bin; "Candidatus" is kept in the name.
With --rank subspecies, bins are split further by "subsp. X" (and "serovar X" with --rank serovar).
--max-bin N caps every bin at N assemblies, keeping complete genomes first, then chromosome,
scaffold and contig level, so dereplication of huge species (Salmonella, E. coli) stays tractable.
--corrections overrides the organism name NCBI reports for named accessions, so a genome deposited
under the wrong species is binned (and dereplicated) with its true relatives.
A summary of bins and sizes is printed; the accession -> bin table is written to binned/bins.tsv.

Bin names and the corrections table are read through mashid's own helpers, so that this script,
apply_label_corrections.py and mashID itself cannot disagree about them; install the package from
the repository root first (``pip install -e .``).
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

try:
    from mashid import MashIDError
    from mashid.names import BIN_RANKS, species_of
    from mashid.names import read_corrections as _read_corrections
except ImportError as exc:  # an installed mashid older than this script, or none at all
    sys.exit(f"error: this script needs the mashid package it ships with, but importing it failed "
             f"({exc}).\nInstall the matching version from the repository root: pip install -e .")

__all__ = ["load_corrections", "main", "species_of"]


def load_corrections(path: Path) -> dict[str, tuple[str, str]]:
    """``{accession: (organism, taxid)}`` from a curated table, exiting cleanly on a bad one."""
    try:
        return _read_corrections(path)
    except MashIDError as exc:
        sys.exit(str(exc))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("report", type=Path, help="assembly_data_report.jsonl from the datasets archive")
    ap.add_argument("data_dir", type=Path, help="ncbi_dataset/data directory holding <accession>/ folders")
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--copy", action="store_true", help="Copy files instead of symlinking")
    ap.add_argument("--rank", choices=list(BIN_RANKS), default="species",
                    help="Bin by species (default), by subspecies, or by subspecies and serovar")
    ap.add_argument("--max-bin", type=int, default=0, metavar="N",
                    help="Keep at most N assemblies per bin, preferring higher assembly levels (0 = no cap)")
    ap.add_argument("--corrections", metavar="FILE.tsv", type=Path, default=None,
                    help="TSV/CSV with columns Accession, Organism and optionally TaxID, overriding "
                         "the organism NCBI reports for those accessions")
    args = ap.parse_args()

    corrections = load_corrections(args.corrections) if args.corrections else {}
    applied: set[str] = set()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    level_rank = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
    records: list[tuple[str, str, str, str, int, list[Path]]] = []  # acc, organism, taxid, bin, level, files
    missing = 0
    with open(args.report) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rec = json.loads(line)
            acc = rec.get("accession") or rec.get("current_accession")
            org = rec.get("organism") or {}
            organism = org.get("organismName") or org.get("organism_name") or "unknown"
            taxid = str(org.get("taxId") or org.get("tax_id") or "")
            if acc in corrections:
                organism, corrected_taxid = corrections[acc]
                taxid = corrected_taxid or taxid
                applied.add(acc)
            info = rec.get("assemblyInfo") or rec.get("assembly_info") or {}
            level = level_rank.get(info.get("assemblyLevel") or info.get("assembly_level") or "", 4)
            files = sorted(p for p in (args.data_dir / acc).glob("*")
                           if re.search(r"\.(fna|fa|fasta)(\.gz)?$", p.name))
            if not files:
                missing += 1
                continue
            records.append((acc, organism, taxid, species_of(organism, args.rank), level, files))

    by_bin: dict[str, list] = {}
    for r in records:
        by_bin.setdefault(r[3], []).append(r)
    dropped = 0
    counts: Counter[str] = Counter()
    with open(args.out_dir / "bins.tsv", "w") as table:
        table.write("Accession\tOrganism\tTaxID\tBin\n")
        for bin_name, recs in by_bin.items():
            recs.sort(key=lambda r: (r[4], r[0]))  # best assembly level first, then accession
            if args.max_bin and len(recs) > args.max_bin:
                dropped += len(recs) - args.max_bin
                recs = recs[: args.max_bin]
            bin_dir = args.out_dir / bin_name
            bin_dir.mkdir(exist_ok=True)
            for acc, organism, taxid, _, _, files in recs:
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
    if dropped:
        print(f"{dropped} assemblies dropped by --max-bin {args.max_bin}", file=sys.stderr)
    if corrections:
        print(f"{len(applied)} of {len(corrections)} label correction(s) applied", file=sys.stderr)
        for acc in sorted(set(corrections) - applied):
            print(f"  unused correction: {acc} not in {args.report.name}", file=sys.stderr)

    total = sum(counts.values())
    print(f"{total} assemblies in {len(counts)} bins; {missing} accession(s) without a sequence file",
          file=sys.stderr)
    for name, n in counts.most_common(15):
        print(f"  {n:6d}  {name}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
