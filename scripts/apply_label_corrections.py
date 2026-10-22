#!/usr/bin/env python3
"""Apply curated label corrections to an already-built database's tables, without rebuilding.

Usage: apply_label_corrections.py corrections.tsv <db>.metadata.tsv [bins.tsv ...]

Rewrites the Organism and TaxID of the named accessions in each target table. Where a table
has a Bin column (bins.tsv) the bin is recomputed from the corrected name, and where it has
a Description column the correction is noted inline so the deposited name stays visible.
Each target is backed up to <name>.bak-<date> unless --no-backup is given; --dry-run reports
what would change and writes nothing.

This fixes reporting for a deployed database. It does not touch the .msh sketches, so the
reference's own FASTA header (the Description shown by mash) keeps the deposited name until
the database is rebuilt with make_mashID_db --metadata corrections.tsv.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import re
import shutil
import sys
from pathlib import Path

NOTE = " [mashID: deposited as {old!r}; corrected to {new!r}, see curation/label_corrections.tsv]"


def bin_of(organism: str) -> str:
    """Mirror scripts/bin_by_species.py species_of() at --rank species."""
    tokens = organism.replace("[", "").replace("]", "").split()
    if not tokens:
        return "unknown"
    if tokens[0] == "Candidatus" and len(tokens) >= 3:
        tokens = ["Candidatus_" + tokens[1]] + tokens[2:]
    species = tokens[1] if len(tokens) > 1 else "sp"
    if species in ("sp.", "sp", "spp.", "bacterium"):
        species = "sp"
    return re.sub(r"[^A-Za-z0-9_.:,\[\]-]", "_", f"{tokens[0]}_{species}")


def normalise_header(fieldnames) -> dict[str, str]:
    """Map lowercase, punctuation-free column names to the actual header names."""
    return {re.sub(r"[^a-z0-9]", "", f.lower()): f for f in (fieldnames or [])}


def read_corrections(path: Path) -> dict[str, tuple[str, str]]:
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="," if path.suffix.lower() == ".csv" else "\t")
        cols = normalise_header(reader.fieldnames)
        acc_col = cols.get("accession") or cols.get("assemblyaccession") or cols.get("acc")
        org_col = (cols.get("organism") or cols.get("organismname") or cols.get("species")
                   or cols.get("name"))
        if not acc_col or not org_col:
            sys.exit(f"{path}: expected columns 'Accession' and 'Organism'; found {reader.fieldnames}")
        tax_col = cols.get("taxid") or cols.get("taxonomyid")
        out = {}
        for row in reader:
            acc = (row.get(acc_col) or "").strip()
            organism = (row.get(org_col) or "").strip()
            if acc and organism:
                out[acc] = (organism, (row.get(tax_col) or "").strip() if tax_col else "")
    return out


def apply_to(target: Path, corrections: dict[str, tuple[str, str]],
             backup: bool, dry_run: bool) -> int:
    lines = target.read_text().splitlines(keepends=True)
    if not lines:
        print(f"{target}: empty, skipped", file=sys.stderr)
        return 0
    header = lines[0].rstrip("\n").split("\t")
    norm = normalise_header(header)
    idx = {key: header.index(name) for key, name in norm.items()}
    acc_i, org_i = idx.get("accession"), idx.get("organism")
    if acc_i is None or org_i is None:
        sys.exit(f"{target}: expected 'Accession' and 'Organism' columns; found {header}")
    tax_i, bin_i, desc_i = idx.get("taxid"), idx.get("bin"), idx.get("description")

    out, changed = [lines[0]], 0
    for ln in lines[1:]:
        fields = ln.rstrip("\n").split("\t")
        acc = fields[acc_i] if acc_i < len(fields) else ""
        if acc in corrections and len(fields) > org_i:
            new_org, new_tax = corrections[acc]
            old_org = fields[org_i]
            if old_org != new_org:
                fields[org_i] = new_org
                if tax_i is not None and new_tax and tax_i < len(fields):
                    fields[tax_i] = new_tax
                if bin_i is not None and bin_i < len(fields):
                    fields[bin_i] = bin_of(new_org)
                if desc_i is not None and desc_i < len(fields):
                    note = NOTE.format(old=old_org, new=new_org)
                    if "[mashID:" not in fields[desc_i]:
                        fields[desc_i] += note
                print(f"  {acc}: {old_org!r} -> {new_org!r}"
                      + (f" (bin {fields[bin_i]})" if bin_i is not None else ""))
                changed += 1
            ln = "\t".join(fields) + "\n"
        out.append(ln)

    if changed and not dry_run:
        if backup:
            stamp = dt.date.today().strftime("%Y%m%d")
            dest = target.with_name(f"{target.name}.bak-{stamp}")
            n = 2
            while dest.exists():            # never clobber an older backup
                dest = target.with_name(f"{target.name}.bak-{stamp}.{n}")
                n += 1
            shutil.copy2(target, dest)
            print(f"  backup: {dest.name}")
        target.write_text("".join(out))
    return changed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("corrections", type=Path)
    ap.add_argument("targets", type=Path, nargs="+", help="<db>.metadata.tsv and/or bins.tsv")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    corrections = read_corrections(args.corrections)
    print(f"{len(corrections)} correction(s) loaded from {args.corrections}")
    total = 0
    for target in args.targets:
        if not target.is_file():
            sys.exit(f"{target}: not found")
        print(f"{target}:")
        n = apply_to(target, corrections, not args.no_backup, args.dry_run)
        print(f"  {n} row(s) {'would change' if args.dry_run else 'changed'}")
        total += n
    if args.dry_run:
        print(f"\ndry run: {total} row(s) would change, nothing written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
