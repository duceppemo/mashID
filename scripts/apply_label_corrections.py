#!/usr/bin/env python3
"""Apply curated label corrections to an already-built database's tables, without rebuilding.

Usage: apply_label_corrections.py corrections.tsv <db>.metadata.tsv [bins.tsv ...]

Rewrites the Organism and TaxID of the named accessions in each target table. Where a table has a
Bin column (bins.tsv) the bin is recomputed from the corrected name at --rank, which must match the
rank the bins were built at; where it has a Description column the correction is noted inline so the
deposited name stays visible; and where it has a Confusable_Organism column (a .confusable.tsv
table) that column is corrected too. Each target is backed up to <name>.bak-<date> unless
--no-backup is given; --dry-run reports what would change and writes nothing.

This fixes reporting for a deployed database. It does not touch the .msh sketches, so the
reference's own FASTA header (the Description shown by mash) keeps the deposited name until
the database is rebuilt with make_mashID_db --metadata corrections.tsv.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    from mashid import MashIDError
    from mashid.fileutil import backup as make_backup
    from mashid.names import BIN_RANKS, normalise_header, species_of
    from mashid.names import read_corrections as _read_corrections
except ImportError as exc:  # an installed mashid older than this script, or none at all
    sys.exit(f"error: this script needs the mashid package it ships with, but importing it failed "
             f"({exc}).\nInstall the matching version from the repository root: pip install -e .")

NOTE = " [mashID: deposited as {old!r}; corrected to {new!r}, see curation/label_corrections.tsv]"


def bin_of(organism: str, rank: str = "species") -> str:
    """The bin name for an organism: mashID's own ``names.species_of``, one implementation."""
    return species_of(organism, rank)


def read_corrections(path: Path) -> dict[str, tuple[str, str]]:
    """``{accession: (organism, taxid)}`` from a curated table, exiting cleanly on a bad one."""
    try:
        return _read_corrections(path)
    except MashIDError as exc:
        sys.exit(str(exc))


def apply_to(target: Path, corrections: dict[str, tuple[str, str]], backup: bool, dry_run: bool,
             rank: str = "species") -> tuple[int, set[str]]:
    """Rewrite the corrected rows of one table. Returns (rows changed, accessions matched)."""
    try:
        text = target.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError as exc:
        sys.exit(f"{target} is not a UTF-8 text table: {exc}")
    # Split on newlines only. str.splitlines() also breaks on \r, \x0b, \x0c, \x1c-\x1e, \x85 and
    # U+2028, any of which can appear inside a Description field and would mangle the row.
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        print(f"{target}: empty, skipped", file=sys.stderr)
        return 0, set()
    header = lines[0].rstrip("\r").split("\t")
    norm = normalise_header(header)
    idx = {key: header.index(name) for key, name in norm.items()}
    acc_i, org_i = idx.get("accession"), idx.get("organism")
    if acc_i is None or org_i is None:
        sys.exit(f"{target}: expected 'Accession' and 'Organism' columns; found {header}")
    tax_i, bin_i, desc_i = idx.get("taxid"), idx.get("bin"), idx.get("description")
    conf_i = idx.get("confusableorganism")

    near_i = idx.get("nearestaccession")
    out, changed, matched = [header], 0, set()
    for ln in lines[1:]:
        fields = ln.rstrip("\r").split("\t")
        acc = fields[acc_i] if acc_i < len(fields) else ""
        row_changed = False
        if acc in corrections:
            matched.add(acc)
            # Fill a ragged row out to the header before writing to it. Skipping the short ones
            # instead reported the correction as unused, and correcting only the fields that
            # happened to be present left a new name beside a stale TaxID.
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))
            new_org, new_tax = corrections[acc]
            old_org = fields[org_i]
            if old_org != new_org:
                fields[org_i] = new_org
                row_changed = True
                if tax_i is not None and new_tax:
                    fields[tax_i] = new_tax
                if bin_i is not None:
                    # The old bin says which rank this table was binned at. Recomputing at the
                    # wrong one leaves this row at a different rank from every other, and the Bin
                    # column is what the next rebuild dereplicates by.
                    expected = bin_of(old_org, rank)
                    if fields[bin_i] and fields[bin_i] != expected:
                        sys.exit(f"{target}: {acc} has Bin {fields[bin_i]!r}, but --rank {rank} "
                                 f"gives {expected!r} for its current name. Pass the --rank the "
                                 f"table was binned at (bin_by_species.py --rank); nothing written.")
                    fields[bin_i] = bin_of(new_org, rank)
                if desc_i is not None:
                    note = NOTE.format(old=old_org, new=new_org)
                    if "[mashID:" not in fields[desc_i]:
                        fields[desc_i] += note
                print(f"  {acc}: {old_org!r} -> {new_org!r}"
                      + (f" (bin {fields[bin_i]})" if bin_i is not None else ""))
        # A .confusable.tsv names the other organism in a second column, and that is the one mashID
        # actually reports; correcting only Organism would leave the note itself stale.
        if conf_i is not None and near_i is not None and max(conf_i, near_i) < len(fields):
            near = fields[near_i]
            if near in corrections and fields[conf_i] != corrections[near][0]:
                matched.add(near)
                print(f"  {acc}: confusable {fields[conf_i]!r} -> {corrections[near][0]!r} "
                      f"(from {near})")
                fields[conf_i] = corrections[near][0]
                row_changed = True
        changed += row_changed
        out.append(fields)

    if changed and not dry_run:
        if backup:
            print(f"  backup: {make_backup(target).name}")
        target.write_text("".join("\t".join(f) + "\n" for f in out), encoding="utf-8")
    return changed, matched


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("corrections", type=Path)
    ap.add_argument("targets", type=Path, nargs="+", help="<db>.metadata.tsv and/or bins.tsv")
    ap.add_argument("--rank", choices=list(BIN_RANKS), default="species",
                    help="Rank the Bin column of a bins.tsv was built at, so a corrected bin keeps "
                         "it (bin_by_species.py --rank). Default: %(default)s")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    corrections = read_corrections(args.corrections)
    print(f"{len(corrections)} correction(s) loaded from {args.corrections}")
    total, matched = 0, set()
    for target in args.targets:
        if not target.is_file():
            sys.exit(f"{target}: not found")
        print(f"{target}:")
        n, hit = apply_to(target, corrections, not args.no_backup, args.dry_run, args.rank)
        print(f"  {n} row(s) {'would change' if args.dry_run else 'changed'}")
        total += n
        matched |= hit
    # A correction that matched nothing is worth saying: the accession may be spelled differently,
    # or this may not be the table it was meant for. But one curation file covers every database,
    # so most of its rows legitimately miss any single target; summarise rather than list each.
    unused = sorted(set(corrections) - matched)
    if unused:
        shown = ", ".join(unused[:4]) + (", ..." if len(unused) > 4 else "")
        print(f"{len(unused)} of {len(corrections)} correction(s) matched no row of any target; "
              f"normal when one curation file covers several databases: {shown}", file=sys.stderr)
    if args.dry_run:
        print(f"\ndry run: {total} row(s) would change, nothing written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
