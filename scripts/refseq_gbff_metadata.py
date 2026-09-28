#!/usr/bin/env python3
"""Extract a mashID metadata table from GenBank flat files (e.g. RefSeq release viral.*.genomic.gbff.gz).

Writes a TSV with Accession (VERSION), Organism, TaxID (db_xref taxon), Length, Description (DEFINITION)
usable with `make_mashID_db --metadata`. Reads .gbff or .gbff.gz, several files allowed.

Usage: refseq_gbff_metadata.py out.tsv file1.gbff.gz [file2.gbff.gz ...]
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from pathlib import Path

_TAXON = re.compile(r'/db_xref="taxon:(\d+)"')


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def records(path: Path):
    """Yield (accession, organism, taxid, length, definition) per GenBank record."""
    acc = organism = taxid = definition = ""
    length = 0
    in_definition = False
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                parts = line.split()
                length = int(parts[2]) if len(parts) > 2 and parts[2].isdigit() else 0
                acc = organism = taxid = definition = ""
                in_definition = False
            elif line.startswith("DEFINITION"):
                definition = line[12:].strip()
                in_definition = True
            elif in_definition and line.startswith(" "):
                definition += " " + line.strip()
            elif line.startswith("VERSION"):
                in_definition = False
                acc = line.split()[1] if len(line.split()) > 1 else ""
            elif line.startswith("  ORGANISM"):
                in_definition = False
                organism = line[12:].strip()
            elif not taxid and "/db_xref=\"taxon:" in line:
                m = _TAXON.search(line)
                if m:
                    taxid = m.group(1)
            elif line.startswith("//"):
                if acc:
                    yield acc, organism, taxid, length, definition
                acc = ""
            else:
                if not line.startswith(" ") or not in_definition:
                    in_definition = False


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__, file=sys.stderr)
        return 2
    out = Path(sys.argv[1])
    n = 0
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Accession", "Organism", "TaxID", "Length", "Description"])
        for f in sys.argv[2:]:
            for acc, organism, taxid, length, definition in records(Path(f)):
                w.writerow([acc, organism, taxid or "NA", length, definition])
                n += 1
    print(f"{n} records written to {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
