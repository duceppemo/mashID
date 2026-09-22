#!/usr/bin/env python3
"""Split a multi-genome fasta (e.g. proGenomes representatives) into one file per genome.

Reads fasta from stdin. Records whose header starts with ``<assembly accession>_<sequence id>``
(``GCA_000005825.2_CP001878.2 Bacillus ...``) are grouped by the assembly accession into
``<out_dir>/<accession>.fna``; other records are grouped by their whole first token.
Also writes ``<out_dir>/accessions.txt`` (one accession per line) and prints counts.

Usage: zcat genomes.fna.gz | split_multifasta_by_genome.py out_dir
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

_ACC = re.compile(rb"^>((?:GC[AF]_\d{9}\.\d+))_\S+")


def genome_of(header: bytes) -> bytes:
    m = _ACC.match(header)
    if m:
        return m.group(1)
    return header[1:].split()[0] if len(header) > 1 else b"unnamed"


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    out_dir = Path(sys.argv[1])
    out_dir.mkdir(parents=True, exist_ok=True)
    seen: set[bytes] = set()
    current: bytes | None = None
    fh = None
    n_seq = 0
    stdin = sys.stdin.buffer
    for line in stdin:
        if line.startswith(b">"):
            n_seq += 1
            acc = genome_of(line)
            if acc != current:
                if fh:
                    fh.close()
                mode = "ab" if acc in seen else "wb"
                seen.add(acc)
                fh = open(out_dir / (acc.decode() + ".fna"), mode, buffering=1 << 20)  # noqa: SIM115
                current = acc
        if fh is None:
            continue  # sequence lines before any header
        fh.write(line)
    if fh:
        fh.close()
    with open(out_dir / "accessions.txt", "w") as af:
        af.write("\n".join(sorted(a.decode() for a in seen)) + "\n")
    print(f"{n_seq} sequences written into {len(seen)} genome files under {out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
