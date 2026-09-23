#!/usr/bin/env python3
"""Check a database built by build_taxon_db.sh against genomes that were left out of it.

The build directory holds every downloaded assembly (ncbi/ncbi_dataset/data/<accession>/) and
binned/bins.tsv lists the ones that went into the database. Assemblies dropped by --max-bin or by
dereplication are labelled genomes the database has never seen: this script screens a sample of them
with mashID and reports how often the reported organism matches the NCBI name, per species.

Usage: validate_db_heldout.py <work_dir> <database.msh> [--per-species N] [--threads T] [--rank species|subspecies]
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import random
import re
import subprocess
import sys
import tempfile
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def species_key(name: str, rank: str) -> str:
    toks = name.replace("[", "").replace("]", "").split()
    if toks and toks[0] == "Candidatus":
        toks = toks[1:]
    key = " ".join(toks[:2])
    if rank == "subspecies":
        m = re.search(r"subsp\. (\S+)", name)
        if m:
            key += f" subsp. {m.group(1)}"
    return key


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("work_dir", type=Path)
    ap.add_argument("database", type=Path)
    ap.add_argument("--per-species", type=int, default=20,
                    help="held-out genomes sampled per species (default 20)")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--rank", choices=["species", "subspecies"], default="species")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, default=None, help="keep mashID output here (default: temporary)")
    args = ap.parse_args()

    report = args.work_dir / "ncbi" / "ncbi_dataset" / "data" / "assembly_data_report.jsonl"
    bins = args.work_dir / "binned" / "bins.tsv"
    in_db = {r["Accession"] for r in read_tsv(bins)}
    # the sidecar knows exactly which accessions are in the database (dereplication drops some)
    sidecar = args.database.with_name(args.database.name[:-4] + ".metadata.tsv")
    if sidecar.is_file():
        in_db = {r["Accession"] for r in read_tsv(sidecar)}

    held: dict[str, list[tuple[str, str]]] = collections.defaultdict(list)
    with open(report) as fh:
        records = [json.loads(line) for line in fh if line.strip()]
    for rec in records:
        acc = rec.get("accession") or rec.get("current_accession")
        org = rec.get("organism") or {}
        name = org.get("organismName") or org.get("organism_name") or ""
        if not acc or acc in in_db or not name:
            continue
        key = species_key(name, args.rank)
        if key.endswith(" sp.") or key.startswith("uncultured"):
            continue
        if list((args.work_dir / "ncbi" / "ncbi_dataset" / "data" / acc).glob("*.fna.gz")):
            held[key].append((acc, name))

    rng = random.Random(args.seed)
    chosen: list[tuple[str, str, str]] = []
    for key, items in sorted(held.items()):
        for acc, name in rng.sample(items, min(args.per_species, len(items))):
            chosen.append((key, acc, name))
    if not chosen:
        print("no held-out genomes found", file=sys.stderr)
        return 1
    print(f"{len(chosen)} held-out genomes from {len(held)} {args.rank} groups", file=sys.stderr)

    tmp = Path(tempfile.mkdtemp(prefix="mashid_heldout_"))
    inp = tmp / "input"
    inp.mkdir()
    truth: dict[str, tuple[str, str]] = {}
    for key, acc, name in chosen:
        src = next((args.work_dir / "ncbi" / "ncbi_dataset" / "data" / acc).glob("*.fna.gz"))
        (inp / f"{acc}.fna.gz").symlink_to(src.resolve())
        truth[acc] = (key, name)
    out = args.out or (tmp / "out")
    cmd = ["mashID", "-i", str(inp), "-o", str(out), "-d", str(args.database), "-t", str(args.threads),
           "-p", str(max(1, args.threads // 4))]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    rows = read_tsv(out / "summary_mashID.tsv")
    per: dict[str, list[int]] = collections.defaultdict(lambda: [0, 0, 0])  # correct, total, noted
    wrong: list[tuple[str, str, str, str]] = []
    for r in rows:
        m = re.match(r"(GC[AF]_\d+\.\d+)", r["Sample"])
        acc = m.group(1) if m else r["Sample"]
        key, name = truth[acc]
        called = species_key(r["Identification"], args.rank)
        ok = called == key
        per[key][1] += 1
        per[key][0] += ok
        per[key][2] += bool(r["Note"])
        if not ok:
            wrong.append((acc, key, r["Identification"], r["Note"]))

    total = sum(v[1] for v in per.values())
    correct = sum(v[0] for v in per.values())
    print(f"\nConcordance at {args.rank} level: {correct}/{total} = {100 * correct / total:.1f}%\n")
    print(f"{'group':45s} {'correct':>8s} {'tested':>7s} {'noted':>6s}")
    for key, (c, t, n) in sorted(per.items()):
        print(f"{key:45s} {c:8d} {t:7d} {n:6d}")
    if wrong:
        print("\nMismatches (accession, NCBI name, mashID call, note):")
        for acc, key, call, note in wrong:
            print(f"  {acc}  {key}  ->  {call}  {note}")
    print(f"\nmashID output: {out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
