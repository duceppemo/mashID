#!/usr/bin/env python3
"""Rebuild a finished run's summary from its per-sample tables and the database's current sidecars.

Usage: reannotate_summary.py RESULTS_DIR [--database DB.msh] [--dry-run]

The summary's ``Identification``, ``TaxID`` and ``Note`` columns are not measurements; they are
looked up from the sidecars when the summary is written. So when a sidecar changes - a reference is
relabelled, or ``make_mashID_db --confusable`` publishes a table the run predates - the summary can
be brought up to date from the hits already recorded in ``<sample>_mashID.tsv``, without screening
the reads again. Everything is produced by mashID's own ``annotate()``, ``summary_row()`` and
``write_*`` functions, so the result is what a fresh run would write; the screening itself cannot
change, and the per-sample tables are left alone.

Unlike the other scripts here this one imports mashid rather than standing alone, which is the
point: reimplementing the annotation would no longer agree with mashID by construction. Install the
package first (``pip install -e .``).

The rewritten files record how they were made under ``reannotated`` in summary_mashID.json. Each is
backed up to <name>.bak-<date> unless --no-backup is given; --dry-run writes nothing.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import shutil
import sys
from pathlib import Path

try:
    from mashid import MashIDError
    from mashid import confusable as confusable_mod
    from mashid.mash import ScreenHit
    from mashid.metadata import read_metadata, sidecar_path
    from mashid.pipeline import (
        JSON_FILENAME,
        MULTIQC_FILENAME,
        RUN_INFO_FILENAME,
        SUMMARY_COLUMNS,
        SUMMARY_FILENAME,
        SampleResult,
        annotate,
        summary_row,
        write_json,
        write_multiqc,
        write_tsv,
    )
    from mashid.samples import Sample
except ImportError as exc:  # an installed mashid older than this script, or none at all
    sys.exit(f"error: this script needs the mashid package it ships with, but importing it failed "
             f"({exc}).\nInstall the matching version from the repository root: pip install -e .")

PER_SAMPLE_SUFFIX = "_mashID.tsv"
REWRITTEN = (SUMMARY_FILENAME, MULTIQC_FILENAME, JSON_FILENAME)


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_hits(path: Path) -> list[ScreenHit]:
    """Rebuild the screen hits recorded in a per-sample table, keeping their recorded order.

    The table is already sorted, already truncated to --n-hits and already past the
    --min-ref-length filter, so re-reading it reproduces the hit list annotate() saw.
    """
    hits = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"Identity", "Shared_Hashes", "Median_Multiplicity", "P_Value", "Accession"} - set(
            reader.fieldnames or [])
        if missing:
            raise MashIDError(f"{path}: missing column(s) {', '.join(sorted(missing))}")
        for row in reader:
            accession = (row.get("Accession") or "").strip()
            if not accession:
                continue
            # query_id only has to carry the accession; accession_from_query_id() parses it back.
            hits.append(ScreenHit(row["Identity"], row["Shared_Hashes"], row["Median_Multiplicity"],
                                  row["P_Value"], f"{accession}.fna", row.get("Description") or ""))
    return hits


def _int_or_none(value: str | None) -> int | None:
    try:
        return int(value)                       # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None


def rebuild(results_dir: Path, database: Path | None = None) -> tuple[list[SampleResult], list[dict],
                                                                      list[str], Path, int]:
    """Re-annotate every sample of a results directory. Returns (results, rows, columns, db, changed)."""
    run_info_path = results_dir / RUN_INFO_FILENAME
    summary_path = results_dir / SUMMARY_FILENAME
    for required in (run_info_path, summary_path):
        if not required.is_file():
            raise MashIDError(f"Not a mashID results directory ({required.name} not found): {results_dir}")
    run_info = json.loads(run_info_path.read_text())
    params = run_info.get("parameters", {})
    database = database or Path(run_info.get("database", {}).get("path", ""))
    if not database.is_file():
        raise MashIDError(f"Database not found: {database or '(none recorded in the run)'}; pass --database")

    metadata = read_metadata(sidecar_path(database))
    confusable = confusable_mod.load(database, None)
    if confusable is None:
        print(f"note: no {confusable_mod.confusable_path(database).name} beside the database, so no "
              f"'not separable' notes; write one with make_mashID_db --confusable", file=sys.stderr)

    # Sequences/Bases and any sample-sheet columns live only in the summary; carry them over.
    with summary_path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        columns = list(reader.fieldnames or SUMMARY_COLUMNS)
        existing = {row["Sample"]: row for row in reader}
    extra_columns = [c for c in columns if c not in SUMMARY_COLUMNS]

    results, changed = [], 0
    for entry in run_info.get("samples", []):
        name = entry["name"]
        per_sample = results_dir / f"{name}{PER_SAMPLE_SUFFIX}"
        if not per_sample.is_file():
            raise MashIDError(f"Per-sample table missing for {name}: {per_sample}")
        old = existing.get(name, {})
        sample = Sample(name=name,
                        files=[Path(f) for f in entry.get("files", [])],
                        seq_type=entry.get("type", "fastq"),
                        sequences=_int_or_none(old.get("Sequences")),
                        bases=_int_or_none(old.get("Bases")),
                        extra={c: old.get(c, "") for c in extra_columns})
        result = SampleResult(sample=sample, hits=read_hits(per_sample))
        result.notes = annotate(result, metadata, params.get("ambiguity_margin", 0.005),
                                params.get("winner_take_all", True), confusable)
        results.append(result)
        if "; ".join(result.notes) != (old.get("Note") or ""):
            changed += 1
    return results, [summary_row(r, metadata) for r in results], columns, database, changed


def backup(path: Path) -> Path:
    stamp = dt.date.today().strftime("%Y%m%d")
    dest = path.with_name(f"{path.name}.bak-{stamp}")
    n = 2
    while dest.exists():                        # never clobber an older backup
        dest = path.with_name(f"{path.name}.bak-{stamp}.{n}")
        n += 1
    shutil.copy2(path, dest)
    return dest


def provenance(database: Path) -> dict:
    record: dict = {
        "at": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        "reason": "Identification, TaxID and Note regenerated from the recorded hits; the run "
                  "predates the sidecars below",
        "screening_unchanged": True,
    }
    for key, path in (("metadata_sidecar", sidecar_path(database)),
                      ("confusable_sidecar", confusable_mod.confusable_path(database))):
        record[key] = {"path": str(path), "md5": md5(path)} if path.is_file() else None
    return record


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("results", metavar="RESULTS_DIR", type=Path,
                    help="A mashID output directory, holding mashID_run.json and the per-sample tables.")
    ap.add_argument("--database", metavar="DB.msh", type=Path, default=None,
                    help="Database to read the sidecars from. Default: the one recorded in the run.")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would change and write nothing.")
    args = ap.parse_args(argv)

    try:
        results, rows, columns, database, changed = rebuild(args.results, args.database)
    except MashIDError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"{len(rows)} sample(s); {changed} with a changed Note; "
          f"{sum(1 for r in rows if r['Note'])} with any note")
    if args.dry_run:
        for row in rows[:3]:
            print(f"  {row['Sample']}: {row['Note'] or '(none)'}")
        print("dry run: nothing written")
        return 0

    for filename in REWRITTEN:
        target = args.results / filename
        if target.is_file() and not args.no_backup:
            print(f"backup: {backup(target).name}")

    write_tsv(args.results / SUMMARY_FILENAME, columns, rows)
    write_multiqc(args.results / MULTIQC_FILENAME, rows)
    json_path = args.results / JSON_FILENAME
    write_json(json_path, results, rows, read_metadata(sidecar_path(database)), database)
    doc = json.loads(json_path.read_text())
    doc["reannotated"] = provenance(database)
    json_path.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"rewrote {', '.join(REWRITTEN)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
