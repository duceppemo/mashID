#!/usr/bin/env python3
"""Rebuild a finished run's summary from its per-sample tables and the current reference names.

Usage: reannotate_summary.py RESULTS_DIR [--database DB.msh] [--dry-run]

The summary's ``Identification``, ``TaxID`` and ``Note`` columns are not measurements; they are
looked up from the sidecars when the summary is written. So when a sidecar changes - a reference is
relabelled, or ``make_mashID_db --confusable`` publishes a table the run predates - the summary can
be brought up to date from the hits already recorded in ``<sample>_mashID.tsv``, without screening
the reads again. Everything is produced by mashID's own ``annotate()``, ``summary_row()`` and
``write_*`` functions, so the result is what a fresh run would write; the screening itself cannot
change, and the per-sample tables are left alone.

The tables read are the ones the run itself used, including any ``--db-metadata`` /
``--db-confusable`` it recorded; ``--database`` reads a given database's own sidecars instead, for a
database that has moved. A run already in step with those tables is left untouched.

Unlike the other scripts here this one imports mashid rather than standing alone, which is the
point: reimplementing the annotation would no longer agree with mashID by construction. Install the
package first (``pip install -e .``).

The rewritten files record how they were made under ``reannotated`` in summary_mashID.json. Each is
backed up to <name>.bak-<date> unless --no-backup is given; --dry-run writes nothing.

The summary being replaced is treated as the authority on which samples must survive: a rebuild that
would drop a row it holds -- an interrupted run, a hand-edited mashID_run.json -- is refused rather
than reported as a success.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import json
import sys
from collections import Counter
from pathlib import Path
from typing import NamedTuple

try:
    from mashid import MashIDError
    from mashid import confusable as confusable_mod
    from mashid.fileutil import backup, md5sum
    from mashid.mash import ScreenHit
    from mashid.metadata import read_metadata, sidecar_path
    from mashid.pipeline import (
        JSON_FILENAME,
        MULTIQC_FILENAME,
        RUN_INFO_FILENAME,
        SUMMARY_COLUMNS,
        SUMMARY_FILENAME,
        SampleResult,
        Settings,
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

# Columns this script recomputes. Everything else in the summary is a measurement and is carried
# over from the file being replaced.
LOOKED_UP_COLUMNS = ("Identification", "TaxID", "Note")
DRY_RUN_SAMPLES = 20


def _read_columns(summary: Path) -> list[str]:
    with summary.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t").fieldnames or [])


def _default(option: str):
    """A mashID option's default, read off ``Settings`` rather than copied.

    A run.json written before an option existed has no entry for it, and re-annotating it under a
    stale copy of the default would disagree with what mashID itself would do today.
    """
    field = next(f for f in dataclasses.fields(Settings) if f.name == option)
    return field.default


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


class Rebuilt(NamedTuple):
    """What a re-annotation produced, and what it was produced from."""

    results: list[SampleResult]
    rows: list[dict]
    columns: list[str]
    database: Path
    changed: list[str]              # samples whose Identification, TaxID or Note would move
    metadata_table: Path            # the sidecars actually read, which --db-metadata can move
    confusable_table: Path | None


def rebuild(results_dir: Path, database: Path | None = None) -> Rebuilt:
    """Re-annotate every sample of a results directory, from the tables the run itself used."""
    database_arg = database
    run_info_path = results_dir / RUN_INFO_FILENAME
    summary_path = results_dir / SUMMARY_FILENAME
    for required in (run_info_path, summary_path):
        if not required.is_file():
            raise MashIDError(f"Not a mashID results directory ({required.name} not found): {results_dir}")
    try:
        run_info = json.loads(run_info_path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise MashIDError(f"{run_info_path} is not readable JSON: {exc}") from exc
    if not isinstance(run_info, dict):
        raise MashIDError(f"{run_info_path} does not hold a run record")
    params = run_info.get("parameters") or {}
    if database is None:
        recorded = (run_info.get("database") or {}).get("path") or ""
        if not recorded:
            raise MashIDError(f"{run_info_path.name} records no database path; pass --database")
        database = Path(recorded)
    if not database.is_file():
        raise MashIDError(f"Database not found: {database}; pass --database")

    # A run made with --db-metadata/--db-confusable was annotated from those tables, not from the
    # ones beside the .msh. Re-annotating it from the sidecars would undo exactly what the override
    # was for -- a corrected sidecar kept outside the database directory, say -- and report the
    # regression as a change. Honour what the run recorded unless --database overrides it.
    overrides: dict[str, Path | None] = {"metadata": None, "confusable": None}
    for key in overrides:
        recorded = params.get(f"db_{key}")
        if not recorded:
            continue
        if database_arg is not None:
            # --database means "read this database's own sidecars", so it deliberately drops the
            # override -- but silently dropping it is how the run's labels would quietly regress.
            print(f"note: the run used --db-{key} {recorded}; --database overrides it, so the "
                  f"names come from the given database's sidecars instead", file=sys.stderr)
            continue
        if not Path(recorded).is_absolute():
            raise MashIDError(f"the run recorded a relative --db-{key} ({recorded}), which cannot "
                              f"be resolved from here; pass --database to re-annotate from a "
                              f"database's own sidecars")
        if not Path(recorded).is_file():
            raise MashIDError(f"the run used --db-{key} {recorded}, which no longer exists; "
                              f"restore it, or pass --database to re-annotate from a "
                              f"database's own sidecars")
        overrides[key] = Path(recorded)
        print(f"note: using the {key} table the run was made with, {recorded}", file=sys.stderr)

    metadata = read_metadata(overrides["metadata"] or sidecar_path(database))
    confusable = confusable_mod.load(database, overrides["confusable"])
    if confusable is None:
        print(f"note: no {confusable_mod.confusable_path(database).name} beside the database, so no "
              f"'not separable' notes; write one with make_mashID_db --confusable", file=sys.stderr)

    # Sequences/Bases and any sample-sheet columns live only in the summary; carry them over.
    with summary_path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        found = list(reader.fieldnames or SUMMARY_COLUMNS)
        if "Sample" not in found:
            raise MashIDError(f"{summary_path} has no Sample column; found {found}")
        rows_read = list(reader)
        existing = {row["Sample"]: row for row in rows_read}
    if len(existing) != len(rows_read):
        # Two rows for one sample: the second would win silently and the first disappear.
        counts = Counter(row["Sample"] for row in rows_read)
        repeated = sorted(name for name, n in counts.items() if n > 1)
        raise MashIDError(f"{summary_path.name} has more than one row for "
                          f"{', '.join(repr(r) for r in repeated[:5])}; rewriting it would keep only "
                          f"one of each. Refusing to write.")
    extra_columns = [c for c in found if c not in SUMMARY_COLUMNS]
    # The column set is mashID's own plus whatever the sample sheet added, exactly as run() builds
    # it -- not the old summary's columns. A summary written before a column existed would otherwise
    # have it dropped again here, and `Note` is the column this script most often has to add.
    columns = list(SUMMARY_COLUMNS) + extra_columns

    # The sample list comes from the run record, but the summary being replaced is the authority on
    # which rows must survive: a run.json with an empty or renamed "samples" key would otherwise
    # truncate the summary to its header and report success.
    entries = run_info.get("samples") or []
    named = {e.get("name") for e in entries if isinstance(e, dict)}
    lost = sorted(set(existing) - {n for n in named if n})
    if lost:
        raise MashIDError(
            f"{run_info_path.name} lists {len(named)} sample(s) but {summary_path.name} has "
            f"{len(existing)} row(s); rewriting it would drop {len(lost)} of them "
            f"({', '.join(repr(n) for n in lost[:5])}{', ...' if len(lost) > 5 else ''}). "
            f"Refusing to write.")

    results, changed = [], []
    for entry in entries:
        name = entry.get("name") if isinstance(entry, dict) else None
        if not name:
            raise MashIDError(f"{run_info_path.name}: a sample entry has no 'name'")
        per_sample = results_dir / f"{name}{PER_SAMPLE_SUFFIX}"
        if not per_sample.is_file():
            raise MashIDError(f"Per-sample table missing for {name}: {per_sample}")
        old = existing.get(name)
        if old is None:
            print(f"warning: {name} is in {run_info_path.name} but not in {summary_path.name}; its "
                  f"Sequences, Bases and any sample-sheet columns will be NA", file=sys.stderr)
            old = {}
        sample = Sample(name=name,
                        files=[Path(f) for f in entry.get("files", [])],
                        seq_type=entry.get("type", "fastq"),
                        sequences=_int_or_none(old.get("Sequences")),
                        bases=_int_or_none(old.get("Bases")),
                        extra={c: old.get(c, "") for c in extra_columns})
        result = SampleResult(sample=sample, hits=read_hits(per_sample))
        result.notes = annotate(result, metadata,
                                params.get("ambiguity_margin", _default("ambiguity_margin")),
                                params.get("winner_take_all", _default("winner_take_all")),
                                confusable)
        results.append(result)

    rows = [summary_row(r, metadata) for r in results]
    for row in rows:
        # A column the old summary never had counts as changed: absent is not the same as agreeing,
        # and treating it as agreement is what made a run predating --confusable look up to date.
        old = existing.get(row["Sample"], {})
        if any(row[c] != (old.get(c) or "") for c in LOOKED_UP_COLUMNS):
            changed.append(row["Sample"])
    return Rebuilt(results, rows, columns, database, changed,
                   overrides["metadata"] or sidecar_path(database),
                   overrides["confusable"] or (confusable_mod.confusable_path(database)
                                               if confusable is not None else None))


def provenance(rebuilt: Rebuilt) -> dict:
    record: dict = {
        "at": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat(),
        "reason": "Identification, TaxID and Note looked up again from the sidecars below, using "
                  "the hits already recorded for each sample",
        "screening_unchanged": True,
        "samples_changed": list(rebuilt.changed),
    }
    for key, path in (("metadata_sidecar", rebuilt.metadata_table),
                      ("confusable_sidecar", rebuilt.confusable_table)):
        record[key] = ({"path": str(path), "md5": md5sum(path)}
                       if path is not None and path.is_file() else None)
    return record


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("results", metavar="RESULTS_DIR", type=Path,
                    help="A mashID output directory, holding mashID_run.json and the per-sample tables.")
    ap.add_argument("--database", metavar="DB.msh", type=Path, default=None,
                    help="Read the names from this database's own sidecars, for a database that has "
                         "moved. Default: the database the run recorded, and the --db-metadata / "
                         "--db-confusable tables it was made with if it used any.")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report what would change and write nothing.")
    args = ap.parse_args(argv)

    try:
        rebuilt = rebuild(args.results, args.database)
    except MashIDError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    rows, changed = rebuilt.rows, rebuilt.changed

    print(f"{len(rows)} sample(s); {len(changed)} with a changed "
          f"{', '.join(LOOKED_UP_COLUMNS)}; {sum(1 for r in rows if r['Note'])} with any note")
    if args.dry_run:
        shown = [r for r in rows if r["Sample"] in changed][:DRY_RUN_SAMPLES]
        for row in shown:
            print(f"  {row['Sample']}: {row['Identification']} ({row['TaxID']}) "
                  f"{row['Note'] or '(no note)'}")
        if len(changed) > len(shown):
            print(f"  ... and {len(changed) - len(shown)} more")
        print("dry run: nothing written")
        return 0

    # Nothing to do is not worth a rewrite: it would stack up another set of backups and stamp the
    # JSON with a re-annotation that changed nothing. The column set counts as something to do,
    # since gaining a Note column is the whole point on a run that predates --confusable, and so
    # does a missing companion file -- "already up to date" must not be said about a summary whose
    # _mqc.tsv or .json is absent.
    header = _read_columns(args.results / SUMMARY_FILENAME)
    missing = [f for f in REWRITTEN if not (args.results / f).is_file()]
    if not changed and header == rebuilt.columns and not missing:
        print(f"already up to date with {rebuilt.metadata_table.name}"
              f"{' and ' + rebuilt.confusable_table.name if rebuilt.confusable_table else ''}; "
              f"nothing written")
        return 0

    for filename in REWRITTEN:
        target = args.results / filename
        if target.is_file() and not args.no_backup:
            print(f"backup: {backup(target).name}")

    write_tsv(args.results / SUMMARY_FILENAME, rebuilt.columns, rows)
    write_multiqc(args.results / MULTIQC_FILENAME, rows)
    json_path = args.results / JSON_FILENAME
    write_json(json_path, rebuilt.results, rows, read_metadata(rebuilt.metadata_table),
               rebuilt.database)
    doc = json.loads(json_path.read_text())
    doc["reannotated"] = provenance(rebuilt)
    json_path.write_text(json.dumps(doc, indent=2) + "\n")
    print(f"rewrote {', '.join(REWRITTEN)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
