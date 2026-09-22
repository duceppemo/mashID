"""The mashID identification pipeline: discover samples, count reads, screen, annotate, report."""

from __future__ import annotations

import csv
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

from mashid.databases import resolve_database
from mashid.mash import ScreenHit, find_mash, info_table, mash_version, screen
from mashid.metadata import (
    NA,
    DbEntry,
    accession_from_query_id,
    entries_from_sketches,
    lookup,
    read_metadata,
    sidecar_path,
    write_metadata,
)
from mashid.samples import Sample, discover_samples
from mashid.seqstats import FastqSubsampler, sample_stats
from mashid.taxonomy import clean_comment, organism_from_comment

log = logging.getLogger(__name__)

NO_HIT = "No significant hit in database"
SUMMARY_FILENAME = "summary_mashID.tsv"

# Annotation heuristics (see annotate()).
# Under winner-take-all a reference keeps only hashes not credited to a better hit. A complete genome
# with >= 95% ANI to the winner (k=21) retains at most ~20% of its sketch, i.e. identity < 0.99, so a
# different organism at >= 0.99 must be supported by its own k-mers.
MIXTURE_MIN_IDENTITY = 0.99
LOW_COVERAGE_MULTIPLICITY = 5  # median multiplicity of the top hit below this suggests low coverage
DEFAULT_MIN_REF_LENGTH = 100_000  # bp; hits to shorter references are ignored (partial records)

SAMPLE_COLUMNS = [
    "Rank", "Identity", "Shared_Hashes", "Median_Multiplicity", "P_Value",
    "Accession", "TaxID", "Identification", "Description",
]
SUMMARY_COLUMNS = [
    "Sample", "Sequences", "Bases", "Identity", "Shared_Hashes", "Median_Multiplicity",
    "P_Value", "Accession", "TaxID", "Identification", "Note",
]


@dataclass
class Settings:
    input: Path
    output: Path
    database: str | Path | None = None  # path, registry name, or None for the default database
    db_metadata: Path | None = None     # override for the <db>.metadata.tsv sidecar
    min_identity: float = 0.9
    max_p_value: float = 0.05
    n_hits: int = 10
    sort_by: str = "identity"  # or "multiplicity"
    threads: int = 1
    parallel: int = 1
    winner_take_all: bool = True
    skip_stats: bool = False
    max_reads: int | None = None  # screen only the first N reads of fastq samples
    ambiguity_margin: float = 0.005
    min_ref_length: int = DEFAULT_MIN_REF_LENGTH  # ignore hits to references shorter than this (0 = keep all)


@dataclass
class SampleResult:
    sample: Sample
    hits: list[ScreenHit]  # sorted, truncated to n_hits
    notes: list[str] = field(default_factory=list)
    n_short_refs: int = 0  # hits dropped because the reference was shorter than min_ref_length

    @property
    def top(self) -> ScreenHit | None:
        return self.hits[0] if self.hits else None


def sort_hits(hits: list[ScreenHit], sort_by: str) -> list[ScreenHit]:
    if sort_by == "multiplicity":
        key = lambda h: (-h.multiplicity_f, -h.identity_f, -h.shared_f, h.p_value_f)  # noqa: E731
    elif sort_by == "identity":
        key = lambda h: (-h.identity_f, -h.shared_f, h.p_value_f, -h.multiplicity_f)  # noqa: E731
    else:
        raise ValueError(f"Unknown sort key: {sort_by}")
    return sorted(hits, key=key)


def drop_short_references(hits: list[ScreenHit], metadata: dict[str, DbEntry] | None,
                          min_length: int) -> tuple[list[ScreenHit], int]:
    """Remove hits whose reference is known (from metadata) to be shorter than ``min_length`` bp."""
    if not metadata or min_length <= 0:
        return hits, 0
    kept: list[ScreenHit] = []
    for hit in hits:
        entry = lookup(metadata, accession_from_query_id(hit.query_id))
        if entry is not None and entry.length is not None and entry.length < min_length:
            log.debug("Ignoring hit to short reference %s (%d bp)", entry.accession, entry.length)
            continue
        kept.append(hit)
    return kept, len(hits) - len(kept)


def identify(hit: ScreenHit, metadata: dict[str, DbEntry] | None) -> tuple[str, str, str]:
    """(accession, organism, taxid) for a hit: from the sidecar when available, else the header."""
    accession = accession_from_query_id(hit.query_id)
    entry = lookup(metadata, accession)
    if entry is not None and entry.organism and entry.organism != NA:
        return accession, entry.organism, entry.taxid or NA
    return accession, organism_from_comment(hit.comment, fallback=accession), NA


def hit_row(rank: int, hit: ScreenHit, metadata: dict[str, DbEntry] | None) -> dict[str, str]:
    accession, organism, taxid = identify(hit, metadata)
    return {
        "Rank": str(rank),
        "Identity": hit.identity,
        "Shared_Hashes": hit.shared_hashes,
        "Median_Multiplicity": hit.median_multiplicity,
        "P_Value": hit.p_value,
        "Accession": accession,
        "TaxID": taxid,
        "Identification": organism,
        "Description": clean_comment(hit.comment),
    }


def _ref_length(hit: ScreenHit, metadata: dict[str, DbEntry] | None) -> int | None:
    entry = lookup(metadata, accession_from_query_id(hit.query_id))
    return entry.length if entry is not None else None


def annotate(result: SampleResult, metadata: dict[str, DbEntry] | None, ambiguity_margin: float,
             winner_take_all: bool = True) -> list[str]:
    """Heuristic warnings for the summary: possible mixture, ambiguous call, low coverage.

    With winner-take-all (the default), hashes shared between references are credited to the best one,
    so a *different* organism that still scores >= MIXTURE_MIN_IDENTITY is supported by its own k-mers and
    is most likely present: "Possible mixture". A different organism below that threshold but within
    ``ambiguity_margin`` of the top hit is reported as "Ambiguous". Without winner-take-all, sister taxa
    all score high, so only the ambiguity check is meaningful.
    """
    notes: list[str] = []
    top = result.top
    if top is None:
        return notes
    _, top_org, _ = identify(top, metadata)
    others: list[tuple[str, ScreenHit]] = []
    for hit in result.hits[1:]:
        _, org, _ = identify(hit, metadata)
        if org != top_org and org not in {o for o, _ in others}:
            others.append((org, hit))  # best hit of each other organism, in rank order

    mixture: list[str] = []
    if winner_take_all:
        top_len = _ref_length(top, metadata)
        labels: list[str] = []
        for org, hit in others:
            if hit.identity_f < MIXTURE_MIN_IDENTITY:
                continue
            mixture.append(org)
            other_len = _ref_length(hit, metadata)
            # A much shorter reference may be a partial assembly whose hashes escape winner-take-all.
            if top_len and other_len and other_len < 0.5 * top_len:
                org = f"{org} (reference only {other_len} bp, may be partial)"
            labels.append(org)
        if labels:
            notes.append("Possible mixture with: " + ", ".join(sorted(labels)))

    if others:
        org, hit = others[0]
        if org not in mixture and top.identity_f - hit.identity_f <= ambiguity_margin:
            notes.append(f"Ambiguous: {org} at identity {hit.identity}")

    if result.sample.seq_type == "fastq" and top.multiplicity_f < LOW_COVERAGE_MULTIPLICITY:
        notes.append(f"Low coverage: median multiplicity {top.median_multiplicity}")
    return notes


def summary_row(result: SampleResult, metadata: dict[str, DbEntry] | None) -> dict[str, str]:
    s = result.sample
    row = {
        "Sample": s.name,
        "Sequences": NA if s.sequences is None else str(s.sequences),
        "Bases": NA if s.bases is None else str(s.bases),
    }
    top = result.top
    if top is None:
        row.update({c: NA for c in SUMMARY_COLUMNS[3:]})
        row["Identification"] = NO_HIT
        row["Note"] = ""
    else:
        h = hit_row(1, top, metadata)
        row.update({c: h[c] for c in SUMMARY_COLUMNS[3:-1]})
        row["Note"] = "; ".join(result.notes)
    return row


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def format_table(columns: list[str], rows: list[dict[str, str]]) -> str:
    widths = {c: max(len(c), *(len(r.get(c, "")) for r in rows)) if rows else len(c) for c in columns}
    lines = ["  ".join(c.ljust(widths[c]) for c in columns)]
    for r in rows:
        lines.append("  ".join(r.get(c, "").ljust(widths[c]) for c in columns))
    return "\n".join(line.rstrip() for line in lines)


def load_metadata(database: Path, override: Path | None, exe: str) -> dict[str, DbEntry]:
    """Reference metadata: ``--db-metadata``, else the sidecar, else built from ``mash info`` on the fly."""
    if override is not None:
        entries = read_metadata(override)
        log.info("Loaded metadata for %d reference(s) from %s", len(entries), override)
        return entries
    path = sidecar_path(database)
    if path.is_file():
        entries = read_metadata(path)
        log.info("Loaded metadata for %d reference(s) from %s", len(entries), path)
        return entries
    log.info("No metadata sidecar found; reading reference names and lengths from the database")
    built, _ = entries_from_sketches(info_table(database, exe))
    entries = {e.accession: e for e in built}
    try:
        write_metadata(path, built)
        log.info("Metadata sidecar written for next time: %s", path)
    except OSError as exc:
        log.debug("Could not write %s (%s); continuing with in-memory metadata", path, exc)
    return entries


def _screen_sample(sample: Sample, settings: Settings, database: Path, metadata: dict[str, DbEntry] | None,
                   threads: int, exe: str) -> SampleResult:
    common = dict(threads=threads, min_identity=settings.min_identity, max_p_value=settings.max_p_value,
                  winner_take_all=settings.winner_take_all, exe=exe)
    if settings.max_reads and sample.seq_type == "fastq":
        subsampler = FastqSubsampler(sample.files, settings.max_reads)
        hits = screen(database, None, stdin_chunks=subsampler.chunks(), **common)
        sample.sequences, sample.bases = subsampler.sequences, subsampler.bases
    else:
        hits = screen(database, sample.files, **common)
    hits, n_short = drop_short_references(hits, metadata, settings.min_ref_length)
    hits = sort_hits(hits, settings.sort_by)[: settings.n_hits]
    return SampleResult(sample=sample, hits=hits, n_short_refs=n_short)


def compute_stats(samples: list[Sample], parallel: int) -> None:
    """Fill in sequences/bases for each sample, using one process per sample."""
    if not samples:
        return
    with ProcessPoolExecutor(max_workers=max(1, parallel)) as pool:
        futures = {pool.submit(sample_stats, s.files, s.seq_type): s for s in samples}
        for fut in as_completed(futures):
            sample = futures[fut]
            try:
                sample.sequences, sample.bases = fut.result()
            except Exception as exc:  # keep going: stats are informational
                log.warning("Could not compute stats for %s: %s", sample.name, exc)


def run(settings: Settings) -> list[SampleResult]:
    exe = find_mash()
    log.info("Using %s (version %s)", exe, mash_version(exe))
    database = resolve_database(settings.database)
    log.info("Database: %s", database)
    metadata = load_metadata(database, settings.db_metadata, exe)

    samples = discover_samples(settings.input)
    log.info("Found %d sample(s) in %s", len(samples), settings.input)
    for s in samples:
        log.debug("  %s (%s): %s", s.name, s.seq_type, ", ".join(f.name for f in s.files))

    settings.output.mkdir(parents=True, exist_ok=True)

    if not settings.skip_stats:
        # With --max-reads, fastq stats come from the reads actually screened.
        to_count = [s for s in samples if not (settings.max_reads and s.seq_type == "fastq")]
        if to_count:
            log.info("Counting sequences and bases in %d sample(s)...", len(to_count))
            compute_stats(to_count, settings.parallel)

    parallel = max(1, min(settings.parallel, len(samples)))
    threads_per_job = max(1, settings.threads // parallel)
    if settings.max_reads:
        log.info("Screening only the first %d reads of each fastq sample", settings.max_reads)
    log.info("Screening %d sample(s), %d at a time with %d thread(s) each...",
             len(samples), parallel, threads_per_job)

    results: dict[str, SampleResult] = {}
    with ThreadPoolExecutor(max_workers=parallel) as pool:
        futures = {pool.submit(_screen_sample, s, settings, database, metadata, threads_per_job, exe): s
                   for s in samples}
        for fut in as_completed(futures):
            result = fut.result()  # raises MashIDError on mash failure
            result.notes = annotate(result, metadata, settings.ambiguity_margin, settings.winner_take_all)
            results[result.sample.name] = result
            top = result.top
            summary = f"{identify(top, metadata)[1]} ({top.identity})" if top else NO_HIT
            if result.notes:
                summary += " [" + "; ".join(result.notes) + "]"
            if result.n_short_refs:
                summary += f" ({result.n_short_refs} hit(s) to references < {settings.min_ref_length} bp ignored)"
            log.info("  %s: %s", result.sample.name, summary)
            rows = [hit_row(i, h, metadata) for i, h in enumerate(result.hits, start=1)]
            write_tsv(settings.output / f"{result.sample.name}_mashID.tsv", SAMPLE_COLUMNS, rows)

    ordered = [results[s.name] for s in samples]
    summary_rows = [summary_row(r, metadata) for r in ordered]
    write_tsv(settings.output / SUMMARY_FILENAME, SUMMARY_COLUMNS, summary_rows)

    print("\nIdentification results:\n", file=sys.stdout)
    print(format_table(SUMMARY_COLUMNS, summary_rows))
    print(f"\nSummary written to {settings.output / SUMMARY_FILENAME}")
    return ordered
