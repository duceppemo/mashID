"""Command-line interface for mashID."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from mashid import MashIDError, __version__
from mashid.databases import DEFAULT_DB_NAME, REGISTRY
from mashid.pipeline import DEFAULT_MIN_REF_LENGTH, Settings, run

log = logging.getLogger("mashid")


def _fraction(value: str) -> float:
    try:
        f = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not a number: {value!r}") from None
    if not 0.0 <= f <= 1.0:
        raise argparse.ArgumentTypeError(f"must be between 0 and 1, got {value}")
    return f


def _positive_int(value: str) -> int:
    try:
        i = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"not an integer: {value!r}") from None
    if i < 1:
        raise argparse.ArgumentTypeError(f"must be at least 1, got {value}")
    return i


def build_parser() -> argparse.ArgumentParser:
    max_cpu = os.cpu_count() or 1
    parser = argparse.ArgumentParser(
        prog="mashID",
        description="Species identification from genome assemblies (fasta) or raw reads (fastq) using Mash.",
    )
    io = parser.add_argument_group("input/output")
    io.add_argument("-i", "--input", metavar="PATH", type=Path, default=None,
                    help="Input directory (searched recursively) with fastq/fasta files, or a single "
                         "fastq/fasta file, gzipped or not. Paired-end files (R1/R2) are screened together.")
    io.add_argument("--sample-sheet", metavar="FILE.tsv", type=Path, default=None,
                    help="Instead of -i: a TSV/CSV with columns 'sample' and 'file' (one row per file, or "
                         "files separated by ';'). Relative paths are resolved from the sheet's directory. "
                         "Other columns are copied into the summary.")
    io.add_argument("-o", "--output", metavar="DIR", type=Path, required=True,
                    help="Output directory (created if needed).")
    io.add_argument("-d", "--database", metavar="FILE.msh|NAME", default=None,
                    help=f"Mash sketch database: a .msh file, or the name of a downloaded pre-built database "
                         f"({', '.join(REGISTRY)}; see mashID_download_db). Default: {DEFAULT_DB_NAME}")
    io.add_argument("--db-metadata", metavar="FILE.tsv", type=Path, default=None,
                    help="Table mapping reference accessions to organism names (and TaxIDs). Default: the "
                         "<database>.metadata.tsv sidecar written by make_mashID_db, if present.")

    flt = parser.add_argument_group("screening")
    flt.add_argument("--identity", metavar="0.9", type=_fraction, default=0.9,
                     help="Minimum identity to report, between 0 and 1. Default: %(default)s")
    flt.add_argument("--p-value", metavar="0.05", type=_fraction, default=0.05,
                     help="Maximum p-value to report. Default: %(default)s")
    flt.add_argument("-n", "--n-hits", metavar="10", type=_positive_int, default=10,
                     help="Number of top hits to report per sample. Default: %(default)s")
    flt.add_argument("-s", "--sort-by", choices=["identity", "multiplicity"], default="identity",
                     help='How to rank hits; determines the "top hit" in the summary. Default: %(default)s')
    flt.add_argument("--no-winner-take-all", action="store_true",
                     help="Disable Mash's winner-take-all strategy (-w). Reports more redundant hits.")
    flt.add_argument("--max-reads", metavar="N", type=_positive_int, default=None,
                     help="Screen only the first N reads of each fastq sample (split across R1/R2), streamed "
                         "to Mash. Much faster on large runs; Sequences/Bases then describe the screened reads.")
    flt.add_argument("--ambiguity-margin", metavar="0.005", type=_fraction, default=0.005,
                     help="Flag a sample as ambiguous when a different organism scores within this identity "
                          "margin of the top hit. Default: %(default)s")
    flt.add_argument("--min-ref-length", metavar="BP", type=int, default=DEFAULT_MIN_REF_LENGTH,
                     help="Ignore hits to references shorter than this many bp (partial records in a database "
                          "otherwise produce false top hits). 0 keeps all hits. Default: %(default)s")
    flt.add_argument("--skip-stats", action="store_true",
                     help="Do not count reads/bases of the input files (faster on very large fastq).")
    flt.add_argument("--fail-on", choices=["none", "no-hit", "note"], default="none",
                     help="Exit with code 2 when any sample has no hit ('no-hit') or has any note or no hit "
                          "('note'). For pipelines. Default: %(default)s")

    perf = parser.add_argument_group("performance")
    perf.add_argument("-t", "--threads", metavar=str(max_cpu), type=_positive_int, default=max_cpu,
                      help="Total number of threads, shared between parallel samples. Default: all (%(default)s)")
    perf.add_argument("-p", "--parallel", metavar="2", type=_positive_int, default=2,
                      help="Number of samples to process concurrently. Default: %(default)s")
    perf.add_argument("-m", "--memory", metavar="GB", help=argparse.SUPPRESS)  # deprecated, ignored

    parser.add_argument("--debug", action="store_true", help="Verbose logging.")
    parser.add_argument("-v", "--version", action="version", version=f"mashID {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.input is None and args.sample_sheet is None:
        parser.error("one of -i/--input or --sample-sheet is required")
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    max_cpu = os.cpu_count() or 1
    threads = args.threads
    if threads > max_cpu:
        log.warning("Requested %d threads but only %d CPU(s) available; using %d", threads, max_cpu, max_cpu)
        threads = max_cpu
    if args.parallel > threads:
        log.warning("--parallel (%d) exceeds --threads (%d); each sample will run with 1 thread",
                    args.parallel, threads)
    if args.memory is not None:
        log.warning("--memory is deprecated and ignored (mashID no longer needs BBMap)")

    settings = Settings(
        input=args.input,
        sample_sheet=args.sample_sheet,
        output=args.output,
        database=args.database,
        db_metadata=args.db_metadata,
        min_identity=args.identity,
        max_p_value=args.p_value,
        n_hits=args.n_hits,
        sort_by=args.sort_by,
        threads=threads,
        parallel=args.parallel,
        winner_take_all=not args.no_winner_take_all,
        skip_stats=args.skip_stats,
        max_reads=args.max_reads,
        ambiguity_margin=args.ambiguity_margin,
        min_ref_length=max(0, args.min_ref_length),
        fail_on=args.fail_on,
        command_line=list(sys.argv if argv is None else ["mashID", *argv]),
    )
    try:
        return run(settings)
    except MashIDError as exc:
        log.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        log.error("Interrupted")
        return 130


if __name__ == "__main__":
    sys.exit(main())
