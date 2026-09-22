"""Build a Mash sketch database (plus metadata sidecar) suitable for mashID from a folder of fasta files."""

from __future__ import annotations

import argparse
import logging
import os
import sys
import tempfile
from pathlib import Path

from mashid import MashIDError, __version__
from mashid.mash import find_mash, info_table, mash_version, sketch, sketch_count
from mashid.metadata import (
    entries_from_sketches,
    read_metadata,
    read_ncbi_assembly_report,
    sidecar_path,
    write_metadata,
)
from mashid.samples import FASTA_EXTENSIONS

log = logging.getLogger("mashid.makedb")

DB_EXTENSIONS = tuple(ext + gz for ext in FASTA_EXTENSIONS for gz in ("", ".gz"))
DEFAULT_MIN_LENGTH = 100_000  # bp; shorter "genomes" are usually partial records that cause false top hits


def list_fasta_files(input_path: Path) -> list[Path]:
    """Fasta files under a directory (recursive), or the paths listed in a text file (one per line)."""
    if input_path.is_file():
        paths = [Path(line.strip()) for line in input_path.read_text().splitlines() if line.strip()]
        missing = [p for p in paths if not p.is_file()]
        if missing:
            raise MashIDError("Listed file(s) not found: " + ", ".join(str(p) for p in missing[:5]))
        return paths
    if not input_path.is_dir():
        raise MashIDError(f"Input path does not exist: {input_path}")
    files: list[Path] = []
    visited: set[Path] = set()
    for root, dirs, filenames in os.walk(input_path, followlinks=True):
        visited.add(Path(root).resolve())
        dirs[:] = sorted(d for d in dirs if not d.startswith(".") and (Path(root) / d).resolve() not in visited)
        for name in sorted(filenames):
            if not name.startswith(".") and name.lower().endswith(DB_EXTENSIONS):
                files.append(Path(root) / name)
    if not files:
        raise MashIDError(
            f"No fasta files found in {input_path}. Accepted extensions: {', '.join(DB_EXTENSIONS)}"
        )
    return files


def load_annotations(metadata: Path | None, assembly_report: Path | None) -> dict[str, tuple[str, str]]:
    """Merge user metadata table and NCBI assembly report into {accession: (organism, taxid)}."""
    annotations: dict[str, tuple[str, str]] = {}
    if assembly_report is not None:
        annotations.update(read_ncbi_assembly_report(assembly_report))
        log.info("Loaded %d organism(s) from NCBI assembly report %s", len(annotations), assembly_report)
    if metadata is not None:
        table = read_metadata(metadata)
        annotations.update({acc: (e.organism, e.taxid) for acc, e in table.items()})
        log.info("Loaded %d organism(s) from %s", len(table), metadata)
    return annotations


def write_sidecar(database: Path, annotations: dict[str, tuple[str, str]], exe: str) -> Path:
    """Write <db>.metadata.tsv from the sketches in the database plus optional annotations."""
    entries, unannotated = entries_from_sketches(info_table(database, exe), annotations)
    path = sidecar_path(database)
    write_metadata(path, entries)
    if annotations and unannotated:
        log.warning("%d of %d reference(s) had no entry in the provided metadata; their names were parsed "
                    "from headers", unannotated, len(entries))
    log.info("Metadata sidecar written: %s (%d references)", path, len(entries))
    return path


def _sketch_files(files: list[Path], output_prefix: Path, threads: int, sketch_size: int, kmer_size: int,
                  exe: str) -> Path:
    fd, list_name = tempfile.mkstemp(prefix="mashID_files_", suffix=".txt")
    list_file = Path(list_name)
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write("\n".join(str(f.resolve()) for f in files) + "\n")
        return sketch(list_file, output_prefix, threads=threads, sketch_size=sketch_size, kmer_size=kmer_size,
                      exe=exe)
    finally:
        list_file.unlink(missing_ok=True)


def make_database(input_path: Path, output_dir: Path, prefix: str, threads: int,
                  sketch_size: int, kmer_size: int,
                  metadata: Path | None = None, assembly_report: Path | None = None,
                  min_length: int = DEFAULT_MIN_LENGTH) -> Path:
    exe = find_mash()
    log.info("Using %s (version %s)", exe, mash_version(exe))
    if not 1 <= kmer_size <= 32:
        raise MashIDError("k-mer size must be between 1 and 32 (inclusive)")
    if sketch_size < 1:
        raise MashIDError("Sketch size must be at least 1")
    if os.sep in prefix or prefix in ("", ".", ".."):
        raise MashIDError(f"Invalid database prefix: {prefix!r}")
    if prefix.endswith(".msh"):
        prefix = prefix[: -len(".msh")]
    annotations = load_annotations(metadata, assembly_report)  # fail early on bad tables

    files = list_fasta_files(input_path)
    log.info("Sketching %d fasta file(s) (k=%d, s=%d) with %d thread(s)...",
             len(files), kmer_size, sketch_size, threads)
    output_dir.mkdir(parents=True, exist_ok=True)
    db = _sketch_files(files, output_dir / prefix, threads, sketch_size, kmer_size, exe)

    # Drop references shorter than --min-length: Mash reports their length, so filter after sketching
    # and re-sketch only when something was removed (the rare case).
    if min_length > 0:
        sketches = info_table(db, exe)
        short = [sk for sk in sketches if sk.length < min_length]
        if short:
            for sk in short:
                log.warning("Excluding %s: only %d bp (< --min-length %d)", Path(sk.query_id).name, sk.length,
                            min_length)
            keep = {sk.query_id for sk in sketches} - {sk.query_id for sk in short}
            kept_files = [f for f in files if str(f.resolve()) in keep]
            if not kept_files:
                db.unlink(missing_ok=True)
                raise MashIDError(f"All references are shorter than --min-length {min_length}; nothing to sketch")
            log.info("Re-sketching %d reference(s) without the %d short one(s)...", len(kept_files), len(short))
            db = _sketch_files(kept_files, output_dir / prefix, threads, sketch_size, kmer_size, exe)

    n = sketch_count(db, exe)
    log.info("Database written: %s (%s sketches)", db, n if n is not None else "unknown number of")
    write_sidecar(db, annotations, exe)
    return db


def annotate_database(database: Path, metadata: Path | None, assembly_report: Path | None) -> Path:
    """Write the metadata sidecar for an existing database without re-sketching."""
    exe = find_mash()
    if not database.is_file():
        raise MashIDError(f"Mash database not found: {database}")
    return write_sidecar(database, load_annotations(metadata, assembly_report), exe)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="make_mashID_db",
        description="Create a Mash sketch database for mashID from genome assemblies, together with a "
                    "<name>.metadata.tsv sidecar mapping each reference to an organism name and TaxID. "
                    "Names come from --metadata / --assembly-report when given, otherwise from the first "
                    "fasta header of each file (NCBI-style headers work).",
    )
    parser.add_argument("-i", "--input", metavar="PATH", type=Path,
                        help="Directory with fasta files (searched recursively; .fna/.fa/.fasta, gzipped or "
                             "not), or a text file listing fasta paths one per line.")
    parser.add_argument("-o", "--output", metavar="DIR", type=Path,
                        help="Output directory for the database.")
    parser.add_argument("-p", "--prefix", metavar="NAME", default="mashID_db",
                        help='Database file name without extension (".msh" is appended). Default: %(default)s')
    parser.add_argument("--metadata", metavar="FILE", type=Path, default=None,
                        help="TSV/CSV with columns Accession, Organism and optionally TaxID, used to name "
                             "references instead of parsing headers.")
    parser.add_argument("--assembly-report", metavar="FILE.jsonl", type=Path, default=None,
                        help="assembly_data_report.jsonl from an NCBI 'datasets download genome' archive; "
                             "provides organism names and TaxIDs for every accession.")
    parser.add_argument("--annotate", metavar="DB.msh", type=Path, default=None,
                        help="Only (re)write the metadata sidecar for an existing database; no sketching. "
                             "-i and -o are then not needed.")
    parser.add_argument("-t", "--threads", metavar="4", type=int, default=4,
                        help="Number of threads. Default: %(default)s")
    parser.add_argument("-s", "--sketch-size", metavar="10000", type=int, default=10000,
                        help="Sketch size (min-hashes per genome). Larger is more sensitive but slower; "
                             "10000 is recommended for species identification. Default: %(default)s")
    parser.add_argument("-k", "--kmer-size", metavar="21", type=int, default=21,
                        help="K-mer size (1-32). Default: %(default)s")
    parser.add_argument("--min-length", metavar="BP", type=int, default=DEFAULT_MIN_LENGTH,
                        help="Exclude references shorter than this many bp (partial records cause false top "
                             "hits). Use 0 to keep everything, e.g. for plasmid or viral databases. "
                             "Default: %(default)s")
    parser.add_argument("--debug", action="store_true", help="Verbose logging.")
    parser.add_argument("-v", "--version", action="version", version=f"make_mashID_db {__version__}")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.debug else logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
    if args.annotate is None and (args.input is None or args.output is None):
        parser.error("-i/--input and -o/--output are required (unless using --annotate)")
    threads = min(max(1, args.threads), os.cpu_count() or 1)
    try:
        if args.annotate is not None:
            annotate_database(args.annotate, args.metadata, args.assembly_report)
        else:
            make_database(args.input, args.output, args.prefix, threads, args.sketch_size, args.kmer_size,
                          metadata=args.metadata, assembly_report=args.assembly_report,
                          min_length=max(0, args.min_length))
    except MashIDError as exc:
        log.error("%s", exc)
        return 1
    except KeyboardInterrupt:
        log.error("Interrupted")
        return 130
    return 0


if __name__ == "__main__":
    sys.exit(main())
