"""Discover input sequence files and group them into samples."""

from __future__ import annotations

import csv
import logging
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from mashid import MashIDError

log = logging.getLogger(__name__)

FASTQ_EXTENSIONS = (".fastq", ".fq")
FASTA_EXTENSIONS = (".fasta", ".fa", ".fna")
SEQ_EXTENSIONS = tuple(
    ext + gz for ext in FASTQ_EXTENSIONS + FASTA_EXTENSIONS for gz in ("", ".gz")
)

# Illumina-style read designations stripped from file stems to obtain the sample name:
#   sample_S10_L001_R1_001, sample_R1, sample_1, sample_L002_2 ...
_READ_SUFFIX = re.compile(r"(?:_S\d+)?(?:_L\d{3})?(?:_R?[12])(?:_001)?$")


@dataclass
class Sample:
    name: str
    files: list[Path]
    seq_type: str  # "fastq" or "fasta"
    sequences: int | None = None  # number of reads (fastq) or contigs (fasta)
    bases: int | None = None
    warnings: list[str] = field(default_factory=list)
    extra: dict[str, str] = field(default_factory=dict)  # extra columns from a sample sheet


def strip_seq_extension(filename: str) -> str | None:
    """Return the file name without its sequence extension, or None if it is not a sequence file."""
    lower = filename.lower()
    for ext in sorted(SEQ_EXTENSIONS, key=len, reverse=True):
        if lower.endswith(ext):
            return filename[: -len(ext)]
    return None


def is_sequence_file(path: Path | str) -> bool:
    return strip_seq_extension(Path(path).name) is not None


def seq_type_of(path: Path | str) -> str:
    lower = Path(path).name.lower()
    if lower.endswith(".gz"):
        lower = lower[:-3]
    if lower.endswith(FASTQ_EXTENSIONS):
        return "fastq"
    if lower.endswith(FASTA_EXTENSIONS):
        return "fasta"
    raise ValueError(f"Not a recognised sequence file: {path}")


def sample_name_from_file(path: Path | str) -> str:
    """Derive a sample name from a file name.

    "MBWGS440_S10_L001_R1_001.fastq.gz" -> "MBWGS440"
    "isolate_A_R2.fq.gz"                -> "isolate_A"
    "GCF_000195955.2.fna"               -> "GCF_000195955.2"
    """
    stem = strip_seq_extension(Path(path).name)
    if stem is None:
        raise ValueError(f"Not a recognised sequence file: {path}")
    name = _READ_SUFFIX.sub("", stem)
    return name or stem


def find_sequence_files(input_path: Path) -> list[Path]:
    """Return sequence files found at ``input_path`` (a file, or a directory searched recursively).

    Paths are absolute but symbolic links are not resolved, so a link's name names the sample."""
    input_path = Path(input_path)
    if input_path.is_symlink() and not input_path.exists():
        raise MashIDError(f"Input is a link to a missing file: {input_path} -> {os.readlink(input_path)}")
    if input_path.is_file():
        if not is_sequence_file(input_path):
            raise MashIDError(
                f"Input file does not have a recognised extension ({', '.join(SEQ_EXTENSIONS)}): {input_path}"
            )
        return [input_path.absolute()]

    if not input_path.is_dir():
        raise MashIDError(f"Input path does not exist: {input_path}")

    files: list[Path] = []
    broken: list[Path] = []
    visited: set[Path] = set()
    for root, dirs, filenames in os.walk(input_path, followlinks=True):
        visited.add(Path(root).resolve())
        # Skip hidden directories and symlink cycles (a directory already visited under another path).
        dirs[:] = sorted(
            d for d in dirs if not d.startswith(".") and (Path(root) / d).resolve() not in visited
        )
        for filename in sorted(filenames):
            if filename.startswith(".") or not is_sequence_file(filename):
                continue
            path = Path(root) / filename
            if not path.exists():  # a symbolic link whose target is gone
                broken.append(path)
                continue
            # Keep the path as found (not resolved): a symbolic link's own name defines the sample name,
            # which is how users rename samples without copying files.
            files.append(path.absolute())
    if broken:
        listed = "\n  ".join(f"{p} -> {os.readlink(p)}" for p in broken[:10])
        more = f"\n  ... and {len(broken) - 10} more" if len(broken) > 10 else ""
        raise MashIDError(f"{len(broken)} sequence file link(s) point to missing files:\n  {listed}{more}")
    if not files:
        raise MashIDError(
            f"No sequence files found in {input_path}. Accepted extensions: {', '.join(SEQ_EXTENSIONS)}"
        )
    return files


def group_samples(files: list[Path]) -> list[Sample]:
    """Group files sharing the same derived sample name (e.g. R1/R2 pairs) into samples."""
    grouped: dict[str, list[Path]] = {}
    for f in files:
        grouped.setdefault(sample_name_from_file(f), []).append(f)

    samples: list[Sample] = []
    for name in sorted(grouped):
        paths = sorted(set(grouped[name]))
        types = {seq_type_of(p) for p in paths}
        if len(types) > 1:
            raise MashIDError(
                f"Sample '{name}' mixes fasta and fastq files; rename them so they do not share a sample name: "
                + ", ".join(str(p) for p in paths)
            )
        basenames = [p.name for p in paths]
        if len(set(basenames)) != len(basenames):
            raise MashIDError(
                f"Sample '{name}' groups files with the same name from different directories, which is "
                "almost certainly two different samples: " + ", ".join(str(p) for p in paths)
                + ". Rename them or use --sample-sheet."
            )
        sample = Sample(name=name, files=paths, seq_type=types.pop())
        if len(paths) > 2:
            msg = f"{len(paths)} files grouped under sample '{name}'; all will be screened together"
            sample.warnings.append(msg)
            log.warning("%s: %s", msg, ", ".join(p.name for p in paths))
        samples.append(sample)
    return samples


def discover_samples(input_path: Path) -> list[Sample]:
    return group_samples(find_sequence_files(input_path))


def read_sample_sheet(path: Path, reserved: set[str] | None = None) -> tuple[list[Sample], list[str]]:
    """Read a TSV/CSV sample sheet: columns ``sample`` and ``file`` (or ``files``), one row per file or
    per sample with files separated by ';' or ','. Relative paths are resolved against the sheet's
    directory. Other columns are kept on each sample (``extra``) and returned as ``extra_columns``;
    a column whose name is in ``reserved`` (an output column) is renamed ``Sheet_<name>`` so it can
    never overwrite a result.
    """
    path = Path(path)
    if not path.is_file():
        raise MashIDError(f"Sample sheet not found: {path}")
    text = path.read_text(encoding="utf-8-sig")
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    if not lines:
        raise MashIDError(f"Sample sheet is empty: {path}")
    delimiter = "\t" if "\t" in lines[0] else ","
    reader = csv.DictReader(lines, delimiter=delimiter)
    cols = {c.strip().lower(): c for c in (reader.fieldnames or [])}
    sample_col = cols.get("sample") or cols.get("sample_id") or cols.get("name")
    file_col = cols.get("file") or cols.get("files") or cols.get("path") or cols.get("fastq")
    if not sample_col or not file_col:
        raise MashIDError(f"{path}: expected columns 'sample' and 'file' (or 'files'); found {reader.fieldnames}")
    extra_columns = [c for c in (reader.fieldnames or []) if c not in (sample_col, file_col)]
    renamed: dict[str, str] = {}
    for c in extra_columns:
        if reserved and c in reserved:
            renamed[c] = f"Sheet_{c}"
            log.warning("Sample sheet column %r clashes with an output column; reported as %r", c, renamed[c])
    extra_columns = [renamed.get(c, c) for c in extra_columns]

    grouped: dict[str, Sample] = {}
    for n, row in enumerate(reader, start=2):
        name = (row.get(sample_col) or "").strip()
        if not name:
            raise MashIDError(f"{path}: line {n} has no sample name")
        if "/" in name or name in (".", ".."):
            raise MashIDError(f"{path}: line {n}: sample name {name!r} cannot contain '/'")
        files = [f.strip() for f in re.split(r"[;,]", row.get(file_col) or "") if f.strip()]
        if not files:
            raise MashIDError(f"{path}: line {n} ({name}) lists no file")
        sample = grouped.get(name)
        if sample is None:
            extra = {renamed.get(c, c): (row.get(c) or "").strip() for c in (reader.fieldnames or [])
                     if c not in (sample_col, file_col)}
            sample = Sample(name=name, files=[], seq_type="", extra=extra)
            grouped[name] = sample
        for f in files:
            fpath = Path(f).expanduser()
            if not fpath.is_absolute():
                fpath = path.parent / fpath
            if not fpath.is_file():
                raise MashIDError(f"{path}: line {n} ({name}): file not found: {fpath}")
            if not is_sequence_file(fpath):
                raise MashIDError(f"{path}: line {n} ({name}): not a recognised sequence file: {fpath}")
            sample.files.append(fpath.resolve())

    samples: list[Sample] = []
    for name, sample in grouped.items():
        types = {seq_type_of(f) for f in sample.files}
        if len(types) > 1:
            raise MashIDError(f"Sample '{name}' mixes fasta and fastq files in the sample sheet")
        sample.seq_type = types.pop()
        sample.files = sorted(set(sample.files))
        samples.append(sample)
    return samples, extra_columns
