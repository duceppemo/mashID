"""Discover input sequence files and group them into samples."""

from __future__ import annotations

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
    """Return sequence files found at ``input_path`` (a file, or a directory searched recursively)."""
    input_path = Path(input_path)
    if input_path.is_file():
        if not is_sequence_file(input_path):
            raise MashIDError(
                f"Input file does not have a recognised extension ({', '.join(SEQ_EXTENSIONS)}): {input_path}"
            )
        return [input_path.resolve()]

    if not input_path.is_dir():
        raise MashIDError(f"Input path does not exist: {input_path}")

    files: list[Path] = []
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
            files.append((Path(root) / filename).resolve())
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
        sample = Sample(name=name, files=paths, seq_type=types.pop())
        if len(paths) > 2:
            msg = f"{len(paths)} files grouped under sample '{name}'; all will be screened together"
            sample.warnings.append(msg)
            log.warning("%s: %s", msg, ", ".join(p.name for p in paths))
        samples.append(sample)
    return samples


def discover_samples(input_path: Path) -> list[Sample]:
    return group_samples(find_sequence_files(input_path))
