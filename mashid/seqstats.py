"""Count sequences and bases in fasta/fastq files (plain or gzipped) without external tools."""

from __future__ import annotations

import io
from pathlib import Path

try:  # python-isal gives a 2-3x faster, drop-in gzip decompressor.
    from isal import igzip as _gzip  # type: ignore
except ImportError:  # pragma: no cover
    import gzip as _gzip  # type: ignore

_BUFFER = 4 * 1024 * 1024
_GZIP_MAGIC = b"\x1f\x8b"


def open_sequence_file(path: Path | str) -> io.BufferedReader:
    """Open a sequence file for binary reading, transparently decompressing gzip (by magic bytes)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == _GZIP_MAGIC:
        return io.BufferedReader(_gzip.open(path, "rb"), _BUFFER)  # type: ignore[arg-type]
    return open(path, "rb", buffering=_BUFFER)


def count_fastq(fh) -> tuple[int, int]:
    reads = 0
    bases = 0
    for i, line in enumerate(fh):
        if i & 3 == 1:  # sequence line of each 4-line record
            reads += 1
            bases += len(line.rstrip(b"\r\n"))
    return reads, bases


def count_fasta(fh) -> tuple[int, int]:
    seqs = 0
    bases = 0
    for line in fh:
        if line.startswith(b">"):
            seqs += 1
        else:
            bases += len(line.rstrip())
    return seqs, bases


def sequence_stats(path: Path | str, seq_type: str) -> tuple[int, int]:
    """Return (number of sequences, total bases) for a fasta or fastq file."""
    counter = count_fastq if seq_type == "fastq" else count_fasta
    with open_sequence_file(path) as fh:
        return counter(fh)


def sample_stats(paths: list[Path], seq_type: str) -> tuple[int, int]:
    """Sum sequence_stats over all files of a sample."""
    n_total = bp_total = 0
    for p in paths:
        n, bp = sequence_stats(p, seq_type)
        n_total += n
        bp_total += bp
    return n_total, bp_total


class FastqSubsampler:
    """Stream the first N reads of one or more fastq files as raw bytes, counting reads and bases.

    Used for ``--max-reads``: the chunks are piped to ``mash screen -`` and ``sequences``/``bases``
    describe exactly what was screened. The per-file quota is ``max_reads`` split evenly.
    """

    def __init__(self, paths: list[Path], max_reads: int, chunk_size: int = _BUFFER) -> None:
        if max_reads < 1:
            raise ValueError("max_reads must be at least 1")
        self.paths = list(paths)
        self.max_reads = max_reads
        self.chunk_size = chunk_size
        self.sequences = 0
        self.bases = 0

    def chunks(self):
        per_file = -(-self.max_reads // len(self.paths))  # ceiling division
        buffer = bytearray()
        for path in self.paths:
            taken = 0
            with open_sequence_file(path) as fh:
                for i, line in enumerate(fh):
                    if i & 3 == 0 and taken >= per_file:
                        break
                    if i & 3 == 1:
                        taken += 1
                        self.sequences += 1
                        self.bases += len(line.rstrip(b"\r\n"))
                    buffer += line
                    if len(buffer) >= self.chunk_size:
                        yield bytes(buffer)
                        buffer.clear()
        if buffer:
            yield bytes(buffer)
