"""Thin wrappers around the Mash executable."""

from __future__ import annotations

import contextlib
import logging
import shutil
import subprocess
import threading
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from mashid import MashIDError

log = logging.getLogger(__name__)


def find_mash() -> str:
    exe = shutil.which("mash")
    if exe is None:
        raise MashIDError(
            "The 'mash' executable was not found on PATH. Install it with: conda install -c bioconda mash"
        )
    return exe


def mash_version(exe: str = "mash") -> str:
    result = subprocess.run([exe, "--version"], capture_output=True, text=True)
    return result.stdout.strip() or result.stderr.strip()


def _fail(cmd: list[str], returncode: int, stderr: str) -> MashIDError:
    lines = stderr.strip().splitlines()
    tail = "\n".join(lines[-10:]) if lines else "(no error output)"
    return MashIDError(f"Command failed (exit {returncode}): {' '.join(cmd)}\n{tail}")


def _run(cmd: list[str]) -> subprocess.CompletedProcess:
    log.debug("Running: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise _fail(cmd, result.returncode, result.stderr)
    return result


def _run_with_stdin(cmd: list[str], chunks: Iterable[bytes]) -> subprocess.CompletedProcess:
    """Run ``cmd`` feeding ``chunks`` to its stdin from a helper thread (avoids pipe deadlocks)."""
    log.debug("Running (stdin): %s", " ".join(cmd))
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Detach stdin from the Popen object: communicate() would otherwise close it under the feeder.
    stdin = proc.stdin
    proc.stdin = None
    assert stdin is not None
    feed_error: list[BaseException] = []

    def feed() -> None:
        try:
            for chunk in chunks:
                stdin.write(chunk)
        except BrokenPipeError:
            pass  # mash exited early; its stderr will explain
        except BaseException as exc:  # noqa: BLE001 - re-raised in the main thread
            feed_error.append(exc)
        finally:
            with contextlib.suppress(OSError):
                stdin.close()

    feeder = threading.Thread(target=feed, daemon=True)
    feeder.start()
    stdout, stderr = proc.communicate()
    feeder.join()
    if feed_error:
        raise feed_error[0]
    if proc.returncode != 0:
        raise _fail(cmd, proc.returncode, stderr.decode(errors="replace"))
    return subprocess.CompletedProcess(cmd, proc.returncode, stdout.decode(errors="replace"),
                                       stderr.decode(errors="replace"))


@dataclass
class ScreenHit:
    """One line of ``mash screen`` output. Numeric fields keep Mash's original text for output."""

    identity: str
    shared_hashes: str  # e.g. "9876/10000"
    median_multiplicity: str
    p_value: str
    query_id: str
    comment: str

    @property
    def identity_f(self) -> float:
        return float(self.identity)

    @property
    def p_value_f(self) -> float:
        return float(self.p_value)

    @property
    def multiplicity_f(self) -> float:
        return float(self.median_multiplicity)

    @property
    def shared_f(self) -> int:
        return int(self.shared_hashes.split("/")[0])


def parse_screen_output(text: str) -> list[ScreenHit]:
    hits: list[ScreenHit] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            log.warning("Skipping malformed mash screen line: %r", line)
            continue
        comment = fields[5] if len(fields) > 5 else ""
        hits.append(ScreenHit(fields[0], fields[1], fields[2], fields[3], fields[4], comment))
    return hits


def screen(
    database: Path,
    files: list[Path] | None,
    threads: int = 1,
    min_identity: float = 0.9,
    max_p_value: float = 0.05,
    winner_take_all: bool = True,
    exe: str = "mash",
    stdin_chunks: Iterable[bytes] | None = None,
) -> list[ScreenHit]:
    """Run ``mash screen`` of ``database`` against read/assembly ``files``, or against ``stdin_chunks``."""
    cmd = [exe, "screen", "-p", str(max(1, threads)), "-i", str(min_identity), "-v", str(max_p_value)]
    if winner_take_all:
        cmd.append("-w")
    cmd.append(str(database))
    if stdin_chunks is not None:
        cmd.append("-")
        result = _run_with_stdin(cmd, stdin_chunks)
    else:
        if not files:
            raise ValueError("screen() needs files or stdin_chunks")
        cmd.extend(str(f) for f in files)
        result = _run(cmd)
    return parse_screen_output(result.stdout)


def sketch(
    list_file: Path,
    output_prefix: Path,
    threads: int = 1,
    sketch_size: int = 10000,
    kmer_size: int = 21,
    exe: str = "mash",
) -> Path:
    """Run ``mash sketch -l`` and return the path of the resulting ``.msh`` file."""
    cmd = [
        exe, "sketch",
        "-p", str(max(1, threads)),
        "-s", str(sketch_size),
        "-k", str(kmer_size),
        "-l", str(list_file),
        "-o", str(output_prefix),
    ]
    _run(cmd)
    return output_prefix.with_name(output_prefix.name + ".msh")


@dataclass
class SketchInfo:
    """One row of ``mash info -t``: the sketch's hash count, sequence length, id and comment."""

    hashes: int
    length: int
    query_id: str
    comment: str


def info_table(database: Path, exe: str = "mash") -> list[SketchInfo]:
    """Every sketch in a database, from ``mash info -t``."""
    result = _run([exe, "info", "-t", str(database)])
    rows: list[SketchInfo] = []
    for line in result.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            continue
        try:
            hashes, length = int(fields[0]), int(fields[1])
        except ValueError:
            hashes, length = 0, 0
        rows.append(SketchInfo(hashes, length, fields[2], fields[3] if len(fields) > 3 else ""))
    return rows


def sketch_count(database: Path, exe: str = "mash") -> int | None:
    """Number of sketches in a database, from ``mash info -H`` (None if unavailable)."""
    try:
        result = _run([exe, "info", "-H", str(database)])
    except MashIDError:
        return None
    for line in result.stdout.splitlines():
        if line.strip().startswith("Sketches:"):
            try:
                return int(line.split(":")[1])
            except ValueError:
                return None
    return None
