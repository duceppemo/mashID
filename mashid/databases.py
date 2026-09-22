"""Registry of pre-built mashID databases and on-demand download with checksum verification."""

from __future__ import annotations

import hashlib
import logging
import os
import shutil
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from mashid import MashIDError, __version__

log = logging.getLogger(__name__)

ENV_DB_DIR = "MASHID_DB_DIR"
DEFAULT_DB_NAME = "mycobacteriaceae"
_CHUNK = 1024 * 1024


@dataclass(frozen=True)
class RemoteDatabase:
    name: str
    filename: str
    url: str
    md5: str
    size: int  # bytes
    description: str
    doi: str


REGISTRY: dict[str, RemoteDatabase] = {
    db.name: db
    for db in (
        RemoteDatabase(
            name="mycobacteriaceae",
            filename="mycobacteriaceae_2025-02-20.msh",
            url="https://ndownloader.figshare.com/files/52609139",
            md5="2fea9a28a2488b07f906e359e37acf6c",
            size=27_733_552,
            description="Mycobacteriaceae (NCBI taxon 1762), all NCBI genomes as of 2025-02-20, "
                        "dereplicated per species at 99.9% identity. k=21, s=1000.",
            doi="10.6084/m9.figshare.28489304.v1",
        ),
        RemoteDatabase(
            name="listeria",
            filename="listeria_2025-02-18.msh",
            url="https://ndownloader.figshare.com/files/52609082",
            md5="b62a2751898f63a7d3bc1091e8583a03",
            size=48_549_192,
            description="Listeria spp., all NCBI genomes as of 2025-02-18, dereplicated.",
            doi="10.6084/m9.figshare.28489262.v1",
        ),
        RemoteDatabase(
            name="progenomes3",
            filename="progenomes3.msh",
            url="https://ndownloader.figshare.com/files/39690601",
            md5="44e174ffc7b388d13e456bf9de92802b",
            size=346_965_392,
            description="proGenomes v3 representative genomes (all bacteria and archaea).",
            doi="10.6084/m9.figshare.22312282.v1",
        ),
        RemoteDatabase(
            name="refseq_bacteria",
            filename="2023-01-19_refseq_bacteria_derep_0.01.msh",
            url="https://ndownloader.figshare.com/files/39690568",
            md5="3d70107aefb96d76396049b1a9298290",
            size=684_479_872,
            description="RefSeq bacteria as of 2023-01-19, dereplicated at 99% identity.",
            doi="10.6084/m9.figshare.22312240.v2",
        ),
    )
}


def db_dir() -> Path:
    """Directory holding downloaded databases: $MASHID_DB_DIR, else $XDG_DATA_HOME/mashID/db."""
    env = os.environ.get(ENV_DB_DIR)
    if env:
        return Path(env).expanduser()
    xdg = os.environ.get("XDG_DATA_HOME")
    base = Path(xdg).expanduser() if xdg else Path.home() / ".local" / "share"
    return base / "mashID" / "db"


def installed_path(name: str, directory: Path | None = None) -> Path:
    return (directory or db_dir()) / REGISTRY[name].filename


def is_installed(name: str, directory: Path | None = None) -> bool:
    return installed_path(name, directory).is_file()


def resolve_database(spec: str | os.PathLike | None, directory: Path | None = None) -> Path:
    """Turn ``-d`` into a path: an existing file, a registry name, or the default database."""
    if spec is None or str(spec) == "":
        if is_installed(DEFAULT_DB_NAME, directory):
            return installed_path(DEFAULT_DB_NAME, directory)
        raise MashIDError(
            f"No database given and the default '{DEFAULT_DB_NAME}' database is not installed. "
            f"Run: mashID_download_db {DEFAULT_DB_NAME}   (or pass -d /path/to/database.msh)"
        )
    path = Path(spec).expanduser()
    if path.is_file():
        return path
    name = str(spec)
    if name in REGISTRY:
        if is_installed(name, directory):
            return installed_path(name, directory)
        raise MashIDError(f"Database '{name}' is not downloaded yet. Run: mashID_download_db {name}")
    raise MashIDError(
        f"Mash database not found: {spec}. Give a path to a .msh file or one of: {', '.join(REGISTRY)}"
    )


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            h.update(chunk)
    return h.hexdigest()


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def download(name: str, directory: Path | None = None, force: bool = False, verify: bool = True) -> Path:
    """Download a registry database to ``directory`` (default: db_dir()), verifying its MD5."""
    if name not in REGISTRY:
        raise MashIDError(f"Unknown database '{name}'. Available: {', '.join(REGISTRY)}")
    remote = REGISTRY[name]
    directory = directory or db_dir()
    dest = directory / remote.filename
    if dest.is_file() and not force:
        if not verify or md5sum(dest) == remote.md5:
            log.info("%s is already installed: %s", name, dest)
            return dest
        log.warning("%s exists but its checksum does not match; re-downloading", dest)

    directory.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(directory).free
    if free < remote.size * 1.1:
        raise MashIDError(f"Not enough free space in {directory}: need {_human(remote.size)}, have {_human(free)}")

    log.info("Downloading %s (%s) from %s", remote.filename, _human(remote.size), remote.url)
    tmp = dest.with_name(dest.name + ".part")
    digest = hashlib.md5()
    received = 0
    request = urllib.request.Request(remote.url, headers={"User-Agent": f"mashID/{__version__}"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response, open(tmp, "wb") as out:
            next_report = 0.1
            while True:
                chunk = response.read(_CHUNK)
                if not chunk:
                    break
                out.write(chunk)
                digest.update(chunk)
                received += len(chunk)
                if remote.size and received / remote.size >= next_report:
                    log.info("  %3.0f%% (%s)", 100 * received / remote.size, _human(received))
                    next_report += 0.1
    except (urllib.error.URLError, OSError) as exc:
        tmp.unlink(missing_ok=True)
        raise MashIDError(f"Download of {remote.url} failed: {exc}") from exc

    if verify and digest.hexdigest() != remote.md5:
        tmp.unlink(missing_ok=True)
        raise MashIDError(
            f"Checksum mismatch for {remote.filename}: expected {remote.md5}, got {digest.hexdigest()}. "
            "The download may be corrupted or the remote file may have changed."
        )
    os.replace(tmp, dest)
    log.info("Installed %s -> %s", name, dest)
    ensure_sidecar(dest)
    return dest


def ensure_sidecar(database: Path) -> Path | None:
    """Write <db>.metadata.tsv (organism names parsed from headers, reference lengths) if missing."""
    from mashid.mash import find_mash, info_table  # local import: avoids a cycle at module load
    from mashid.metadata import entries_from_sketches, sidecar_path, write_metadata

    path = sidecar_path(database)
    if path.is_file():
        return path
    try:
        exe = find_mash()
        entries, _ = entries_from_sketches(info_table(database, exe))
        write_metadata(path, entries)
    except (MashIDError, OSError) as exc:
        log.warning("Could not write metadata sidecar %s: %s", path, exc)
        return None
    log.info("Metadata sidecar written: %s (%d references)", path, len(entries))
    return path


def list_rows(directory: Path | None = None) -> list[dict[str, str]]:
    rows = []
    for name, db in REGISTRY.items():
        installed = is_installed(name, directory)
        rows.append({
            "Name": name + (" (default)" if name == DEFAULT_DB_NAME else ""),
            "Installed": "yes" if installed else "no",
            "Size": _human(db.size),
            "Description": db.description,
        })
    return rows
