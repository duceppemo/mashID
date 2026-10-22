"""Small file helpers used in more than one place: checksums and non-clobbering backups."""

from __future__ import annotations

import datetime as dt
import hashlib
import shutil
from pathlib import Path

_CHUNK = 1 << 20


def sidecar_for(database: Path, suffix: str) -> Path:
    """``/x/db.msh`` + ``.metadata.tsv`` -> ``/x/db.metadata.tsv``; a name without ``.msh`` keeps it.

    One derivation for every sidecar, so that the writer, the reader and the downloader cannot
    disagree about where a database's companion files live.
    """
    database = Path(database)
    stem = database.name[:-4] if database.name.endswith(".msh") else database.name
    return database.with_name(stem + suffix)


def md5sum(path: Path) -> str:
    """MD5 of a file, read in chunks so that a 284 MB database does not land in memory."""
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def backup(path: Path) -> Path:
    """Copy ``path`` to ``<name>.bak-<date>``, never overwriting an existing backup.

    Running a correction twice in one day must not destroy the copy of the original that the first
    run made, so a taken name gets a ``.2``, ``.3``, ... suffix.
    """
    path = Path(path)
    stamp = dt.date.today().strftime("%Y%m%d")
    dest = path.with_name(f"{path.name}.bak-{stamp}")
    n = 2
    while dest.exists():
        dest = path.with_name(f"{path.name}.bak-{stamp}.{n}")
        n += 1
    shutil.copy2(path, dest)
    return dest
