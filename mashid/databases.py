"""Registry of pre-built mashID databases and on-demand download with checksum verification."""

from __future__ import annotations

import glob
import hashlib
import http.client
import logging
import os
import secrets
import shutil
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from mashid import MashIDError, __version__
from mashid.fileutil import md5sum

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
    # Optional metadata sidecar (<stem>.metadata.tsv) published with the database.
    metadata_url: str | None = None
    metadata_md5: str | None = None
    metadata_size: int = 0
    # Optional confusable-reference sidecar (<stem>.confusable.tsv), from make_mashID_db
    # --confusable. Costs O(n^2) in references to compute, so publishing it with the database
    # saves every user that work. Without it mashID simply does not report the note.
    confusable_url: str | None = None
    confusable_md5: str | None = None
    confusable_size: int = 0


REGISTRY: dict[str, RemoteDatabase] = {
    db.name: db
    for db in (
        RemoteDatabase(
            name="mycobacteriaceae",
            filename="mycobacteriaceae_2026-10-22.msh",
            url="https://ndownloader.figshare.com/files/70524381",
            md5="a26156de99e6ee2f868a27d097fef702",
            size=299_559_160,
            description="Mycobacteriaceae (NCBI taxon 1762): 16,433 GenBank assemblies as of 2026-10-22, "
                        "dereplicated per species at 99.9% identity to 3728 references, with NCBI TaxIDs. "
                        "k=21, s=10000.",
            doi="10.6084/m9.figshare.33965176.v4",
            # v3 is a rebuild with the curated label corrections applied before binning, so the two
            # mislabelled references v2 renamed are gone: binned with their true relatives they were
            # near-duplicates and dereplication dropped them.
            metadata_url="https://ndownloader.figshare.com/files/70524378",
            metadata_md5="7f49515304fcf8afe77bc4a080109c75",
            metadata_size=1_114_036,
            confusable_url="https://ndownloader.figshare.com/files/70524375",
            confusable_md5="64b31294b3f7b9c6c9ef3a05a62b6b33",
            confusable_size=129_266,
        ),
        RemoteDatabase(
            name="mycobacteriaceae-2025",
            filename="mycobacteriaceae_2025-02-20.msh",
            url="https://ndownloader.figshare.com/files/52609139",
            md5="2fea9a28a2488b07f906e359e37acf6c",
            size=27_733_552,
            description="Previous Mycobacteriaceae database (2025-02-20): 3309 references, k=21, s=1000, no "
                        "TaxIDs. Kept for reproducibility; prefer 'mycobacteriaceae'.",
            doi="10.6084/m9.figshare.28489304.v1",
        ),
        RemoteDatabase(
            name="listeria",
            filename="listeria_2026-09-22.msh",
            url="https://ndownloader.figshare.com/files/69241489",
            md5="4d0f8f53d71dbcea526508f1b246ab2f",
            size=115_771_056,
            description="Listeria (NCBI taxon 1637): 7453 RefSeq assemblies as of 2026-09-22, dereplicated per "
                        "species at 99.9% identity to 1441 references (42 species/subspecies), with NCBI TaxIDs. "
                        "k=21, s=10000.",
            doi="10.6084/m9.figshare.33968620.v1",
            metadata_url="https://ndownloader.figshare.com/files/69241486",
            metadata_md5="ea50b80cf146e2027375c08e3ccdc94d",
            metadata_size=404_165,
        ),
        RemoteDatabase(
            name="listeria-2025",
            filename="listeria_2025-02-18.msh",
            url="https://ndownloader.figshare.com/files/52609082",
            md5="b62a2751898f63a7d3bc1091e8583a03",
            size=48_549_192,
            description="Previous Listeria database (2025-02-18), no TaxIDs. Kept for reproducibility; prefer "
                        "'listeria'.",
            doi="10.6084/m9.figshare.28489262.v1",
        ),
        RemoteDatabase(
            name="salmonella",
            filename="salmonella_2026-10-22.msh",
            url="https://ndownloader.figshare.com/files/70524918",
            md5="98576e2fac7f26333c7af6b25cef5b90",
            size=115_370_968,
            description="Salmonella (NCBI taxon 590): 21,543 RefSeq assemblies as of 2026-10-22, binned by "
                        "subspecies, dereplicated to 1436 references covering S. bongori and all six S. enterica "
                        "subspecies, with NCBI TaxIDs. k=21, s=10000. 29 mislabelled references were corrected "
                        "before binning. Species/subspecies level; serovars in names are those of the nearest "
                        "reference, not a serotype call.",
            doi="10.6084/m9.figshare.33968917.v3",
            metadata_url="https://ndownloader.figshare.com/files/70524915",
            metadata_md5="c1872d8a0be34c604f6ad1ec298bda32",
            metadata_size=425_546,
            confusable_url="https://ndownloader.figshare.com/files/70524912",
            confusable_md5="0b4e8427628a87e3e7ffbd2a7f3ec4b6",
            confusable_size=116_911,
        ),
        RemoteDatabase(
            name="brucella",
            filename="brucella_2026-10-22.msh",
            url="https://ndownloader.figshare.com/files/70524939",
            md5="64a85eb2a226e55f95bc71f57aa071bd",
            size=31_728_432,
            description="Brucella (NCBI taxon 234, including former Ochrobactrum): 2022 GenBank assemblies as of "
                        "2026-10-22, dereplicated at 99.97% identity to 395 references (29 species) with NCBI "
                        "TaxIDs. k=21, s=10000. 3 mislabelled references were corrected before binning. "
                        "Classical species are nearest-reference calls; biovars not resolved.",
            doi="10.6084/m9.figshare.33970168.v3",
            metadata_url="https://ndownloader.figshare.com/files/70524936",
            metadata_md5="c8759d5b7c31b30dac1dcee8bc0cc686",
            metadata_size=103_110,
            confusable_url="https://ndownloader.figshare.com/files/70524933",
            confusable_md5="e7b72a4060a8ced167dc5a1d81eb150d",
            confusable_size=137_285,
        ),
        RemoteDatabase(
            name="escherichia",
            filename="escherichia_shigella_2026-10-22.msh",
            url="https://ndownloader.figshare.com/files/70524879",
            md5="15c8890547a6171350095fb2af02c31a",
            size=362_190_880,
            description="Escherichia and Shigella (NCBI taxa 561, 620): 58,198 RefSeq assemblies as of "
                        "2026-10-22, dereplicated to 4510 references with NCBI TaxIDs. k=21, s=10000. Two "
                        "mislabelled Shigella references were corrected before binning. Species level; Shigella "
                        "species other than sonnei/flexneri are not reliably separable.",
            doi="10.6084/m9.figshare.33969760.v3",
            metadata_url="https://ndownloader.figshare.com/files/70524876",
            metadata_md5="a88cd194435a742f1b6055bcd043cfbf",
            metadata_size=1_095_360,
            confusable_url="https://ndownloader.figshare.com/files/70524873",
            confusable_md5="8273957a447c83bd67d6dc74cf1fc968",
            confusable_size=174_956,
        ),
        RemoteDatabase(
            name="refseq_viral",
            filename="refseq_viral_r237_2026-09-28.msh",
            url="https://ndownloader.figshare.com/files/69440466",
            md5="7394b54d24840d343ac77e527a53482d",
            size=304_371_760,
            description="All 19,625 viral genomes and segments of NCBI RefSeq release 237 (6427 phages), ICTV "
                        "species names and TaxIDs. k=21, s=2000. Virus identification and prophage or viral "
                        "contamination detection; no reference-length filter is applied.",
            doi="10.6084/m9.figshare.34018548.v1",
            metadata_url="https://ndownloader.figshare.com/files/69440463",
            metadata_md5="eaff27b0ff31105d6b2619247cfde6be",
            metadata_size=3_956_486,
        ),
        RemoteDatabase(
            name="progenomes4",
            filename="progenomes4_2026-09-22.msh",
            url="https://ndownloader.figshare.com/files/69242068",
            md5="daa010ff78267f365aef8d9a05eb9b02",
            size=536_227_432,
            description="proGenomes 4 representative genomes: 32,887 species-cluster representatives across "
                        "4384 genera of bacteria and archaea, current NCBI names and TaxIDs. k=21, s=2000. "
                        "For species-level screening of anything.",
            doi="10.6084/m9.figshare.33968806.v2",
            # v2 is a sidecar-only update: two Salmonella references regained the serovar formula
            # that was being truncated at its first comma, and Source_File now names the proGenomes
            # archive the references came from rather than the split directory the build deleted.
            # The sketch is unchanged and keeps its original URL.
            metadata_url="https://ndownloader.figshare.com/files/70525026",
            metadata_md5="d83eb894f4b896b25cdfbed1bd38b6a3",
            metadata_size=8_269_377,
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
            description="RefSeq bacteria as of 2023-01-19, dereplicated at 99% identity: 82,557 references "
                        "across 32,750 organisms, with current NCBI names and TaxIDs for 82,547 of them. "
                        "k=21, s=1000. Three years old; prefer 'progenomes4' for current names at species level.",
            doi="10.6084/m9.figshare.22312240.v3",
            # The sketch is from 2023 and unchanged. v3 adds the metadata sidecar it never had, so its
            # names and TaxIDs now come from NCBI rather than from each reference's fasta header, and
            # the curated label corrections reach it. The ten references withdrawn from RefSeq since
            # 2023 keep their header names and have no TaxID.
            metadata_url="https://ndownloader.figshare.com/files/70524954",
            metadata_md5="a6e545c7942a645fb2841d06cc3baab5",
            metadata_size=19_351_078,
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
        return path.absolute()  # never hand Mash a relative name that could look like an option
    name = str(spec)
    if name in REGISTRY:
        if is_installed(name, directory):
            return installed_path(name, directory)
        raise MashIDError(f"Database '{name}' is not downloaded yet. Run: mashID_download_db {name}")
    raise MashIDError(
        f"Mash database not found: {spec}. Give a path to a .msh file or one of: {', '.join(REGISTRY)}"
    )


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


_PART_PREFIX = ".part-"
_PART_ATTEMPTS = 8
# Only sweep temporary files old enough that no live download could own one. A shorter window would
# delete the partial file of a download running in another process, which is what unique names fixed.
_STALE_PART_AGE = 24 * 3600


def _sweep_stale_parts(dest: Path) -> None:
    """Remove this destination's abandoned temporary files, left by a killed download.

    Cleanup on the way out only covers exceptions; SIGKILL from a scheduler's walltime limit leaves
    the partial file behind, and for a 284 MB sketch that is worth reclaiming.
    """
    cutoff = time.time() - _STALE_PART_AGE
    # glob.escape: a destination named my[old].msh would otherwise read as a character class and
    # sweep another destination's files instead of its own.
    pattern = glob.escape(dest.name) + _PART_PREFIX + "*"
    for leftover in dest.parent.glob(pattern):
        try:
            if leftover.is_file() and leftover.stat().st_mtime < cutoff:
                leftover.unlink()
                log.info("Removed an abandoned partial download: %s", leftover.name)
        except OSError as exc:  # someone else's file, or a racing cleanup: not our problem
            log.debug("Could not remove %s (%s)", leftover, exc)


def _open_new_part(dest: Path) -> tuple[int, Path]:
    """Create a uniquely named temporary file beside ``dest``; returns its fd and path.

    ``O_CREAT | O_EXCL`` with mode 0o666 lets the kernel apply the user's umask as the file is
    created. That is unique like ``tempfile.mkstemp`` but without mkstemp's fixed 0600, so no chmod
    is needed afterwards -- and a chmod was both a regression on filesystems that refuse one
    (vfat, some bind mounts), where it threw away a complete verified download, and a way to ignore
    the umask a shared installation may rely on.
    """
    for _ in range(_PART_ATTEMPTS):
        path = dest.with_name(dest.name + _PART_PREFIX + secrets.token_hex(6))
        try:
            return os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666), path
        except FileExistsError:
            continue
    raise MashIDError(f"Could not create a temporary file for {dest.name} in {dest.parent}")


def _fetch(url: str, dest: Path, size: int, md5: str | None, verify: bool) -> None:
    """Stream ``url`` to ``dest`` through a temporary file, verifying the MD5 when given.

    The temporary name is unique per call, so two processes installing the same database at once --
    two cluster jobs that each fetch it on first use -- cannot write into one another's partial
    file. Sharing one ``<name>.part`` meant the second truncated what the first was still writing,
    and whichever finished first published a file with a hole of zeros in it under a checksum
    computed from its own complete stream.
    """
    log.info("Downloading %s (%s) from %s", dest.name, _human(size), url)
    dest.parent.mkdir(parents=True, exist_ok=True)
    _sweep_stale_parts(dest)
    fd, tmp = _open_new_part(dest)
    digest = hashlib.md5()
    received = 0
    request = urllib.request.Request(url, headers={"User-Agent": f"mashID/{__version__}"})
    try:
        # os.fdopen first, so the descriptor is owned by a context manager before anything that can
        # raise: evaluated after urlopen() in one `with`, a refused connection leaked it every time.
        with os.fdopen(fd, "wb") as out, urllib.request.urlopen(request, timeout=60) as response:
            next_report = 0.1
            while True:
                chunk = response.read(_CHUNK)
                if not chunk:
                    break
                out.write(chunk)
                digest.update(chunk)
                received += len(chunk)
                if size and received / size >= next_report:
                    log.info("  %3.0f%% (%s)", 100 * received / size, _human(received))
                    next_report += 0.1
        if verify and md5 and digest.hexdigest() != md5:
            raise MashIDError(
                f"Checksum mismatch for {dest.name}: expected {md5}, got {digest.hexdigest()}. "
                "The download may be corrupted or the remote file may have changed."
            )
        # Atomic, so an interrupted download never leaves a half-written database in place.
        os.replace(tmp, dest)
    # http.client.IncompleteRead -- a response cut short mid-stream -- is neither a URLError nor an
    # OSError, and a 284 MB transfer is long enough for that to happen.
    except (urllib.error.URLError, http.client.HTTPException, OSError) as exc:
        tmp.unlink(missing_ok=True)
        raise MashIDError(f"Download of {url} failed: {exc}") from exc
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def registry_name_of(database: Path) -> str | None:
    """The registry name a path belongs to, or None for a self-built database."""
    return next((n for n, r in REGISTRY.items() if r.filename == Path(database).name), None)


def refresh_advice(database: Path, rebuild_flag: str) -> str:
    """How to bring a stale sidecar back in step with its database.

    For a pre-built database the answer is to fetch the published sidecar again, never to rebuild
    it: ``make_mashID_db --annotate`` without ``--metadata`` takes every name from the sketch
    headers, which throws away the curated label corrections and all of the TaxIDs.
    """
    name = registry_name_of(database)
    if name is not None:
        return f"mashID_download_db {name}"
    return f"make_mashID_db {rebuild_flag} {Path(database).name}"


def _up_to_date(path: Path, md5: str | None, verify: bool) -> bool:
    if not path.is_file():
        return False
    if not verify or md5 is None:
        return True
    try:
        return md5sum(path) == md5
    except OSError as exc:
        # Unreadable, so it cannot be shown to match; re-downloading is the recoverable answer.
        log.warning("Could not read %s to check it (%s); treating it as out of date", path.name, exc)
        return False


def _sidecar_up_to_date(path: Path, md5: str | None) -> bool:
    """Like ``_up_to_date`` but always checksums, whatever ``--no-verify`` says.

    --no-verify is there to skip hashing a 284 MB sketch; a sidecar is a megabyte, and skipping it
    would leave a corrected label table beside the sketch it contradicts -- or, worse, install a new
    confusable table beside the deposited names it was computed to replace.

    The download of a sidecar is verified for the same reason (see ``download``): a decision that
    always hashes, paired with a fetch that does not, could never converge.
    """
    return _up_to_date(path, md5, verify=True)


def download(name: str, directory: Path | None = None, force: bool = False, verify: bool = True) -> Path:
    """Download a registry database (and its sidecars, when published) to ``directory``
    (default: db_dir()), verifying MD5 checksums. Existing, matching files are kept.

    A sidecar whose checksum differs from the registry's is replaced, so local edits to an installed
    ``.metadata.tsv`` -- the curated label corrections, for instance -- do not survive a later
    download of that database. Re-apply them, or keep the corrected copy elsewhere and point
    ``--db-metadata`` at it.
    """
    if name not in REGISTRY:
        raise MashIDError(f"Unknown database '{name}'. Available: {', '.join(REGISTRY)}")
    from mashid.confusable import confusable_path  # local imports: keep module load light
    from mashid.metadata import sidecar_path

    remote = REGISTRY[name]
    directory = directory or db_dir()
    dest = directory / remote.filename
    # One derivation of each sidecar path, shared with the readers (metadata.sidecar_path and
    # confusable.confusable_path), so that the two cannot disagree about where a sidecar lives.
    sidecar = sidecar_path(dest)
    confusable = confusable_path(dest)

    need_db = force or not _up_to_date(dest, remote.md5, verify)
    need_meta = remote.metadata_url is not None and (
        force or not _sidecar_up_to_date(sidecar, remote.metadata_md5))
    need_confusable = remote.confusable_url is not None and (
        force or not _sidecar_up_to_date(confusable, remote.confusable_md5))
    for path, needed in ((dest, need_db), (sidecar, need_meta), (confusable, need_confusable)):
        if path.is_file() and needed and not force:
            log.warning("%s exists but its checksum does not match; re-downloading", path.name)
    # Before the early return, not only inside _fetch: a download killed outright leaves its
    # partial file behind, and the next run of this command usually finds everything installed and
    # never reaches _fetch, so the 284 MB would sit there for good.
    if directory.is_dir():
        for path in (dest, sidecar, confusable):
            _sweep_stale_parts(path)
    if not need_db and not need_meta and not need_confusable:
        log.info("%s is already installed: %s", name, dest)
        return dest

    directory.mkdir(parents=True, exist_ok=True)
    needed = ((remote.size if need_db else 0) + (remote.metadata_size if need_meta else 0)
              + (remote.confusable_size if need_confusable else 0))
    free = shutil.disk_usage(directory).free
    if free < needed * 1.1:
        raise MashIDError(f"Not enough free space in {directory}: need {_human(needed)}, have {_human(free)}")

    if need_db:
        _fetch(remote.url, dest, remote.size, remote.md5, verify)
        log.info("Installed %s -> %s", name, dest)
    # The sidecars are fetched with verify=True whatever --no-verify says, to match the decision
    # above. Hashing to decide a file is stale and then not hashing what replaces it would reinstall
    # the same mismatching bytes on every run, for ever, and report success each time.
    if need_meta:
        assert remote.metadata_url is not None
        _fetch(remote.metadata_url, sidecar, remote.metadata_size, remote.metadata_md5, verify=True)
        log.info("Installed metadata -> %s", sidecar)
    if need_confusable:
        assert remote.confusable_url is not None
        _fetch(remote.confusable_url, confusable, remote.confusable_size,
               remote.confusable_md5, verify=True)
        log.info("Installed inseparable-organism pairs -> %s", confusable)
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
