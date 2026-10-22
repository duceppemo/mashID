import hashlib
import http.server
import threading
from functools import partial
from pathlib import Path

import pytest

from mashid import MashIDError
from mashid.databases import (
    DEFAULT_DB_NAME,
    REGISTRY,
    RemoteDatabase,
    db_dir,
    download,
    list_rows,
    resolve_database,
)
from mashid.download_cli import main as download_main


def test_registry_is_well_formed():
    assert DEFAULT_DB_NAME in REGISTRY
    for name, db in REGISTRY.items():
        assert db.name == name and db.filename.endswith(".msh")
        assert db.url.startswith("https://") and len(db.md5) == 32 and db.size > 0
        if db.metadata_url:
            assert db.metadata_url.startswith("https://") and len(db.metadata_md5 or "") == 32 and db.metadata_size > 0
    assert REGISTRY[DEFAULT_DB_NAME].metadata_url  # the default database ships its TaxID sidecar


def test_db_dir_env(monkeypatch, tmp_path):
    monkeypatch.setenv("MASHID_DB_DIR", str(tmp_path / "dbs"))
    assert db_dir() == tmp_path / "dbs"
    monkeypatch.delenv("MASHID_DB_DIR")
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg"))
    assert db_dir() == tmp_path / "xdg" / "mashID" / "db"


def test_resolve_database(monkeypatch, tmp_path):
    monkeypatch.setenv("MASHID_DB_DIR", str(tmp_path))
    with pytest.raises(MashIDError, match="mashID_download_db mycobacteriaceae"):
        resolve_database(None)
    with pytest.raises(MashIDError, match="not downloaded yet"):
        resolve_database("listeria")
    with pytest.raises(MashIDError, match="not found"):
        resolve_database(tmp_path / "nope.msh")
    f = tmp_path / "custom.msh"
    f.write_bytes(b"x")
    assert resolve_database(f) == f
    installed = tmp_path / REGISTRY[DEFAULT_DB_NAME].filename
    installed.write_bytes(b"x")
    assert resolve_database(None) == installed
    assert resolve_database(DEFAULT_DB_NAME) == installed


@pytest.fixture
def http_file(tmp_path):
    """Serve a small fake database over HTTP and patch the registry to point at it."""
    payload = b"MASH" + bytes(range(256)) * 100
    served = tmp_path / "served"
    served.mkdir()
    (served / "fake.msh").write_bytes(payload)
    (served / "fake.metadata.tsv").write_text("Accession\tOrganism\tTaxID\nX1\tGenusa one\t7\n")
    handler = partial(http.server.SimpleHTTPRequestHandler, directory=str(served))
    handler.log_message = lambda *a, **k: None  # type: ignore[assignment]
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    url = f"http://127.0.0.1:{server.server_address[1]}/fake.msh"
    yield url, payload, hashlib.md5(payload).hexdigest()
    server.shutdown()


def test_download_verifies_checksum(monkeypatch, tmp_path, http_file):
    url, payload, md5 = http_file
    good = RemoteDatabase("fake", "fake.msh", url, md5, len(payload), "test", "doi")
    bad = RemoteDatabase("bad", "bad.msh", url, "0" * 32, len(payload), "test", "doi")
    monkeypatch.setitem(REGISTRY, "fake", good)
    monkeypatch.setitem(REGISTRY, "bad", bad)
    dest_dir = tmp_path / "dbs"

    path = download("fake", dest_dir)
    assert path == dest_dir / "fake.msh" and path.read_bytes() == payload
    assert not list(dest_dir.glob("fake.msh.part-*"))
    # second call is a no-op (already installed, checksum matches)
    assert download("fake", dest_dir) == path
    # corrupted local copy is re-downloaded
    path.write_bytes(b"corrupt")
    assert download("fake", dest_dir).read_bytes() == payload

    with pytest.raises(MashIDError, match="Checksum mismatch"):
        download("bad", dest_dir)
    assert not (dest_dir / "bad.msh").exists() and not list(dest_dir.glob("bad.msh.part-*"))
    with pytest.raises(MashIDError, match="Unknown database"):
        download("nope", dest_dir)


def test_download_cli(monkeypatch, tmp_path, http_file, capsys):
    url, payload, md5 = http_file
    monkeypatch.setitem(REGISTRY, "fake", RemoteDatabase("fake", "fake.msh", url, md5, len(payload), "t", "d"))
    monkeypatch.setenv("MASHID_DB_DIR", str(tmp_path / "dbs"))
    assert download_main([]) == 0
    out = capsys.readouterr().out
    assert "mycobacteriaceae (default)" in out and "Installed" in out
    assert download_main(["fake"]) == 0
    assert (tmp_path / "dbs" / "fake.msh").is_file()
    assert [r["Installed"] for r in list_rows() if r["Name"].startswith("fake")] == ["yes"]
    assert download_main(["nope"]) == 1


def test_download_fetches_the_metadata_sidecar(monkeypatch, tmp_path, http_file):
    url, payload, md5 = http_file
    meta_url = url.replace("fake.msh", "fake.metadata.tsv")
    meta_bytes = b"Accession\tOrganism\tTaxID\nX1\tGenusa one\t7\n"
    remote = RemoteDatabase("withmeta", "fake.msh", url, md5, len(payload), "t", "d",
                            metadata_url=meta_url, metadata_md5=hashlib.md5(meta_bytes).hexdigest(),
                            metadata_size=len(meta_bytes))
    monkeypatch.setitem(REGISTRY, "withmeta", remote)
    dest_dir = tmp_path / "dbs"
    download("withmeta", dest_dir)
    assert (dest_dir / "fake.metadata.tsv").read_bytes() == meta_bytes  # published sidecar, not a generated one
    # a stale local sidecar is replaced on the next call, the database itself is kept
    (dest_dir / "fake.metadata.tsv").write_text("stale")
    download("withmeta", dest_dir)
    assert (dest_dir / "fake.metadata.tsv").read_bytes() == meta_bytes
    bad = RemoteDatabase("badmeta", "bad.msh", url, md5, len(payload), "t", "d",
                         metadata_url=meta_url, metadata_md5="0" * 32, metadata_size=len(meta_bytes))
    monkeypatch.setitem(REGISTRY, "badmeta", bad)
    with pytest.raises(MashIDError, match="Checksum mismatch for bad.metadata.tsv"):
        download("badmeta", dest_dir)


def test_download_fetches_the_confusable_sidecar(monkeypatch, tmp_path, http_file):
    url, payload, md5 = http_file
    conf_bytes = (b"Accession\tOrganism\tConfusable_Organism\tDistance\tNearest_Accession\n"
                  b"X1\tGenusa one\tGenusa two\t0.001000\tX2\n")
    (tmp_path / "served" / "fake.confusable.tsv").write_bytes(conf_bytes)
    conf_url = url.replace("fake.msh", "fake.confusable.tsv")
    remote = RemoteDatabase("withconf", "fake.msh", url, md5, len(payload), "t", "d",
                            confusable_url=conf_url,
                            confusable_md5=hashlib.md5(conf_bytes).hexdigest(),
                            confusable_size=len(conf_bytes))
    monkeypatch.setitem(REGISTRY, "withconf", remote)
    dest_dir = tmp_path / "dbs"
    download("withconf", dest_dir)
    assert (dest_dir / "fake.confusable.tsv").read_bytes() == conf_bytes
    # a stale local copy is replaced on the next call
    (dest_dir / "fake.confusable.tsv").write_text("stale")
    download("withconf", dest_dir)
    assert (dest_dir / "fake.confusable.tsv").read_bytes() == conf_bytes
    bad = RemoteDatabase("badconf", "bad.msh", url, md5, len(payload), "t", "d",
                         confusable_url=conf_url, confusable_md5="0" * 32,
                         confusable_size=len(conf_bytes))
    monkeypatch.setitem(REGISTRY, "badconf", bad)
    with pytest.raises(MashIDError, match="Checksum mismatch for bad.confusable.tsv"):
        download("badconf", dest_dir)


def test_refresh_advice_points_at_the_right_command():
    """make_mashID_db --annotate without --metadata rebuilds a sidecar from the sketch headers,
    which discards curated labels and every TaxID. Never advise it for a pre-built database."""
    from mashid.databases import REGISTRY as R
    from mashid.databases import refresh_advice, registry_name_of
    known = Path("/db") / R["mycobacteriaceae"].filename
    assert registry_name_of(known) == "mycobacteriaceae"
    assert refresh_advice(known, "--annotate") == "mashID_download_db mycobacteriaceae"
    assert refresh_advice(known, "--confusable") == "mashID_download_db mycobacteriaceae"
    mine = Path("/db/my_own.msh")
    assert registry_name_of(mine) is None
    assert refresh_advice(mine, "--annotate") == "make_mashID_db --annotate my_own.msh"


def test_no_verify_still_refreshes_the_sidecars(monkeypatch, tmp_path, http_file):
    """--no-verify skips hashing a 284 MB sketch. It must not leave a new confusable table beside
    the deposited names it was computed to replace."""
    url, payload, md5 = http_file
    meta_bytes = b"Accession\tOrganism\tTaxID\nX1\tGenusa one\t7\n"
    conf_bytes = (b"Accession\tOrganism\tConfusable_Organism\tDistance\tNearest_Accession\n"
                  b"X1\tGenusa one\tGenusa two\t0.001000\tX2\n")
    (tmp_path / "served" / "fake.confusable.tsv").write_bytes(conf_bytes)
    remote = RemoteDatabase("nv", "fake.msh", url, md5, len(payload), "t", "d",
                            metadata_url=url.replace("fake.msh", "fake.metadata.tsv"),
                            metadata_md5=hashlib.md5(meta_bytes).hexdigest(),
                            metadata_size=len(meta_bytes),
                            confusable_url=url.replace("fake.msh", "fake.confusable.tsv"),
                            confusable_md5=hashlib.md5(conf_bytes).hexdigest(),
                            confusable_size=len(conf_bytes))
    monkeypatch.setitem(REGISTRY, "nv", remote)
    dest_dir = tmp_path / "dbs"
    dest_dir.mkdir()
    (dest_dir / "fake.msh").write_bytes(payload)
    (dest_dir / "fake.metadata.tsv").write_text("Accession\tOrganism\nX1\tDeposited name\n")
    (dest_dir / "fake.confusable.tsv").write_text("stale")
    download("nv", dest_dir, verify=False)
    assert (dest_dir / "fake.metadata.tsv").read_bytes() == meta_bytes
    assert (dest_dir / "fake.confusable.tsv").read_bytes() == conf_bytes


def test_a_sidecar_whose_registry_md5_is_wrong_fails_instead_of_looping(monkeypatch, tmp_path,
                                                                        http_file):
    """The decision always hashes the sidecar. If the fetch did not, a registry MD5 that disagrees
    with what the URL serves would re-download and reinstall the same bytes on every single run,
    reporting success each time and never converging."""
    url, payload, md5 = http_file
    meta_bytes = b"Accession\tOrganism\tTaxID\nX1\tGenusa one\t7\n"
    remote = RemoteDatabase("wrongmd5", "fake.msh", url, md5, len(payload), "t", "d",
                            metadata_url=url.replace("fake.msh", "fake.metadata.tsv"),
                            metadata_md5="0" * 32, metadata_size=len(meta_bytes))
    monkeypatch.setitem(REGISTRY, "wrongmd5", remote)
    dest_dir = tmp_path / "dbs"
    for verify in (True, False):
        with pytest.raises(MashIDError, match="Checksum mismatch for fake.metadata.tsv"):
            download("wrongmd5", dest_dir, verify=verify)
    assert not list(dest_dir.glob("*.part-*"))


def test_an_unreadable_sidecar_is_treated_as_stale(monkeypatch, tmp_path, http_file, caplog):
    """_up_to_date now hashes the sidecar whatever --no-verify says, so it must survive a file it
    cannot read rather than raising OSError out of a download."""
    url, payload, md5 = http_file
    meta_bytes = b"Accession\tOrganism\tTaxID\nX1\tGenusa one\t7\n"
    remote = RemoteDatabase("unread", "fake.msh", url, md5, len(payload), "t", "d",
                            metadata_url=url.replace("fake.msh", "fake.metadata.tsv"),
                            metadata_md5=hashlib.md5(meta_bytes).hexdigest(),
                            metadata_size=len(meta_bytes))
    monkeypatch.setitem(REGISTRY, "unread", remote)
    dest_dir = tmp_path / "dbs"
    dest_dir.mkdir()
    sidecar = dest_dir / "fake.metadata.tsv"
    sidecar.write_text("unreadable")
    sidecar.chmod(0o000)
    try:
        with caplog.at_level("WARNING"):
            download("unread", dest_dir, verify=False)
    finally:
        sidecar.chmod(0o644)
    assert "Could not read fake.metadata.tsv" in caplog.text
    assert sidecar.read_bytes() == meta_bytes


def test_two_downloads_at_once_do_not_write_into_one_partial_file(monkeypatch, tmp_path, http_file):
    """Two cluster jobs installing the same database each fetched through one <name>.part: the
    second truncated what the first was still writing, and whichever finished first published a
    file with a hole of zeros in it under a checksum taken from its own complete stream."""
    import threading

    url, payload, md5 = http_file
    monkeypatch.setitem(REGISTRY, "race", RemoteDatabase("race", "race.msh", url, md5, len(payload),
                                                         "t", "d"))
    dest_dir = tmp_path / "dbs"
    dest_dir.mkdir()
    started = threading.Barrier(2)
    errors: list[BaseException] = []

    def go():
        try:
            started.wait(timeout=10)
            download("race", dest_dir)
        except BaseException as exc:            # noqa: BLE001 - recorded and re-raised below
            errors.append(exc)

    threads = [threading.Thread(target=go) for _ in range(2)]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=30)
    assert not errors, errors
    assert (dest_dir / "race.msh").read_bytes() == payload
    assert not list(dest_dir.glob("*.part*"))   # each call cleans up its own temporary file


def test_installed_files_are_readable_by_everyone(tmp_path, http_file, monkeypatch):
    """mkstemp creates 0600 and os.replace carries the mode over, which would install a database
    only its downloader can read -- and the shared MASHID_DB_DIR case is the whole point."""
    import stat

    url, payload, md5 = http_file
    meta_bytes = b"Accession\tOrganism\tTaxID\nX1\tGenusa one\t7\n"
    monkeypatch.setitem(REGISTRY, "modes", RemoteDatabase(
        "modes", "fake.msh", url, md5, len(payload), "t", "d",
        metadata_url=url.replace("fake.msh", "fake.metadata.tsv"),
        metadata_md5=hashlib.md5(meta_bytes).hexdigest(), metadata_size=len(meta_bytes)))
    dest_dir = tmp_path / "dbs"
    download("modes", dest_dir)
    for name in ("fake.msh", "fake.metadata.tsv"):
        mode = stat.S_IMODE((dest_dir / name).stat().st_mode)
        assert mode & stat.S_IRGRP and mode & stat.S_IROTH, f"{name} is {oct(mode)}"


def test_a_failed_download_does_not_leak_a_descriptor(tmp_path, monkeypatch):
    """mkstemp's fd is only adopted by os.fdopen. Evaluated after urlopen() in one `with`, a
    refused connection leaked it on every attempt."""
    import os
    import urllib.error
    import urllib.request

    from mashid.databases import _fetch

    def refuse(*a, **k):
        raise urllib.error.URLError("connection refused")

    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    before = len(os.listdir("/proc/self/fd"))
    for _ in range(5):
        with pytest.raises(MashIDError, match="failed"):
            _fetch("http://localhost:1/x.msh", tmp_path / "x.msh", 10, None, True)
    assert len(os.listdir("/proc/self/fd")) == before
    assert not list(tmp_path.glob("*.part-*"))


def test_an_abandoned_partial_download_is_swept_but_a_live_one_is_not(tmp_path, http_file):
    """A scheduler's walltime SIGKILL leaves the temp file behind; a download running right now in
    another process must keep its own."""
    import os

    url, payload, md5 = http_file
    dest_dir = tmp_path / "dbs"
    dest_dir.mkdir()
    old = dest_dir / "fake.msh.part-oldone"
    old.write_bytes(b"abandoned")
    os.utime(old, (0, 0))
    fresh = dest_dir / "fake.msh.part-livenow"
    fresh.write_bytes(b"in flight elsewhere")
    other = dest_dir / "another.msh.part-oldone"
    other.write_bytes(b"not ours")
    os.utime(other, (0, 0))

    from mashid.databases import _fetch
    _fetch(url, dest_dir / "fake.msh", len(payload), md5, True)
    assert not old.exists()                 # a day-old leftover of this destination goes
    assert fresh.exists()                   # one that could still be in use stays
    assert other.exists()                   # another destination's is none of our business


def test_a_truncated_transfer_is_a_clean_error(monkeypatch, tmp_path, http_file):
    """http.client.IncompleteRead is neither a URLError nor an OSError, so it used to traceback."""
    import http.client
    import urllib.request

    url, payload, md5 = http_file
    monkeypatch.setitem(REGISTRY, "cut", RemoteDatabase("cut", "cut.msh", url, md5, len(payload),
                                                        "t", "d"))

    class Cut:
        def read(self, _n):
            raise http.client.IncompleteRead(b"half")

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **k: Cut())
    dest_dir = tmp_path / "dbs"
    with pytest.raises(MashIDError, match="failed"):
        download("cut", dest_dir)
    assert not list(dest_dir.glob("*.part-*"))        # no orphaned partial download left behind


def test_download_skips_a_confusable_sidecar_that_is_not_published(monkeypatch, tmp_path, http_file):
    """Every registry entry today has no confusable_url; nothing must be fetched or demanded."""
    url, payload, md5 = http_file
    remote = RemoteDatabase("noconf", "fake.msh", url, md5, len(payload), "t", "d")
    monkeypatch.setitem(REGISTRY, "noconf", remote)
    dest_dir = tmp_path / "dbs"
    download("noconf", dest_dir)
    assert (dest_dir / "fake.msh").is_file()
    assert not (dest_dir / "fake.confusable.tsv").exists()
    download("noconf", dest_dir)             # second call is a no-op, not a re-download
