import hashlib
import http.server
import threading
from functools import partial

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
    assert not (dest_dir / "fake.msh.part").exists()
    # second call is a no-op (already installed, checksum matches)
    assert download("fake", dest_dir) == path
    # corrupted local copy is re-downloaded
    path.write_bytes(b"corrupt")
    assert download("fake", dest_dir).read_bytes() == payload

    with pytest.raises(MashIDError, match="Checksum mismatch"):
        download("bad", dest_dir)
    assert not (dest_dir / "bad.msh").exists() and not (dest_dir / "bad.msh.part").exists()
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
