"""Tests for database-level confusability: inseparable references of different organisms."""

import pytest

from mashid import MashIDError
from mashid.confusable import (
    NOTE_PREFIX,
    Neighbour,
    confusable_path,
    distinguishable,
    note_for,
    read,
    write,
)
from mashid.metadata import DbEntry


def test_confusable_path(tmp_path):
    assert confusable_path(tmp_path / "db.msh") == tmp_path / "db.confusable.tsv"
    assert confusable_path(tmp_path / "db") == tmp_path / "db.confusable.tsv"


# --- which name pairs are a genuine conflict ----------------------------------------------------

@pytest.mark.parametrize("a,b", [
    ("Mycobacterium tuberculosis variant bovis", "Mycobacterium tuberculosis"),
    ("Mycolicibacterium goodii", "Mycolicibacterium wolinskyi"),
    ("Mycobacterium canetti", "Mycobacterium tuberculosis"),
    ("Mycobacterium avium subsp. paratuberculosis", "Mycobacterium avium subsp. hominissuis"),
    ("Mycobacterium intracellulare", "Mycobacterium paraintracellulare"),
])
def test_different_organisms_conflict(a, b):
    assert distinguishable(a, b) and distinguishable(b, a)


@pytest.mark.parametrize("a,b", [
    ("Mycobacterium abscessus", "Mycobacterium abscessus"),               # the same name
    ("Mycobacterium abscessus", "Mycobacteroides abscessus"),             # renamed genus
    ("Mycobacteroides chelonae", "[Mycobacterium] chelonae"),             # renamed genus, bracketed
    ("Mycobacterium pinniadriaticum", "Mycobacterium sp. P3"),            # unnamed species
    ("Mycobacterium sp.", "Mycobacterium fortuitum"),
    ("Mycolicibacterium sp. D3", "Mycolicibacterium parafortuitum"),
    ("Mycobacterium bovis", "NA"),
    ("Mycobacterium bovis", "unknown"),
    ("Mycobacterium bovis", ""),
    ("Lausannevirus", "Mycobacterium bovis"),                             # not a binomial
])
def test_same_or_unusable_names_do_not_conflict(a, b):
    assert not distinguishable(a, b) and not distinguishable(b, a)


def test_candidatus_is_read_past_the_prefix():
    assert not distinguishable("Candidatus Mycobacterium abscessus", "Mycobacteroides abscessus")
    assert distinguishable("Candidatus Mycobacterium alpha", "Mycobacterium beta")


# --- the table round-trips ----------------------------------------------------------------------

def _metadata():
    return {"GCA_1.1": DbEntry("GCA_1.1", "Genusa one", "1"),
            "GCA_2.1": DbEntry("GCA_2.1", "Genusa two", "2")}


def test_write_read_roundtrip_and_ordering(tmp_path):
    p = tmp_path / "db.confusable.tsv"
    neighbours = {"GCA_1.1": [Neighbour("Genusa two", 0.004, "GCA_2.1"),
                              Neighbour("Genusa three", 0.001, "GCA_3.1")]}
    write(p, neighbours, _metadata())
    rows = p.read_text().splitlines()
    assert rows[0].split("\t")[0] == "Accession"
    assert rows[1].split("\t")[2] == "Genusa three"      # nearest pair sorts first in the file
    back = read(p)
    assert [n.organism for n in back["GCA_1.1"]] == ["Genusa three", "Genusa two"]
    assert back["GCA_1.1"][0].distance == pytest.approx(0.001)
    assert back["GCA_1.1"][0].accession == "GCA_3.1"


def test_read_rejects_a_table_without_the_columns(tmp_path):
    p = tmp_path / "db.confusable.tsv"
    p.write_text("foo\tbar\nx\ty\n")
    with pytest.raises(MashIDError, match="missing column"):
        read(p)
    with pytest.raises(MashIDError, match="not found"):
        read(tmp_path / "nope.tsv")


def test_read_skips_unusable_rows(tmp_path):
    p = tmp_path / "db.confusable.tsv"
    p.write_text("Accession\tOrganism\tConfusable_Organism\tDistance\tNearest_Accession\n"
                 "\tGenusa one\tGenusa two\t0.001\tGCA_2.1\n"          # no accession
                 "GCA_1.1\tGenusa one\t\t0.001\tGCA_2.1\n"             # no organism
                 "GCA_1.1\tGenusa one\tGenusa two\tnotanumber\tGCA_2.1\n"
                 "GCA_1.1\tGenusa one\tGenusa three\t0.002\tGCA_3.1\n")
    assert [n.organism for n in read(p)["GCA_1.1"]] == ["Genusa three"]


# --- the note -----------------------------------------------------------------------------------

def test_note_lists_other_organisms_nearest_first():
    table = {"GCA_1.1": [Neighbour("Genusa two", 0.001, "GCA_2.1"),
                         Neighbour("Genusa three", 0.004, "GCA_3.1")]}
    note = note_for("GCA_1.1", "Genusa one", table)
    assert note == NOTE_PREFIX + "Genusa two, Genusa three"


def test_note_caps_the_list_and_counts_the_rest():
    table = {"GCA_1.1": [Neighbour(f"Genusa {n}", 0.001 * i, f"GCA_{i}.1")
                         for i, n in enumerate(("two", "three", "four", "five", "six"), 2)]}
    note = note_for("GCA_1.1", "Genusa one", table)
    assert note.endswith("and 2 more")
    assert note.count(",") == 2


def test_note_is_none_when_nothing_to_say():
    table = {"GCA_1.1": [Neighbour("Genusa two", 0.001, "GCA_2.1")]}
    assert note_for("GCA_9.9", "Genusa one", table) is None       # reference not in the table
    assert note_for("GCA_1.1", "Genusa one", None) is None        # no table loaded
    assert note_for("GCA_1.1", "Genusa one", {}) is None
    # the top hit's own organism never counts as a conflict with itself
    assert note_for("GCA_1.1", "Genusa two", {"GCA_1.1": [Neighbour("Genusa two", 0.001, "x")]}) is None


def test_load_finds_the_sidecar_or_returns_none(tmp_path):
    from mashid.confusable import load
    db = tmp_path / "db.msh"
    db.write_bytes(b"not really a sketch")
    assert load(db, None) is None                      # no sidecar: the check is simply off
    write(confusable_path(db), {"GCA_1.1": [Neighbour("Genusa two", 0.001, "GCA_2.1")]}, _metadata())
    assert list(load(db, None)) == ["GCA_1.1"]
    other = tmp_path / "elsewhere.tsv"
    write(other, {"GCA_9.9": [Neighbour("Genusa two", 0.001, "GCA_2.1")]}, _metadata())
    assert list(load(db, other)) == ["GCA_9.9"]        # explicit override wins


def test_load_warns_when_the_table_predates_the_database(tmp_path, caplog):
    import os

    from mashid.confusable import load
    db = tmp_path / "db.msh"
    db.write_bytes(b"not really a sketch")
    table = confusable_path(db)
    write(table, {"GCA_1.1": [Neighbour("Genusa two", 0.001, "GCA_2.1")]}, _metadata())
    os.utime(table, (1, 1))                            # table older than the database
    with caplog.at_level("WARNING"):
        assert load(db, None) is not None
    assert "older than the database" in caplog.text


def test_an_empty_table_is_not_an_error(tmp_path):
    from mashid.confusable import load
    db = tmp_path / "db.msh"
    db.write_bytes(b"not really a sketch")
    write(confusable_path(db), {}, _metadata())
    table = load(db, None)
    assert table == {}
    assert note_for("GCA_1.1", "Genusa one", table) is None


def test_confusable_distance_must_be_positive(tmp_path):
    from mashid.makedb import confusable_database
    db = tmp_path / "db.msh"
    db.write_bytes(b"not really a sketch")
    write(confusable_path(db), {}, _metadata())
    with pytest.raises(MashIDError, match="greater than 0"):
        confusable_database(db, None, 1, 0.0)


# --- mash dist streaming ------------------------------------------------------------------------

def _fake_mash(tmp_path, stdout_lines, exit_code=0, stderr=""):
    """An executable stub standing in for the mash binary."""
    import os
    import sys as _sys
    stub = tmp_path / "fake_mash.py"
    body = "\n".join(f"print({line!r})" for line in stdout_lines)
    stub.write_text(f"#!{_sys.executable}\nimport sys\n{body}\n"
                    f"sys.stderr.write({stderr!r})\nsys.exit({exit_code})\n")
    os.chmod(stub, 0o755)
    return str(stub)


def test_dist_pairs_parses_and_skips_unusable_lines(tmp_path):
    from mashid.mash import dist_pairs
    db = tmp_path / "db.msh"
    db.write_bytes(b"x")
    exe = _fake_mash(tmp_path, [
        "/db/GCA_1.1.fna\t/db/GCA_2.1.fna\t0.001\t0\t900/1000",
        "short\tline",                                        # too few fields
        "/db/GCA_1.1.fna\t/db/GCA_3.1.fna\tnotanumber\t0\t1",  # unparseable distance
        "/db/GCA_2.1.fna\t/db/GCA_3.1.fna\t0.004\t0\t800/1000",
    ])
    assert list(dist_pairs(db, db, 1, 0.005, exe)) == [
        ("/db/GCA_1.1.fna", "/db/GCA_2.1.fna", 0.001),
        ("/db/GCA_2.1.fna", "/db/GCA_3.1.fna", 0.004),
    ]


def test_dist_pairs_reports_a_failed_command(tmp_path):
    from mashid.mash import dist_pairs
    db = tmp_path / "db.msh"
    db.write_bytes(b"x")
    exe = _fake_mash(tmp_path, [], exit_code=1, stderr="ERROR: could not open db\n")
    with pytest.raises(MashIDError, match="could not open db"):
        list(dist_pairs(db, db, 1, None, exe))


def test_dist_pairs_stays_quiet_when_the_caller_stops_early(tmp_path):
    """Closing the pipe kills Mash with SIGPIPE; that is ours to cause, not a failure to report."""
    from mashid.mash import dist_pairs
    db = tmp_path / "db.msh"
    db.write_bytes(b"x")
    exe = _fake_mash(tmp_path, [f"/db/GCA_{i}.1.fna\t/db/GCA_9.1.fna\t0.001\t0\t1/1"
                                for i in range(2000)])
    gen = dist_pairs(db, db, 1, 0.005, exe)
    assert next(gen)[2] == 0.001
    gen.close()                     # must not raise


# --- compute ------------------------------------------------------------------------------------

def _compute_with(monkeypatch, tmp_path, pairs, metadata):
    from mashid import confusable as mod
    db = tmp_path / "db.msh"
    db.write_bytes(b"x")

    def fake(database, other, threads, max_distance, exe):
        yield from pairs

    monkeypatch.setattr(mod, "dist_pairs", fake)
    return mod.compute(db, metadata)


def test_compute_keeps_the_nearest_reference_of_each_other_organism(monkeypatch, tmp_path):
    meta = {"GCA_1.1": DbEntry("GCA_1.1", "Genusa one"),
            "GCA_2.1": DbEntry("GCA_2.1", "Genusa two"),
            "GCA_3.1": DbEntry("GCA_3.1", "Genusa two"),
            "GCA_4.1": DbEntry("GCA_4.1", "Genusa four")}
    pairs = [("/db/GCA_1.1.fna", "/db/GCA_2.1.fna", 0.004),
             ("/db/GCA_1.1.fna", "/db/GCA_3.1.fna", 0.001),   # same organism, nearer: wins
             ("/db/GCA_1.1.fna", "/db/GCA_4.1.fna", 0.003)]
    out = _compute_with(monkeypatch, tmp_path, pairs, meta)
    assert [(n.organism, n.distance, n.accession) for n in out["GCA_1.1"]] == [
        ("Genusa two", 0.001, "GCA_3.1"), ("Genusa four", 0.003, "GCA_4.1")]
    # the relation is recorded from both sides
    assert [n.organism for n in out["GCA_2.1"]] == ["Genusa one"]
    assert out["GCA_3.1"][0].distance == pytest.approx(0.001)


def test_compute_skips_self_pairs_same_names_and_synonyms(monkeypatch, tmp_path):
    meta = {"GCA_1.1": DbEntry("GCA_1.1", "Genusa one"),
            "GCA_2.1": DbEntry("GCA_2.1", "Genusa one"),      # same organism
            "GCA_3.1": DbEntry("GCA_3.1", "Genusb one"),      # genus synonym
            "GCA_4.1": DbEntry("GCA_4.1", "Genusa sp."),      # unnamed species
            "GCA_5.1": DbEntry("GCA_5.1", "Genusa five")}
    pairs = [("/db/GCA_1.1.fna", "/db/GCA_1.1.fna", 0.0),
             ("/db/GCA_1.1.fna", "/db/GCA_2.1.fna", 0.001),
             ("/db/GCA_1.1.fna", "/db/GCA_3.1.fna", 0.001),
             ("/db/GCA_1.1.fna", "/db/GCA_4.1.fna", 0.001),
             ("/db/GCA_1.1.fna", "/db/GCA_5.1.fna", 0.002)]
    out = _compute_with(monkeypatch, tmp_path, pairs, meta)
    assert list(out) == ["GCA_1.1", "GCA_5.1"]
    assert [n.organism for n in out["GCA_1.1"]] == ["Genusa five"]


def test_compute_warns_about_references_missing_from_the_sidecar(monkeypatch, tmp_path, caplog):
    meta = {"GCA_1.1": DbEntry("GCA_1.1", "Genusa one")}
    pairs = [("/db/GCA_1.1.fna", "/db/GCA_unknown.1.fna", 0.001)]
    with caplog.at_level("WARNING"):
        out = _compute_with(monkeypatch, tmp_path, pairs, meta)
    assert out == {}
    assert "no organism in the metadata sidecar" in caplog.text


def test_compute_needs_a_database_and_metadata(tmp_path):
    from mashid.confusable import compute
    db = tmp_path / "db.msh"
    with pytest.raises(MashIDError, match="not found"):
        compute(db, {"GCA_1.1": DbEntry("GCA_1.1", "Genusa one")})
    db.write_bytes(b"x")
    with pytest.raises(MashIDError, match="No reference metadata"):
        compute(db, {})


def test_confusable_database_writes_the_sidecar(monkeypatch, tmp_path):
    from mashid import confusable as mod
    from mashid import makedb
    db = tmp_path / "db.msh"
    db.write_bytes(b"x")
    (tmp_path / "db.metadata.tsv").write_text(
        "Accession\tOrganism\tTaxID\nGCA_1.1\tGenusa one\t1\nGCA_2.1\tGenusa two\t2\n")
    monkeypatch.setattr(makedb, "find_mash", lambda: "mash")
    monkeypatch.setattr(mod, "compute",
                        lambda *a, **k: {"GCA_1.1": [Neighbour("Genusa two", 0.001, "GCA_2.1")]})
    out = makedb.confusable_database(db, None, 1, 0.005)
    assert out == confusable_path(db)
    assert read(out)["GCA_1.1"][0].organism == "Genusa two"


def test_confusable_database_error_paths(monkeypatch, tmp_path):
    from mashid import confusable as mod
    from mashid import makedb
    monkeypatch.setattr(makedb, "find_mash", lambda: "mash")
    db = tmp_path / "db.msh"

    with pytest.raises(MashIDError, match="Mash database not found"):
        makedb.confusable_database(db, None, 1, 0.005)

    db.write_bytes(b"x")
    with pytest.raises(MashIDError, match="No metadata sidecar"):
        makedb.confusable_database(db, None, 1, 0.005)

    (tmp_path / "db.metadata.tsv").write_text("Accession\tOrganism\tTaxID\nGCA_1.1\tGenusa one\t1\n")
    monkeypatch.setattr(mod, "compute", lambda *a, **k: {})

    def unwritable(*_a, **_k):
        raise OSError("read-only file system")

    monkeypatch.setattr(mod, "write", unwritable)
    with pytest.raises(MashIDError, match="Could not write"):
        makedb.confusable_database(db, None, 1, 0.005)


def test_note_omits_organisms_already_reported():
    table = {"GCA_1.1": [Neighbour("Genusa two", 0.001, "GCA_2.1"),
                         Neighbour("Genusa three", 0.002, "GCA_3.1")]}
    note = note_for("GCA_1.1", "Genusa one", table, exclude={"Genusa two"})
    assert note == NOTE_PREFIX + "Genusa three"
    assert note_for("GCA_1.1", "Genusa one", table,
                    exclude={"Genusa two", "Genusa three"}) is None
