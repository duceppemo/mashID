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


def test_note_omits_organisms_already_reported():
    table = {"GCA_1.1": [Neighbour("Genusa two", 0.001, "GCA_2.1"),
                         Neighbour("Genusa three", 0.002, "GCA_3.1")]}
    note = note_for("GCA_1.1", "Genusa one", table, exclude={"Genusa two"})
    assert note == NOTE_PREFIX + "Genusa three"
    assert note_for("GCA_1.1", "Genusa one", table,
                    exclude={"Genusa two", "Genusa three"}) is None
