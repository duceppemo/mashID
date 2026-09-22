"""The shipped example dataset must keep producing its documented result."""

import csv
import shutil
from pathlib import Path

import pytest

from mashid.cli import main as mashid_main

EXAMPLE = Path(__file__).resolve().parent.parent / "example"

pytestmark = pytest.mark.skipif(shutil.which("mash") is None, reason="mash not installed")


def _read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_example_matches_expected(tmp_path):
    out = tmp_path / "out"
    sidecar = EXAMPLE / "db" / "example.metadata.tsv"
    sidecar_before = sidecar.read_bytes()
    rc = mashid_main(["-i", str(EXAMPLE / "reads"), "-o", str(out), "-d", str(EXAMPLE / "db" / "example.msh"),
                      "-t", "2", "-p", "2"])
    assert rc == 0
    got = _read_tsv(out / "summary_mashID.tsv")
    expected = _read_tsv(EXAMPLE / "expected" / "summary_mashID.tsv")
    assert got == expected
    for row in expected:
        assert _read_tsv(out / f"{row['Sample']}_mashID.tsv") == \
            _read_tsv(EXAMPLE / "expected" / f"{row['Sample']}_mashID.tsv")
    # the example database ships with its sidecar, so nothing is rewritten inside the repository
    assert sidecar.read_bytes() == sidecar_before


def test_example_illustrates_each_case():
    rows = {r["Sample"]: r for r in _read_tsv(EXAMPLE / "expected" / "summary_mashID.tsv")}
    assert rows["alpha"]["Identification"] == "Exemplaria alpha" and rows["alpha"]["Note"] == ""
    assert rows["delta_assembly"]["Identification"] == "Fictivibrio delta" and rows["delta_assembly"]["Sequences"] == "2"
    assert rows["mixed_beta_delta"]["Note"].startswith("Possible mixture with: ")
    assert rows["unknown"]["Identification"] == "No significant hit in database"
