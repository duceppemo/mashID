"""Tests for the curated label-correction helpers in scripts/."""

import datetime as dt
import importlib.util
import json
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"

SIDECAR_HEADER = "Accession\tOrganism\tTaxID\tLength\tHashes\tSource_File\tDescription\n"
CORRECTION = ("Accession\tOrganism\tTaxID\tDeposited_Organism\tReason\n"
              "GCA_1.1\tGenusb two\t222\tGenusa one\tswapped in the collection\n")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


binner = _load("bin_by_species")
applier = _load("apply_label_corrections")


# --- reading the corrections table -------------------------------------------------------------

@pytest.mark.parametrize("read", [binner.load_corrections, applier.read_corrections])
def test_extra_columns_are_ignored(tmp_path, read):
    p = tmp_path / "corrections.tsv"
    p.write_text(CORRECTION)
    assert read(p) == {"GCA_1.1": ("Genusb two", "222")}


@pytest.mark.parametrize("read", [binner.load_corrections, applier.read_corrections])
def test_csv_and_flexible_headers(tmp_path, read):
    p = tmp_path / "corrections.csv"
    p.write_text("assembly accession,organism name\nGCA_1.1,Genusb two\n")
    assert read(p) == {"GCA_1.1": ("Genusb two", "")}


@pytest.mark.parametrize("read", [binner.load_corrections, applier.read_corrections])
def test_rows_without_accession_or_organism_are_skipped(tmp_path, read):
    p = tmp_path / "corrections.tsv"
    p.write_text("Accession\tOrganism\n\tGenusb two\nGCA_2.1\t\nGCA_1.1\tGenusb two\n")
    assert read(p) == {"GCA_1.1": ("Genusb two", "")}


@pytest.mark.parametrize("read", [binner.load_corrections, applier.read_corrections])
def test_missing_columns_is_an_error(tmp_path, read):
    p = tmp_path / "corrections.tsv"
    p.write_text("foo\tbar\nx\ty\n")
    with pytest.raises(SystemExit, match="expected columns"):
        read(p)


# --- the two bin-name implementations must not drift apart -------------------------------------

@pytest.mark.parametrize("organism", [
    "Mycobacterium tuberculosis variant bovis",
    "Mycolicibacterium wolinskyi",
    "Mycobacterium sp.",
    "Mycobacterium sp",
    "[Mycobacterium] manitobense",
    "Candidatus Mycobacterium something",
    "Mycobacterium avium subsp. paratuberculosis",
    "Escherichia coli O157:H7",
    "Unnamed",
    "",
])
def test_bin_of_mirrors_species_of(organism):
    assert applier.bin_of(organism) == binner.species_of(organism, "species")


# --- rewriting an already-built database's tables -----------------------------------------------

def _sidecar(tmp_path, organism="Genusa one", taxid="111"):
    p = tmp_path / "db.metadata.tsv"
    p.write_text(SIDECAR_HEADER
                 + f"GCA_1.1\t{organism}\t{taxid}\t4000000\t1000\t/x/GCA_1.1.fna\tcontig one\n"
                 + "GCA_9.9\tGenusc nine\t999\t4000000\t1000\t/x/GCA_9.9.fna\tcontig nine\n")
    return p


def _corrections(tmp_path):
    p = tmp_path / "corrections.tsv"
    p.write_text(CORRECTION)
    return p


def test_sidecar_organism_taxid_and_description(tmp_path):
    p = _sidecar(tmp_path)
    assert applier.apply_to(p, applier.read_corrections(_corrections(tmp_path)), True, False) == 1
    rows = [ln.split("\t") for ln in p.read_text().splitlines()[1:]]
    assert rows[0][1] == "Genusb two" and rows[0][2] == "222"
    assert "deposited as 'Genusa one'" in rows[0][6]
    assert rows[1][1] == "Genusc nine" and rows[1][2] == "999"      # untouched
    stamp = dt.date.today().strftime("%Y%m%d")
    assert (tmp_path / f"db.metadata.tsv.bak-{stamp}").is_file()


def test_bins_table_bin_is_recomputed(tmp_path):
    p = tmp_path / "bins.tsv"
    p.write_text("Accession\tOrganism\tTaxID\tBin\nGCA_1.1\tGenusa one\t111\tGenusa_one\n")
    assert applier.apply_to(p, applier.read_corrections(_corrections(tmp_path)), False, False) == 1
    assert p.read_text().splitlines()[1].split("\t")[3] == "Genusb_two"


def test_idempotent(tmp_path):
    p = _sidecar(tmp_path)
    corrections = applier.read_corrections(_corrections(tmp_path))
    assert applier.apply_to(p, corrections, False, False) == 1
    first = p.read_text()
    assert applier.apply_to(p, corrections, False, False) == 0
    assert p.read_text() == first                                   # no second note appended


def test_dry_run_changes_nothing(tmp_path):
    p = _sidecar(tmp_path)
    before = p.read_text()
    assert applier.apply_to(p, applier.read_corrections(_corrections(tmp_path)), True, True) == 1
    assert p.read_text() == before
    assert not list(tmp_path.glob("*.bak-*"))


def test_existing_backup_is_never_clobbered(tmp_path):
    p = _sidecar(tmp_path)
    corrections = applier.read_corrections(_corrections(tmp_path))
    stamp = dt.date.today().strftime("%Y%m%d")
    older = p.with_name(f"{p.name}.bak-{stamp}")
    older.write_text("precious\n")
    applier.apply_to(p, corrections, True, False)
    assert older.read_text() == "precious\n"
    assert p.with_name(f"{p.name}.bak-{stamp}.2").is_file()


def test_missing_columns_in_target_is_an_error(tmp_path):
    p = tmp_path / "db.metadata.tsv"
    p.write_text("foo\tbar\nx\ty\n")
    with pytest.raises(SystemExit, match="expected 'Accession' and 'Organism'"):
        applier.apply_to(p, {"GCA_1.1": ("Genusb two", "222")}, False, False)


# --- binning honours the corrections -----------------------------------------------------------

def test_binner_uses_corrected_name(tmp_path, monkeypatch):
    data = tmp_path / "data"
    for acc in ("GCA_1.1", "GCA_9.9"):
        (data / acc).mkdir(parents=True)
        (data / acc / f"{acc}_genomic.fna").write_text(">x\nACGT\n")
    report = tmp_path / "report.jsonl"
    report.write_text("".join(
        json.dumps({"accession": acc, "organism": {"organismName": "Genusa one", "taxId": "111"},
                    "assemblyInfo": {"assemblyLevel": "Contig"}}) + "\n"
        for acc in ("GCA_1.1", "GCA_9.9")))
    out = tmp_path / "binned"
    monkeypatch.setattr(sys, "argv", ["bin_by_species.py", str(report), str(data), str(out),
                                      "--corrections", str(_corrections(tmp_path))])
    assert binner.main() == 0
    bins = {ln.split("\t")[0]: ln.split("\t")[1:] for ln in (out / "bins.tsv").read_text().splitlines()[1:]}
    assert bins["GCA_1.1"] == ["Genusb two", "222", "Genusb_two"]
    assert bins["GCA_9.9"] == ["Genusa one", "111", "Genusa_one"]
    assert (out / "Genusb_two" / "GCA_1.1_genomic.fna").exists()
    assert (out / "Genusa_one" / "GCA_9.9_genomic.fna").exists()
