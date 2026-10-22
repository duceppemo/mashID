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
    assert applier.apply_to(p, applier.read_corrections(_corrections(tmp_path)), True, False) == (1, {"GCA_1.1"})
    rows = [ln.split("\t") for ln in p.read_text().splitlines()[1:]]
    assert rows[0][1] == "Genusb two" and rows[0][2] == "222"
    assert "deposited as 'Genusa one'" in rows[0][6]
    assert rows[1][1] == "Genusc nine" and rows[1][2] == "999"      # untouched
    stamp = dt.date.today().strftime("%Y%m%d")
    assert (tmp_path / f"db.metadata.tsv.bak-{stamp}").is_file()


def test_bins_table_bin_is_recomputed(tmp_path):
    p = tmp_path / "bins.tsv"
    p.write_text("Accession\tOrganism\tTaxID\tBin\nGCA_1.1\tGenusa one\t111\tGenusa_one\n")
    assert applier.apply_to(p, applier.read_corrections(_corrections(tmp_path)), False, False) == (1, {"GCA_1.1"})
    assert p.read_text().splitlines()[1].split("\t")[3] == "Genusb_two"


def test_idempotent(tmp_path):
    p = _sidecar(tmp_path)
    corrections = applier.read_corrections(_corrections(tmp_path))
    assert applier.apply_to(p, corrections, False, False)[0] == 1
    first = p.read_text()
    assert applier.apply_to(p, corrections, False, False)[0] == 0
    assert p.read_text() == first                                   # no second note appended


def test_dry_run_changes_nothing(tmp_path):
    p = _sidecar(tmp_path)
    before = p.read_text()
    assert applier.apply_to(p, applier.read_corrections(_corrections(tmp_path)), True, True) == (1, {"GCA_1.1"})
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


def test_a_ragged_row_is_filled_out_rather_than_skipped(tmp_path):
    """It used to raise IndexError. Skipping such a row instead would report the correction as
    unused and, where the row reached Organism but not TaxID, leave a new name on a stale TaxID."""
    p = tmp_path / "bins.tsv"
    p.write_text("Accession\tOrganism\tTaxID\tBin\n"
                 "GCA_1.1\tGenusa one\n"                 # short of TaxID and Bin
                 "GCA_9.9\tGenusc nine\n")               # ragged and uncorrected: left alone
    changed, matched = applier.apply_to(p, applier.read_corrections(_corrections(tmp_path)),
                                        False, False)
    assert changed == 1 and matched == {"GCA_1.1"}
    rows = [ln.split("\t") for ln in p.read_text().splitlines()[1:]]
    assert rows[0] == ["GCA_1.1", "Genusb two", "222", "Genusb_two"]
    assert rows[1] == ["GCA_9.9", "Genusc nine"]


def test_a_correction_matching_a_ragged_row_is_not_called_unused(tmp_path, capsys):
    p = tmp_path / "db.metadata.tsv"
    p.write_text(SIDECAR_HEADER + "GCA_1.1\tGenusa one\n")
    applier.main([str(_corrections(tmp_path)), str(p), "--no-backup"])
    assert "unused correction" not in capsys.readouterr().err


def _serovar_bins(tmp_path, name):
    p = tmp_path / name
    p.write_text("Accession\tOrganism\tTaxID\tBin\n"
                 "GCA_1.1\tSalmonella enterica subsp. enterica serovar Dublin\t28901"
                 "\tSalmonella_enterica_subsp_enterica_serovar_Dublin\n")
    return p


def _serovar_correction(tmp_path):
    p = tmp_path / "serovar_corrections.tsv"
    p.write_text("Accession\tOrganism\n"
                 "GCA_1.1\tSalmonella enterica subsp. enterica serovar Typhimurium\n")
    return p


def test_bin_is_recomputed_at_the_rank_the_bins_were_built_at(tmp_path):
    p = _serovar_bins(tmp_path, "bins.tsv")
    applier.apply_to(p, applier.read_corrections(_serovar_correction(tmp_path)), False, False,
                     "serovar")
    assert p.read_text().splitlines()[1].split("\t")[3] == (
        "Salmonella_enterica_subsp_enterica_serovar_Typhimurium")


def test_the_wrong_rank_is_an_error_rather_than_a_flattened_bin(tmp_path):
    """The old bin says which rank the table was binned at, so getting --rank wrong is detectable;
    silently writing a species-level bin among serovar-level ones misdirects the next rebuild."""
    p = _serovar_bins(tmp_path, "bins.tsv")
    before = p.read_text()
    with pytest.raises(SystemExit, match="Pass the --rank the table was binned at"):
        applier.apply_to(p, applier.read_corrections(_serovar_correction(tmp_path)), False, False)
    assert p.read_text() == before


def test_a_confusable_table_has_both_organism_columns_corrected(tmp_path):
    """Confusable_Organism is the column the note is built from, so leaving it stale would keep
    reporting the deposited name."""
    p = tmp_path / "db.confusable.tsv"
    p.write_text("Accession\tOrganism\tConfusable_Organism\tDistance\tNearest_Accession\n"
                 "GCA_1.1\tGenusa one\tGenusc nine\t0.001\tGCA_9.9\n"
                 "GCA_9.9\tGenusc nine\tGenusa one\t0.001\tGCA_1.1\n")
    changed, matched = applier.apply_to(p, applier.read_corrections(_corrections(tmp_path)),
                                        False, False)
    assert changed == 2 and matched == {"GCA_1.1"}
    rows = [ln.split("\t") for ln in p.read_text().splitlines()[1:]]
    assert rows[0][1] == "Genusb two"       # its own Organism
    assert rows[1][2] == "Genusb two"       # and where another row names it


def test_a_crlf_table_is_rewritten_without_a_stray_carriage_return(tmp_path):
    p = tmp_path / "db.metadata.tsv"
    p.write_bytes(SIDECAR_HEADER.replace("\n", "\r\n").encode()
                  + b"GCA_1.1\tGenusa one\t111\t4000000\t1000\t/x/a.fna\tcontig one\r\n")
    applier.apply_to(p, applier.read_corrections(_corrections(tmp_path)), False, False)
    text = p.read_text()
    assert "\r" not in text
    assert text.splitlines()[1].split("\t")[6].endswith("]")      # the note ends the field


def test_a_missing_corrections_file_exits_cleanly(tmp_path):
    with pytest.raises(SystemExit, match="Corrections file not found"):
        applier.read_corrections(tmp_path / "nope.tsv")
    with pytest.raises(SystemExit, match="Corrections file not found"):
        binner.load_corrections(tmp_path / "nope.tsv")


def test_a_correction_that_matches_nothing_is_reported(tmp_path, capsys):
    p = tmp_path / "db.metadata.tsv"
    p.write_text(SIDECAR_HEADER + "GCA_7.7\tGenusc nine\t999\t4000000\t1000\t/x/c.fna\td\n")
    assert applier.main([str(_corrections(tmp_path)), str(p), "--dry-run"]) == 0
    assert "unused correction: GCA_1.1" in capsys.readouterr().err


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
