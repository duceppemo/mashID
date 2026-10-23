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
    # the RefSeq accession of the same assembly is registered alongside the row's own
    assert read(p) == {"GCA_1.1": ("Genusb two", "222"), "GCF_1.1": ("Genusb two", "222")}


@pytest.mark.parametrize("read", [binner.load_corrections, applier.read_corrections])
def test_csv_and_flexible_headers(tmp_path, read):
    p = tmp_path / "corrections.csv"
    p.write_text("assembly accession,organism name\nGCA_1.1,Genusb two\n")
    assert read(p) == {"GCA_1.1": ("Genusb two", ""), "GCF_1.1": ("Genusb two", "")}


@pytest.mark.parametrize("read", [binner.load_corrections, applier.read_corrections])
def test_rows_without_accession_or_organism_are_skipped(tmp_path, read):
    p = tmp_path / "corrections.tsv"
    p.write_text("Accession\tOrganism\n\tGenusb two\nGCA_2.1\t\nGCA_1.1\tGenusb two\n")
    assert read(p) == {"GCA_1.1": ("Genusb two", ""), "GCF_1.1": ("Genusb two", "")}


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
    assert "matched no row" not in capsys.readouterr().err


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
    err = capsys.readouterr().err
    assert "1 of 1 correction(s) matched no row" in err and "GCA_1.1" in err


def test_a_genbank_row_corrects_the_refseq_accession_of_the_same_assembly(tmp_path, capsys):
    """One assembly, two accessions: a database built from RefSeq carries GCF_1.1 where the
    curation file, written from the GenBank record, says GCA_1.1. The correction still applies,
    and the pair counts as one correction rather than two."""
    p = tmp_path / "db.metadata.tsv"
    p.write_text(SIDECAR_HEADER + "GCF_1.1\tGenusa one\t111\t4000000\t1000\t/x/GCF_1.1.fna\tcontig one\n")
    assert applier.main([str(_corrections(tmp_path)), str(p)]) == 0
    row = p.read_text().splitlines()[1].split("\t")
    assert row[1] == "Genusb two" and row[2] == "222"
    out = capsys.readouterr()
    assert "1 correction(s) loaded" in out.out and "matched no row" not in out.err


def test_the_version_is_not_varied_when_pairing_accessions(tmp_path, capsys):
    """GCF_1.2 is a different assembly from GCA_1.1, so nothing is corrected."""
    p = tmp_path / "db.metadata.tsv"
    p.write_text(SIDECAR_HEADER + "GCF_1.2\tGenusa one\t111\t4000000\t1000\t/x/GCF_1.2.fna\tcontig one\n")
    assert applier.main([str(_corrections(tmp_path)), str(p), "--dry-run"]) == 0
    assert "1 of 1 correction(s) matched no row" in capsys.readouterr().err


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


def test_the_metadata_build_path_pairs_accessions_like_the_corrections_reader(tmp_path):
    """build_taxon_db.sh bins with --corrections and then names with --metadata, from one file.

    Both readers have to agree about GenBank/RefSeq pairing. When only the binning side paired
    them, a row written for GCA_1.1 binned GCF_1.1 under its corrected species and then let the
    assembly report's deposited name into the sidecar, with nothing logged: the report registers
    both accessions, so the reference did not count as unannotated either.
    """
    import json

    from mashid.makedb import load_annotations
    from mashid.metadata import entries_from_sketches

    # Accessions have to be well formed here: the sidecar writer recovers one from the sketch path.
    corrections = tmp_path / "corrections.tsv"
    corrections.write_text("Accession\tOrganism\tTaxID\n"
                           "GCA_000000001.1\tGenusb two\t222\n")
    report = tmp_path / "assembly_data_report.jsonl"
    report.write_text(json.dumps({
        "accession": "GCF_000000001.1", "pairedAccession": "GCA_000000001.1",
        "organism": {"organismName": "Genusa one", "taxId": 111},
    }) + "\n")
    annotations = load_annotations(corrections, report)  # the row is written for GCA_000000001.1
    assert annotations["GCF_000000001.1"] == ("Genusb two", "222"), "the curated name must win for both"
    assert annotations["GCA_000000001.1"] == ("Genusb two", "222")

    from mashid.mash import SketchInfo
    entries, unannotated = entries_from_sketches(
        [SketchInfo(1000, 4_000_000, "/derep/Genusb_two/GCF_000000001.1_genomic.fna.gz",
                    "NZ_1 Genusa one strain x")],
        annotations)
    assert (entries[0].organism, entries[0].taxid) == ("Genusb two", "222") and unannotated == 0


def test_binner_counts_corrections_per_assembly_not_per_key(tmp_path, capsys, monkeypatch):
    """A paired mapping holds two keys per assembly; counting keys announced twice as many
    corrections as the file has and called every applied row's twin unmatched."""
    import json

    data = tmp_path / "data" / "GCF_1.1"
    data.mkdir(parents=True)
    (data / "GCF_1.1_genomic.fna").write_text(">x\nACGT\n")
    report = tmp_path / "report.jsonl"
    report.write_text(json.dumps({
        "accession": "GCF_1.1", "organism": {"organismName": "Genusa one", "taxId": 111},
        "assemblyInfo": {"assemblyLevel": "Complete Genome"},
    }) + "\n")
    monkeypatch.setattr("sys.argv", ["bin_by_species.py", str(report), str(tmp_path / "data"),
                                     str(tmp_path / "out"), "--corrections", str(_corrections(tmp_path))])
    assert binner.main() == 0
    err = capsys.readouterr().err
    assert "1 of 1 label correction(s) applied" in err
    assert "match no accession" not in err


def test_bins_table_appears_only_once_it_is_complete(tmp_path, monkeypatch):
    """build_taxon_db.sh skips binning when bins.tsv exists, so a half-written one must not be
    left behind under that name."""
    import json

    data = tmp_path / "data" / "GCF_1.1"
    data.mkdir(parents=True)
    (data / "GCF_1.1_genomic.fna").write_text(">x\nACGT\n")
    report = tmp_path / "report.jsonl"
    report.write_text(json.dumps({"accession": "GCF_1.1",
                                  "organism": {"organismName": "Genusa one", "taxId": 111}}) + "\n")
    out = tmp_path / "out"
    real_replace = __import__("os").replace
    seen = {}

    def spy(src, dst):
        seen["partial_existed"] = Path(src).is_file()
        seen["final_existed_before"] = Path(dst).exists()
        return real_replace(src, dst)

    monkeypatch.setattr(binner.os, "replace", spy)
    monkeypatch.setattr("sys.argv", ["bin_by_species.py", str(report), str(tmp_path / "data"), str(out)])
    assert binner.main() == 0
    assert seen == {"partial_existed": True, "final_existed_before": False}
    assert (out / "bins.tsv").is_file() and not (out / "bins.tsv.partial").exists()
