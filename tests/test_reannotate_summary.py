"""Tests for scripts/reannotate_summary.py."""

import importlib.util
import json
from pathlib import Path

import pytest

from mashid import MashIDError
from mashid.confusable import NOTE_PREFIX, Neighbour, confusable_path
from mashid.confusable import write as write_confusable
from mashid.metadata import DbEntry, sidecar_path
from mashid.metadata import write_metadata as write_sidecar

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def _load():
    spec = importlib.util.spec_from_file_location("reannotate_summary", SCRIPTS / "reannotate_summary.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


script = _load()

PER_SAMPLE = (
    "Rank\tIdentity\tShared_Hashes\tMedian_Multiplicity\tP_Value\tAccession\tTaxID\tIdentification\tDescription\n"
    "1\t0.9999\t9900/10000\t50\t0\tGCA_1.1\t111\tOld name\tNZ_1 Old name strain x\n"
    "2\t0.9000\t1100/10000\t40\t0\tGCA_2.1\t222\tGenusa two\tNZ_2 Genusa two strain y\n"
)
SUMMARY = (
    "Sample\tSequences\tBases\tIdentity\tShared_Hashes\tMedian_Multiplicity\tEst_Depth\tP_Value"
    "\tAccession\tTaxID\tIdentification\tNote\n"
    "s1\t1000\t400000\t0.9999\t9900/10000\t50\t0.1\t0\tGCA_1.1\t111\tOld name\t\n"
)


def _results_dir(tmp_path, *, confusable=True, summary=SUMMARY, per_sample=PER_SAMPLE):
    """A finished results directory whose database sidecars have since been corrected."""
    db = tmp_path / "db.msh"
    db.write_bytes(b"not really a sketch")
    write_sidecar(sidecar_path(db), [
        DbEntry("GCA_1.1", "Genusa one", "1765", 4_000_000, 10000, "/x/GCA_1.1.fna", "d"),
        DbEntry("GCA_2.1", "Genusa two", "222", 4_000_000, 10000, "/x/GCA_2.1.fna", "d"),
    ])
    if confusable:
        write_confusable(confusable_path(db), {"GCA_1.1": [Neighbour("Genusa two", 0.001, "GCA_2.1")]}, {})

    out = tmp_path / "out"
    out.mkdir()
    (out / "s1_mashID.tsv").write_text(per_sample)
    (out / "summary_mashID.tsv").write_text(summary)
    (out / "summary_mashID_mqc.tsv").write_text("# id: mashid\nSample\tNote\ns1\t\n")
    (out / "summary_mashID.json").write_text(json.dumps({"samples": []}) + "\n")
    (out / "mashID_run.json").write_text(json.dumps({
        "database": {"path": str(db)},
        "parameters": {"ambiguity_margin": 0.005, "winner_take_all": True},
        "samples": [{"name": "s1", "type": "fastq", "files": ["/r/s1_R1.fastq.gz"]}],
    }) + "\n")
    return out, db


def test_rebuild_picks_up_the_corrected_label_and_the_new_note(tmp_path):
    out, _ = _results_dir(tmp_path)
    results, rows, columns, database, changed = script.rebuild(out)
    assert len(rows) == 1 and changed == 1
    row = rows[0]
    assert row["Identification"] == "Genusa one"       # relabelled since the run
    assert row["TaxID"] == "1765"
    assert row["Note"] == NOTE_PREFIX + "Genusa two"
    # measurements are carried through untouched
    assert row["Identity"] == "0.9999" and row["Shared_Hashes"] == "9900/10000"
    assert row["Sequences"] == "1000" and row["Bases"] == "400000"
    assert columns[0] == "Sample"
    assert results[0].sample.seq_type == "fastq"


def test_est_depth_is_recomputed_from_the_sidecar_length(tmp_path):
    out, _ = _results_dir(tmp_path)
    _, rows, _, _, _ = script.rebuild(out)
    assert rows[0]["Est_Depth"] == "0.1"               # 400000 bases / 4 Mb reference


def test_writes_all_three_summary_files_and_records_provenance(tmp_path):
    out, db = _results_dir(tmp_path)
    assert script.main([str(out)]) == 0
    assert "Genusa one" in (out / "summary_mashID.tsv").read_text()
    assert NOTE_PREFIX in (out / "summary_mashID_mqc.tsv").read_text()
    doc = json.loads((out / "summary_mashID.json").read_text())
    assert doc["samples"][0]["top_hit"]["Identification"] == "Genusa one"
    assert doc["samples"][0]["notes"] == [NOTE_PREFIX + "Genusa two"]
    record = doc["reannotated"]
    assert record["screening_unchanged"] is True
    assert record["metadata_sidecar"]["path"] == str(sidecar_path(db))
    assert len(record["metadata_sidecar"]["md5"]) == 32
    assert record["confusable_sidecar"]["path"] == str(confusable_path(db))


def test_backs_up_what_it_overwrites_without_clobbering(tmp_path):
    out, _ = _results_dir(tmp_path)
    before = (out / "summary_mashID.tsv").read_text()
    script.main([str(out)])
    backups = sorted(p.name for p in out.glob("summary_mashID.tsv.bak-*"))
    assert len(backups) == 1
    assert (out / backups[0]).read_text() == before
    script.main([str(out)])                            # second run keeps the first backup
    assert len(list(out.glob("summary_mashID.tsv.bak-*"))) == 2
    assert (out / backups[0]).read_text() == before


def test_no_backup_and_dry_run(tmp_path):
    out, _ = _results_dir(tmp_path)
    before = (out / "summary_mashID.tsv").read_text()
    assert script.main([str(out), "--dry-run"]) == 0
    assert (out / "summary_mashID.tsv").read_text() == before
    assert not list(out.glob("*.bak-*"))
    assert script.main([str(out), "--no-backup"]) == 0
    assert not list(out.glob("*.bak-*"))
    assert (out / "summary_mashID.tsv").read_text() != before


def test_works_without_a_confusable_sidecar(tmp_path, capsys):
    out, _ = _results_dir(tmp_path, confusable=False)
    _, rows, _, _, _ = script.rebuild(out)
    assert rows[0]["Identification"] == "Genusa one"   # the relabel still lands
    assert rows[0]["Note"] == ""                       # but no database-wide note
    assert "make_mashID_db --confusable" in capsys.readouterr().err


def test_database_override(tmp_path):
    out, db = _results_dir(tmp_path)
    moved = tmp_path / "moved.msh"
    db.rename(moved)
    sidecar_path(db).rename(sidecar_path(moved))
    confusable_path(db).rename(confusable_path(moved))
    with pytest.raises(MashIDError, match="Database not found"):
        script.rebuild(out)
    _, rows, _, used, _ = script.rebuild(out, moved)
    assert used == moved and rows[0]["Identification"] == "Genusa one"


def test_rejects_a_directory_that_is_not_a_results_directory(tmp_path):
    (tmp_path / "empty").mkdir()
    with pytest.raises(MashIDError, match="Not a mashID results directory"):
        script.rebuild(tmp_path / "empty")
    out, _ = _results_dir(tmp_path)
    (out / "s1_mashID.tsv").unlink()
    with pytest.raises(MashIDError, match="Per-sample table missing for s1"):
        script.rebuild(out)


def test_rejects_a_per_sample_table_missing_columns(tmp_path):
    out, _ = _results_dir(tmp_path, per_sample="Rank\tIdentity\nn1\t0.9\n")
    with pytest.raises(MashIDError, match="missing column"):
        script.rebuild(out)


def test_main_reports_an_error_without_a_traceback(tmp_path, capsys):
    (tmp_path / "empty").mkdir()
    assert script.main([str(tmp_path / "empty")]) == 1
    assert "Not a mashID results directory" in capsys.readouterr().err


def test_read_hits_preserves_order_and_skips_blank_accessions(tmp_path):
    p = tmp_path / "s_mashID.tsv"
    p.write_text(PER_SAMPLE + "3\t0.8\t1/10000\t1\t0\t\t\t\t\n")
    hits = script.read_hits(p)
    assert [h.identity for h in hits] == ["0.9999", "0.9000"]
    assert hits[0].query_id == "GCA_1.1.fna"
