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
    rebuilt = script.rebuild(out)
    results, rows, columns, changed = (rebuilt.results, rebuilt.rows,
                                       rebuilt.columns, rebuilt.changed)
    assert len(rows) == 1 and changed == ["s1"]
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
    rows = script.rebuild(out).rows
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
    out, db = _results_dir(tmp_path)
    before = (out / "summary_mashID.tsv").read_text()
    script.main([str(out)])
    backups = sorted(p.name for p in out.glob("summary_mashID.tsv.bak-*"))
    assert len(backups) == 1
    assert (out / backups[0]).read_text() == before
    # Relabel again so the next run has something to do, and check it keeps the first backup.
    write_sidecar(sidecar_path(db), [DbEntry("GCA_1.1", "Genusa three", "1765", 4_000_000, 10000),
                                     DbEntry("GCA_2.1", "Genusa two", "222", 4_000_000, 10000)])
    script.main([str(out)])
    assert len(list(out.glob("summary_mashID.tsv.bak-*"))) == 2
    assert (out / backups[0]).read_text() == before


def test_a_run_already_up_to_date_is_left_alone(tmp_path, capsys):
    """A second identical run used to stack up another set of backups and stamp the JSON with a
    re-annotation that changed nothing."""
    out, _ = _results_dir(tmp_path)
    assert script.main([str(out)]) == 0
    after = (out / "summary_mashID.tsv").read_text()
    stamped = json.loads((out / "summary_mashID.json").read_text())["reannotated"]["at"]
    capsys.readouterr()
    assert script.main([str(out)]) == 0
    assert "nothing written" in capsys.readouterr().out
    assert (out / "summary_mashID.tsv").read_text() == after
    assert len(list(out.glob("summary_mashID.tsv.bak-*"))) == 1
    assert json.loads((out / "summary_mashID.json").read_text())["reannotated"]["at"] == stamped


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
    rows = script.rebuild(out).rows
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
    rebuilt = script.rebuild(out, moved)
    rows, used = rebuilt.rows, rebuilt.database
    assert used == moved and rows[0]["Identification"] == "Genusa one"


def test_uses_the_sidecar_tables_the_run_was_made_with(tmp_path, capsys):
    """A run made with --db-metadata was annotated from that table. Re-annotating it from the
    sidecars beside the .msh would undo exactly what the override was for, and call it a change."""
    out, db = _results_dir(tmp_path)
    curated = tmp_path / "curated.metadata.tsv"
    write_sidecar(curated, [DbEntry("GCA_1.1", "Genusa curated", "9999", 4_000_000, 10000),
                            DbEntry("GCA_2.1", "Genusa two", "222", 4_000_000, 10000)])
    _rewrite_run_json(out, parameters={"ambiguity_margin": 0.005, "winner_take_all": True,
                                       "db_metadata": str(curated)})
    rebuilt = script.rebuild(out)
    assert rebuilt.metadata_table == curated
    assert rebuilt.rows[0]["Identification"] == "Genusa curated"
    assert rebuilt.rows[0]["TaxID"] == "9999"
    assert "using the metadata table the run was made with" in capsys.readouterr().err
    # --database is the way to say "re-annotate from this database's own sidecars instead"
    assert script.rebuild(out, db).rows[0]["Identification"] == "Genusa one"


def test_a_recorded_override_that_has_gone_missing_is_an_error(tmp_path):
    out, _ = _results_dir(tmp_path)
    _rewrite_run_json(out, parameters={"db_confusable": str(tmp_path / "gone.tsv")})
    with pytest.raises(MashIDError, match="which no longer exists"):
        script.rebuild(out)


def test_the_json_records_which_tables_were_used_and_what_moved(tmp_path):
    out, db = _results_dir(tmp_path)
    assert script.main([str(out), "--no-backup"]) == 0
    record = json.loads((out / "summary_mashID.json").read_text())["reannotated"]
    assert record["samples_changed"] == ["s1"]
    assert record["metadata_sidecar"]["path"] == str(sidecar_path(db))
    assert record["confusable_sidecar"]["path"] == str(confusable_path(db))
    assert "predates" not in record["reason"]          # it did not necessarily predate anything


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


def _rewrite_run_json(out, **changes):
    path = out / "mashID_run.json"
    doc = json.loads(path.read_text())
    doc.update(changes)
    path.write_text(json.dumps(doc) + "\n")


def test_refuses_to_write_when_the_run_lists_no_samples(tmp_path):
    """The destructive case: rebuilding from an empty sample list truncated the summary to its
    header and exited 0, which with --no-backup cannot be undone."""
    out, _ = _results_dir(tmp_path)
    before = (out / "summary_mashID.tsv").read_text()
    for samples in ([], None):
        _rewrite_run_json(out, samples=samples)
        with pytest.raises(MashIDError, match="would drop 1 of them"):
            script.rebuild(out)
        assert (out / "summary_mashID.tsv").read_text() == before


def test_refuses_to_drop_a_row_the_run_no_longer_lists(tmp_path):
    out, _ = _results_dir(tmp_path)
    (out / "summary_mashID.tsv").write_text(
        (out / "summary_mashID.tsv").read_text()
        + "s2\t1\t1\t0.9\t1/10000\t1\t0.1\t0\tGCA_2.1\t222\tGenusa two\t\n")
    with pytest.raises(MashIDError, match=r"would drop 1 of them \('s2'\)"):
        script.rebuild(out)


def test_warns_about_a_sample_the_summary_does_not_have(tmp_path, capsys):
    out, _ = _results_dir(tmp_path)
    (out / "s2_mashID.tsv").write_text(PER_SAMPLE)
    doc = json.loads((out / "mashID_run.json").read_text())
    doc["samples"].append({"name": "s2", "type": "fastq", "files": []})
    (out / "mashID_run.json").write_text(json.dumps(doc) + "\n")
    rows = script.rebuild(out).rows
    assert [r["Sample"] for r in rows] == ["s1", "s2"]
    row = next(r for r in rows if r["Sample"] == "s2")
    assert row["Sequences"] == "NA" and row["Bases"] == "NA"
    assert "s2 is in mashID_run.json but not in summary_mashID.tsv" in capsys.readouterr().err


@pytest.mark.parametrize("doc,message", [
    ("not json at all", "is not readable JSON"),
    ('["a list"]', "does not hold a run record"),
])
def test_rejects_an_unreadable_run_record(tmp_path, doc, message):
    out, _ = _results_dir(tmp_path)
    (out / "mashID_run.json").write_text(doc)
    with pytest.raises(MashIDError, match=message):
        script.rebuild(out)


def test_rejects_a_run_record_without_a_database(tmp_path):
    out, _ = _results_dir(tmp_path)
    _rewrite_run_json(out, database=None)
    with pytest.raises(MashIDError, match="records no database path"):
        script.rebuild(out)


def test_rejects_a_sample_entry_without_a_name(tmp_path):
    out, _ = _results_dir(tmp_path)
    (out / "summary_mashID.tsv").write_text(SUMMARY.splitlines()[0] + "\n")   # header only
    _rewrite_run_json(out, samples=[{"type": "fastq"}])
    with pytest.raises(MashIDError, match="has no 'name'"):
        script.rebuild(out)


def test_rejects_a_summary_without_a_sample_column(tmp_path):
    out, _ = _results_dir(tmp_path, summary="Foo\tBar\nx\ty\n")
    with pytest.raises(MashIDError, match="no Sample column"):
        script.rebuild(out)


def test_a_column_the_old_summary_lacks_is_added_and_counted(tmp_path, capsys):
    """The headline use case: a run that predates make_mashID_db --confusable has no Note column.
    Taking the column set from the old file dropped Note from the rewrite, and treating an absent
    column as agreement reported the run as already up to date."""
    header, row = SUMMARY.splitlines()
    no_note = "\t".join(header.split("\t")[:-1]) + "\n" + "\t".join(row.split("\t")[:-1]) + "\n"
    out, _ = _results_dir(tmp_path, summary=no_note)
    rebuilt = script.rebuild(out)
    assert "Note" in rebuilt.columns and rebuilt.changed == ["s1"]
    assert script.main([str(out), "--no-backup"]) == 0
    written = (out / "summary_mashID.tsv").read_text().splitlines()
    assert written[0].split("\t")[-1] == "Note"
    assert written[1].endswith(NOTE_PREFIX + "Genusa two")
    # and the two files that were rewritten together agree with each other
    assert NOTE_PREFIX in (out / "summary_mashID_mqc.tsv").read_text()


def test_sample_sheet_columns_are_kept_but_mashid_columns_are_not_taken_from_the_old_file(tmp_path):
    out, _ = _results_dir(tmp_path, summary=SUMMARY.replace("\tNote\n", "\tNote\tExpected\n")
                                                 .replace("\tOld name\t\n", "\tOld name\t\tMbovis\n"))
    rebuilt = script.rebuild(out)
    rows, columns = rebuilt.rows, rebuilt.columns
    assert columns[-1] == "Expected" and rows[0]["Expected"] == "Mbovis"


def test_refuses_a_summary_with_two_rows_for_one_sample(tmp_path):
    """Keyed by name, the second row silently won and the first vanished."""
    out, _ = _results_dir(tmp_path, summary=SUMMARY + SUMMARY.splitlines()[1] + "\n")
    with pytest.raises(MashIDError, match="more than one row for 's1'"):
        script.rebuild(out)


def test_dry_run_reports_an_identification_change_not_only_a_note(tmp_path, capsys):
    """The run's summary says 'Old name'; the sidecar now says 'Genusa one' with a new TaxID."""
    out, _ = _results_dir(tmp_path, confusable=False)
    assert script.main([str(out), "--dry-run"]) == 0
    printed = capsys.readouterr().out
    assert "1 with a changed Identification, TaxID, Note" in printed
    assert "s1: Genusa one (1765)" in printed


def test_an_unchanged_run_reports_nothing_to_do(tmp_path, capsys):
    out, _ = _results_dir(tmp_path)
    script.main([str(out)])                            # bring it up to date
    capsys.readouterr()
    assert script.main([str(out), "--dry-run"]) == 0
    assert "0 with a changed" in capsys.readouterr().out


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
