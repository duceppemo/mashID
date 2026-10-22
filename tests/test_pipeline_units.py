from mashid.mash import ScreenHit, parse_screen_output
from mashid.metadata import accession_from_query_id
from mashid.pipeline import format_table, sort_hits

RAW = (
    "0.99914\t9821/10000\t1\t0\t/db/GCF_000463275.1.fna\t[1092 seqs] NZ_AKYV01000001.1 Mycobacterium tuberculosis variant bovis BCG str. Sweden Contig_1 [...]\n"
    "0.905133\t1233/10000\t5\t1e-10\t/db/GCF_002703865.1.fna.gz\tNZ_NAZK01000001.1 Mycobacterium bovis strain X\n"
    "\n"
)


def test_parse_and_sort():
    hits = parse_screen_output(RAW)
    assert len(hits) == 2
    assert hits[0].shared_f == 9821 and hits[1].p_value_f == 1e-10
    assert [h.identity for h in sort_hits(hits[::-1], "identity")] == ["0.99914", "0.905133"]
    assert [h.identity for h in sort_hits(hits, "multiplicity")] == ["0.905133", "0.99914"]


def test_parse_tolerates_missing_comment():
    hits = parse_screen_output("0.9\t1/10\t1\t0\tq\n")
    assert hits[0].comment == ""


def test_accession_from_query_id():
    assert accession_from_query_id("/db/GCF_000463275.1.fna") == "GCF_000463275.1"
    assert accession_from_query_id("/db/GCF_000463275.1.fna.gz") == "GCF_000463275.1"
    assert accession_from_query_id("reads_sketch") == "reads_sketch"


def test_format_table():
    out = format_table(["A", "Bee"], [{"A": "1", "Bee": "x"}, {"A": "22", "Bee": ""}])
    assert out.splitlines() == ["A   Bee", "1   x", "22"]


def test_accession_from_ncbi_datasets_filename():
    assert accession_from_query_id("/x/GCF_000195955.2_ASM19595v2_genomic.fna.gz") == "GCF_000195955.2"
    assert accession_from_query_id("GCA_000000001.1.fa") == "GCA_000000001.1"


def _hit(identity, mult, acc, comment):
    return ScreenHit(str(identity), "9000/10000", str(mult), "0", f"/db/{acc}.fna", comment)


def _result(hits, seq_type="fastq"):
    from mashid.pipeline import SampleResult
    from mashid.samples import Sample
    return SampleResult(sample=Sample(name="s", files=[], seq_type=seq_type), hits=hits)


def test_annotate_ambiguous_mixture_low_coverage():
    from mashid.pipeline import annotate
    top = _hit(0.999, 50, "A1", "NZ_1 Genusa one strain x")
    same_org = _hit(0.998, 50, "A2", "NZ_2 Genusa one strain y")
    close_other = _hit(0.996, 50, "B1", "NZ_3 Genusb two strain z")
    far_other = _hit(0.975, 40, "C1", "NZ_4 Genusc three strain w")
    hits = [top, same_org, close_other, far_other]
    # winner-take-all: another organism at >= 0.99 identity is a mixture, not an ambiguity;
    # Genusc at 0.975 is what a sister strain looks like after winner-take-all and is not reported
    assert annotate(_result(hits), None, 0.005, winner_take_all=True) == \
        ["Possible mixture with: Genusb two"]
    # a much shorter reference is flagged as possibly partial
    from mashid.metadata import DbEntry
    meta = {"A1": DbEntry("A1", "Genusa one", length=4_000_000), "B1": DbEntry("B1", "Genusb two", length=900_000)}
    assert annotate(_result(hits), meta, 0.005, winner_take_all=True) == \
        ["Possible mixture with: Genusb two (reference only 900000 bp, may be partial)"]
    # without winner-take-all only the ambiguity check applies
    assert annotate(_result(hits), None, 0.005, winner_take_all=False) == \
        ["Ambiguous: Genusb two at identity 0.996"]
    # ambiguity below the mixture threshold
    low = [_hit(0.96, 50, "A1", "NZ_1 Genusa one"), _hit(0.958, 50, "B1", "NZ_3 Genusb two")]
    assert annotate(_result(low), None, 0.005) == ["Ambiguous: Genusb two at identity 0.958"]
    # a clearly separated second organism gives no note
    assert annotate(_result([top, _hit(0.95, 50, "B1", "NZ_3 Genusb two")]), None, 0.005) == []
    # low coverage only applies to reads
    assert annotate(_result([_hit(0.99, 2, "A1", "NZ_1 Genusa one")]), None, 0.005) == \
        ["Low coverage: median multiplicity 2"]
    assert annotate(_result([_hit(0.99, 1, "A1", "NZ_1 Genusa one")], seq_type="fasta"), None, 0.005) == []
    assert annotate(_result([]), None, 0.005) == []


def test_annotate_reports_inseparable_organisms_from_the_database():
    from mashid.confusable import NOTE_PREFIX, Neighbour
    from mashid.pipeline import annotate
    top = _hit(0.999, 50, "A1", "NZ_1 Genusa one strain x")
    sister = _hit(0.90, 50, "B1", "NZ_3 Genusb two strain z")  # crushed by winner-take-all
    table = {"A1": [Neighbour("Genusb two", 0.001, "B1"), Neighbour("Genusd four", 0.004, "D1")]}

    # the margin cannot see the sister taxon, but the database knows the two cannot be separated
    assert annotate(_result([top, sister]), None, 0.005, True, None) == []
    assert annotate(_result([top, sister]), None, 0.005, True, table) == \
        [NOTE_PREFIX + "Genusb two, Genusd four"]
    # a reference with no inseparable neighbour is not annotated
    assert annotate(_result([top, sister]), None, 0.005, True, {"Z9": []}) == []
    # an organism already named by the mixture check is not repeated
    mixture = _hit(0.995, 50, "B1", "NZ_3 Genusb two strain z")
    assert annotate(_result([top, mixture]), None, 0.005, True, table) == \
        ["Possible mixture with: Genusb two", NOTE_PREFIX + "Genusd four"]


def test_database_info_records_the_sidecars_actually_used(tmp_path):
    """--db-metadata / --db-confusable drive Identification, TaxID and Note, so recording whatever
    happens to sit beside the .msh would assert the wrong provenance."""
    from mashid.confusable import confusable_path
    from mashid.metadata import sidecar_path
    from mashid.pipeline import database_info
    db = tmp_path / "db.msh"
    db.write_bytes(b"not really a sketch")
    sidecar_path(db).write_text("Accession\tOrganism\nA1\tGenusa one\n")
    confusable_path(db).write_text("Accession\tOrganism\tConfusable_Organism\tDistance\n")

    beside = database_info(db, None, md5=False)
    assert beside["metadata_file"] == str(sidecar_path(db))
    assert beside["confusable_file"] == str(confusable_path(db))

    elsewhere = tmp_path / "other.tsv"
    other_conf = tmp_path / "other.confusable.tsv"
    used = database_info(db, None, md5=False, metadata_file=elsewhere, confusable_file=other_conf)
    assert used["metadata_file"] == str(elsewhere)
    assert used["confusable_file"] == str(other_conf)

    # and nothing is claimed when neither exists
    bare = tmp_path / "bare.msh"
    bare.write_bytes(b"x")
    empty = database_info(bare, None, md5=False)
    assert empty["metadata_file"] is None and empty["confusable_file"] is None


def test_run_info_records_the_sidecar_overrides(tmp_path, monkeypatch):
    import json

    from mashid import pipeline
    from mashid.pipeline import Settings, write_run_info
    monkeypatch.setattr(pipeline, "mash_version", lambda exe: "2.3")
    db = tmp_path / "db.msh"
    db.write_bytes(b"not really a sketch")
    settings = Settings(input=tmp_path, output=tmp_path, database=db,
                        db_metadata=tmp_path / "m.tsv", db_confusable=tmp_path / "c.tsv")
    path = tmp_path / "mashID_run.json"
    write_run_info(path, settings, db, None, [], "mash", 0.0, 0)
    params = json.loads(path.read_text())["parameters"]
    assert params["db_metadata"] == str(tmp_path / "m.tsv")
    assert params["db_confusable"] == str(tmp_path / "c.tsv")


def test_run_info_records_the_overrides_absolute(tmp_path, monkeypatch):
    """scripts/reannotate_summary.py reads these back from wherever it happens to run, so a bare
    'curated.tsv' would resolve against that directory: a different file, or none."""
    import json
    from pathlib import Path

    from mashid import pipeline
    from mashid.pipeline import Settings, write_run_info
    monkeypatch.setattr(pipeline, "mash_version", lambda exe: "2.3")
    db = tmp_path / "db.msh"
    db.write_bytes(b"not really a sketch")
    monkeypatch.chdir(tmp_path)
    settings = Settings(input=tmp_path, output=tmp_path, database=db,
                        db_metadata=Path("curated.tsv"), db_confusable=Path("conf.tsv"))
    path = tmp_path / "mashID_run.json"
    write_run_info(path, settings, db, None, [], "mash", 0.0, 0)
    params = json.loads(path.read_text())["parameters"]
    for key in ("db_metadata", "db_confusable"):
        assert Path(params[key]).is_absolute(), params[key]


def test_database_wide_note_does_not_trip_fail_on_note():
    from mashid.confusable import NOTE_PREFIX
    from mashid.pipeline import actionable_notes, exit_code_for
    separable = _result([_hit(0.999, 50, "A1", "NZ_1 Genusa one")])
    separable.notes = [NOTE_PREFIX + "Genusb two"]
    assert actionable_notes(separable.notes) == []
    assert exit_code_for([separable], "note") == 0
    separable.notes.append("Low coverage: median multiplicity 2")
    assert exit_code_for([separable], "note") == 2
    assert exit_code_for([_result([])], "note") == 2        # no hit still fails


def test_identify_prefers_metadata():
    from mashid.metadata import DbEntry
    from mashid.pipeline import identify
    hit = _hit(0.99, 1, "GCF_000000001.1", "contig_1 no organism here")
    assert identify(hit, None) == ("GCF_000000001.1", "GCF_000000001.1", "NA")
    meta = {"GCF_000000001.1": DbEntry("GCF_000000001.1", "Genusa one", "1234")}
    assert identify(hit, meta) == ("GCF_000000001.1", "Genusa one", "1234")
    # version-less match
    meta = {"GCF_000000001.2": DbEntry("GCF_000000001.2", "Genusa one", "1234")}
    assert identify(hit, meta)[1] == "Genusa one"


def test_drop_short_references():
    from mashid.metadata import DbEntry
    from mashid.pipeline import drop_short_references
    hits = [_hit(0.9999, 30, "GCF_000000009.1", "tiny"), _hit(0.999, 30, "GCF_000000001.1", "full")]
    meta = {"GCF_000000009.1": DbEntry("GCF_000000009.1", "Genusx tiny", length=900),
            "GCF_000000001.1": DbEntry("GCF_000000001.1", "Genusa one", length=4_000_000)}
    kept, dropped = drop_short_references(hits, meta, 100_000)
    assert [h.query_id for h in kept] == ["/db/GCF_000000001.1.fna"] and dropped == 1
    assert drop_short_references(hits, meta, 0) == (hits, 0)
    assert drop_short_references(hits, None, 100_000) == (hits, 0)
    # unknown length: kept
    assert drop_short_references(hits, {"GCF_000000009.1": DbEntry("GCF_000000009.1", "x")}, 100_000)[1] == 0


def test_effective_min_ref_length():
    from mashid.metadata import DbEntry
    from mashid.pipeline import DEFAULT_MIN_REF_LENGTH, effective_min_ref_length
    bacteria = {f"A{i}": DbEntry(f"A{i}", "Genusa one", length=4_000_000) for i in range(5)}
    viruses = {f"V{i}": DbEntry(f"V{i}", "Virus x", length=30_000) for i in range(5)}
    assert effective_min_ref_length(None, bacteria) == DEFAULT_MIN_REF_LENGTH
    assert effective_min_ref_length(None, viruses) == 0
    assert effective_min_ref_length(None, {**bacteria, **viruses, "V9": DbEntry("V9", "v", length=1000)}) == 0
    assert effective_min_ref_length(None, None) == DEFAULT_MIN_REF_LENGTH
    assert effective_min_ref_length(50_000, viruses) == 50_000  # explicit values are respected
    assert effective_min_ref_length(-5, viruses) == 0
