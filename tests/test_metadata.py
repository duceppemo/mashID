import json

import pytest

from mashid import MashIDError
from mashid.metadata import DbEntry, read_metadata, read_ncbi_assembly_report, sidecar_path, write_metadata


def test_sidecar_path(tmp_path):
    assert sidecar_path(tmp_path / "foo.msh") == tmp_path / "foo.metadata.tsv"
    assert sidecar_path(tmp_path / "foo") == tmp_path / "foo.metadata.tsv"


def test_roundtrip(tmp_path):
    p = tmp_path / "db.metadata.tsv"
    entries = [DbEntry("GCF_1.1", "Genusa one", "10", 4_000_000, 1000, "/x/GCF_1.1.fna", "desc"),
               DbEntry("GCF_2.1", "Genusb two")]
    write_metadata(p, entries)
    back = read_metadata(p)
    assert back["GCF_1.1"].organism == "Genusa one" and back["GCF_1.1"].taxid == "10"
    assert back["GCF_1.1"].length == 4_000_000 and back["GCF_1.1"].hashes == 1000
    assert back["GCF_1.1"].description == "desc" and back["GCF_2.1"].taxid == "NA"
    assert back["GCF_2.1"].length is None and back["GCF_2.1"].hashes is None


def test_read_user_table_flexible_headers(tmp_path):
    p = tmp_path / "meta.csv"
    p.write_text("﻿assembly accession,organism name,taxid\nGCF_1.1,Genusa one,10\n,skipped,\n")
    back = read_metadata(p)
    assert list(back) == ["GCF_1.1"] and back["GCF_1.1"].taxid == "10"
    p.write_text("foo\tbar\nx\ty\n")
    with pytest.raises(MashIDError, match="expected columns"):
        read_metadata(p)
    with pytest.raises(MashIDError, match="not found"):
        read_metadata(tmp_path / "nope.tsv")


def test_read_ncbi_assembly_report(tmp_path):
    p = tmp_path / "assembly_data_report.jsonl"
    rows = [
        {"accession": "GCF_000195955.2", "pairedAccession": "GCA_000195955.2",
         "organism": {"organismName": "Mycobacterium tuberculosis H37Rv", "taxId": 83332}},  # strain text dropped
        {"accession": "GCF_000000002.1", "organism": {"organismName": "Genusb two"}},
    ]
    p.write_text("\n".join(json.dumps(r) for r in rows) + "\n\n")
    back = read_ncbi_assembly_report(p)
    assert back["GCF_000195955.2"] == ("Mycobacterium tuberculosis", "83332")
    assert back["GCA_000195955.2"] == ("Mycobacterium tuberculosis", "83332")
    assert back["GCF_000000002.1"] == ("Genusb two", "NA")
    # `datasets summary ... --as-json-lines` uses snake_case keys
    p.write_text(json.dumps({"accession": "GCF_000000005.1", "paired_accession": "GCA_000000005.1",
                             "organism": {"organism_name": "Salmonella bongori", "tax_id": 54736}}) + "\n")
    back = read_ncbi_assembly_report(p)
    assert back["GCF_000000005.1"] == ("Salmonella bongori", "54736") and "GCA_000000005.1" in back
    p.write_text("{not json\n")
    with pytest.raises(MashIDError, match="not valid JSON"):
        read_ncbi_assembly_report(p)


def test_entries_from_sketches():
    from mashid.mash import SketchInfo
    from mashid.metadata import entries_from_sketches
    sketches = [SketchInfo(1000, 5_000_000, "/db/GCF_1.1.fna", "NZ_1 Genusa one strain x"),
                SketchInfo(947, 900, "/db/GCF_2.1.fna", "NZ_2 Genusb two")]
    entries, unannotated = entries_from_sketches(sketches, {"GCF_1.1": ("Curated one", "11")})
    assert unannotated == 1
    assert entries[0].organism == "Curated one" and entries[0].taxid == "11" and entries[0].length == 5_000_000
    assert entries[1].organism == "Genusb two" and entries[1].hashes == 947 and entries[1].description == "NZ_2 Genusb two"


def test_normalise_organism_name():
    from mashid.metadata import normalise_organism_name as n
    assert n("Mycobacterium tuberculosis variant bovis BCG str. Sweden") == "Mycobacterium tuberculosis variant bovis"
    assert n("Mycobacterium tuberculosis H37Rv") == "Mycobacterium tuberculosis"
    assert n("Mycobacteroides abscessus subsp. massiliense str. GO 06") == "Mycobacteroides abscessus subsp. massiliense"
    assert n("Mycobacterium sp. JS623") == "Mycobacterium sp. JS623"
    assert n("uncultured Mycobacterium sp.") == "Mycobacterium sp."  # MAG prefix dropped
    assert n("Bacterium X") == "Bacterium X"  # no binomial: unchanged
    assert n("") == "NA" and n("NA") == "NA"


def test_empty_or_binary_metadata_is_a_clean_error(tmp_path):
    p = tmp_path / "empty.tsv"
    p.write_text("")
    with pytest.raises(MashIDError, match="is empty"):
        read_metadata(p)
    p.write_bytes(b"\xff\xfe\x00garbage")
    with pytest.raises(MashIDError, match="not a UTF-8"):
        read_metadata(p)
