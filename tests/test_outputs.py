"""Sample sheets, JSON / MultiQC / run-info outputs, --fail-on and Est_Depth (need mash)."""

import csv
import gzip
import json
import random
import shutil
from pathlib import Path

import pytest

from mashid import MashIDError
from mashid.cli import main as mashid_main
from mashid.samples import read_sample_sheet

pytestmark = pytest.mark.skipif(shutil.which("mash") is None, reason="mash not installed")

EXAMPLE = Path(__file__).resolve().parent.parent / "example"
DB = EXAMPLE / "db" / "example.msh"


def _read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_sample_sheet_parsing(tmp_path):
    (tmp_path / "a_R1.fq").write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "a_R2.fq").write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "b.fa").write_text(">c\nACGT\n")
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text("sample\tfile\texpected\nS1\ta_R1.fq\tGenusa\nS1\ta_R2.fq\tGenusa\nS2\tb.fa\tGenusb\n")
    samples, extra = read_sample_sheet(sheet)
    assert extra == ["expected"]
    assert [s.name for s in samples] == ["S1", "S2"]
    assert len(samples[0].files) == 2 and samples[0].seq_type == "fastq" and samples[0].extra == {"expected": "Genusa"}
    assert samples[1].seq_type == "fasta"

    sheet.write_text("sample,files\nS1,a_R1.fq;a_R2.fq\n")  # csv, files joined with ';'
    samples, extra = read_sample_sheet(sheet)
    assert len(samples[0].files) == 2 and extra == []

    sheet.write_text("sample\tfile\nS1\tmissing.fq\n")
    with pytest.raises(MashIDError, match="file not found"):
        read_sample_sheet(sheet)
    sheet.write_text("foo\tbar\nx\ty\n")
    with pytest.raises(MashIDError, match="expected columns"):
        read_sample_sheet(sheet)
    sheet.write_text("sample\tfile\nS1\ta_R1.fq\nS1\tb.fa\n")
    with pytest.raises(MashIDError, match="mixes"):
        read_sample_sheet(sheet)


def test_sample_sheet_run_and_outputs(tmp_path):
    reads = EXAMPLE / "reads"
    sheet = tmp_path / "sheet.tsv"
    sheet.write_text(
        "sample\tfile\tsite\n"
        f"alphaX\t{reads / 'alpha_S1_L001_R1_001.fastq.gz'}\tLab1\n"
        f"alphaX\t{reads / 'alpha_S1_L001_R2_001.fastq.gz'}\tLab1\n"
        f"nothing\t{reads / 'unknown.fasta'}\tLab2\n"
    )
    out = tmp_path / "out"
    rc = mashid_main(["--sample-sheet", str(sheet), "-o", str(out), "-d", str(DB), "-t", "2"])
    assert rc == 0
    rows = {r["Sample"]: r for r in _read_tsv(out / "summary_mashID.tsv")}
    assert list(rows) == ["alphaX", "nothing"]
    assert rows["alphaX"]["site"] == "Lab1" and rows["alphaX"]["Identification"] == "Exemplaria alpha"
    # Est_Depth: 8000 reads * 150 bp / 120 kb = 10.0
    assert rows["alphaX"]["Est_Depth"] == "10.0" and rows["nothing"]["Est_Depth"] == "NA"

    doc = json.loads((out / "summary_mashID.json").read_text())
    assert doc["samples"][0]["sample"] == "alphaX" and doc["samples"][0]["top_hit"]["Accession"] == "GCF_000000001.1"
    assert doc["samples"][0]["metadata"] == {"site": "Lab1"} and doc["samples"][1]["top_hit"] is None
    assert doc["samples"][0]["hits"][0]["Rank"] == "1"

    mqc = (out / "summary_mashID_mqc.tsv").read_text().splitlines()
    assert mqc[0] == "# id: mashid" and "plot_type: table" in "\n".join(mqc[:8])
    table = [ln for ln in mqc if not ln.startswith("#")]
    assert table[0].split("\t")[0] == "Sample" and table[1].startswith("alphaX\tExemplaria alpha\t")

    info = json.loads((out / "mashID_run.json").read_text())
    assert info["mashid_version"] and info["mash_version"].startswith("2.")
    assert info["database"]["path"] == str(DB) and len(info["database"]["md5"]) == 32
    assert info["database"]["references"] == 3 and info["parameters"]["threads"] == 2
    assert [s["name"] for s in info["samples"]] == ["alphaX", "nothing"] and info["exit_code"] == 0
    assert info["command_line"][0] == "mashID" and info["duration_seconds"] >= 0


def test_fail_on(tmp_path):
    out = tmp_path / "out"
    args = ["-i", str(EXAMPLE / "reads"), "-o", str(out), "-d", str(DB), "-t", "2"]
    assert mashid_main(args) == 0
    assert mashid_main([*args, "--fail-on", "no-hit"]) == 2     # 'unknown' has no hit
    assert mashid_main([*args, "--fail-on", "note"]) == 2       # the mixture has a note
    solo = tmp_path / "solo"
    solo.mkdir()
    shutil.copy(EXAMPLE / "reads" / "delta_assembly.fasta.gz", solo)
    assert mashid_main(["-i", str(solo), "-o", str(out), "-d", str(DB), "--fail-on", "note", "-t", "1"]) == 0
    info = json.loads((out / "mashID_run.json").read_text())
    assert info["exit_code"] == 0


def test_neither_input_nor_sheet_is_an_error(capsys):
    with pytest.raises(SystemExit):
        mashid_main(["-o", "x"])
    assert "sample-sheet" in capsys.readouterr().err


def test_db_check_reports_issues(tmp_path):
    from mashid.makedb import check_database
    from mashid.makedb import main as makedb_main
    rng = random.Random(9)
    genomes = tmp_path / "g"
    genomes.mkdir()
    seq = "".join(rng.choice("ACGT") for _ in range(120_000))
    (genomes / "GCF_000000001.1.fna").write_text(">NZ_1 Genusa one strain x\n" + seq + "\n")
    (genomes / "GCF_000000002.1.fna").write_text(">NZ_2 Genusb one strain y\n" + seq[::-1] + "\n")  # synonym-like
    (genomes / "GCF_000000003.1.fna").write_text(">contig_1\n" + seq[:900] + "\n")  # short, no name
    with gzip.open(genomes / "GCF_000000004.1.fna.gz", "wt") as fh:
        fh.write(">NZ_4 Genusc two\n" + seq[5000:] + "\n")
    assert makedb_main(["-i", str(genomes), "-o", str(tmp_path), "-p", "chk", "-s", "500", "--min-length", "0"]) == 0
    lines = check_database(tmp_path / "chk.msh")
    text = "\n".join(lines)
    assert "sketches: 4" in text and "sketch size: 500" in text
    assert "1 reference(s) shorter than 100000 bp" in text and "GCF_000000003.1" in text
    assert "without a usable organism name" in text
    assert "TaxIDs: 0/4" in text
    assert "one: Genusa, Genusb" in text
    assert "issue(s) found" in text
    assert makedb_main(["--check", str(tmp_path / "chk.msh")]) == 0
    assert makedb_main(["--check", str(tmp_path / "nope.msh")]) == 1
