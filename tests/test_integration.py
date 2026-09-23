"""End-to-end tests that need the ``mash`` executable (skipped otherwise)."""

import csv
import gzip
import random
import shutil
import subprocess
from pathlib import Path

import pytest

from mashid.cli import main as mashid_main
from mashid.makedb import main as makedb_main

pytestmark = pytest.mark.skipif(shutil.which("mash") is None, reason="mash not installed")

GENOME_LEN = 120_000
ORGANISMS = {
    "GCF_000000001.1": "NZ_TEST01.1 Genusa speciesone strain A chromosome, complete genome",
    "GCF_000000002.1": "NZ_TEST02.1 Genusb speciestwo subsp. three strain B chromosome",
    "GCF_000000003.1": "NZ_TEST03.1 Genusc sp. XY99 chromosome",
}


def _random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    root = tmp_path_factory.mktemp("db")
    rng = random.Random(1)
    genomes = {}
    fasta_dir = root / "fasta"
    fasta_dir.mkdir()
    for i, (acc, header) in enumerate(ORGANISMS.items()):
        seq = _random_seq(rng, GENOME_LEN)
        genomes[acc] = seq
        content = f">{header}\n" + "\n".join(seq[j:j + 80] for j in range(0, len(seq), 80)) + "\n"
        if i == 1:  # exercise gzipped database input
            with gzip.open(fasta_dir / f"{acc}.fna.gz", "wt") as fh:
                fh.write(content)
        else:
            (fasta_dir / f"{acc}.fna").write_text(content)
    meta = root / "meta.tsv"
    meta.write_text("Accession\tOrganism\tTaxID\nGCF_000000001.1\tGenusa speciesone (curated)\t111\n")
    rc = makedb_main(["-i", str(fasta_dir), "-o", str(root), "-p", "test_db.msh", "-s", "2000", "-t", "2",
                      "--metadata", str(meta), "--min-length", "0"])
    assert rc == 0
    db_path = root / "test_db.msh"
    assert db_path.is_file()
    sidecar = root / "test_db.metadata.tsv"
    assert sidecar.is_file()
    rows = {r["Accession"]: r for r in _read_tsv(sidecar)}
    assert rows["GCF_000000001.1"]["Organism"] == "Genusa speciesone (curated)" and rows["GCF_000000001.1"]["TaxID"] == "111"
    assert rows["GCF_000000002.1"]["Organism"] == "Genusb speciestwo subsp. three" and rows["GCF_000000002.1"]["TaxID"] == "NA"
    assert rows["GCF_000000001.1"]["Length"] == str(GENOME_LEN) and rows["GCF_000000001.1"]["Hashes"] == "2000"
    return db_path, genomes


def _write_reads(path: Path, genome: str, rng: random.Random, n: int = 8000, length: int = 150,
                 gz: bool = True) -> None:
    opener = gzip.open if gz else open
    with opener(path, "wt") as fh:
        for i in range(n):
            start = rng.randrange(0, len(genome) - length)
            fh.write(f"@read{i}\n{genome[start:start + length]}\n+\n{'I' * length}\n")


def _read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_identify_reads_pairs_and_assembly(db, tmp_path):
    db_path, genomes = db
    rng = random.Random(2)
    inp = tmp_path / "input"
    inp.mkdir()
    # paired-end reads from genome 1 (two files, screened together)
    _write_reads(inp / "sampleA_S1_L001_R1_001.fastq.gz", genomes["GCF_000000001.1"], rng, n=5000)
    _write_reads(inp / "sampleA_S1_L001_R2_001.fastq.gz", genomes["GCF_000000001.1"], rng, n=5000)
    # single-end plain fastq from genome 2, in a sub-folder
    (inp / "sub").mkdir()
    _write_reads(inp / "sub" / "sampleB.fastq", genomes["GCF_000000002.1"], rng, gz=False)
    # assembly of genome 3
    (inp / "sampleC.fasta").write_text(">contig1\n" + genomes["GCF_000000003.1"] + "\n")
    # random sequence matching nothing
    (inp / "sampleD.fa").write_text(">contig1\n" + _random_seq(rng, 20_000) + "\n")

    out = tmp_path / "out"
    rc = mashid_main(["-i", str(inp), "-o", str(out), "-d", str(db_path), "-t", "2", "-p", "2"])
    assert rc == 0

    summary = {r["Sample"]: r for r in _read_tsv(out / "summary_mashID.tsv")}
    assert list(summary) == ["sampleA", "sampleB", "sampleC", "sampleD"]
    assert summary["sampleA"]["Identification"] == "Genusa speciesone (curated)"
    assert summary["sampleA"]["Accession"] == "GCF_000000001.1" and summary["sampleA"]["TaxID"] == "111"
    assert summary["sampleB"]["TaxID"] == "NA"
    assert all(summary[s]["Note"] == "" for s in summary)
    assert summary["sampleA"]["Sequences"] == "10000" and summary["sampleA"]["Bases"] == str(10000 * 150)
    assert summary["sampleB"]["Identification"] == "Genusb speciestwo subsp. three"
    assert summary["sampleC"]["Identification"] == "Genusc sp. XY99"
    assert summary["sampleC"]["Sequences"] == "1" and summary["sampleC"]["Bases"] == str(GENOME_LEN)
    assert summary["sampleD"]["Identification"] == "No significant hit in database"
    assert summary["sampleD"]["Identity"] == "NA"
    for s in ("sampleA", "sampleB", "sampleC"):
        assert float(summary[s]["Identity"]) > 0.95

    per_sample = _read_tsv(out / "sampleA_mashID.tsv")
    assert per_sample[0]["Rank"] == "1" and per_sample[0]["Identification"] == "Genusa speciesone (curated)"
    assert "Description" in per_sample[0]
    assert (out / "sampleD_mashID.tsv").read_text().count("\n") == 1  # header only
    assert not (out / "tmp").exists()


def test_skip_stats_and_single_file(db, tmp_path):
    db_path, genomes = db
    reads = tmp_path / "solo_R1.fq.gz"
    _write_reads(reads, genomes["GCF_000000003.1"], random.Random(3), n=6000)
    out = tmp_path / "out"
    rc = mashid_main(["-i", str(reads), "-o", str(out), "-d", str(db_path), "--skip-stats", "-t", "1", "-p", "1"])
    assert rc == 0
    row = _read_tsv(out / "summary_mashID.tsv")[0]
    assert row["Sample"] == "solo" and row["Sequences"] == "NA" and row["Identification"] == "Genusc sp. XY99"


def test_bad_database_is_reported(db, tmp_path):
    reads = tmp_path / "x.fq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    (tmp_path / "bad.msh").write_text("not a sketch")
    assert mashid_main(["-i", str(reads), "-o", str(tmp_path / "o"), "-d", str(tmp_path / "bad.msh")]) == 1
    assert mashid_main(["-i", str(reads), "-o", str(tmp_path / "o"), "-d", str(tmp_path / "missing.msh")]) == 1


def test_makedb_validation(tmp_path):
    (tmp_path / "g.fna").write_text(">x\nACGT\n")
    assert makedb_main(["-i", str(tmp_path), "-o", str(tmp_path / "db"), "-k", "40"]) == 1
    assert makedb_main(["-i", str(tmp_path / "nope"), "-o", str(tmp_path / "db")]) == 1
    assert makedb_main(["-i", str(tmp_path), "-o", str(tmp_path / "db"), "-p", "a/b"]) == 1


def test_max_reads_streams_a_subset(db, tmp_path):
    db_path, genomes = db
    inp = tmp_path / "input"
    inp.mkdir()
    rng = random.Random(4)
    _write_reads(inp / "sub_R1.fastq.gz", genomes["GCF_000000002.1"], rng, n=5000)
    _write_reads(inp / "sub_R2.fastq.gz", genomes["GCF_000000002.1"], rng, n=5000)
    out = tmp_path / "out"
    rc = mashid_main(["-i", str(inp), "-o", str(out), "-d", str(db_path), "--max-reads", "3001", "-t", "2"])
    assert rc == 0
    row = _read_tsv(out / "summary_mashID.tsv")[0]
    # 3001 reads split over two files -> 1501 per file
    assert row["Sequences"] == "3002" and row["Bases"] == str(3002 * 150)
    assert row["Identification"] == "Genusb speciestwo subsp. three"


def test_notes_flag_mixture_and_ambiguity(db, tmp_path):
    db_path, genomes = db
    inp = tmp_path / "input"
    inp.mkdir()
    rng = random.Random(5)
    # a mixture: reads from genome 1 and genome 3 in one sample
    mixed = inp / "mixed.fastq"
    _write_reads(mixed, genomes["GCF_000000001.1"], rng, n=6000, gz=False)
    reads3 = tmp_path / "g3.fastq"
    _write_reads(reads3, genomes["GCF_000000003.1"], rng, n=6000, gz=False)
    with open(mixed, "a") as dst, open(reads3) as src:
        dst.write(src.read())
    out = tmp_path / "out"
    rc = mashid_main(["-i", str(inp), "-o", str(out), "-d", str(db_path), "-t", "2"])
    assert rc == 0
    row = _read_tsv(out / "summary_mashID.tsv")[0]
    assert "Possible mixture with:" in row["Note"]

    # ambiguity: a second reference nearly identical to genome 1 but named differently
    twin_dir = tmp_path / "twin"
    twin_dir.mkdir()
    g1 = genomes["GCF_000000001.1"]
    twin = g1[:30_000] + ("T" if g1[30_000] != "T" else "A") + g1[30_001:]
    for acc, header, seq in [("GCF_000000001.1", ORGANISMS["GCF_000000001.1"], g1),
                             ("GCF_000000009.1", "NZ_TEST09.1 Genusa twin strain Q", twin)]:
        (twin_dir / f"{acc}.fna").write_text(f">{header}\n{seq}\n")
    assert makedb_main(["-i", str(twin_dir), "-o", str(twin_dir), "-p", "twin", "-s", "2000", "--min-length", "0"]) == 0
    inp2 = tmp_path / "input2"
    inp2.mkdir()
    _write_reads(inp2 / "one.fastq.gz", g1, rng, n=6000)
    out2 = tmp_path / "out2"
    rc = mashid_main(["-i", str(inp2), "-o", str(out2), "-d", str(twin_dir / "twin.msh"), "-t", "2",
                      "--no-winner-take-all"])
    assert rc == 0
    row = _read_tsv(out2 / "summary_mashID.tsv")[0]
    assert row["Note"].startswith("Ambiguous: Genusa ")


def test_low_coverage_note(db, tmp_path):
    db_path, genomes = db
    reads = tmp_path / "shallow.fastq"
    _write_reads(reads, genomes["GCF_000000001.1"], random.Random(6), n=1000, gz=False)  # ~1.25x
    out = tmp_path / "out"
    assert mashid_main(["-i", str(reads), "-o", str(out), "-d", str(db_path), "-t", "1"]) == 0
    row = _read_tsv(out / "summary_mashID.tsv")[0]
    assert row["Identification"] == "Genusa speciesone (curated)" and row["Note"].startswith("Low coverage")


def test_annotate_existing_database(db, tmp_path):
    db_path, _ = db
    report = tmp_path / "assembly_data_report.jsonl"
    report.write_text('{"accession": "GCF_000000003.1", "organism": {"organismName": "Genusc curated", "taxId": 333}}\n')
    assert makedb_main(["--annotate", str(db_path), "--assembly-report", str(report)]) == 0
    rows = {r["Accession"]: r for r in _read_tsv(db_path.with_name("test_db.metadata.tsv"))}
    assert rows["GCF_000000003.1"]["Organism"] == "Genusc curated" and rows["GCF_000000003.1"]["TaxID"] == "333"
    assert rows["GCF_000000001.1"]["TaxID"] == "NA"  # no longer curated: rebuilt from headers
    # restore the sidecar used by the other tests
    meta = db_path.parent / "meta.tsv"
    assert makedb_main(["--annotate", str(db_path), "--metadata", str(meta)]) == 0


def test_short_reference_is_ignored_by_default_and_excluded_at_build(db, tmp_path):
    db_path, genomes = db
    g1 = genomes["GCF_000000001.1"]
    # a 900 bp fragment of genome 1 labelled as another organism: fully "contained" in any genome-1 sample
    frag_dir = tmp_path / "frag"
    frag_dir.mkdir()
    (frag_dir / "GCF_000000001.1.fna").write_text(f">{ORGANISMS['GCF_000000001.1']}\n{g1}\n")
    (frag_dir / "GCF_000000008.1.fna").write_text(">NZ_TEST08.1 Genusz fragment strain F\n" + g1[1000:1900] + "\n")

    # 1. kept in the database when --min-length 0 ...
    assert makedb_main(["-i", str(frag_dir), "-o", str(frag_dir), "-p", "withfrag", "-s", "500", "--min-length", "0"]) == 0
    # sample: genome 1 with ~1% mutations everywhere except the fragment region, so the full reference
    # scores < 1.0 while the tiny fragment is perfectly contained (what a partial record does on real data)
    rng = random.Random(7)
    mutated = list(g1)
    for i in range(0, len(mutated), 100):
        if not 1000 <= i < 1900:
            mutated[i] = "A" if mutated[i] != "A" else "C"
    reads = tmp_path / "s.fastq.gz"
    _write_reads(reads, "".join(mutated), rng, n=6000)
    out = tmp_path / "out_keep"
    assert mashid_main(["-i", str(reads), "-o", str(out), "-d", str(frag_dir / "withfrag.msh"), "--min-ref-length", "0",
                        "-t", "1"]) == 0
    assert _read_tsv(out / "summary_mashID.tsv")[0]["Identification"] == "Genusz fragment"  # the false call
    # ... but ignored by mashID by default
    out = tmp_path / "out_default"
    assert mashid_main(["-i", str(reads), "-o", str(out), "-d", str(frag_dir / "withfrag.msh"), "-t", "1"]) == 0
    assert _read_tsv(out / "summary_mashID.tsv")[0]["Identification"] == "Genusa speciesone"

    # 2. excluded at build time with the default --min-length
    assert makedb_main(["-i", str(frag_dir), "-o", str(frag_dir), "-p", "clean", "-s", "500"]) == 0
    rows = _read_tsv(frag_dir / "clean.metadata.tsv")
    assert [r["Accession"] for r in rows] == ["GCF_000000001.1"]
    # 3. everything too short: error
    only_frag = tmp_path / "onlyfrag"
    only_frag.mkdir()
    (only_frag / "GCF_000000008.1.fna").write_text(">NZ_TEST08.1 Genusz fragment\n" + g1[:900] + "\n")
    assert makedb_main(["-i", str(only_frag), "-o", str(only_frag), "-p", "none", "-s", "500"]) == 1


def test_metadata_is_built_on_the_fly_without_sidecar(db, tmp_path):
    db_path, genomes = db
    bare = tmp_path / "bare"
    bare.mkdir()
    shutil.copy(db_path, bare / "bare.msh")
    reads = tmp_path / "s.fastq.gz"
    _write_reads(reads, genomes["GCF_000000002.1"], random.Random(8), n=6000)
    out = tmp_path / "out"
    assert mashid_main(["-i", str(reads), "-o", str(out), "-d", str(bare / "bare.msh"), "-t", "1"]) == 0
    assert (bare / "bare.metadata.tsv").is_file()  # written for next time
    row = _read_tsv(out / "summary_mashID.tsv")[0]
    assert row["Identification"] == "Genusb speciestwo subsp. three" and row["TaxID"] == "NA"


def test_non_utf8_reference_header_does_not_crash(tmp_path):
    """Mash echoes headers verbatim; a Latin-1 byte must not raise UnicodeDecodeError."""
    from mashid.mash import find_mash, info_table
    from mashid.metadata import entries_from_sketches
    fasta = tmp_path / "latin1.fna"
    fasta.write_bytes(b">NZ_X1.1 Bacillus subtilis strain caf\xe9 chromosome\n" + b"ACGT" * 300 + b"\n")
    exe = find_mash()
    subprocess.run([exe, "sketch", "-o", str(tmp_path / "latin1"), str(fasta)], check=True, capture_output=True)
    rows = info_table(tmp_path / "latin1.msh", exe)
    assert len(rows) == 1 and "Bacillus subtilis" in rows[0].comment
    entries, _ = entries_from_sketches(rows)
    assert entries[0].organism == "Bacillus subtilis"


def test_sample_named_summary_is_rejected(db, tmp_path):
    db_path, genomes = db
    inp = tmp_path / "in"
    inp.mkdir()
    (inp / "summary.fasta").write_text(">c\n" + genomes["GCF_000000001.1"][:5000] + "\n")
    assert mashid_main(["-i", str(inp), "-o", str(tmp_path / "out"), "-d", str(db_path), "-t", "1"]) == 1


def test_makedb_list_file_relative_paths_and_bad_prefix(db, tmp_path):
    db_path, genomes = db
    gdir = tmp_path / "g"
    gdir.mkdir()
    (gdir / "GCF_000000001.1.fna").write_text(">NZ_1 Genusa one\n" + genomes["GCF_000000001.1"] + "\n")
    listing = tmp_path / "list.txt"
    listing.write_text("# relative to this file\ng/GCF_000000001.1.fna\n")
    assert makedb_main(["-i", str(listing), "-o", str(tmp_path / "db"), "-p", "rel", "-s", "500"]) == 0
    assert (tmp_path / "db" / "rel.msh").is_file()
    assert makedb_main(["-i", str(listing), "-o", str(tmp_path / "db"), "--prefix=-bad", "-s", "500"]) == 1
