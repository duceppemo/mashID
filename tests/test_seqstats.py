import gzip

from mashid.seqstats import sample_stats, sequence_stats

FASTQ = "@r1\nACGT\n+\nIIII\n@r2\nACGTACGT\n+\nIIIIIIII\n"
FASTA = ">c1 desc\nACGT\nAC\n>c2\n\nACGTACGT\n"


def test_fastq_plain_and_gz(tmp_path):
    (tmp_path / "a.fastq").write_text(FASTQ)
    with gzip.open(tmp_path / "a.fastq.gz", "wt") as fh:
        fh.write(FASTQ)
    assert sequence_stats(tmp_path / "a.fastq", "fastq") == (2, 12)
    assert sequence_stats(tmp_path / "a.fastq.gz", "fastq") == (2, 12)


def test_fasta_plain_and_gz(tmp_path):
    (tmp_path / "a.fna").write_text(FASTA)
    with gzip.open(tmp_path / "a.fna.gz", "wt") as fh:
        fh.write(FASTA)
    assert sequence_stats(tmp_path / "a.fna", "fasta") == (2, 14)
    assert sequence_stats(tmp_path / "a.fna.gz", "fasta") == (2, 14)


def test_crlf_and_gzip_detected_by_content(tmp_path):
    (tmp_path / "a.fastq").write_bytes(FASTQ.replace("\n", "\r\n").encode())
    assert sequence_stats(tmp_path / "a.fastq", "fastq") == (2, 12)
    # gzipped content with a misleading (non-.gz) name still works
    with gzip.open(tmp_path / "b.fastq", "wt") as fh:
        fh.write(FASTQ)
    assert sequence_stats(tmp_path / "b.fastq", "fastq") == (2, 12)


def test_sample_stats_sums_pairs(tmp_path):
    (tmp_path / "r1.fq").write_text(FASTQ)
    (tmp_path / "r2.fq").write_text(FASTQ)
    assert sample_stats([tmp_path / "r1.fq", tmp_path / "r2.fq"], "fastq") == (4, 24)
