from pathlib import Path

import pytest

from mashid import MashIDError
from mashid.samples import (
    discover_samples,
    find_sequence_files,
    sample_name_from_file,
    seq_type_of,
    strip_seq_extension,
)


@pytest.mark.parametrize(
    "filename, expected",
    [
        ("MBWGS440_S10_L001_R1_001.fastq.gz", "MBWGS440"),
        ("MBWGS440_S10_L001_R2_001.fastq.gz", "MBWGS440"),
        ("isolate_A_R1.fq.gz", "isolate_A"),
        ("isolate_A_1.fastq", "isolate_A"),
        ("isolate_A_2.fastq", "isolate_A"),
        ("sample.fastq", "sample"),
        ("GCF_000195955.2.fna", "GCF_000195955.2"),
        ("assembly.FASTA.GZ", "assembly"),
        ("barcode01_pass.fastq.gz", "barcode01_pass"),
    ],
)
def test_sample_name_from_file(filename, expected):
    assert sample_name_from_file(filename) == expected


def test_strip_seq_extension():
    assert strip_seq_extension("a.fastq.gz") == "a"
    assert strip_seq_extension("a.fa") == "a"
    assert strip_seq_extension("a.txt") is None
    assert strip_seq_extension("a.gz") is None


def test_seq_type_of():
    assert seq_type_of("x.fq.gz") == "fastq"
    assert seq_type_of("x.fna") == "fasta"
    with pytest.raises(ValueError):
        seq_type_of("x.txt")


def _touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(">x\nACGT\n" if "fa" in path.suffixes[0] else "@r\nACGT\n+\nIIII\n")
    return path


def test_discover_groups_pairs_and_recurses(tmp_path):
    _touch(tmp_path / "S1_S1_L001_R1_001.fastq.gz")
    _touch(tmp_path / "S1_S1_L001_R2_001.fastq.gz")
    _touch(tmp_path / "sub" / "asm.fasta")
    _touch(tmp_path / "notes.txt")
    _touch(tmp_path / ".hidden.fastq")

    samples = discover_samples(tmp_path)
    assert [s.name for s in samples] == ["S1", "asm"]
    assert len(samples[0].files) == 2 and samples[0].seq_type == "fastq"
    assert len(samples[1].files) == 1 and samples[1].seq_type == "fasta"


def test_discover_single_file(tmp_path):
    f = _touch(tmp_path / "one_R1.fq")
    samples = discover_samples(f)
    assert len(samples) == 1 and samples[0].name == "one" and samples[0].files == [f.absolute()]


def test_symlink_name_defines_the_sample(tmp_path):
    real = _touch(tmp_path / "data" / "GCF_000012005.1_ASM1200v1_genomic.fna.gz")
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "isolate42.fna.gz").symlink_to(real)
    samples = discover_samples(tmp_path / "in")
    assert [s.name for s in samples] == ["isolate42"]
    assert samples[0].files[0].read_bytes() == real.read_bytes()  # still readable through the link


def test_more_than_two_files_warns_but_groups(tmp_path):
    for lane in ("L001", "L002"):
        for r in ("R1", "R2"):
            _touch(tmp_path / f"S1_S1_{lane}_{r}_001.fastq")
    samples = discover_samples(tmp_path)
    assert len(samples) == 1 and len(samples[0].files) == 4 and samples[0].warnings


def test_mixed_types_is_an_error(tmp_path):
    _touch(tmp_path / "S1.fastq")
    _touch(tmp_path / "S1.fasta")
    with pytest.raises(MashIDError, match="mixes"):
        discover_samples(tmp_path)


def test_errors(tmp_path):
    with pytest.raises(MashIDError, match="No sequence files"):
        find_sequence_files(tmp_path)
    with pytest.raises(MashIDError, match="does not exist"):
        find_sequence_files(tmp_path / "missing")
    _touch(tmp_path / "notes.txt")
    with pytest.raises(MashIDError, match="recognised extension"):
        find_sequence_files(tmp_path / "notes.txt")


def test_dangling_symlinks_are_reported(tmp_path):
    _touch(tmp_path / "good.fastq")
    (tmp_path / "gone_R1.fastq.gz").symlink_to(tmp_path / "missing" / "gone_R1.fastq.gz")
    with pytest.raises(MashIDError, match=r"1 sequence file link\(s\) point to missing files"):
        find_sequence_files(tmp_path)
    with pytest.raises(MashIDError, match="link to a missing file"):
        find_sequence_files(tmp_path / "gone_R1.fastq.gz")


def test_symlink_cycle_does_not_loop(tmp_path):
    _touch(tmp_path / "a" / "S1.fastq")
    (tmp_path / "a" / "loop").symlink_to(tmp_path, target_is_directory=True)
    samples = discover_samples(tmp_path)
    assert [s.name for s in samples] == ["S1"] and len(samples[0].files) == 1
