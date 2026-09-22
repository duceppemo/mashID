#!/usr/bin/env python3
"""Regenerate the example dataset (deterministic; needs mash and mashID installed).

Three fictional organisms with random 120 kb genomes are sketched into example/db/example.msh, and four
samples are derived from them:

  reads/alpha_S1_L001_R{1,2}_001.fastq.gz   paired-end reads of Exemplaria alpha
  reads/mixed_beta_delta.fastq.gz           single-end reads, half E. beta subsp. gamma, half Fictivibrio delta
                                            (a mixed culture: mashID reports the mixture in the Note column)
  reads/delta_assembly.fasta.gz             two-contig assembly of Fictivibrio delta
  reads/unknown.fasta                       a random contig matching nothing

expected/summary_mashID.tsv holds the summary mashID produces for them.
"""

from __future__ import annotations

import gzip
import random
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
GENOME_LEN = 120_000
READ_LEN = 150
QUAL = "I" * READ_LEN  # constant qualities keep the gzipped files small
ORGANISMS = {
    "GCF_000000001.1": "NZ_EXAMPL01000001.1 Exemplaria alpha strain EX-1 chromosome, complete genome",
    "GCF_000000002.1": "NZ_EXAMPL02000001.1 Exemplaria beta subsp. gamma strain EX-2 chromosome, complete genome",
    "GCF_000000003.1": "NZ_EXAMPL03000001.1 Fictivibrio delta strain EX-3 chromosome, complete genome",
}


def random_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            fh.write("\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")


def write_reads(fh, genome: str, rng: random.Random, n: int, prefix: str) -> None:
    for i in range(n):
        start = rng.randrange(0, len(genome) - READ_LEN)
        read = genome[start:start + READ_LEN]
        if rng.random() < 0.5:  # reverse strand
            read = read.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        fh.write(f"@{prefix}_{i}\n{read}\n+\n{QUAL}\n")


def main() -> int:
    if shutil.which("mash") is None or shutil.which("make_mashID_db") is None:
        print("mash and make_mashID_db must be on PATH", file=sys.stderr)
        return 1
    rng = random.Random(2026)
    genomes = {acc: random_seq(rng, GENOME_LEN) for acc in ORGANISMS}

    genome_dir = HERE / "genomes"
    reads_dir = HERE / "reads"
    db_dir = HERE / "db"
    for d in (genome_dir, reads_dir, db_dir):
        shutil.rmtree(d, ignore_errors=True)
        d.mkdir()

    for acc, header in ORGANISMS.items():
        write_fasta(genome_dir / f"{acc}.fna", [(header, genomes[acc])])

    # Sketch with relative file names so the database is portable (make_mashID_db stores absolute paths).
    list_file = genome_dir / "list.txt"
    list_file.write_text("".join(f"{acc}.fna\n" for acc in ORGANISMS))
    cmd = ["mash", "sketch", "-p", "1", "-s", "1000", "-k", "21", "-l", "list.txt", "-o", "../db/example"]
    subprocess.run(cmd, cwd=genome_dir, check=True, capture_output=True)
    list_file.unlink()
    subprocess.run(["make_mashID_db", "--annotate", str(db_dir / "example.msh")], check=True, capture_output=True)

    g1, g2, g3 = (genomes[a] for a in ORGANISMS)
    with gzip.open(reads_dir / "alpha_S1_L001_R1_001.fastq.gz", "wt") as fh:
        write_reads(fh, g1, rng, 4000, "alpha_r1")
    with gzip.open(reads_dir / "alpha_S1_L001_R2_001.fastq.gz", "wt") as fh:
        write_reads(fh, g1, rng, 4000, "alpha_r2")
    with gzip.open(reads_dir / "mixed_beta_delta.fastq.gz", "wt") as fh:
        write_reads(fh, g2, rng, 5500, "beta")
        write_reads(fh, g3, rng, 5500, "delta")
    write_fasta(reads_dir / "delta_assembly.fasta.gz",
                [("contig_1 length=70000", g3[:70_000]), ("contig_2 length=50000", g3[70_000:])])
    write_fasta(reads_dir / "unknown.fasta", [("contig_1", random_seq(rng, 30_000))])
    shutil.rmtree(genome_dir)  # the sketch database is what the example needs; genomes are reproducible
    print("Example regenerated. Now run:  bash example/run_example.sh --update-expected")
    return 0


if __name__ == "__main__":
    sys.exit(main())
