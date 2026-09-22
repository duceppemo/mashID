# Example dataset

A small, self-contained dataset to check that mashID works, and to see each kind of result. The
organisms are fictional: their genomes are random sequences, so nothing here resembles a real species.

```bash
bash example/run_example.sh
# or, by hand:
mashID -i example/reads -o example_out -d example/db/example.msh
```

| Sample | Files | What it shows |
| --- | --- | --- |
| `alpha` | `alpha_S1_L001_R1_001.fastq.gz`, `..._R2_001.fastq.gz` | Paired-end reads, grouped into one sample and screened together. |
| `delta_assembly` | `delta_assembly.fasta.gz` | An assembly (two contigs): multiplicity 1, identity ≈ 1. |
| `mixed_beta_delta` | `mixed_beta_delta.fastq.gz` | Reads from two organisms: the `Note` column reports the mixture. |
| `unknown` | `unknown.fasta` | A sequence matching nothing in the database. |

`db/example.msh` holds sketches (k=21, s=1000) of three 120 kb genomes and `db/example.metadata.tsv`
their names. `expected/` holds the tables mashID produces; `run_example.sh` compares them.

`make_example.py` regenerates everything deterministically (needs `mash` and mashID installed):

```bash
python example/make_example.py && bash example/run_example.sh --update-expected
```
