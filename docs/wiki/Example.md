# Example

The repository ships a small dataset in
[`example/`](https://github.com/duceppemo/mashID/tree/master/example) that runs offline in a few
seconds and shows every kind of result mashID produces. The organisms are fictional: their genomes are
random 120 kb sequences, so nothing resembles a real species.

```bash
git clone https://github.com/duceppemo/mashID && cd mashID   # if not already cloned
bash example/run_example.sh
```

`run_example.sh` runs mashID on `example/reads` against `example/db/example.msh` and compares the
summary with `example/expected/summary_mashID.tsv`, printing `OK` on a match. To run it by hand:

```bash
mashID -i example/reads -o example_out -d example/db/example.msh
```

## The database

`example/db/example.msh` holds sketches (k=21, s=1000) of three genomes, and
`example/db/example.metadata.tsv` their names, parsed from the fasta headers:

| Accession | Organism |
| --- | --- |
| GCF_000000001.1 | *Exemplaria alpha* |
| GCF_000000002.1 | *Exemplaria beta* subsp. *gamma* |
| GCF_000000003.1 | *Fictivibrio delta* |

## The samples and what they show

| Sample | Files | Result |
| --- | --- | --- |
| `alpha` | `alpha_S1_L001_R1_001.fastq.gz`, `alpha_S1_L001_R2_001.fastq.gz` | Paired-end reads. Both files share the sample name once the Illumina suffixes are stripped, so they are screened together: 8000 reads, 999/1000 hashes of *E. alpha* found, median multiplicity 8 and an estimated depth of 10× (1.2 Mb of reads over a 120 kb reference). |
| `mixed_beta_delta` | `mixed_beta_delta.fastq.gz` | Half the reads come from *E. beta* subsp. *gamma*, half from *F. delta*. The top hit is one of them; the other keeps 996/1000 of its hashes under winner-take-all, so the `Note` reads `Possible mixture with: Fictivibrio delta`. |
| `delta_assembly` | `delta_assembly.fasta.gz` | A two-contig assembly of *F. delta*: 2 sequences, identity ≈ 1, multiplicity 1 as always for assemblies. |
| `unknown` | `unknown.fasta` | A random sequence: nothing reaches the identity threshold, so `No significant hit in database` and `NA` elsewhere. |

## Expected summary

```
Sample            Sequences  Bases    Identity  Shared_Hashes  Median_Multiplicity  Est_Depth  P_Value  Accession        TaxID  Identification                  Note
alpha             8000       1200000  0.999952  999/1000       8                    10.0       0        GCF_000000001.1  NA     Exemplaria alpha
delta_assembly    2          120000   0.999952  999/1000       1                    NA         0        GCF_000000003.1  NA     Fictivibrio delta
mixed_beta_delta  11000      1650000  0.999857  997/1000       6                    13.8       0        GCF_000000002.1  NA     Exemplaria beta subsp. gamma    Possible mixture with: Fictivibrio delta
unknown           1          30000    NA        NA             NA                   NA         NA       NA               NA     No significant hit in database
```

The per-sample table of the mixed sample shows both organisms with almost all of their hashes present,
which is what distinguishes a mixture from two similar references competing for the same k-mers
(see [Outputs](Outputs#the-note-column)):

```
Rank  Identity  Shared_Hashes  Accession        Identification
1     0.999857  997/1000       GCF_000000002.1  Exemplaria beta subsp. gamma
2     0.999809  996/1000       GCF_000000003.1  Fictivibrio delta
```

`TaxID` is `NA` because the names were parsed from headers; a database built with `--metadata` or
`--assembly-report` fills it in (see [Databases](Databases#the-metadata-sidecar)).

## Things to try

```bash
# Rank by multiplicity instead of identity
mashID -i example/reads -o out -d example/db/example.msh -s multiplicity

# Screen only the first 2000 reads of each fastq sample: still identified, lower multiplicity
mashID -i example/reads -o out -d example/db/example.msh --max-reads 2000

# Without winner-take-all the mixture note is not available
mashID -i example/reads -o out -d example/db/example.msh --no-winner-take-all
```

## Regenerating the example

`example/make_example.py` rebuilds the genomes, database, reads and expected tables deterministically
from a fixed seed (needs `mash` and mashID installed):

```bash
python example/make_example.py && bash example/run_example.sh --update-expected
```

A test in the repository (`tests/test_example.py`) runs the example in CI and fails if its results
drift from `example/expected`.
