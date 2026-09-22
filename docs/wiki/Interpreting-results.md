# Interpreting results

mashID reports, for each sample, the reference whose sketch is best contained in the sample, plus
notes. This page shows what typical situations look like in `summary_mashID.tsv` and the per-sample
table, and what to do about them. Identity values below are for a database with sketch size 1000
or more and k=21; exact numbers vary with the database.

## A clean isolate

```
Sample  Identity  Shared_Hashes  Median_Multiplicity  Est_Depth  Identification            Note
S1      0.9997    994/1000       122                  95.3       Mycobacterium bovis
```

Identity above 0.99, nearly all hashes shared, a depth consistent with the run, and no note. The
per-sample table shows the same organism in the first rows; other organisms, if any, at identities
far below the top hit. This is a species-level identification you can act on.

## Contamination or a mixed culture

```
Sample  Identity  Shared_Hashes  Median_Multiplicity  Identification         Note
S2      0.9995    995/1000       60                   Escherichia coli       Possible mixture with: Staphylococcus aureus
```

Per-sample table:

```
Rank  Identity  Shared_Hashes  Median_Multiplicity  Identification
1     0.9995    995/1000       60                   Escherichia coli
2     0.9991    992/1000       9                    Staphylococcus aureus
```

Two organisms both have almost all of their hashes present. With winner-take-all, shared k-mers go to
the best reference, so the second organism can only score this high with k-mers of its own: it is
in the sample. The multiplicities give the relative abundance, here roughly six to one. Whether this
is a contaminated culture, a polymicrobial sample or index hopping is for you to decide; mashID only
says both are there. `-s multiplicity` ranks by abundance instead of identity when the dominant
organism matters more than the best match.

## Two names for one organism

```
Note: Possible mixture with: Mycobacteroides abscessus
Identification: Mycobacterium abscessus subsp. abscessus
```

If the "mixture" partner is the same species under another genus name, the database lists the
organism under synonyms and mashID counts them as different. Run `make_mashID_db --check` on the
database: it lists species that appear under several genera. Harmonise the names with a
`--metadata` table or rebuild with TaxIDs from NCBI (`--assembly-report`).

## A partial reference on top

```
Identification: Mycobacterium tuberculosis    Shared_Hashes: 943/947
```

A shared-hash denominator below the database's sketch size, or a reference length far below a genome
in the per-sample `Description`, means the top hit is a fragment: a 1 kb record is entirely contained
in any related sample. mashID ignores references shorter than `--min-ref-length` (100 kb) when the
database has a metadata sidecar, which it creates on first use, so this should not happen unless
`--min-ref-length 0` was given. `make_mashID_db --check` lists such records.

## Ambiguous between close species

```
Identity  Identification                  Note
0.9862    Mycobacterium intracellulare    Ambiguous: Mycobacterium paraintracellulare at identity 0.9851
```

The best hits of two organisms are within `--ambiguity-margin` of each other and neither clears the
mixture threshold: the database cannot separate them for this sample. Mash resolves species, and
these two are at the species boundary. Report the complex or group, or use a method with more
resolution (ANI against complete genomes, MLST, SNP typing).

## The *Mycobacterium tuberculosis* complex and other very close groups

MTBC members differ by well under 0.1% of their genome, and Mash identities among them differ in the
fourth decimal. Which member comes first is decided by which reference happens to share one more hash;
`Identification` will alternate between *M. tuberculosis* and *M. tuberculosis* variant *bovis* from
sample to sample. Treat the result as "MTBC" and assign the variant or lineage with a SNP-based tool.
The same applies to *Bacillus cereus* group members, *Escherichia coli* / *Shigella*, and other
groups defined below the species level.

## A novel species or the wrong database

```
Identity  Shared_Hashes  Identification            Note
0.9312    218/1000       Mycobacterium sp. XY-12
```

Identity between 0.90 and 0.97 with a small fraction of shared hashes means the sample is related to
the top hit but is not it: a species absent from the database, or a database for the wrong group.
Screen against a broad database (`refseq_bacteria` or `progenomes3`) to find the genus, then a
dedicated one if it exists. Below 0.90 nothing is reported at all:

```
Identification: No significant hit in database
```

Check that the input is what you think it is (host reads, adapters only, an empty file) and try the
broad databases.

## Low coverage

```
Median_Multiplicity  Est_Depth  Note
2                    2.4        Low coverage: median multiplicity 2
```

The identification can still be right, but with two observations per k-mer a contaminant at a few
percent is invisible and the identity estimate is noisy. `Est_Depth` says how deep the sample is if
it is the reported organism; a large gap between depth and multiplicity (many bases, low
multiplicity) suggests that much of the sample is something else, such as host DNA.

## A quick decision table

| Observation | Read it as |
| --- | --- |
| Identity ≥ 0.99, no note | Species-level identification. |
| `Possible mixture with: X` | X is present too; check multiplicities for proportions. |
| `Ambiguous: X` | Database cannot separate the two; report the group. |
| Identity 0.90–0.97, few shared hashes | Related organism not in the database. |
| `No significant hit` | Not in this database, or not the expected kind of input. |
| `Low coverage` | Plausible but weakly supported; sequence deeper for contaminant detection. |
