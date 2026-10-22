# Curated label corrections

`label_corrections.tsv` holds reference labels that are wrong at the source (NCBI) and must be
overridden. `Organism`/`TaxID` are the corrected values; `Deposited_Organism`/`Deposited_TaxID`
record what NCBI says, so every correction is reversible and self-documenting. `Reason` carries
the evidence. Only `Accession`, `Organism` and `TaxID` are read — the rest are ignored by both
consumers below.

## On rebuild

`make_mashID_db` applies `--metadata` *after* the NCBI assembly report, so entries here win over
the deposited name. Pass it to the binner too, otherwise a mislabelled genome is dereplicated
against the wrong species and can survive as a near-duplicate of the cluster it really belongs to:

```bash
bin_by_species.py assembly_data_report.jsonl ncbi_dataset/data binned/ \
    --corrections curation/label_corrections.tsv
make_mashID_db -i binned/ -o /path/to/db -p my_database \
    --assembly-report assembly_data_report.jsonl \
    --metadata curation/label_corrections.tsv
```

(`-o` is the output directory and `-p` the file name without `.msh`.)

## On an already-built database

To fix reporting without rebuilding, rewrite the tables in place (`--dry-run` first):

```bash
scripts/apply_label_corrections.py curation/label_corrections.tsv \
    ~/.local/share/mashID/db/<db>.metadata.tsv <build_dir>/bins.tsv
```

It is idempotent and backs each target up. Pass `--rank` if the `bins.tsv` was built below species
(`bin_by_species.py --rank subspecies|serovar`), or the corrected row's bin is written at species
level while every other row keeps its rank, and the next rebuild dereplicates against the wrong
cluster. A `<db>.confusable.tsv` can be given as a target too: it names organisms in two columns and
both are corrected, since the second is the one the `Note` is built from.

This does not touch the `.msh`, so the reference's own FASTA header — the `Description` column mash
reports — keeps the deposited name until a rebuild. `Identification` and `TaxID` are correct either
way. A correction that matches no row of any target is reported on stderr rather than passing
silently.

## Upstream status

Every correction here has been reported to GenBank/RefSeq (`gb-admin@ncbi.nlm.nih.gov`, cc
`genomes@ncbi.nlm.nih.gov`) with the evidence from its `Reason` column, so none needs filing again:

- **2026-10-22** — `GCA_013697125.1` and `GCA_025822535.1` (mycobacteriaceae).
- **2026-10-22** — the 3 *Brucella* and 15 *Salmonella enterica* records from the sweep below, with
  the MLST-9 profiles and the SISTR table, plus the note that `match_status` already separates them
  and that nine of the *Salmonella* records are internally inconsistent (the deposited subspecies
  contradicts the deposited serovar, which is a subsp. I serovar). The *Shigella* pair and
  `GCF_002266055.1` were explicitly listed as *not* reported, since only negative evidence exists
  for them.

NCBI maintains the ANI flag but not the records themselves, so the organism names change only if the
submitting groups agree; a report is not a fix.

Check whether anything has changed upstream before trusting the `Checked` date in a `Reason`:

```bash
while read -r a; do
  curl -s "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$a/dataset_report" \
    | python3 -c "import json,sys; r=json.load(sys.stdin)['reports'][0]; \
        print(r['accession'], '|', r['organism']['organism_name'], '|', \
              (r.get('average_nucleotide_identity') or {}).get('match_status'))"
done < <(cut -f1 curation/label_corrections.tsv | tail -n +2)
```

If a name is corrected at NCBI, its row here becomes a no-op rather than a conflict — a rebuild
takes the corrected name from the assembly report and `apply_label_corrections.py` reports zero
rows changed — so the row can simply be dropped.

## When to add a row

Only with positive evidence of the true identity, not merely an NCBI ANI flag. Within the
*M. tuberculosis* complex every member sits above 98% ANI and several species have no type
assembly, so `taxonomy_check_status: Inconclusive` and `below_threshold_mismatch` are normal
there and are not grounds for relabelling on their own. Of 3,534 references in the
mycobacteriaceae database, 38 carried a hard `mismatch`; all but the two listed here were
explained by *Mycobacterium*/*Mycolicibacterium* reclassification, reciprocal synonym pairs, or
complex members lacking a type strain. Automated ANI filtering at build time would do more harm
than good.

`below_threshold_match` is a different signal and not a mislabel: the best match *is* the declared
species but below the species threshold, which usually means a divergent or undescribed lineage.
`GCA_001187505.1` (*M. goodii* X7B) is an example — 88.75% ANI to the *M. goodii* type, 0.102 mash
distance versus 0.009-0.012 for the six genuine *goodii*. It is left labelled as deposited, since
*goodii* is the closest named species and no better label exists.

## The 2026-10-22 sweep of every pre-built database

Each database's `--confusable` table was computed and then read the way the Databases page
describes, as a mislabel detector: for every reference, how far it sits from its own namesakes
against how far from the nearest genuinely different organism. A correct label is nearest its own
species. The candidates that fell out were then put to a group-specific test, because within all
four groups that produced candidates ANI is uninformative — the classical *Brucella* species are a
single genomospecies, *Shigella* is polyphyletic inside *E. coli*, and NCBI's ANI check resolves
*Salmonella* only to species.

| database | references | pairs ≤ 0.005 | corrections found |
| --- | --- | --- | --- |
| `listeria` | 1441 | 0 | none |
| `progenomes4` | 32887 | 0 | none (one representative per species cluster) |
| `refseq_viral` | 19625 | 765 over 421 refs | none; see below |
| `brucella` | 387 | 1772 over 147 refs | 3 |
| `salmonella` | 1442 | 803 over 267 refs | 15 |
| `escherichia` | 4503 | 2356 over 1631 refs | 1 |

The decisive tests were MLST-9 for *Brucella* (`mlst --scheme brucella`), the Achtman *E. coli*
scheme for *Shigella*, and SISTR cgMLST for *Salmonella*, each run against controls drawn from the
same database so that the test was shown to separate the species before it was trusted. MLST also
settled which member of a near-identical pair was wrong, which distance alone cannot do.

**Where the bad labels concentrate.** All 87 non-*enterica* subspecies labels in the salmonella
database were checked with SISTR, plus 30 *enterica* as a control: 15 of the 102 conclusive calls
are contradicted, and the rate tracks how rare the label is — `arizonae` 7 of 10, `houtenae` 2 of 3,
`indica` 1 of 5, `diarizonae` 3 of 29, `salamae` 2 of 23, `enterica` 0 of 30. A rare subspecies name
has few correct examples to dilute a submission error.

**What `match_status` is worth.** For those 102 references, `subspecies_match` was never wrong (0 of
13), so it is a sound positive signal. The converse is not: 55 of the 70 conclusive labels carrying
`species_match` or `derived_species_match` are nonetheless correct, so demanding
`subspecies_match` at build time would discard 55 good subspecies labels to catch 15 bad ones.
Treat anything short of `subspecies_match` as unverified and check it — for *Salmonella*, SISTR over
the subspecies-labelled references is cheap — rather than as grounds to downgrade the label.

**Deliberately not corrected.** `GCF_003977345.1` (*S. boydii*) and `GCF_003977375.1` (*S. sonnei*)
are 0.000096 apart and share ST405. *S. sonnei* is dominated by ST152 (8 of 14 sampled) and no
*sonnei* control is ST405, so the *sonnei* label is the implausible one; but *S. boydii* spans 11 STs
in 14 genomes, so ST405 cannot be positively assigned and neither row meets the bar above.
`GCF_002266055.1` (*arizonae*) matched only 85 cgMLST alleles, too few to call.

**`refseq_viral` has no mislabels.** Its 20 distance-zero pairs are 7 RefSeq `AC_`/`NC_` curated
duplicates of one virus under an old and a current ICTV name (`AC_000007.1` "Human adenovirus 2" and
`NC_001405.1` "Human mastadenovirus C") and 13 phage records named per isolate rather than per
species (`NC_024360.1` "Listeria phage LMSP-25" and `NC_042048.1` "Listeria phage LMTA-34"). Both
are naming granularity, not error, and neither is fixed by a correction row.
