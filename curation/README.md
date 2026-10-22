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

It is idempotent and backs each target up. This does not touch the `.msh`, so the reference's own
FASTA header — the `Description` column mash reports — keeps the deposited name until a rebuild.
`Identification` and `TaxID` are correct either way.

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
