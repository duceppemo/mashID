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

### The deployed September sidecars, repaired 2026-10-22

The databases built on 2026-09-22 were annotated before `_tokens` stopped splitting a serovar
designation at its commas, so their sidecars carry names truncated at the first comma — most
consequentially `serovar 1` and `serovar 4` for the monophasic Typhimurium formulae `1,4,[5],12:i:-`
and `4,[5],12:i:-`. A rebuild is not needed to repair that, because the sidecar keeps the
reference's own fasta header verbatim in `Description`: parsing that column again with the fixed
code restores the full name. Done on 2026-10-22 for the deployed copies and their build-directory
twins (backups `*.metadata.tsv.bak-20261022`): 44 names in salmonella, 2 in progenomes4, none in
brucella, escherichia, listeria, mycobacteriaceae or refseq_viral — the last of those stores its
names verbatim from the GenBank `ORGANISM` line and never passes them through the parser.

The `bins.tsv` of those builds needed the corrections but not the repair, and it is worth being
clear why: `species_of` derives the `Bin` column with its own literal regex, deliberately kept
independent of the comment parser so bin names stay byte-stable across versions, and `bins.tsv`
stores the raw NCBI name with its strain text rather than a parsed one. The comma bug could not
reach either column. The corrections were applied on 2026-10-22 (backups `bins.tsv.bak-20261022`):
29 rows in salmonella, 3 in brucella, 2 in escherichia, none in listeria or mycobacteriaceae. A
binning table covers every genome that was binned, not just the references that survived
dereplication, which is why salmonella changes 29 rows here against 16 in its sidecar.

That twin match changes nothing anywhere else, which was checked rather than assumed. Of the 36
rows, 31 are written as GCF and 5 as GCA, because each was recorded with the accession as it
appeared in the database where the mislabel was found — so every other database already matched its
corrections under its own prefix. Audited against the full input sets a rebuild would read, not
just the binned subsets: 2,022 *Brucella*, 21,543 *Salmonella*, 58,198 *Escherichia*/*Shigella*,
16,433 *Mycobacteriaceae* and 32,887 proGenomes 4 assemblies, 131,083 in all, plus every built
table. **Zero** are reachable only via the twin, and none is the same assembly at a different
version. So rebuilding those databases for the twin match would reproduce what is already on disk;
`refseq_bacteria` was the only database it could reach, and the only one never curated against,
having been built in 2023 from RefSeq while the *M. goodii* genome was examined in the
GenBank-based mycobacteriaceae set.

`refseq_bacteria`, which had no sidecar at all, was given one on the same day (82,547 of 82,557
references annotated from NCBI Datasets) and five of the corrections fall inside it. One of those
five is the reason `read_corrections` now registers each row under both accessions of its assembly:
the curation file records the *M. goodii* genome as `GCA_025822535.1` because that is the record we
examined, while `refseq_bacteria` carries it as `GCF_025822535.1`, and the correction used to miss
it in silence. The two accessions are one assembly in GenBank and in RefSeq, so a curated name for
either holds for both. The version is not varied, which would be a different assembly.

The `derep/` tree keeps each genome in a directory named for the bin it was binned as, so a
mislabelled reference stays filed under the wrong species until its file is moved there too — which
the tables alone do not do. Moved on 2026-10-22: 16 genomes in salmonella, 3 in brucella, 1 in
escherichia and the 2 mycobacteriaceae ones corrected the day before, 22 in all, and the
`Source_File` column of each sidecar was repointed at the new location. Only a row whose recorded
path no longer resolved was repointed, and only to a file of the same accession, so the column
cannot be made to name a different genome; afterwards every `Source_File` in all five taxon
databases resolves. Two things the move does not fix, by design:

- `mash info` on the `.msh` still shows the path the genome had when it was sketched. The sketch is
  self-contained and identification never reads these files, so this is a provenance record, not a
  broken reference.
- the destination bin was dereplicated without this genome, so the tree now holds a member that was
  never compared against its new neighbours, and the per-bin `*.log` files still describe the
  original run. The 2026-10-22 rebuild is what resolves that: there the corrected genomes were
  binned with their true relatives from the start, and the near-duplicates among them were dropped.

All four binning tables of the 2026-10-22 rebuilds report zero changes against the same curation
file, and every corrected genome in their `binned/` and `derep/` trees is already in the right
directory, which is the check that `bin_by_species.py --corrections` really did apply all 36 at
build time. Re-running on the repaired September tables also changes nothing, and every row's `Bin` still
agrees with `species_of` of its own `Organism` at that table's rank — the guard against correcting
one row at a rank different from every other.

Two independent checks agree on which rows were truncated. progenomes4 was re-derived from its own
`assembly_data_report.jsonl` instead of from the headers and gave the same two rows. For salmonella,
1,424 of the 1,426 accessions shared with the 2026-10-22 rebuild — which was annotated from a fresh
assembly report — match name for name after the repair; the two that differ are the two curated
corrections, which `apply_label_corrections.py` then applies. Only a change where the stored name is
a prefix of the re-derived one followed by a comma was accepted: a curated annotation is *entitled*
to differ from the fasta header (12 of the 387 brucella rows are `Ochrobactrum anthropi` headers annotated
as `Brucella anthropi`), and overwriting those from the header would undo real curation.

## Upstream status

Corrections are reported to GenBank/RefSeq (`gb-admin@ncbi.nlm.nih.gov`, cc
`genomes@ncbi.nlm.nih.gov`) with the evidence from the `Reason` column. **Sent:**

- **2026-10-22** — `GCA_013697125.1` and `GCA_025822535.1` (mycobacteriaceae).
- **2026-10-22** — the 3 *Brucella* and 15 *Salmonella enterica* records from the sweep below, with
  the MLST-9 profiles and the SISTR table, plus the note that `match_status` already separates them
  and that nine of the *Salmonella* records are internally inconsistent (the deposited subspecies
  contradicts the deposited serovar, which is a subsp. I serovar). The *Shigella* pair and
  `GCF_002266055.1` were explicitly listed as *not* reported, since only negative evidence exists
  for them.

- **2026-10-23** — the 15 records found after that report went out, all of them on the rebuilds of
  2026-10-22: 14 *Salmonella* (11 of the twelve `diarizonae serovar 38:[k]:-` records in the CFSAN
  2016-10-24 submission, `GCF_026132465.1` and `GCF_026132575.1` from the Sichuan 2022-11-10
  submission, and `GCF_016903075.1`, the one deposited under the common label) and `GCF_002290505.1`
  (*Shigella*). Organised by submission rather than by record, because that turned out to be the
  unit the errors occur in, and including the eleven same-submission controls that type exactly as
  deposited — evidence that the discrepancy is not an artefact of our method or of one submitter's
  workflow.

That accounts for 35 of the 36 rows: 2 mycobacteriaceae, 3 *Brucella* and 15 *Salmonella* in the
first report, 14 *Salmonella* and `GCF_002290505.1` in the second. **Nothing is awaiting a send**,
with one loose end worth checking rather than assuming: `GCF_002290495.1`, the other half of the
*Shigella* pair, is named in the 2026-10-23 report only as the cross-referenced partner of
`GCF_002290505.1` ("which we reported to you separately"), while the entry above for 2026-10-22
records that the pair was *excluded* from the first report for resting on negative evidence. Both
cannot be true. If it was never reported in its own right, it is the one record still to send, and
the evidence for it is now positive rather than negative: its partner types as a deposited *S.
sonnei* ST58, and the two cannot both be correctly named.

Keep this section honest: a row in `label_corrections.tsv` is not evidence that anyone upstream has
been told about it, so a new row means this section is out of date until it names that row's report.

The text of the second report is kept at
[`reports/2026-10-23_ncbi_salmonella_shigella.txt`](reports/2026-10-23_ncbi_salmonella_shigella.txt),
with the sender's contact address omitted, so what was claimed upstream stays checkable against
what the evidence actually supports. The first report predates this directory.

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
mycobacteriaceae database, 46 carried a hard `mismatch`; all but the two listed here were
explained by *Mycobacterium*/*Mycolicibacterium* reclassification, reciprocal synonym pairs, or
complex members lacking a type strain. Automated ANI filtering at build time would do more harm
than good.

`below_threshold_match` is a different signal and not a mislabel: the best match *is* the declared
species but below the species threshold, which usually means a divergent or undescribed lineage.
`GCA_001187505.1` (*M. goodii* X7B) is an example — 88.75% ANI to the *M. goodii* type, 0.102 mash
distance versus 0.009-0.012 for the six genuine *goodii*. It is left labelled as deposited, since
*goodii* is the closest named species and no better label exists.

## The 2026-10-22 sweep of the six databases built from NCBI assemblies

Six of the pre-built databases, the ones whose references carry assembly accessions and NCBI names.
`mycobacteriaceae` is not in the table because its two corrections were found earlier, by the same
method; `mycobacteriaceae-2025`, `listeria-2025`, `progenomes3` and `refseq_bacteria` are not,
because at the time they had no metadata sidecar to audit.

Each database's `--confusable` table was computed and then read the way the Databases page
describes, as a mislabel detector: for every reference, how far it sits from its own namesakes
against how far from the nearest genuinely different organism. A correct label is nearest its own
species. The candidates that fell out were then put to a group-specific test, because within all
four groups that produced candidates ANI is uninformative — the classical *Brucella* species are a
single genomospecies, *Shigella* is polyphyletic inside *E. coli*, and NCBI's ANI check resolves
*Salmonella* only to species.

The reference counts are those of the databases as they then stood; the pair counts are from the
confusable tables the sweep computed for them, which were working files and are not kept. The last
column is what the sweep itself found, by auditing those references. Rows reached later — by typing
whole submission batches, and by re-auditing each rebuild — brought the curated totals to 29
*Salmonella* and 2 *Shigella*, as the sections below describe; a sweep of surviving references is
where this started, not the whole of it.

| database | references | pairs ≤ 0.005 | the sweep found |
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

**Where the bad labels concentrate.** All 87 references of the salmonella database not labelled
subsp. *enterica* — 70 carrying another subspecies name, 17 *S. bongori* — were checked with SISTR,
plus 30 *enterica* as a control: 15 of the 102 conclusive calls
are contradicted, and the rate tracks how rare the label is — `arizonae` 7 of 10, `houtenae` 2 of 3,
`indica` 1 of 5, `diarizonae` 3 of 29, `salamae` 2 of 23, `enterica` 0 of 30. A rare subspecies name
has few correct examples to dilute a submission error.

The 2026-10-22 rebuild pulled 410 assemblies newer than the audited database, and the two of them
carrying a rare subspecies name were **both** wrong: `GCF_026132465.1` (deposited *arizonae*, SISTR
subsp. *enterica* serovar Cerro) and `GCF_001832525.1` (deposited *diarizonae* 38:[k]:-, SISTR subsp.
*enterica* serovar Inverness). Each shares its serovar and a near-adjacent accession with a record
already corrected here, so they arrive in submission batches. Re-audit the rare subspecies labels on
every rebuild; an audit of the previous database does not cover the new assemblies.

**Type the whole submission batch, not the surviving references.** Auditing survivors re-discovers
one batch a record at a time, because every correction changes which genomes dereplication keeps:
two of the rows here were found exactly that way, one on each of two later salmonella rebuilds.
Typing both submissions exhaustively — all 27 assemblies making an explicit non-*enterica*
subspecies claim, where only 6 were references in the audited database and 3 of those 6 type as
deposited — ended it and produced 11 further corrections, with no inconclusive calls. The result is sharper than "batches are
unreliable": in the CFSAN 2016-10-24 submission **all twelve** `subsp. diarizonae serovar 38:[k]:-`
records are subsp. *enterica* serovar Inverness, while **all eleven** other non-*enterica* labels in
that same submission (`60:r:e,n,x` → IIIb 60:r:z, `16:z10:e,n` → Redlands, `59:[k]:z35` →
IIIb 38:(k):z35) type exactly as deposited. The error attaches to one serovar designation, not to
the submitter or the date. All four `subsp. arizonae` records in the Sichuan Agricultural University
2022-11-10 submission are subsp. *enterica* (three Cerro, one Apeyeme). None of the 27 carried
`subspecies_match`, which is why a single submission could carry a dozen copies of one wrong label.

**The audit's own blind spot: a common label on a rare subspecies.** Screening rare labels is right
because that is where the *rate* is highest, but it cannot find the reverse error, and 30 subsp.
*enterica* controls out of roughly 3,000 references is not coverage of the common class. Once the
other corrections had collapsed their pairs, exactly one sub-threshold cross-name pair was left in
the salmonella database, and it was this: `GCF_016903075.1`, deposited subsp. *enterica*, is subsp.
*salamae* serovar II 58:d:z8 on 330 alleles — 0.000542 from `GCF_002034565.1`, deposited *and*
independently confirmed as *salamae* II 58:d:z6. Two *salamae* genomes, one named wrongly. The
lesson is the method rather than the record: after correcting the rare labels, re-read the
confusable table, because collapsing the easy pairs is what makes the remaining one visible.

With that row applied the salmonella database converges: **its confusable table holds no cross-name
pair below the 0.001 dereplication threshold at all**, the tightest being serovar-level distinctions
inside subsp. *enterica* (Typhimurium against Typhimurium var. 5-, Yovokome against Manhattan) at
0.00103, above the threshold and mere naming granularity. That emptiness is the stopping condition
worth aiming for, because a collapsed pair is exactly what conceals the next error: while eighteen
such pairs cluttered the table the *salamae* record was indistinguishable from them. Reaching it
took four passes over this database — the cost of *discovering* the 29 corrections, not of applying
them; a rebuild from scratch now needs one pass.

**Mislabels come in submission batches, so check the batch, not the record.** The same thing happened
in `escherichia` on the same rebuild, which is why that row says 2: correcting `GCF_002290495.1` from
*S. flexneri* to *S. sonnei* collapsed its 0.000074 cross-name pair, and the rebuilt database
immediately showed a 0.000036 pair in its place, against `GCF_002290505.1` — the consecutive
accession, same submitter, same release date, same wrong name. One correction had exposed the next.
After adding a row, list its submitter and release date from the assembly report and type the whole
batch; here the other ten members of the JCVI September 2017 *Shigella* batch were spread over nine
unrelated STs and needed nothing, but that is a result, not an assumption. The cheap signal that a
correction is incomplete is a sub-dereplication-threshold cross-name pair surviving a rebuild: the
corrected reference should have no near namesake left under a different name.

**That signal has a benign class, so read it before acting on it.** The rebuilt mycobacteriaceae
database has nine such pairs and not one is a submission error. The tightest is distance **zero**:
`GCA_062013771.1` (*M. intracellulare*) and `GCA_055794805.1` (*M. paraintracellulare*), different
biosamples from different submitters, yet NCBI gives each a `species_match` to its own declared
name at 99.11% and 99.47% — which can only happen if the two type assemblies are themselves within
the species threshold of each other. That is a synonym conflict in the taxonomy, not a wrong label,
and the same goes for *M. timonense*/*M. bouchedurhonense* at 0.000041, the
*marinum*/*pseudoshottsii*/*liflandii* trio, and *M. mungi* against *M. tuberculosis* variant
*africanum*. The distinguishing question is not how close the pair is but whether each member
matches *its own* type: in the *Shigella* case the suspect's nearest same-name reference was 150
times further away than its nearest differently-named one, and MLST put it in a sequence type
attested only among the other name. Where both names are genuinely defensible, the `--confusable`
table doing its job — warning that the call could be either — is the correct outcome, not a
correction row.

Verification of the mycobacteriaceae rebuild, for the record: 3,728 references against 3,534
deployed, 200 added and 6 dropped. Both corrected genomes disappeared entirely rather than being
renamed, because once binned with their true relatives they were near-duplicates and dereplication
dropped them — the deployed database had the right *names* after `apply_label_corrections.py`, but
still carried the redundant sequences, and only a rebuild fixes the FASTA headers that script
cannot touch. Scored against the same assembly report, the deployed reference set carries 46 hard
`mismatch` records and the rebuilt set 44; the difference is exactly those two accessions, and none
of the 200 new references brought a hard mismatch with it.

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
`GCF_002290485.1` (*S. boydii*) and `GCF_002317875.1` (*S. sonnei*) share ST877 and all seven
alleles, but sit 0.003349 apart and each is nearer an *E. coli* reference than the other, so a
shared seven-locus profile at that distance is ordinary *Shigella*/*E. coli* polyphyly rather than
the signature of a mislabel — contrast the 0.000036 of the pair above. Two more members of that
batch carry STs typical of *E. coli* rather than their deposited species (`GCF_002317065.1`
*S. sonnei* ST10, `GCF_002317935.1` *S. boydii* ST48), which is negative evidence only and does not
meet the bar.

**`refseq_viral` has no mislabels.** Its 20 distance-zero pairs are 7 RefSeq `AC_`/`NC_` curated
duplicates of one virus under an old and a current ICTV name (`AC_000007.1` "Human adenovirus 2" and
`NC_001405.1` "Human mastadenovirus C") and 13 phage records named per isolate rather than per
species (`NC_024360.1` "Listeria phage LMSP-25" and `NC_042048.1` "Listeria phage LMTA-34"). Both
are naming granularity, not error, and neither is fixed by a correction row.
