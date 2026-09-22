# Pipelines

mashID is a single command with predictable outputs, which makes it easy to embed in a workflow.

## Inputs

Give a directory (`-i`) and let file names define samples, or give a
[sample sheet](Usage#sample-sheets) when the pipeline already knows the sample names:

```bash
mashID --sample-sheet samples.tsv -o mashID -d mycobacteriaceae -t 16 -p 4
```

Extra sample-sheet columns are copied into the summary and the JSON, so an expected species from the
LIMS ends up next to the identification.

## Exit codes

| Code | Meaning |
| --- | --- |
| 0 | Completed. |
| 1 | Error: bad input, missing database, Mash failure. |
| 2 | `--fail-on no-hit`: a sample had no hit; `--fail-on note`: a sample had a note or no hit. |
| 130 | Interrupted. |

`--fail-on` lets a pipeline stop, or route samples to review, without parsing tables.

## Outputs to consume

- `summary_mashID.json`: everything, including all reported hits per sample and any sample-sheet
  metadata. With [jq](https://jqlang.github.io/jq/):

  ```bash
  jq -r '.samples[] | [.sample, .top_hit.Identification // "none", (.notes | join("; "))] | @tsv' mashID/summary_mashID.json
  ```

- `summary_mashID.tsv`: one line per sample, stable column order, `NA` for missing values.
- `summary_mashID_mqc.tsv`: picked up automatically by MultiQC, one row per sample with
  identification, identity, depth and notes:

  ```bash
  multiqc results/    # finds every *_mqc.tsv under results/
  ```

- `mashID_run.json`: provenance for reports and audits (versions, parameters, database MD5, inputs).

## Nextflow

```groovy
process MASHID {
    tag "$batch"
    conda 'bioconda::mash=2.3 conda-forge::python=3.12'   // plus: pip install mashID, or bioconda::mashid once available
    publishDir "${params.outdir}/mashID", mode: 'copy'

    input:
    tuple val(batch), path(reads), path(db), path(db_meta)

    output:
    path "mashID/summary_mashID.tsv", emit: summary
    path "mashID/summary_mashID.json", emit: json
    path "mashID/summary_mashID_mqc.tsv", emit: multiqc
    path "mashID/*_mashID.tsv", emit: per_sample
    path "mashID/mashID_run.json", emit: run_info

    script:
    """
    mashID -i . -o mashID -d ${db} -t ${task.cpus} -p 2 --max-reads 500000
    """
}
```

Stage the database's `.metadata.tsv` next to the `.msh` so mashID does not need to regenerate it, or
point `--db-metadata` at it. In a `--sample-sheet` design, write the sheet in the process script from
the tuple's sample names.

## Snakemake

```python
rule mashid:
    input:
        reads=lambda wc: SAMPLES[wc.batch],
        db=config["mashid_db"],
    output:
        summary="mashID/{batch}/summary_mashID.tsv",
        json="mashID/{batch}/summary_mashID.json",
        mqc="mashID/{batch}/summary_mashID_mqc.tsv",
    threads: 8
    conda: "envs/mashid.yaml"      # python>=3.10, mash=2.3, pip: mashID
    shell:
        "mashID -i {input.reads} -o mashID/{wildcards.batch} -d {input.db} -t {threads} -p 2"
```

## Containers

Once the [bioconda recipe](https://github.com/bioconda/bioconda-recipes/pull/69474) is merged, a
Biocontainer image is built automatically and can be used with Docker or Apptainer:

```bash
docker run --rm -v $PWD:/data quay.io/biocontainers/mashid:<tag> mashID -i /data/reads -o /data/out -d /data/db.msh
```

Databases live outside the container; mount them and pass the path with `-d`, or set
`MASHID_DB_DIR` to a mounted directory that holds downloaded databases.

## Throughput

Screening reads every base, so runtime scales with input size. `--max-reads 200000` to `500000`
gives the same identification in a fraction of the time on large runs (see [Usage](Usage#quick-screening-of-large-runs)).
`-p` runs samples concurrently and `-t` is the total thread budget; `-t 32 -p 4` gives each sample
eight Mash threads.
