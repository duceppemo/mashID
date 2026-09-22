# Development

```bash
git clone https://github.com/duceppemo/mashID
cd mashID
conda env create -f environment.yml && conda activate mashID
pip install -e ".[test]"
pre-commit install          # optional: ruff on every commit
ruff check .
pytest
```

The integration tests build a small synthetic database and are skipped when `mash` is not on `PATH`.

## Layout

| Path | Content |
| --- | --- |
| `mashid/cli.py`, `makedb.py`, `download_cli.py` | The three command-line entry points. |
| `mashid/pipeline.py` | Sample discovery → stats → screening → annotation → tables. |
| `mashid/mash.py` | Wrappers around `mash screen`, `mash sketch`, `mash info`. |
| `mashid/samples.py` | Input discovery and sample naming. |
| `mashid/seqstats.py` | Read/base counting and the `--max-reads` streamer. |
| `mashid/taxonomy.py` | Organism name parsing from fasta headers. |
| `mashid/metadata.py` | The metadata sidecar. |
| `mashid/databases.py` | Registry and download of pre-built databases. |
| `docs/wiki/` | Source of this wiki, published by `.github/workflows/wiki.yml`. |
| `recipe/` | Bioconda recipe. |

## Continuous integration

`.github/workflows/ci.yml` runs on every push and pull request: environment from `environment.yml`
on Linux with Python 3.10 and 3.12 and on macOS with Python 3.12, `ruff check`, `pytest`, and the
CLIs' `--help`. `tests/test_taxonomy_corpus.py` checks the header parser against 3309 real NCBI
headers, and `tests/test_example.py` keeps the shipped example's results exact.

## Releasing

1. Update the version in `mashid/__init__.py`, `pyproject.toml` and `CITATION.cff` (a test checks they agree) and add a section to `CHANGELOG.md`. Run `pytest` and check its exit status before tagging.
2. Commit, then tag: `git tag -a vX.Y.Z -m "mashID X.Y.Z" && git push origin master vX.Y.Z`.
3. `.github/workflows/release.yml` builds the wheel and creates the GitHub release with the
   changelog section as notes. Zenodo archives the release automatically and mints a version DOI
   under the concept DOI [10.5281/zenodo.22888109](https://doi.org/10.5281/zenodo.22888109); put the
   new version DOI in `CITATION.cff` (`doi:` and `identifiers`) in the next commit.
4. Update `recipe/meta.yaml`: `version`, `sha256` of the new tag's tarball
   (`curl -sL https://github.com/duceppemo/mashID/archive/refs/tags/vX.Y.Z.tar.gz | sha256sum`), and
   `number: 0`. Copy `recipe/` to `bioconda-recipes/recipes/mashid` and open a pull request; once its
   checks pass, comment `@BiocondaBot please add label`.

## Editing this wiki

Edit the Markdown files in `docs/wiki/` and push to `master`. The `wiki.yml` workflow copies them to
the wiki repository. Page names are file names (`Databases.md` → `Databases`); `_Sidebar.md` is the
navigation.
