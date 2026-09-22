# Contributing

Thanks for helping improve mashID. Bug reports, database recipes and documentation fixes are as
welcome as code.

## Reporting a bug

Open an [issue](https://github.com/duceppemo/mashID/issues/new/choose) with the exact command, the
log (`--debug`), and the versions of mashID and Mash. If a database misidentifies a sample, the
per-sample `_mashID.tsv` and the sample's origin are the most useful details.

## Changing code

```bash
git clone https://github.com/duceppemo/mashID && cd mashID
conda env create -f environment.yml && conda activate mashID
pip install -e ".[test]"
pre-commit install
```

- Keep `ruff check .` and `pytest` green; add a test for every fix or feature.
- No new runtime dependencies without discussion: mashID is deliberately stdlib-only.
- User-visible changes go in `CHANGELOG.md`; options and outputs are documented in `docs/wiki/`.
- Open the pull request against `master`. CI runs lint and tests on Python 3.10 and 3.12.

## Sharing a database

Pre-built databases are hosted on Figshare and listed in `mashid/databases.py` with their size and
MD5. To propose one, open an issue with the Figshare link, how it was built (source, dereplication,
`-k`/`-s`) and its intended scope.
