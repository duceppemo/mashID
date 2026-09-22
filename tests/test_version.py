"""pyproject.toml, mashid.__version__ and CITATION.cff must agree."""

import re
from pathlib import Path

import mashid

ROOT = Path(__file__).resolve().parent.parent


def test_versions_agree():
    pyproject = (ROOT / "pyproject.toml").read_text()
    assert re.search(r'^version = "([^"]+)"$', pyproject, re.M).group(1) == mashid.__version__
    citation = (ROOT / "CITATION.cff").read_text()
    assert re.search(r"^version: (\S+)$", citation, re.M).group(1) == mashid.__version__
