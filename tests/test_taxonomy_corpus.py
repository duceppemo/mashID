"""Regression corpus: 3309 real NCBI headers from the 2025-02-20 Mycobacteriaceae database."""

import gzip
import re
from pathlib import Path

from mashid.taxonomy import organism_from_comment

CORPUS = Path(__file__).resolve().parent / "data" / "ncbi_headers.txt.gz"
BINOMIAL = re.compile(r"^(Candidatus )?\[?[A-Z][a-z]+\]? (?:[a-z][a-z-]+|sp\.)")


def test_every_ncbi_header_yields_a_binomial():
    with gzip.open(CORPUS, "rt") as fh:
        headers = [ln.rstrip("\n") for ln in fh if ln.strip() and not ln.startswith("#")]
    assert len(headers) == 3309
    failures = [(h, organism_from_comment(h)) for h in headers if not BINOMIAL.match(organism_from_comment(h))]
    assert failures == []
    organisms = {organism_from_comment(h) for h in headers}
    assert 600 <= len(organisms) <= 700  # 624 at the time of writing; a big change means a parser regression
    assert "Mycobacterium tuberculosis variant bovis" in organisms
    assert "Mycobacteroides abscessus subsp. massiliense" in organisms
    assert "Candidatus Mycobacterium methanotrophicum" in organisms
