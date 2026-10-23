"""Organism names read as taxonomic claims, plus the tables they arrive in.

Shared by the package and the scripts in ``scripts/`` so that one answer to "do these two names
mean different organisms?" serves screening, database building and curation alike. There were three
divergent copies of this logic before, and they did not agree.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Iterable
from pathlib import Path

from mashid import MashIDError

# An unnamed reference makes no claim to a species, so it cannot contradict one.
PLACEHOLDER_SPECIES = frozenset({"sp.", "sp", "spp.", "bacterium", "unknown", "na"})

# NCBI spells the same rank several ways; normalise before comparing.
RANK_SYNONYMS = {
    "ssp.": "subsp.", "ssp": "subsp.", "subsp": "subsp.",
    "serotype": "serovar", "sv.": "serovar",
    "pv.": "pathovar",
    "var.": "variant",
    "strain": "str.",
}

# "subsp. X", "serovar X" and a bare strain designation all subdivide the species above them: a
# sample that is the finer one is also, correctly, the coarser one, so the two names do not disagree
# and there is nothing to warn about.
#
# "variant" is the exception, because NCBI uses it where it formerly used species: M. bovis,
# M. africanum, M. microti and M. caprae are now "Mycobacterium tuberculosis variant <x>". Those are
# names a clinician still reads as different organisms -- and telling an M. bovis sample that this
# database cannot separate it from M. tuberculosis is the most useful thing the confusable table has
# to say about the Mycobacteriaceae.
SEPARATE_ORGANISM_RANKS = frozenset({"variant"})

# Deliberately *not* PLACEHOLDER_SPECIES: species_of() names directories in databases that are
# already built, so adding "unknown" here would rename existing bins. Keep the two sets apart.
_BIN_PLACEHOLDERS = ("sp.", "sp", "spp.", "bacterium")

BIN_RANKS = ("species", "subspecies", "serovar")


def normalise_header(fieldnames: list[str] | None) -> dict[str, str]:
    """Map lowercase, punctuation-free column names to the actual header names."""
    return {re.sub(r"[^a-z0-9]", "", f.lower()): f for f in (fieldnames or [])}


def open_table(path: Path) -> csv.DictReader:
    """A ``DictReader`` over a TSV or CSV, tolerating a BOM and sniffing the delimiter."""
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8-sig")
    except UnicodeDecodeError as exc:
        raise MashIDError(f"{path} is not a UTF-8 text table: {exc}") from exc
    # Split on newlines only. str.splitlines() also breaks on \r, \x0b, \x0c, \x1c-\x1e, \x85 and
    # U+2028, any of which can reach a Description field from a fasta header, where they would end
    # the row early and turn its tail into a phantom reference.
    lines = [line.rstrip("\r") for line in text.split("\n")]
    if lines and lines[-1] == "":
        lines.pop()
    if not lines or not lines[0].strip():
        raise MashIDError(f"{path} is empty")
    delimiter = "\t" if "\t" in lines[0] else ","
    return csv.DictReader(lines, delimiter=delimiter)


def taxon_tokens(organism: str) -> tuple[str, ...]:
    """An organism name as a comparable chain of taxonomic tokens, coarsest first.

    Brackets (NCBI's ``[Mycobacterium] chelonae``), a ``Candidatus`` prefix, letter case and the
    several spellings of each rank word are normalised away, so only real differences remain.
    """
    tokens = organism.replace("[", "").replace("]", "").split()
    if tokens[:1] == ["Candidatus"]:
        tokens = tokens[1:]
    return tuple(RANK_SYNONYMS.get(t.casefold(), t.casefold()) for t in tokens)


def binomial(organism: str) -> tuple[str, str] | None:
    """``(genus, species epithet)`` as written, or None when the name claims no species.

    Brackets and a ``Candidatus`` prefix are stripped, but letter case is kept so that the pair can
    be shown back to the user; compare with ``taxon_tokens()``, which normalises for comparison.
    """
    tokens = organism.replace("[", "").replace("]", "").split()
    if tokens[:1] == ["Candidatus"]:
        tokens = tokens[1:]
    if len(tokens) < 2:
        return None
    if tokens[0].casefold() in PLACEHOLDER_SPECIES or tokens[1].casefold() in PLACEHOLDER_SPECIES:
        return None
    return tokens[0], tokens[1]


def distinguishable(a: str, b: str) -> bool:
    """True when ``a`` and ``b`` are claims about genuinely different organisms.

    False when the two names make the same claim; when either claims no species (an unnamed
    ``Genus sp.`` reference contradicts nothing); when they differ only by a renamed genus
    (*Mycobacterium* / *Mycobacteroides abscessus*); and when one is the other at a finer rank
    (*Mycobacterium avium* and *M. avium* subsp. *hominissuis*), which is a redundant statement
    rather than a conflict. ``SEPARATE_ORGANISM_RANKS`` records the one rank word that escapes the
    last rule.

    The asymmetry is deliberate: a name against its own ancestor says nothing, but *any* two named
    things below one species disagree with each other and are worth reporting, whether they are
    divided at a formal rank (*abscessus* subsp. *abscessus* against subsp. *massiliense*) or only
    by strain (*E. coli* K-12 against O157:H7, where "not separable" is the useful answer).
    """
    if binomial(a) is None or binomial(b) is None:
        return False
    path_a, path_b = taxon_tokens(a), taxon_tokens(b)
    if path_a[1] == path_b[1] and path_a[0] != path_b[0]:
        # The same epithet under two genera is a renamed genus; judge on what follows it.
        path_a, path_b = path_a[1:], path_b[1:]
    if path_a == path_b:
        return False
    shorter, longer = sorted((path_a, path_b), key=len)
    if longer[: len(shorter)] == shorter:
        return longer[len(shorter)] in SEPARATE_ORGANISM_RANKS
    return True


def species_of(organism: str, rank: str = "species") -> str:
    """A filesystem-safe bin name for an organism, at ``species``, ``subspecies`` or ``serovar``.

    This names directories and the ``Bin`` column of a built database, so its output must stay
    byte-stable across versions; it is deliberately more literal than ``taxon_tokens()``.
    """
    tokens = organism.replace("[", "").replace("]", "").split()
    if not tokens:
        return "unknown"
    if tokens[0] == "Candidatus" and len(tokens) >= 3:
        tokens = ["Candidatus_" + tokens[1]] + tokens[2:]
    genus = tokens[0]
    species = tokens[1] if len(tokens) > 1 else "sp"
    if species in _BIN_PLACEHOLDERS:
        species = "sp"
    name = f"{genus}_{species}"
    if rank in ("subspecies", "serovar"):
        m = re.search(r"subsp\.? (\S+)", organism)
        if m:
            name += f"_subsp_{m.group(1)}"
    if rank == "serovar":
        m = re.search(r"serovar (\S+)", organism)
        if m:
            name += f"_serovar_{m.group(1)}"
    return re.sub(r"[^A-Za-z0-9_.:,\[\]-]", "_", name)


def paired_accession(accession: str) -> str | None:
    """The same assembly's other accession: GCA_000195955.2 <-> GCF_000195955.2, else None."""
    prefix = accession[:4]
    if prefix == "GCA_":
        return "GCF_" + accession[4:]
    if prefix == "GCF_":
        return "GCA_" + accession[4:]
    return None


def with_paired_accessions(mapping: dict[str, tuple[str, str]]) -> dict[str, tuple[str, str]]:
    """Register every entry under both accessions of its assembly, keeping explicit keys.

    One assembly has a GenBank and a RefSeq accession, so a curated organism name given for either
    holds for both, and one curation file can name references in a database built from either
    source. Anything that reads a curated table has to agree about this: pairing in one reader and
    not another gives a build that bins a reference as corrected and then names it as deposited.
    """
    out = dict(mapping)
    for accession, value in mapping.items():
        twin = paired_accession(accession)
        if twin:
            out.setdefault(twin, value)
    return out


def assemblies(accessions: Iterable[str]) -> set[str]:
    """Collapse each GenBank/RefSeq pair to one assembly, for counting and reporting.

    A paired mapping holds two keys per assembly, which would otherwise double a "corrections
    loaded" count and report one unmatched row twice, once under each accession.
    """
    return {min(a, paired_accession(a) or a) for a in accessions}


def read_corrections(path: Path) -> dict[str, tuple[str, str]]:
    """Read a TSV/CSV of curated labels into ``{accession: (organism, taxid)}``.

    Needs columns Accession and Organism; TaxID is optional, and every other column is ignored so
    that the same file can carry the evidence for each correction.
    """
    path = Path(path)
    if not path.is_file():
        raise MashIDError(f"Corrections file not found: {path}")
    reader = open_table(path)
    cols = normalise_header(reader.fieldnames)
    acc_col = cols.get("accession") or cols.get("assemblyaccession") or cols.get("acc")
    org_col = (cols.get("organism") or cols.get("organismname") or cols.get("species")
               or cols.get("name"))
    if not acc_col or not org_col:
        raise MashIDError(f"{path}: expected columns 'Accession' and 'Organism' (optionally "
                          f"'TaxID'); found {reader.fieldnames}")
    tax_col = cols.get("taxid") or cols.get("taxonomyid")
    out: dict[str, tuple[str, str]] = {}
    for row in reader:
        acc = (row.get(acc_col) or "").strip()
        organism = (row.get(org_col) or "").strip()
        if acc and organism:
            out[acc] = (organism, (row.get(tax_col) or "").strip() if tax_col else "")
    # Every row also answers for the assembly's other accession; an explicit row always wins over
    # the twin it is registered under. The version is never varied: that would be a different
    # assembly.
    return with_paired_accessions(out)
