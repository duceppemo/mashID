"""Extract an organism name from a Mash sketch comment (usually an NCBI fasta header)."""

from __future__ import annotations

import re

# Mash prefixes the comment of multi-sequence files with "[N seqs] " and suffixes it with " [...]".
_SEQ_COUNT_PREFIX = re.compile(r"^\[\d+\s+seqs?\]\s*")
_TRAILING_ELLIPSIS = re.compile(r"\s*\[\.\.\.\]\s*$")

# A genus is a capitalised Latin word, optionally in brackets (e.g. "[Clostridium]").
_GENUS = re.compile(r"^\[?[A-Z][a-z]{2,}\]?$")
# A species epithet is lowercase (possibly hyphenated) or an abbreviation such as "sp." / "spp.".
_SPECIES = re.compile(r"^(?:[a-z][a-z-]+|spp?\.?)$")

# Capitalised English words that commonly precede a lowercase word in fasta headers but are not genera.
_NOT_A_GENUS = {
    "Complete", "Whole", "Draft", "Genome", "Chromosome", "Plasmid", "Strain", "Contig",
    "Scaffold", "Unnamed", "Sequence", "Isolate", "Sample", "Unknown", "Uncultured", "Bacterium",
    "Assembly", "Reference", "Partial", "Linear", "Circular", "Segment", "Clone", "Node",
    "Human", "Mouse", "Bovine", "Marine", "Soil", "Gut", "Metagenome", "Bacteria", "Archaea",
    "Virus", "Phage", "Environmental", "Synthetic", "Hypothetical",
}

# Infraspecific rank markers -> label used in the reported name.
_RANK_LABEL = {
    "subsp.": "subsp.", "subsp": "subsp.", "subspecies": "subsp.", "ssp.": "subsp.", "ssp": "subsp.",
    "variant": "variant", "var.": "var.", "var": "var.",
    "serovar": "serovar", "sv.": "serovar", "serotype": "serotype",
    "biovar": "biovar", "bv.": "biovar",
    "pathovar": "pv.", "pv.": "pv.",
    "genomovar": "genomovar", "morphovar": "morphovar",
}

# Tokens after which infraspecific information is no longer expected.
_STOP_WORDS = {
    "strain", "str.", "isolate", "chromosome", "plasmid", "contig", "scaffold", "genome",
    "sequence", "complete", "whole", "draft", "dna", "segment", "clone", "node", "unnamed",
}


def clean_comment(comment: str) -> str:
    """Remove the "[N seqs]" prefix and "[...]" suffix that Mash adds to comments."""
    text = _SEQ_COUNT_PREFIX.sub("", comment or "")
    return _TRAILING_ELLIPSIS.sub("", text).strip()


def _tokens(text: str) -> list[str]:
    # Commas and semicolons carry no taxonomic meaning; drop them so "441," becomes "441".
    return [t for t in text.replace(",", " ").replace(";", " ").split() if t]


def _find_binomial(tokens: list[str]) -> tuple[int, str, str] | None:
    """Return (index of genus token, genus, species) or None."""
    for i, tok in enumerate(tokens[:-1]):
        nxt = tokens[i + 1]
        if tok == "Candidatus" and i + 2 < len(tokens):
            if _GENUS.match(nxt) and _SPECIES.match(tokens[i + 2]):
                return i, f"Candidatus {nxt}", tokens[i + 2]
            continue
        if tok in _NOT_A_GENUS or not _GENUS.match(tok):
            continue
        if _SPECIES.match(nxt):
            return i, tok, nxt
    return None


def organism_from_comment(comment: str, fallback: str = "") -> str:
    """Build a short organism name (genus, species, infraspecific ranks) from a sketch comment.

    Examples of accepted inputs::

        "[153 seqs] NZ_LZJQ01000001.1 Mycobacterium mantenii strain E2660 contig_1, ... [...]"
        "NC_002945.4 Mycobacterium tuberculosis variant bovis AF2122/97 chromosome, complete genome"
        "NZ_FVSX01000015.1 Mycobacteroides abscessus subsp. massiliense strain 441, ..."

    When no binomial can be recognised, ``fallback`` is returned if given, otherwise the cleaned
    comment (or "unknown" if the comment is empty). This function never raises on odd input.
    """
    text = clean_comment(comment)
    tokens = _tokens(text)
    found = _find_binomial(tokens)
    if found is None:
        return fallback or text or "unknown"

    genus_idx, genus, species = found
    n_consumed = 3 if genus.startswith("Candidatus ") else 2
    rest = tokens[genus_idx + n_consumed:]
    parts = [genus, species]

    # "Mycobacterium sp. JS623": keep the designation that follows "sp.".
    if re.fullmatch(r"spp?\.?", species) and rest and rest[0].lower() not in _STOP_WORDS \
            and rest[0].lower() not in _RANK_LABEL:
        parts.append(rest[0])
        rest = rest[1:]

    # Collect infraspecific ranks until a stop word is met: "subsp. enterica serovar Typhimurium".
    i = 0
    while i < len(rest) - 1:
        tok = rest[i]
        low = tok.lower()
        if low in _STOP_WORDS:
            break
        if low in _RANK_LABEL:
            epithet = rest[i + 1]
            if epithet.lower() not in _STOP_WORDS and epithet.lower() not in _RANK_LABEL:
                parts.extend([_RANK_LABEL[low], epithet])
                i += 2
                continue
        i += 1

    return " ".join(parts)
