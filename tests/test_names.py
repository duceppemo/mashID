"""Tests for mashid.names: when two organism names are a real conflict, and the shared readers."""

import importlib.util
from pathlib import Path

import pytest

from mashid import MashIDError
from mashid.names import (
    binomial,
    distinguishable,
    normalise_header,
    open_table,
    read_corrections,
    species_of,
    taxon_tokens,
)

SCRIPTS = Path(__file__).resolve().parent.parent / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# --- names read as a chain of ranks -------------------------------------------------------------

@pytest.mark.parametrize("a,b", [
    # Two species.
    ("Mycolicibacterium goodii", "Mycolicibacterium wolinskyi"),
    ("Mycobacterium canetti", "Mycobacterium tuberculosis"),
    ("Mycobacterium intracellulare", "Mycobacterium paraintracellulare"),
    # NCBI's "variant" stands where it used to say species, so these really are different
    # organisms -- and this is the pair the whole feature exists for.
    ("Mycobacterium tuberculosis", "Mycobacterium tuberculosis variant bovis"),
    ("Mycobacterium tuberculosis variant bovis", "Mycobacterium tuberculosis variant africanum"),
    # Siblings below a species disagree with each other even though each agrees with the parent.
    ("Mycobacteroides abscessus subsp. abscessus", "Mycobacteroides abscessus subsp. massiliense"),
    ("Mycobacterium avium subsp. paratuberculosis", "Mycobacterium avium subsp. hominissuis"),
    ("Salmonella enterica subsp. enterica serovar Dublin",
     "Salmonella enterica subsp. enterica serovar Typhimurium"),
])
def test_different_organisms_conflict(a, b):
    assert distinguishable(a, b) and distinguishable(b, a)


@pytest.mark.parametrize("a,b", [
    # One name is the other at a finer rank: a sample that is the subspecies is also the species,
    # so pairing the two states nothing. This was two fifths of the shipped mycobacteriaceae table.
    ("Mycobacterium avium", "Mycobacterium avium subsp. hominissuis"),
    ("Mycobacteroides abscessus", "Mycobacteroides abscessus subsp. abscessus"),
    ("Mycobacterium intracellulare", "Mycobacterium intracellulare subsp. chimaera"),
    ("Salmonella enterica", "Salmonella enterica subsp. enterica serovar Dublin"),
    ("Escherichia coli", "Escherichia coli O157:H7"),             # a strain designation
    ("Mycobacterium avium", "Mycobacterium avium ssp. avium"),    # the other spelling of subsp.
    # Notation rather than taxonomy.
    ("Mycobacterium abscessus", "Mycobacterium abscessus"),
    ("[Clostridium] difficile", "Clostridium difficile"),
    ("Candidatus Mycobacterium lepraemurium", "Mycobacterium lepraemurium"),
    ("Mycobacterium Abscessus", "Mycobacterium abscessus"),
    # A renamed genus is one organism under two spellings.
    ("Mycobacterium abscessus", "Mycobacteroides abscessus"),
    ("Mycobacteroides chelonae", "[Mycobacterium] chelonae"),
    ("Mycobacterium abscessus subsp. bolletii", "Mycobacteroides abscessus subsp. bolletii"),
    # No species claimed, so nothing to contradict.
    ("Mycobacterium pinniadriaticum", "Mycobacterium sp. P3"),
    ("Mycobacterium sp.", "Mycobacterium fortuitum"),
    ("Mycobacterium bovis", "NA"),
    ("Mycobacterium bovis", "unknown"),
    ("Mycobacterium bovis", ""),
    ("Lausannevirus", "Mycobacterium bovis"),
])
def test_same_or_unusable_names_do_not_conflict(a, b):
    assert not distinguishable(a, b) and not distinguishable(b, a)


def test_two_strains_of_one_species_still_disagree_with_each_other():
    """Only a name against its own ancestor is vacuous. Two things named below one species are a
    real distinction, formal rank or not -- and for E. coli that is the distinction that matters."""
    assert distinguishable("Escherichia coli K-12", "Escherichia coli O157:H7")
    assert distinguishable("Escherichia coli str. K-12", "Escherichia coli strain O157:H7")
    assert not distinguishable("Escherichia coli strain K-12", "Escherichia coli str. K-12")
    assert not distinguishable("Escherichia coli", "Escherichia coli str. K-12")


def test_a_renamed_genus_is_still_judged_on_what_follows_the_epithet():
    # Strip the genus difference and the ranks below still decide.
    assert distinguishable("Mycobacterium abscessus subsp. abscessus",
                           "Mycobacteroides abscessus subsp. massiliense")
    assert not distinguishable("Mycobacterium abscessus subsp. abscessus",
                               "Mycobacteroides abscessus")


def test_taxon_tokens_normalises_notation_and_rank_spellings():
    assert taxon_tokens("[Mycobacterium] Chelonae") == ("mycobacterium", "chelonae")
    assert taxon_tokens("Candidatus Mycobacterium x") == ("mycobacterium", "x")
    assert taxon_tokens("M. avium ssp. avium") == taxon_tokens("M. avium subsp. avium")
    assert taxon_tokens("Genus species var. x") == taxon_tokens("Genus species variant x")
    assert taxon_tokens("") == ()


def test_binomial_keeps_the_case_it_was_given():
    assert binomial("[Mycobacterium] chelonae") == ("Mycobacterium", "chelonae")
    assert binomial("Candidatus Mycobacterium lepraemurium") == ("Mycobacterium", "lepraemurium")
    assert binomial("Mycobacterium sp.") is None
    assert binomial("Lausannevirus") is None
    assert binomial("") is None


# --- bin names: one implementation, shared with the curation scripts ----------------------------

PARITY_NAMES = [
    "Mycobacterium tuberculosis variant bovis",
    "Mycolicibacterium wolinskyi",
    "Mycobacterium sp.",
    "Mycobacterium sp",
    "[Mycobacterium] manitobense",
    "Candidatus Mycobacterium something",
    "Mycobacterium avium subsp. paratuberculosis",
    "Salmonella enterica subsp. enterica serovar Dublin",
    "Escherichia coli O157:H7",
    "Unnamed",
    "",
]


@pytest.mark.parametrize("rank", ["species", "subspecies", "serovar"])
@pytest.mark.parametrize("organism", PARITY_NAMES)
def test_bin_names_agree_across_the_scripts_at_every_rank(organism, rank):
    """The old divergence only showed below species, which is where salmonella is binned."""
    binner = _load("bin_by_species")
    applier = _load("apply_label_corrections")
    expected = species_of(organism, rank)
    assert binner.species_of(organism, rank) == expected
    assert applier.bin_of(organism, rank) == expected


def test_species_of_keeps_the_rank_it_is_asked_for():
    name = "Salmonella enterica subsp. enterica serovar Dublin"
    assert species_of(name) == "Salmonella_enterica"
    assert species_of(name, "subspecies") == "Salmonella_enterica_subsp_enterica"
    assert species_of(name, "serovar") == "Salmonella_enterica_subsp_enterica_serovar_Dublin"


def test_species_of_does_not_treat_unknown_as_a_placeholder():
    """Unlike distinguishable(): renaming existing bins would be a breaking change."""
    assert species_of("Genus unknown") == "Genus_unknown"
    assert species_of("Genus sp.") == "Genus_sp"


# --- the shared table readers -------------------------------------------------------------------

def test_normalise_header_folds_case_and_punctuation():
    assert normalise_header(["Assembly Accession", "organism_name", "TaxID"]) == {
        "assemblyaccession": "Assembly Accession", "organismname": "organism_name", "taxid": "TaxID"}
    assert normalise_header(None) == {}


def test_open_table_handles_a_bom_and_sniffs_the_delimiter(tmp_path):
    tsv = tmp_path / "a.tsv"
    tsv.write_bytes("﻿Accession\tOrganism\nGCA_1.1\tGenusa one\n".encode())
    assert next(iter(open_table(tsv)))["Accession"] == "GCA_1.1"
    csv_path = tmp_path / "a.csv"
    csv_path.write_text("Accession,Organism\nGCA_1.1,Genusa one\n")
    assert next(iter(open_table(csv_path)))["Organism"] == "Genusa one"


@pytest.mark.parametrize("char", ["\x0b", "\x0c", "\x1c", "\x1e", "\x85", " "])
def test_open_table_only_breaks_rows_on_a_newline(tmp_path, char):
    """str.splitlines() breaks on all of these, and any can reach a Description from a fasta
    header: the row would end early and its tail would read as a second, phantom reference."""
    p = tmp_path / "a.tsv"
    p.write_text(f"Accession\tOrganism\tDescription\nGCA_1.1\tGenusa one\tCP1{char}odd\n")
    rows = list(open_table(p))
    assert len(rows) == 1
    assert rows[0]["Description"] == f"CP1{char}odd"


def test_open_table_tolerates_crlf_and_a_missing_final_newline(tmp_path):
    p = tmp_path / "a.tsv"
    p.write_text("Accession\tOrganism\r\nGCA_1.1\tGenusa one\r\nGCA_2.1\tGenusa two")
    assert [r["Organism"] for r in open_table(p)] == ["Genusa one", "Genusa two"]


def test_open_table_rejects_bad_bytes_and_empty_files(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_bytes(b"Accession\tOrganism\nGCA_1.1\tCaf\xe9 one\n")
    with pytest.raises(MashIDError, match="not a UTF-8 text table"):
        open_table(bad)
    empty = tmp_path / "empty.tsv"
    empty.write_text("")
    with pytest.raises(MashIDError, match="is empty"):
        open_table(empty)


def test_read_corrections_reports_a_missing_file(tmp_path):
    with pytest.raises(MashIDError, match="Corrections file not found"):
        read_corrections(tmp_path / "nope.tsv")


def test_read_corrections_accepts_a_bom_and_extra_columns(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_bytes("﻿Accession\tOrganism\tTaxID\tReason\nGCA_1.1\tGenusb two\t222\twhy\n".encode())
    assert read_corrections(p) == {"GCA_1.1": ("Genusb two", "222"),
                                   "GCF_1.1": ("Genusb two", "222")}  # the assembly's RefSeq twin
