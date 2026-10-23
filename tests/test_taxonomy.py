import pytest

from mashid.taxonomy import clean_comment, organism_from_comment


@pytest.mark.parametrize(
    "comment, expected",
    [
        ("[153 seqs] NZ_LZJQ01000001.1 Mycobacterium mantenii strain E2660 contig_1, whole genome shotgun sequence [...]",
         "Mycobacterium mantenii"),
        ("NZ_CP060409.1 Mycolicibacterium fortuitum strain W4 chromosome", "Mycolicibacterium fortuitum"),
        ("NC_002945.4 Mycobacterium tuberculosis variant bovis AF2122/97 chromosome, complete genome",
         "Mycobacterium tuberculosis variant bovis"),
        ("[15 seqs] NZ_FVSX01000015.1 Mycobacteroides abscessus subsp. massiliense strain 441, whole genome shotgun sequence [...]",
         "Mycobacteroides abscessus subsp. massiliense"),
        ("NC_015848.1 Mycobacterium canettii CIPT 140010059, complete sequence", "Mycobacterium canettii"),
        ("NZ_X.1 Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome",
         "Salmonella enterica subsp. enterica serovar Typhimurium"),
        ("NZ_X.1 Mycobacterium sp. JS623 chromosome, complete genome", "Mycobacterium sp. JS623"),
        ("NZ_X.1 Candidatus Mycobacterium methanotrophicum isolate X", "Candidatus Mycobacterium methanotrophicum"),
        ("NZ_X.1 [Clostridium] difficile strain 630", "[Clostridium] difficile"),
        ("NZ_X.1 Escherichia coli O157:H7 str. Sakai DNA, complete genome", "Escherichia coli"),
        # organism first, no accession (custom databases)
        ("Mycobacterium bovis AF2122/97", "Mycobacterium bovis"),
        # "variant" after a stop word is ignored
        ("NZ_X.1 Mycobacterium bovis strain X plasmid variant 2", "Mycobacterium bovis"),
        # trailing "sub" in unrelated words must not be treated as a subspecies marker
        ("NZ_X.1 Mycobacterium bovis substrain BCG Pasteur", "Mycobacterium bovis"),
        # a comma inside a serovar designation is part of the name
        ("NZ_X.1 Salmonella enterica subsp. arizonae serovar 18:z4,z23:- str. CVM N23850",
         "Salmonella enterica subsp. arizonae serovar 18:z4,z23:-"),
        ("NZ_X.1 Salmonella enterica subsp. enterica serovar 1,4,[5],12:i:- chromosome, complete genome",
         "Salmonella enterica subsp. enterica serovar 1,4,[5],12:i:-"),
        # ... but a header that omits the space after the comma is still an ordinary word break,
        # which proGenomes 3 does ("ORS278,complete sequence")
        ("114615.SAMEA3138227.CU234118 Bradyrhizobium sp. ORS278,complete sequence",
         "Bradyrhizobium sp. ORS278"),
        ("NZ_X.1 Mycobacterium mantenii strain E2660 chromosome,complete genome", "Mycobacterium mantenii"),
        # a space-less comma is cut at the stop word, so a formula before it keeps its own commas
        ("NZ_X.1 Salmonella enterica subsp. enterica serovar 1,4,[5],12:i:-,complete genome",
         "Salmonella enterica subsp. enterica serovar 1,4,[5],12:i:-"),
        # semicolons separate words in a header exactly as commas do
        ("114615.SAMEA3138227.CU234118 Bradyrhizobium sp. ORS278;complete sequence",
         "Bradyrhizobium sp. ORS278"),
        # a doubled or mixed separator must not be left hanging on the kept part
        ("114615.SAMEA3138227.CU234118 Bradyrhizobium sp. ORS278,,complete sequence",
         "Bradyrhizobium sp. ORS278"),
        ("114615.SAMEA3138227.CU234118 Bradyrhizobium sp. ORS278,;complete sequence",
         "Bradyrhizobium sp. ORS278"),
        # a stop word opening the token still stops the name, rather than becoming part of it
        ("NZ_X.1 Mycobacterium bovis genome,ORS278,complete", "Mycobacterium bovis"),
    ],
)
def test_organism_from_comment(comment, expected):
    assert organism_from_comment(comment) == expected


def test_no_binomial_falls_back():
    assert organism_from_comment("contig_1 length=1234", fallback="GCF_1") == "GCF_1"
    assert organism_from_comment("contig_1 length=1234") == "contig_1 length=1234"
    assert organism_from_comment("") == "unknown"
    assert organism_from_comment("", fallback="ACC") == "ACC"


def test_never_raises_on_odd_input():
    for odd in ["[", "]", "[3 seqs]", "Mycobacterium", "Mycobacterium ", "  ", "A b", "Complete genome"]:
        assert isinstance(organism_from_comment(odd), str)


def test_clean_comment():
    assert clean_comment("[3 seqs] NZ_1 Genus species [...]") == "NZ_1 Genus species"
    assert clean_comment("NZ_1 Genus species") == "NZ_1 Genus species"


def test_english_words_are_not_genera():
    assert organism_from_comment("Human gut metagenome assembly", fallback="X") == "X"
    assert organism_from_comment("Marine sediment bacterium", fallback="X") == "X"
