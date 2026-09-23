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
