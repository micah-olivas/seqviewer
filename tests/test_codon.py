import pytest

from seqviewer.codon import CODON_TABLE, translate, translate_codon


def test_table_is_complete():
    assert len(CODON_TABLE) == 64


def test_standard_code():
    assert translate_codon("ATG") == "M"
    assert translate_codon("TGG") == "W"
    assert translate_codon("TAA") == translate_codon("TAG") == translate_codon("TGA") == "*"


def test_case_and_rna_are_accepted():
    assert translate_codon("atg") == "M"
    assert translate_codon("AUG") == "M"


@pytest.mark.parametrize("codon,residue", [
    ("ACN", "T"),   # every ACx codon is threonine
    ("CGN", "R"),
    ("GTN", "V"),
    ("CCN", "P"),
])
def test_fourfold_degenerate_codons_resolve(codon, residue):
    """A gap becomes N before translation, so these are routine, not exotic."""
    assert translate_codon(codon) == residue


@pytest.mark.parametrize("codon,residue", [
    ("RAY", "B"),   # aspartate or asparagine
    ("SAR", "Z"),   # glutamate or glutamine
    ("MTA", "J"),   # leucine or isoleucine
])
def test_recognized_ambiguous_residues(codon, residue):
    assert translate_codon(codon) == residue


@pytest.mark.parametrize("codon", ["NNN", "AT", "ATGC", "AT!", ""])
def test_unresolvable_codons_are_X(codon):
    assert translate_codon(codon) == "X"


def test_translate_ignores_trailing_partial_codon():
    assert translate("ATGATG") == "MM"
    assert translate("ATGATGA") == "MM"
    assert translate("AT") == ""


def test_translate_matches_biopython_across_iupac():
    """The reason this module exists is to drop the Biopython dependency."""
    Seq = pytest.importorskip("Bio.Seq").Seq
    from itertools import product

    alphabet = "ACGTUNRYSWKMBDHV"
    for codon in ("".join(p) for p in product(alphabet, repeat=3)):
        assert translate_codon(codon) == str(Seq(codon).translate()), codon
