"""Standard genetic code translation, without a Biopython dependency.

Handles IUPAC ambiguity the way Biopython does: an ambiguous codon translates
to a single residue when every base it could stand for gives the same residue
(``ACN`` is threonine however the third position resolves), to one of the
ambiguous residue codes when the possibilities form a recognized pair, and to
``X`` otherwise.  This matters here because gaps in a consensus are written as
``N`` before translation, so ambiguous codons are routine rather than exotic.

Verified against ``Bio.Seq.Seq.translate()`` over all 4,096 codons drawn from
the full IUPAC alphabet.
"""

from __future__ import annotations

__all__ = ["translate_codon", "translate"]

_BASES = "TCAG"
_RESIDUES = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)

#: Unambiguous codon -> residue, standard genetic code (NCBI table 1).
CODON_TABLE = {
    a + b + c: _RESIDUES[i]
    for i, (a, b, c) in enumerate(
        (a, b, c) for a in _BASES for b in _BASES for c in _BASES
    )
}

#: IUPAC nucleotide code -> the unambiguous bases it stands for.
IUPAC_BASES = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

#: Residue pairs that have their own single-letter ambiguity code.
_AMBIGUOUS_RESIDUES = {
    frozenset("DN"): "B",   # aspartate or asparagine
    frozenset("EQ"): "Z",   # glutamate or glutamine
    frozenset("LI"): "J",   # leucine or isoleucine
}


def translate_codon(codon: str) -> str:
    """Translate a single three-base codon to a one-letter residue.

    Returns ``X`` for anything that cannot be resolved to one residue or to a
    recognized ambiguous residue code, including codons of the wrong length or
    containing characters outside the IUPAC alphabet.
    """
    codon = codon.upper().replace("U", "T")
    residue = CODON_TABLE.get(codon)
    if residue is not None:
        return residue
    if len(codon) != 3:
        return "X"

    try:
        options = [IUPAC_BASES[base] for base in codon]
    except KeyError:
        return "X"

    possible = {
        CODON_TABLE[a + b + c]
        for a in options[0] for b in options[1] for c in options[2]
    }
    if len(possible) == 1:
        return possible.pop()
    return _AMBIGUOUS_RESIDUES.get(frozenset(possible), "X")


def translate(seq: str) -> str:
    """Translate a nucleotide sequence in frame 0, ignoring any trailing partial codon."""
    end = len(seq) - len(seq) % 3
    return "".join(translate_codon(seq[i:i + 3]) for i in range(0, end, 3))
