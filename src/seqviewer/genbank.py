"""Read an annotated construct file into a :class:`~seqviewer.construct.Reference`.

GenBank, SnapGene ``.dna``, and ApE ``.gb`` files all carry the annotations a
pileup wants to draw over its reference bar.  Biopython does the parsing; this
module's work is turning its objects into the package's own flat ``Feature``
records, which means three things Biopython will not do for you:

* deciding whether a multi-part location genuinely crosses the origin, which
  ``len(location.parts) > 1`` gets wrong for any spliced feature,
* picking one label out of a handful of competing qualifiers, and
* recovering the colour a human chose in SnapGene or ApE, which each tool
  stores in a different place and one of them hides inside ``/note``.

Biopython is imported inside the loader, so importing this module — or the
package — never requires it.  Only annotated formats do.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .construct import Feature, Reference

__all__ = ["load_reference", "DEFAULT_SKIP_TYPES", "feature_spans"]

#: Suffixes Biopython can read, and the format name it wants for each.
_ANNOTATED = {
    ".dna": "snapgene",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".genbank": "genbank",
    ".ape": "genbank",
}

#: Types dropped before anything is drawn.
#:
#: ``source`` spans the whole record, so it neighbours everything and costs a
#: lane to say nothing.  ``primer_bind`` is the flood risk: a SnapGene-annotated
#: vector carries many of them, each a few dozen pixels wide and each with a
#: label far wider than itself.  ``gene`` is handled separately — it is dropped
#: only where a CDS already covers the same span, because a record that
#: annotates genes and no CDS still needs its open reading frames.
DEFAULT_SKIP_TYPES = frozenset({"source", "primer_bind"})

#: Qualifiers searched in order for a feature's display name.  This follows
#: dna-features-viewer's ``label_fields``, minus ``source`` (that is the
#: organism, not a name) and plus ``standard_name``.
_LABEL_QUALIFIERS = (
    "label", "gene", "product", "standard_name", "locus_tag", "note",
)

#: SnapGene round-trips its own styling through ``/note``, so a note reading
#: "color: #ffd281; direction: BOTH" is markup rather than a name.  The same
#: applies to its multi-segment description.
_NOTE_IS_MARKUP = re.compile(
    r"^\s*(color\s*:|direction\s*:|This\s.*feature\shas\s\d+\ssegments)",
    re.IGNORECASE,
)

#: ``/note="color: #ffd281; direction: BOTH"`` — SnapGene's convention.
_NOTE_COLOR = re.compile(r"\bcolor\s*:\s*(#[0-9a-fA-F]{3,6}|[a-zA-Z]+)")

_HEX = re.compile(r"^#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})$")

#: ApE writes colour names as well as hexes into its qualifiers.
_NAMED_COLORS = {
    "cyan": "#00ffff", "green": "#00ff00", "red": "#ff0000",
    "blue": "#0000ff", "yellow": "#ffff00", "magenta": "#ff00ff",
    "orange": "#ffa500", "pink": "#ffc0cb", "purple": "#800080",
    "white": "#ffffff", "black": "#000000", "gray50": "#7f7f7f",
    "cornflower blue": "#6495ed",
}

#: Longest label kept; a whole abstract pasted into /note is not a name.
_MAX_LABEL = 40


def _normalise_color(raw: str) -> Optional[str]:
    """Return *raw* as ``#rrggbb``/``#rgb``, or None if it isn't a colour."""
    value = raw.strip().lower()
    if _HEX.match(value):
        return value
    return _NAMED_COLORS.get(value)


def _feature_color(qualifiers: Dict[str, Sequence[str]],
                   strand: Optional[int]) -> Optional[str]:
    """Recover the colour a human chose, from wherever their tool put it.

    ``/color`` is the plain convention, ApE uses ``/ApEinfo_fwdcolor`` and
    ``/ApEinfo_revcolor``, and SnapGene hides its value inside ``/note``.
    """
    direct = ["color", "Color", "colour"]
    if strand == -1:
        direct = ["ApEinfo_revcolor"] + direct + ["ApEinfo_fwdcolor"]
    else:
        direct = ["ApEinfo_fwdcolor"] + direct + ["ApEinfo_revcolor"]
    for key in direct:
        for value in qualifiers.get(key, ()):
            found = _normalise_color(value)
            if found:
                return found
    for note in qualifiers.get("note", ()):
        match = _NOTE_COLOR.search(note)
        if match:
            found = _normalise_color(match.group(1))
            if found:
                return found
    return None


def _feature_label(qualifiers: Dict[str, Sequence[str]], fallback: str) -> str:
    """Pick a display name, skipping notes that are really styling markup."""
    for key in _LABEL_QUALIFIERS:
        for value in qualifiers.get(key, ()):
            text = value.strip()
            if not text:
                continue
            if key == "note" and _NOTE_IS_MARKUP.match(text):
                continue
            text = " ".join(text.split())
            if len(text) > _MAX_LABEL:
                text = text[: _MAX_LABEL - 1].rstrip() + "…"
            return text
    return fallback


def _wraps_origin(location, seq_len: int, circular: bool) -> bool:
    """Decide whether *location* genuinely crosses base 1.

    ``len(location.parts) > 1`` is not that test: a spliced CDS such as
    ``join(20..25,35..45)`` has two parts and stays put.  A feature crosses the
    origin only on a circular construct and only when its pieces reach both ends
    of the sequence.
    """
    parts = list(location.parts)
    if len(parts) < 2 or not circular:
        return False
    spans = sorted((int(p.start), int(p.end)) for p in parts)
    return spans[0][0] == 0 and spans[-1][1] == seq_len


def _wrapped_bounds(location) -> Tuple[int, int]:
    """Return ``(start, end)`` for an origin-crossing location.

    ``start > end`` on purpose.  Biopython's ``CompoundLocation.start``/``.end``
    give the min and max over the parts, which for a wrapping feature is the
    whole sequence — the hull says nothing about the real extent.  Encoding the
    5' piece's start and the 3' piece's end keeps it, and a reader that sees
    ``start > end`` knows to draw two spans.

    Minus-strand parts arrive in biological rather than reference order, so the
    pieces are sorted by position first.
    """
    spans = sorted((int(p.start), int(p.end)) for p in location.parts)
    return spans[-1][0], spans[0][1]


def feature_spans(feature: Feature, ref_len: int) -> List[Tuple[int, int]]:
    """Return the drawable ``[start, end)`` spans of *feature*, clipped to *ref_len*.

    One span normally; two when the feature crosses the origin, which is encoded
    as ``start > end``.  Features lying wholly outside the reference return no
    spans at all, which is how a view whose groups have different reference
    lengths drops annotations that do not apply.
    """
    if feature.wraps_origin or feature.start > feature.end:
        pieces = [(feature.start, ref_len), (0, feature.end)]
    else:
        pieces = [(feature.start, feature.end)]
    out = []
    for start, end in pieces:
        start, end = max(0, start), min(ref_len, end)
        if end > start:
            out.append((start, end))
    return out


def _drop_redundant_genes(features: List[Feature]) -> List[Feature]:
    """Drop ``gene`` features a CDS already covers.

    Plasmid files routinely annotate both, and drawing both puts two glyphs on
    one open reading frame.  A ``gene`` with no matching CDS is kept, because a
    record that annotates genes and nothing else still has reading frames worth
    showing.
    """
    cds_spans = {(f.start, f.end) for f in features if f.type == "CDS"}
    if not cds_spans:
        return features
    return [
        f for f in features
        if not (f.type == "gene" and (f.start, f.end) in cds_spans)
    ]


def _read_fasta(path: Path) -> Reference:
    """Read the first record of a FASTA file, which carries no annotations."""
    name, chunks = path.stem, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if chunks:
                break
            name = line[1:].strip().split()[0] or path.stem
        else:
            chunks.append(line.strip())
    return Reference(seq="".join(chunks).upper(), name=name)


def load_reference(path, skip_types: Optional[Sequence[str]] = None) -> Reference:
    """Load *path* as a :class:`Reference`, with its annotations.

    FASTA is read directly.  GenBank, SnapGene, and ApE files go through
    Biopython, which is imported here rather than at module scope so that an
    unannotated workflow never needs it installed.

    *skip_types* replaces :data:`DEFAULT_SKIP_TYPES`; pass ``()`` to keep
    everything the file declares.
    """
    path = Path(path)
    fmt = _ANNOTATED.get(path.suffix.lower())
    if fmt is None:
        return _read_fasta(path)

    from Bio import SeqIO  # only annotated formats need biopython

    record = next(SeqIO.parse(str(path), fmt))
    seq = str(record.seq).upper()
    topology = record.annotations.get("topology", "linear")
    circular = topology == "circular"
    skip = frozenset(DEFAULT_SKIP_TYPES if skip_types is None else skip_types)

    features: List[Feature] = []
    for raw in record.features:
        if raw.type in skip:
            continue
        location = raw.location
        strand = location.strand
        if _wraps_origin(location, len(seq), circular):
            start, end = _wrapped_bounds(location)
            wraps = True
        else:
            start, end = int(location.start), int(location.end)
            wraps = False
        features.append(Feature(
            type=raw.type,
            start=start,
            end=end,
            strand=strand,
            label=_feature_label(raw.qualifiers, raw.type),
            wraps_origin=wraps,
            color=_feature_color(raw.qualifiers, strand),
        ))

    return Reference(
        seq=seq,
        name=record.name if record.name and record.name != "unknown" else path.stem,
        topology=topology,
        features=_drop_redundant_genes(features),
    )
