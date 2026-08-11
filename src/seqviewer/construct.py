"""What a sequencing construct is, independent of how it is viewed.

``Reference`` and ``Feature`` describe a plasmid or an amplicon-with-insert and
are the substrate every viewer in this package shares.  Nothing here knows about
reads, grids, or HTML: a construct is a sequence, a topology, and what is
annotated on it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

__all__ = ["Feature", "Reference"]


@dataclass(frozen=True)
class Feature:
    """One annotated region on a reference, kept as a plain record.

    ``[start, end)`` is a half-open span on the reference as stored.  ``strand``
    is 1, -1, or None.  ``wraps_origin`` marks a feature that crosses base 1 of a
    circular construct; such a feature is stored with ``start > end``, because
    the span from its 5' piece's start to its 3' piece's end is the only pair
    that preserves its real extent.  Read those with
    :func:`seqviewer.genbank.feature_spans`, which returns the one or two
    drawable spans and clips them to a reference length.

    ``color`` is the colour a human chose in SnapGene or ApE, when the file
    carried one.  It is a presentation hint rather than data: a renderer is free
    to adjust it for contrast, and to fall back to its own palette when None.
    """

    type: str
    start: int
    end: int
    strand: Optional[int] = None
    label: Optional[str] = None
    wraps_origin: bool = False
    color: Optional[str] = None

    def __len__(self) -> int:
        """Bases spanned; 0 for an origin-crossing feature, which has no single span.

        Use :func:`seqviewer.genbank.feature_spans` to size a wrapping feature —
        ``end - start`` is negative for one and this clamps it to zero.
        """
        return max(0, self.end - self.start)


@dataclass
class Reference:
    """A construct sequence and what is annotated on it."""

    seq: str
    name: str = ""
    topology: str = "linear"
    features: List[Feature] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.seq)

    @property
    def is_circular(self) -> bool:
        return self.topology == "circular"

    def find(self, feature_type: str) -> Optional[Feature]:
        """Return the first feature of *feature_type*, or None."""
        for feature in self.features:
            if feature.type == feature_type:
                return feature
        return None

    def flank_lengths(self, insert_type: str = "insert") -> Optional[Tuple[int, int]]:
        """Derive ``(5' length, 3' length)`` from the annotated insert.

        Returns None when there is no such feature, when it wraps the origin, or
        when it spans the whole reference and so leaves no flanks to label.
        """
        insert = self.find(insert_type)
        if insert is None or insert.wraps_origin:
            return None
        five, three = insert.start, len(self.seq) - insert.end
        if five <= 0 and three <= 0:
            return None
        return (max(0, five), max(0, three))
