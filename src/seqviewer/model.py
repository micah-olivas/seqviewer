"""The data a construct view is built from.

Two layers live here.  ``Reference`` and ``Feature`` describe a construct — a
plasmid or an amplicon-with-insert — and are the substrate both this package's
viewers are meant to share.  ``PileupGroup`` and ``PileupView`` describe one
rendered pileup page and are what :func:`seqviewer.render` consumes.

The pileup's older vector/insert split survives as ``PileupView.flanks``, a
plain ``(5' length, 3' length)`` pair.  :meth:`PileupView.from_reference`
derives that pair from a reference's annotated insert, which is the seam along
which the hard-coded three-region model gives way to arbitrary features.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

__all__ = [
    "Feature", "Reference", "PileupGroup", "PileupView", "Theme", "Cell", "Row",
]

#: One position of one read: the base observed, and whether it matched the
#: reference there.  A base of ``"-"`` is a gap or an uncovered position.
Cell = Tuple[str, bool]

#: One read across the full width of the reference.
Row = Sequence[Cell]


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


@dataclass(frozen=True)
class Theme:
    """Names the page uses to bridge into a host application's theme.

    A page written with the defaults is self-contained.  An application that
    already stores a light/dark preference passes its own names so the pileup
    follows the same setting as the rest of its output.
    """

    storage_key: str = "seqviewer-theme"
    css_prefix: str = "cv"
    style_id: str = "cv-theme-bridge"
    script_id: str = "cv-theme-sync"


@dataclass
class PileupGroup:
    """One reference and the reads shown against it.

    A page carries one group per subpopulation, ordered by the caller.  ``rows``
    is the grid: one entry per read, each as wide as ``ref_seq``.
    """

    name: str
    ref_seq: str
    rows: List[Row] = field(default_factory=list)
    n_reads: int = 0
    fraction: float = 0.0
    status: str = ""
    highlighted: bool = False

    def __post_init__(self) -> None:
        for i, row in enumerate(self.rows):
            if len(row) != len(self.ref_seq):
                raise ValueError(
                    f"group {self.name!r}: row {i} is {len(row)} cells wide but "
                    f"the reference is {len(self.ref_seq)} bases"
                )


@dataclass
class PileupView:
    """Everything one pileup page shows.

    ``highlight_ids`` are group names called out in the page header, and
    ``highlight_label`` is the word used to introduce them — "Recoverable" in a
    streak-out report, something else elsewhere.  Groups whose ``highlighted``
    flag is set get a star next to their name.

    ``features`` are drawn as an annotation track over the reference bar.
    ``ref_len`` is the length they are stated against, which is what lets a page
    whose groups have different reference lengths drop the annotations that do
    not apply rather than drawing them off the end.

    ``flanks`` remains the focus region: the one feature whose edges are worth
    dashed boundary lines through the full height of the pileup, and the window
    the amino-acid track translates.  It is None when no feature was singled
    out, and then the page draws neither — an annotated plasmid with no
    designated insert gets its features and nothing else.
    """

    title: str
    groups: List[PileupGroup] = field(default_factory=list)
    total_reads: int = 0
    top_fraction: float = 0.0
    highlight_ids: List[str] = field(default_factory=list)
    highlight_label: str = "Highlighted"
    flanks: Optional[Tuple[int, int]] = None
    features: List[Feature] = field(default_factory=list)
    ref_len: Optional[int] = None
    theme: Theme = field(default_factory=Theme)

    @classmethod
    def from_reference(
        cls,
        title: str,
        reference: Reference,
        groups: Sequence[PileupGroup],
        insert_type: str = "insert",
        **kwargs,
    ) -> "PileupView":
        """Build a view from *reference*: its features, length, and focus region.

        The focus region comes from a feature of *insert_type*; when the
        reference has none, ``flanks`` is None and the page draws no boundary
        lines and no translation.
        """
        return cls(
            title=title,
            groups=list(groups),
            flanks=reference.flank_lengths(insert_type),
            features=list(reference.features),
            ref_len=len(reference),
            **kwargs,
        )
