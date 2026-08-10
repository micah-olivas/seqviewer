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
    is 1, -1, or None.  ``wraps_origin`` marks a feature stored in several parts
    because it crosses base 1 of a circular construct, which means its real
    extent is not ``[start, end)`` and callers should refuse it rather than read
    it as one stretch.
    """

    type: str
    start: int
    end: int
    strand: Optional[int] = None
    label: Optional[str] = None
    wraps_origin: bool = False

    def __len__(self) -> int:
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
    """

    title: str
    groups: List[PileupGroup] = field(default_factory=list)
    total_reads: int = 0
    top_fraction: float = 0.0
    highlight_ids: List[str] = field(default_factory=list)
    highlight_label: str = "Highlighted"
    flanks: Optional[Tuple[int, int]] = None
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
        """Build a view whose flanks come from *reference*'s annotated insert."""
        return cls(
            title=title,
            groups=list(groups),
            flanks=reference.flank_lengths(insert_type),
            **kwargs,
        )
