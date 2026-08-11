"""What one rendered pileup page shows.

These are the types :func:`seqviewer.render` consumes, and the only ones in the
package specific to the pileup viewer.  The pileup's older vector/insert split
survives as ``PileupView.flanks``, a plain ``(5' length, 3' length)`` pair;
:meth:`PileupView.from_reference` derives that pair from a reference's annotated
insert, which is the seam along which the hard-coded three-region model gives way
to arbitrary features.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

from .construct import Reference
from .grid import Row
from .theme import Theme

__all__ = ["PileupGroup", "PileupView"]


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

    ``flanks`` remains the focus region: the one span whose edges are worth
    dashed boundary lines through the full height of the pileup.  It is None
    when no feature was singled out, and then the page draws no boundaries.

    ``translate`` decides whether that region also gets amino-acid rows.  It is
    a separate question from whether a focus region exists — a vector's payload
    may be a span worth marking and not a reading frame worth translating — so
    set it False for a construct where a protein readout would be meaningless.
    Nothing is translated when there is no focus region to translate.
    """

    title: str
    groups: List[PileupGroup] = field(default_factory=list)
    total_reads: int = 0
    highlight_ids: List[str] = field(default_factory=list)
    highlight_label: str = "Highlighted"
    flanks: Optional[Tuple[int, int]] = None
    features: List[Feature] = field(default_factory=list)
    ref_len: Optional[int] = None
    translate: bool = True
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
