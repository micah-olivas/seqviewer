"""The primitives a pileup grid is made of.

These live apart from both the construct substrate and the pileup view types
because they are the seam between them: :mod:`seqviewer.align` produces rows and
:mod:`seqviewer.pileup` consumes them, and neither needs to know about the other.
They are plain type aliases, so importing this module costs nothing.
"""

from __future__ import annotations

from typing import Sequence, Tuple

__all__ = ["Cell", "Row"]

#: One position of one read: the base observed, and whether it matched the
#: reference there.  A base of ``"-"`` is a gap or an uncovered position.
Cell = Tuple[str, bool]

#: One read across the full width of the reference.
Row = Sequence[Cell]
