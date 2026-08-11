"""Lightweight viewers for sequencing constructs.

The pileup viewer takes a grid of ``(base, is_match)`` cells and writes a
self-contained HTML page.  Input can come in at three levels, each costing only
what it needs:

============================  ==========================  ========================
Level                         Entry point                 Requires
============================  ==========================  ========================
a grid you already have       :func:`render`              nothing
a SAM or BAM                  :func:`grid_from_alignment` pysam
raw reads                     :func:`grid_from_reads`     pysam, minimap2, samtools
============================  ==========================  ========================

The alignment helpers live in :mod:`seqviewer.align` and are imported
lazily, so importing this package never requires pysam::

    from seqviewer import PileupGroup, PileupView, render

    view = PileupView(title="pUC19-WT", groups=[group], flanks=(100, 100))
    Path("pileup.html").write_text(render(view))

Out of scope, deliberately: sequence editing, restriction-site and REBASE
calculation, primer design, ORF finding, implementing aligners, and chromatogram
viewing.
"""

from __future__ import annotations

from .construct import Feature, Reference
from .grid import Cell, Row
from .pileup import PileupGroup, PileupView
from .render import render
from .theme import Theme

__version__ = "0.1.0"

__all__ = [
    "Cell",
    "Feature",
    "PileupGroup",
    "PileupView",
    "Read",
    "Reference",
    "Row",
    "Theme",
    "grid_from_alignment",
    "grid_from_reads",
    "reads_from_alignment",
    "render",
]

_LAZY = {
    "Read": "align",
    "grid_from_alignment": "align",
    "grid_from_reads": "align",
    "reads_from_alignment": "align",
}


def __getattr__(name: str):
    """Import the alignment helpers on first use, so pysam stays optional."""
    module = _LAZY.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from importlib import import_module

    return getattr(import_module(f".{module}", __name__), name)


def __dir__():
    return sorted(__all__)
