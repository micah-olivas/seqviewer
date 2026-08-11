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

The same view renders two ways.  :func:`render` draws every read; for a page
that fits on a screen, :class:`~seqviewer.summary.SummaryView` reduces the reads
to coverage and called variants and :func:`render_summary` draws that as an
annotated map — features, a lollipop per variant, and a coverage profile::

    from seqviewer import SummaryView, render_summary

    Path("summary.html").write_text(render_summary(SummaryView.from_view(view)))

Out of scope, deliberately: sequence editing, restriction-site and REBASE
calculation, primer design, ORF finding, implementing aligners, and chromatogram
viewing.
"""

from __future__ import annotations

from .construct import Feature, Reference
from .grid import Cell, Row
from .pileup import PileupGroup, PileupView
from .render import render
from .render_summary import render_summary
from .summary import GroupSummary, SummaryView, Variant, summarize_group
from .theme import Theme

__version__ = "0.1.0"

__all__ = [
    "Cell",
    "Feature",
    "GroupSummary",
    "PileupGroup",
    "PileupView",
    "Read",
    "Reference",
    "Row",
    "SummaryView",
    "Theme",
    "Variant",
    "grid_from_alignment",
    "grid_from_reads",
    "reads_from_alignment",
    "render",
    "render_summary",
    "summarize_group",
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
