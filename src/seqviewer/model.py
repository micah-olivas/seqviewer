"""Compatibility re-export of the types that used to live here.

The three layers this module held have moved to the modules that own them:

============================  ===========================================
Module                        Holds
============================  ===========================================
:mod:`seqviewer.construct`    ``Reference``, ``Feature`` — the substrate
:mod:`seqviewer.grid`         ``Cell``, ``Row`` — pileup grid primitives
:mod:`seqviewer.theme`        ``Theme`` — the host-application bridge
:mod:`seqviewer.pileup`       ``PileupGroup``, ``PileupView`` — one page
============================  ===========================================

Import from those directly, or from the package root.  This module stays so that
``from .model import ...`` keeps working; it can go once nothing imports it.
"""

from __future__ import annotations

from .construct import Feature, Reference
from .grid import Cell, Row
from .pileup import PileupGroup, PileupView
from .theme import Theme

__all__ = [
    "Cell", "Feature", "PileupGroup", "PileupView", "Reference", "Row", "Theme",
]
