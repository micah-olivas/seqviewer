"""Where a sample came from on its plate, drawn small enough for a corner.

A demultiplexed run is one page per well, and the page names the well in its
title if it names it at all.  A plate map places it: a grid of the plate's
wells with the current one filled, in the corner of both pages, so a reader
moving between wells keeps their bearings without reading a label.

The map is an SVG built here and styled by the page that carries it, since the
two pages name their colour tokens differently.  This module emits structure
and class names only.
"""
from __future__ import annotations

import re
import string
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

__all__ = ["FORMATS", "Well", "parse_well", "plate_svg"]

#: Rows and columns of each plate format the map can draw.
FORMATS: Dict[int, Tuple[int, int]] = {96: (8, 12), 384: (16, 24)}

#: How wide one well's cell is, per format, in the drawing's own units.  384
#: wells are packed twice as tight so both maps come out about the same size.
_PITCH = {96: 7.0, 384: 3.6}

#: Margin around the grid inside the plate outline.
_INSET = 3.0

#: How much of the plate's top-left corner is cut off.  A real plate is
#: chamfered there so it can only be loaded one way round, and the same cut on
#: the map is what tells a reader which corner A1 is without a label.
_CHAMFER = 5.0

#: Height of the title line over the plate, and of the label line under it.
#: Both hold 8px type: 7px measured legible only at full size, and the map is
#: drawn at whatever width its corner allows.
_TITLE_H = 12.0
_LABEL_H = 11.0

_WELL_RE = re.compile(r"^([A-Za-z])0*([1-9][0-9]?)$")


@dataclass(frozen=True)
class Well:
    """A position on a plate, zero-indexed, with the format it is read in."""

    row: int
    col: int
    plate: int

    @property
    def label(self) -> str:
        """The name as it is written on the plate: ``A1``, ``P24``."""
        return f"{string.ascii_uppercase[self.row]}{self.col + 1}"


def parse_well(text: str, plate: Optional[int] = None) -> Well:
    """Read ``A1``, ``a01`` or ``P24`` into a :class:`Well`.

    The format is *plate* when given, and otherwise the smallest one the
    position fits: a well past row H or column 12 is on a 384-well plate, and
    anything else is read as 96.  A position that fits neither is an error,
    as is a format the map cannot draw.
    """
    match = _WELL_RE.match(text.strip())
    if not match:
        raise ValueError(f"{text!r} is not a well: expected a row letter and a "
                         "column number, such as A1 or P24")
    row = string.ascii_uppercase.index(match.group(1).upper())
    col = int(match.group(2)) - 1

    if plate is not None and plate not in FORMATS:
        raise ValueError(f"no {plate}-well plate; one of "
                         f"{', '.join(str(k) for k in sorted(FORMATS))}")
    candidates = [plate] if plate is not None else sorted(FORMATS)
    for fmt in candidates:
        rows, cols = FORMATS[fmt]
        if row < rows and col < cols:
            return Well(row, col, fmt)
    rows, cols = FORMATS[plate] if plate is not None else FORMATS[max(FORMATS)]
    raise ValueError(f"{text.strip().upper()} is off the "
                     f"{plate or max(FORMATS)}-well plate, which runs "
                     f"A1 to {string.ascii_uppercase[rows - 1]}{cols}")


def plate_svg(well: Well, prefix: str = "sv") -> str:
    """The plate as an inline SVG, with *well* filled.

    Every well is drawn, since a 384 map with only its edges shown reads as a
    96, and the filled one is drawn last so it is never under a neighbour.
    The plate's own name and the well's are in a ``<title>`` for hover and
    for a screen reader; the label under the grid repeats the well.
    """
    rows, cols = FORMATS[well.plate]
    pitch = _PITCH[well.plate]
    radius = pitch * 0.32
    width = cols * pitch + 2 * _INSET
    plate_h = rows * pitch + 2 * _INSET
    top = _TITLE_H
    height = top + plate_h + _LABEL_H

    c = _CHAMFER
    outline = (f"M{c:.1f},{top:.1f} H{width:.1f} V{top + plate_h:.1f} H0 "
               f"V{top + c:.1f} Z")

    parts = [
        f'<svg class="{prefix}-plate" viewBox="0 0 {width:.1f} {height:.1f}" '
        f'width="{width:.0f}" height="{height:.0f}" role="img" '
        f'aria-label="Well {well.label} of a {well.plate}-well plate">',
        f"<title>Well {well.label} of {well.plate}</title>",
        f'<text class="{prefix}-plate-title" x="0" y="{top - 3:.1f}">'
        f"{well.plate}-well plate</text>",
        f'<path class="{prefix}-plate-frame" d="{outline}" />',
    ]

    def centre(r: int, k: int) -> Tuple[float, float]:
        return (_INSET + (k + 0.5) * pitch, top + _INSET + (r + 0.5) * pitch)

    wells = []
    for r in range(rows):
        for k in range(cols):
            if (r, k) == (well.row, well.col):
                continue
            x, y = centre(r, k)
            wells.append(f"M{x:.1f},{y:.1f} m{-radius:.2f},0 "
                         f"a{radius:.2f},{radius:.2f} 0 1,0 {2 * radius:.2f},0 "
                         f"a{radius:.2f},{radius:.2f} 0 1,0 {-2 * radius:.2f},0")
    parts.append(f'<path class="{prefix}-plate-well" d="{" ".join(wells)}" />')

    x, y = centre(well.row, well.col)
    parts.append(f'<circle class="{prefix}-plate-here" cx="{x:.1f}" cy="{y:.1f}" '
                 f'r="{radius * 1.45:.2f}" />')
    parts.append(f'<text class="{prefix}-plate-label" x="{width:.1f}" '
                 f'y="{top + plate_h + _LABEL_H - 2:.1f}" text-anchor="end">'
                 f"{well.label}</text>")
    parts.append("</svg>")
    return "".join(parts)
