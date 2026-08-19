"""Draw a bounded window of an alignment at base resolution.

The pileup draws a whole reference as one pixel per base or two, which is the
right scale for seeing a subpopulation as a block and the wrong one for reading a
codon.  This module draws the opposite: a window of sixty to a hundred and twenty
columns, wide enough that every base is a letter, with each codon bracketed under
the three bases it is translated from.

Why a window rather than a zoom of the whole reference: a canvas is capped near
32,767 pixels a side and the pileup multiplies by ``devicePixelRatio``, so at
retina density and twelve pixels a base the drawing fails silently somewhere
around 1,400 bases.  Bounding the window sidesteps the ceiling entirely — this
draws the same size box whatever the reference length — and, being SVG, it has no
pixel ceiling to begin with.

The register rule is the whole point of the module, so it is stated once here and
tested directly: a codon of the reading frame starting at *frame_start* occupies
exactly ``3 * CELL_W`` pixels, beginning at a base whose offset from
*frame_start* is divisible by three.  Nothing rounds, and no minimum width is
imposed on a codon — imposing one is what put the pileup's own amino-acid track
out of register with the bases it describes.

Everything here is pure geometry over a grid, so it is testable without a
browser, and the emitted markup uses only constructs that behave identically in
Chrome and WebKit: no ``dominant-baseline`` (long broken in Safari, and ``dy`` in
ems is exact everywhere), no foreignObject, and no CSS the two engines disagree
about.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

from .codon import translate_codon
from .grid import Row

__all__ = [
    "CELL_W",
    "DEFAULT_COLUMNS",
    "MAX_READ_ROWS",
    "Window",
    "window_bounds",
    "window_css",
    "window_svg",
]

#: Pixels per base.  Wide enough for a monospace capital at 11px with air around
#: it, which is what makes this a sequence rather than a heat map.
CELL_W = 12.0

#: Height of one row of letters, and the type size inside it.
ROW_H = 15.0
FONT_SIZE = 11.0

#: Columns in a window by default.  Ninety is about as much sequence as reads as
#: sequence rather than as a wall; it is also three codons short of a line of
#: thirty-one, which keeps a whole window on one screen at any zoom.
DEFAULT_COLUMNS = 90

#: Reads drawn before the window says how many it left out.  A base-resolution
#: window is for reading bases, not for counting reads.
MAX_READ_ROWS = 16

#: Height of the band holding the codon brackets and their residues.
CODON_H = 24.0

#: Space above the ruler and below the last row.
PAD_TOP = 14.0


@dataclass
class Window:
    """A laid-out base-resolution window.

    ``start`` and ``end`` are the half-open reference span drawn.  ``rows_shown``
    and ``rows_hidden`` account for every read the group had, so a caller can say
    what was left out rather than implying it saw everything.
    """

    start: int
    end: int
    width: float
    height: float
    svg: str
    rows_shown: int = 0
    rows_hidden: int = 0
    residues: int = 0

    @property
    def columns(self) -> int:
        return self.end - self.start


def window_bounds(
    center: int, ref_len: int, columns: int = DEFAULT_COLUMNS,
) -> Tuple[int, int]:
    """A half-open span of *columns* bases around *center*, clamped to the reference.

    Clamping shifts the window rather than shrinking it, so a variant near either
    end still gets a full window of context instead of half of one.
    """
    columns = max(1, min(columns, ref_len))
    start = center - columns // 2
    start = max(0, min(start, ref_len - columns))
    return start, start + columns


def _e(text: str) -> str:
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _consensus(rows: Sequence[Row], ref_seq: str, start: int, end: int) -> List[str]:
    """Majority called base per column, ``-`` where deletions win, ``""`` if uncovered."""
    out = []
    for i in range(start, end):
        counts: Counter = Counter()
        deletions = 0
        for row in rows:
            base = row[i][0]
            if base == "-":
                deletions += 1
            else:
                counts[base.upper()] += 1
        if not counts and not deletions:
            out.append("")
        elif not counts or deletions > max(counts.values()):
            out.append("-")
        else:
            out.append(counts.most_common(1)[0][0])
    return out


def _interesting_first(
    rows: Sequence[Row], ref_seq: str, start: int, end: int,
) -> List[int]:
    """Read indices, those disagreeing inside the window first.

    An inspector is opened because something is wrong at this position, so the
    reads that show it belong at the top rather than wherever the caller's own
    ordering left them.
    """
    def score(index: int) -> tuple:
        row = rows[index]
        disagreements = 0
        covered = 0
        for i in range(start, end):
            base = row[i][0]
            if base == "-":
                continue
            covered += 1
            if base.upper() != ref_seq[i]:
                disagreements += 1
        # Most disagreements first; among equals, the read covering most of the
        # window; ties by original order so the result is deterministic.
        return (-disagreements, -covered, index)

    return sorted(range(len(rows)), key=score)


def _letter(x: float, y: float, char: str, cls: str) -> str:
    """One monospace glyph, centred in its cell.

    ``dy="0.35em"`` rather than ``dominant-baseline``: the latter was unsupported
    in WebKit for years and the two engines still disagree about ``central``
    versus ``middle``, while a dy in ems is exact in both.
    """
    return (
        f'<text class="{cls}" x="{x:.1f}" y="{y:.1f}" dy="0.35em" '
        f'text-anchor="middle">{_e(char)}</text>'
    )


def _base_class(base: str, matches: bool, prefix: str) -> str:
    if matches:
        return f"{prefix}-same"
    return f"{prefix}-b {prefix}-{base.lower()}"


def window_svg(
    ref_seq: str,
    rows: Sequence[Row],
    start: int,
    end: int,
    frame: Optional[Tuple[int, int]] = None,
    prefix: str = "svz",
    max_read_rows: int = MAX_READ_ROWS,
    label: str = "",
) -> Window:
    """Draw ``[start, end)`` of an alignment at one letter per base.

    Args:
        ref_seq: The reference, whose window is drawn as letters.
        rows: The grid, one row per read, each as wide as *ref_seq*.
        start: First reference position drawn, 0-based.
        end: One past the last position drawn.
        frame: The reading frame as ``[start, end)`` on the reference.  Codons are
            laid out from its start, so a window opened anywhere still brackets
            whole codons rather than the window's own first three bases.
        prefix: Class-name prefix, so two windows on one page cannot collide.
        max_read_rows: Reads drawn before the rest are counted instead.
        label: Optional caption drawn above the ruler.

    Returns:
        A :class:`Window` carrying the SVG and what it had to leave out.
    """
    ref = ref_seq.upper()
    start = max(0, start)
    end = min(len(ref), end)
    columns = end - start
    if columns <= 0:
        return Window(start=start, end=start, width=0, height=0, svg="")

    order = _interesting_first(rows, ref, start, end)
    shown = order[:max_read_rows]
    consensus = _consensus(rows, ref, start, end)

    has_frame = frame is not None and frame[1] - frame[0] >= 3
    codon_band = CODON_H if has_frame else 0.0

    width = columns * CELL_W
    ruler_y = PAD_TOP
    ref_y = ruler_y + ROW_H
    cons_y = ref_y + ROW_H
    codon_y = cons_y + ROW_H
    reads_y = codon_y + codon_band
    height = reads_y + len(shown) * ROW_H + 4

    parts = [
        f'<svg class="{prefix}" viewBox="0 0 {width:.0f} {height:.0f}" '
        f'width="{width:.0f}" height="{height:.0f}" '
        f'preserveAspectRatio="xMinYMin meet" role="img" '
        f'aria-label="{_e(label or "sequence window")}">'
    ]

    # --- alternating column shading, so the eye can hold a column ---
    for column in range(columns):
        if (start + column) % 10 < 5:
            continue
        parts.append(
            f'<rect class="{prefix}-stripe" x="{column * CELL_W:.1f}" '
            f'y="{ruler_y:.1f}" width="{CELL_W:.1f}" '
            f'height="{height - ruler_y - 4:.1f}" />'
        )

    # --- ruler: a number every ten bases, on the base it labels ---
    for column in range(columns):
        position = start + column
        if (position + 1) % 10:
            continue
        x = column * CELL_W + CELL_W / 2
        parts.append(
            f'<text class="{prefix}-pos" x="{x:.1f}" y="{ruler_y - 3:.1f}" '
            f'text-anchor="middle">{position + 1}</text>'
        )

    # --- reference and consensus letters ---
    for column in range(columns):
        i = start + column
        x = column * CELL_W + CELL_W / 2
        parts.append(_letter(x, ref_y + ROW_H / 2, ref[i], f"{prefix}-ref"))

        called = consensus[column]
        if not called:
            continue
        if called == "-":
            parts.append(
                f'<rect class="{prefix}-del" x="{column * CELL_W + 2:.1f}" '
                f'y="{cons_y + ROW_H / 2 - 1:.1f}" '
                f'width="{CELL_W - 4:.1f}" height="2" />'
            )
        else:
            parts.append(_letter(
                x, cons_y + ROW_H / 2, called,
                _base_class(called, called == ref[i], prefix),
            ))

    # --- codon brackets and residues, laid out from the frame, not the window ---
    residues = 0
    if has_frame:
        frame_start, frame_end = frame
        bracket = 3 * CELL_W
        # The first codon of the frame that reaches into the window.  Floor
        # division on a negative numerator rounds toward minus infinity, which is
        # what picks the codon a window starting mid-codon sits inside.
        index = max(0, (start - frame_start) // 3)
        while True:
            codon_start = frame_start + index * 3
            if codon_start >= end or codon_start + 3 > frame_end:
                break
            x0 = (codon_start - start) * CELL_W
            centre = x0 + bracket / 2
            ref_codon = ref[codon_start:codon_start + 3]
            seen = "".join(
                consensus[p - start] if start <= p < end else ref[p]
                for p in range(codon_start, codon_start + 3)
            )
            seen = "".join(c if c and c != "-" else "N" for c in seen)
            was, now = translate_codon(ref_codon), translate_codon(seen)
            state = (
                f"{prefix}-aa-diff" if was != now else f"{prefix}-aa-same"
            )
            title = f"codon {index + 1}: {ref_codon} → {seen}, {was} → {now}"
            bracket_d = (
                f"M{x0 + 1:.1f},{codon_y + 3:.1f} "
                f"V{codon_y + 7:.1f} H{x0 + bracket - 1:.1f} "
                f"V{codon_y + 3:.1f}"
            )
            glyph = f"<g><title>{_e(title)}</title>"
            glyph += f'<path class="{prefix}-bracket" d="{bracket_d}" />'
            # A codon straddling an edge keeps its bracket, which is cut off and
            # so reads as "continues past here", but only gets its residue when
            # the letter would actually land inside the window.
            if 0 <= centre <= width:
                glyph += _letter(centre, codon_y + 15, now or "?", state)
                residues += 1
            parts.append(glyph + "</g>")
            index += 1

    # --- reads ---
    for slot, read_index in enumerate(shown):
        row = rows[read_index]
        y = reads_y + slot * ROW_H
        # A "-" inside the read's own covered span is a deletion and gets a dash;
        # one outside it is absence of coverage and is left blank.  Without the
        # distinction the two look identical, which is exactly the ambiguity the
        # reducer exists to resolve.
        first = last = None
        for i, (base, _) in enumerate(row):
            if base != "-":
                if first is None:
                    first = i
                last = i
        for column in range(columns):
            i = start + column
            base = row[i][0]
            x = column * CELL_W + CELL_W / 2
            if base == "-":
                if first is not None and first < i < last:
                    parts.append(
                        f'<rect class="{prefix}-rdel" '
                        f'x="{column * CELL_W + 3:.1f}" '
                        f'y="{y + ROW_H / 2 - 1:.1f}" '
                        f'width="{CELL_W - 6:.1f}" height="1.5" />'
                    )
                continue
            upper = base.upper()
            parts.append(_letter(
                x, y + ROW_H / 2, upper,
                _base_class(upper, upper == ref[i], prefix),
            ))

    parts.append("</svg>")

    return Window(
        start=start, end=end, width=width, height=height,
        svg="".join(parts), rows_shown=len(shown),
        rows_hidden=len(rows) - len(shown), residues=residues,
    )


def window_css(prefix: str = "svz", token_prefix: str = "cv") -> str:
    """The rules a window's classes need, reading the page's own colour tokens.

    Kept as a function rather than a static block because the token prefix is the
    host application's to choose, the same way :class:`~seqviewer.theme.Theme`
    lets it name the theme bridge.
    """
    t = token_prefix
    return f"""
svg.{prefix} {{ display: block; max-width: 100%; height: auto; }}
.{prefix} text {{
    font-family: var(--mono);
    font-size: {FONT_SIZE:.0f}px;
    /* A monospace advance is what lets the layout be measured in Python; a
       fallback with a different advance would shift letters off their cells. */
    font-variant-ligatures: none;
    -webkit-font-smoothing: antialiased;
}}
.{prefix}-stripe {{ fill: var(--{t}-grid); opacity: 0.35; }}
.{prefix}-pos {{ font-size: 9px; fill: var(--muted); font-family: var(--mono); }}
.{prefix}-ref {{ fill: var(--{t}-tick-label); font-weight: 700; }}
.{prefix}-same {{ fill: var(--muted); }}
.{prefix}-b {{ font-weight: 700; }}
.{prefix}-a {{ fill: var(--{t}-a); }}
.{prefix}-t {{ fill: var(--{t}-t); }}
.{prefix}-c {{ fill: var(--{t}-c); }}
.{prefix}-g {{ fill: var(--{t}-g); }}
.{prefix}-n {{ fill: var(--muted); }}
.{prefix}-del {{ fill: var(--{t}-bad); }}
.{prefix}-rdel {{ fill: var(--{t}-bad); opacity: 0.65; }}
.{prefix}-bracket {{
    fill: none; stroke: var(--{t}-stem); stroke-width: 1.2;
}}
.{prefix}-aa-same {{ fill: var(--{t}-neutral); font-weight: 700; }}
.{prefix}-aa-diff {{ fill: var(--{t}-bad); font-weight: 700; }}
"""
