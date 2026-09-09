"""Render a summarized pileup as a self-contained HTML page.

Where the pileup page draws every read as a row of pixels, this one draws what
the reads add up to: an annotated map of the construct, and under it one compact
band per subpopulation — a reference ribbon, a lollipop per called variant, and
a coverage profile — followed by the variants as a table.  A page that took a
scroll and a squint becomes a few centimetres you can take in at once.

The drawing is **SVG generated in Python**, not canvas, for the same reasons
:mod:`seqviewer.annotate` is: a summary is a few dozen outlined shapes with text
in them rather than a pixel matrix.  SVG buys crisp strokes at any pixel ratio,
real text, native tooltips from ``<title>``, and per-theme fills from CSS.  It
also means the page carries no JavaScript data payload at all, so free-text
feature labels never reach a ``<script>`` block — a class of escaping bug this
page simply does not have.

The feature track is not drawn here.  :mod:`seqviewer.annotate` owns feature
geometry for the whole package; this module calls it with ``max_lanes`` turned
down, which is what that parameter is for.
"""

from __future__ import annotations

import math

import html as _html
from dataclasses import replace
from typing import List, Optional, Sequence, Tuple

from .annotate import TrackPlan, plan_track, track_style, track_svg
from .summary import GroupSummary, SummaryView, Variant
from .zoom import window_bounds, window_css, window_svg

# The pileup owns the package's colour table.  Imported rather than copied so
# the two views cannot drift apart; if it moves, this fails loudly at import
# time, which is the failure worth having.
from .render import _PALETTE as _PILEUP_PALETTE

__all__ = ["render_summary"]

#: Nominal drawing width.  Everything is laid out in this coordinate space and
#: the SVG is scaled to its container by its ``viewBox``, so the page is
#: responsive without measuring anything in the browser.
WIDTH = 680.0

#: Narrowest the drawing is allowed to be scaled to before its container scrolls
#: instead.  A viewBox scales type along with geometry, so without a floor a
#: narrow window shrinks the whole band into illegibility.
MIN_DRAW_WIDTH = 420

RULER_H = 18
RIBBON_H = 11
DEPTH_H = 72

#: Height of the disagreement track.
MISMATCH_H = 46

#: Height of the strip that marks called positions.  Flat, and above the track
#: rather than over it: a mark drawn to its variant frequency put a tall glyph on
#: a position whose disagreement the bars already state, so an annotated column
#: read as a worse one.
CALLS_H = 11

#: Width of the left gutter, in the same units as the drawing.  It holds the
#: rate labels, which over the bars were drawn on top of the noise they were
#: labelling, and the name of each row.
GUTTER = 98.0

#: Names for the three rows, in the order they are stacked.
ROW_LABELS = ("Mismatches", "Features", "Reads")

#: Where each of the gutter's two columns ends, as an x in the drawing: the
#: row's name, then the rate ticks nearer the plot.
LABEL_X = -34.0
TICK_X = -6.0

#: How far below the top of a track a label centred on its line has to sit
#: before it stops hanging over the row above.  Half its own type size.
TICK_INSET = 4.5

#: Rates the disagreement track rules a line at.
MISMATCH_MARKS = (0.10, 0.50)

#: Rate at or above which a column is drawn as a finding rather than as noise.
#: Below it the bar is muted: every position disagrees a little, and drawing all
#: of it in the alert colour spends the reader's attention on the noise floor.
MISMATCH_ALERT = 0.10

#: Separation between the reference ribbon and the coverage profile, so a group
#: at full depth does not read as one thick bar.
RIBBON_GAP = 3


#: Heads closer together than this are lifted onto another tier rather than
#: drawn on top of each other, up to this many tiers.


#: Base-resolution windows drawn per group, and their size.  Bounded because each
#: window is one SVG text element per base per row: a window is cheap, forty are
#: not, and a summary that carries forty has stopped summarising.
MAX_INSPECTORS = 6
INSPECTOR_COLUMNS = 96
INSPECTOR_READS = 16

#: Consequence -> the tier its glyph and chip are drawn in.
_TIER_OF = {
    "frameshift": "bad",
    "nonsense": "bad",
    "inframe_indel": "warn",
    "missense": "warn",
    "silent": "ok",
    "noncoding": "neutral",
    "": "neutral",
}

#: Verdict -> the words a group's chip carries.
_VERDICT_TEXT = {
    "clean": "No variants called",
    "variant": "Variants called",
    "noncoding": "Outside the reading frame",
    "silent": "Silent only",
    "missense": "Missense",
    "inframe_indel": "In-frame indel",
    "nonsense": "Premature stop",
    "frameshift": "Frameshift",
}

#: Colours this view needs that the pileup's table has no use for.
_EXTRA_PALETTE = {
    "light": {
        "depth": "#d3d9e0",
        "depth-edge": "#9aa3ad",
        "stem": "#8b95a1",
        "ribbon": "#c8ccd0",
        "focus": "#aab3bd",
        "grid": "#e5e7eb",
        "mm-bed": "#f4f6f8",
        "mm-quiet": "#b9c2cb",
        "mm": "#c22f2f",
        "ok": "#0f7a52",
        "warn": "#a1650a",
        "bad": "#c22f2f",
        "neutral": "#5b6672",
        "ok-bg": "#e8f6ef",
        "warn-bg": "#fdf3e2",
        "bad-bg": "#fdecec",
        "neutral-bg": "#eef1f4",
    },
    "dark": {
        "depth": "#333f55",
        "depth-edge": "#6b7688",
        "stem": "#6b7688",
        "ribbon": "#4a5568",
        "focus": "#7a8698",
        "grid": "#334155",
        "mm-bed": "#1c2536",
        "mm-quiet": "#4a5568",
        "mm": "#f0655f",
        "ok": "#35c493",
        "warn": "#e0a355",
        "bad": "#f0655f",
        "neutral": "#93a1b0",
        "ok-bg": "#14312a",
        "warn-bg": "#33280f",
        "bad-bg": "#35191c",
        "neutral-bg": "#202b45",
    },
}


def _palette() -> dict:
    """The pileup's colours plus this view's, per theme."""
    return {
        theme: {**_PILEUP_PALETTE[theme], **_EXTRA_PALETTE[theme]}
        for theme in ("light", "dark")
    }


def _e(text: str) -> str:
    """Escape for both HTML body text and SVG attribute values."""
    return _html.escape(str(text), quote=True)


def _nice_step(ref_len: int, target_ticks: int = 10) -> int:
    """A round number of bases between ruler ticks, near *target_ticks* of them."""
    if ref_len <= target_ticks:
        return 1
    rough = ref_len / target_ticks
    magnitude = 10 ** max(0, len(str(int(rough))) - 1)
    for multiple in (1, 2, 2.5, 5, 10):
        step = magnitude * multiple
        if step >= rough:
            return int(step)
    return int(magnitude * 10)


def _x(pos: float, cell_w: float) -> float:
    return pos * cell_w


# --------------------------------------------------------------------------
# The map: ruler plus the annotation track, drawn once for the whole page
# --------------------------------------------------------------------------

def _ruler_parts(ref_len: int, cell_w: float, height: int) -> List[str]:
    """Tick marks and base labels along the top of the map."""
    parts = [
        f'<line class="sv-axis" x1="0" y1="{height}" '
        f'x2="{WIDTH:.1f}" y2="{height}" />'
    ]
    step = _nice_step(ref_len)
    # The end label is anchored to the right edge, so a tick label landing
    # under it is dropped rather than drawn through it.  Its width is estimated
    # from the monospace advance the stylesheet sets.
    end_text = f"{ref_len} bp"
    end_left = WIDTH - len(end_text) * 6.0
    pos = step
    while pos < ref_len:
        x = _x(pos, cell_w)
        parts.append(
            f'<line class="sv-tick" x1="{x:.1f}" y1="{height - 4}" '
            f'x2="{x:.1f}" y2="{height}" />'
        )
        if x + len(str(pos)) * 3.0 < end_left:
            parts.append(
                f'<text class="sv-tick-label" x="{x:.1f}" y="{height - 7}" '
                f'text-anchor="middle">{pos}</text>'
            )
        pos += step
    parts.append(
        f'<text class="sv-tick-label sv-tick-end" x="{WIDTH:.1f}" '
        f'y="{height - 7}" text-anchor="end">{ref_len} bp</text>'
    )
    return parts


def _focus_parts(focus: Optional[Tuple[int, int]], cell_w: float,
                 height: float) -> List[str]:
    """Dashed boundaries of the reading frame, drawn through a band's height."""
    if focus is None:
        return []
    parts = []
    for edge in focus:
        x = _x(edge, cell_w)
        parts.append(
            f'<line class="sv-focus-edge" x1="{x:.1f}" y1="0" '
            f'x2="{x:.1f}" y2="{height:.1f}" />'
        )
    return parts


# --------------------------------------------------------------------------
# One group's band: disagreement above the reference, reads below
# --------------------------------------------------------------------------

def _depth_parts(group: GroupSummary, cell_w: float, ceiling: int,
                 top: float, peaks: Optional[Sequence[float]] = None) -> List[str]:
    """A filled coverage profile under the reference.

    Each pixel column reports the *thinnest* coverage it spans, not the mean, so
    a dropout narrower than one pixel still shows as a notch.  On a page whose
    job is to be trusted at a glance, a coverage hole that averages away is the
    failure worth avoiding.

    The height is logarithmic.  A run where one region draws many times the
    reads of the rest is common -- primer dimer and truncated product both do
    it -- and on a linear axis that one column sets the ceiling and presses the
    rest of the profile flat, so an ordinary dropout stops being visible at all.
    """
    if not group.depth or not ceiling:
        return []
    base = top + DEPTH_H
    scale = math.log10(1.0 + ceiling)
    columns = min(int(WIDTH), group.ref_len)
    points = []
    for column in range(columns):
        start = group.ref_len * column // columns
        end = max(start + 1, group.ref_len * (column + 1) // columns)
        value = min(group.depth[start:end])
        x = WIDTH * column / columns
        share = math.log10(1.0 + value) / scale if scale else 0.0
        y = base - DEPTH_H * min(1.0, share)
        points.append(f"{x:.1f},{y:.1f}")
    # Close the polygon along the baseline so it fills.
    outline = " ".join(points)
    parts = [
        f'<polygon class="sv-depth" points="0,{base} {outline} '
        f'{WIDTH:.1f},{base} " />',
        f'<polyline class="sv-depth-edge" points="{outline}" />',
    ]
    # The reads that disagree, as the top slice of their own column.  The
    # column's height is logarithmic so that one deep region does not press the
    # rest flat; the split within a column is proportional, so a column half
    # red is a position half of whose reads disagree.
    if peaks:
        strokes = []
        for column, peak in enumerate(peaks):
            if peak < MISMATCH_ALERT:
                continue
            start = group.ref_len * column // len(peaks)
            end = max(start + 1, group.ref_len * (column + 1) // len(peaks))
            value = min(group.depth[start:end])
            if not value:
                continue
            drawn = DEPTH_H * min(1.0, math.log10(1.0 + value) / scale)
            top_y = base - drawn
            x = WIDTH * column / len(peaks) + WIDTH / len(peaks) / 2.0
            strokes.append(f"M{x:.1f},{top_y + drawn * peak:.1f}V{top_y:.1f}")
        if strokes:
            parts.append(
                f'<path class="sv-depth-hot" d="{"".join(strokes)}" '
                f'stroke-width="{max(1.0, WIDTH / len(peaks)):.2f}" />'
            )
    parts.append(_row_label(ROW_LABELS[2], top, DEPTH_H))
    return parts


def _row_label(text: str, top: float, height: float) -> str:
    """The row's name, right-aligned in the gutter and centred on the row."""
    return (
        f'<text class="sv-row-label" x="{LABEL_X:.0f}" '
        f'y="{top + height / 2:.1f}" '
        f'text-anchor="end" dominant-baseline="middle">{_e(text)}</text>'
    )


def _mismatch_height(fraction: float, height: float) -> float:
    """Return the drawn height of *fraction*.

    Proportional, so a rate reads off the axis directly and the noise floor
    stays near the baseline where it belongs.  A log scale lifts a one-percent
    column to a third of the height, which fills the track with sequencing error
    and leaves nothing to distinguish the columns worth acting on.
    """
    if fraction <= 0.0:
        return 0.0
    return height * min(1.0, fraction)


def _mismatch_peaks(group: GroupSummary, columns: int) -> List[float]:
    """Return the worst disagreement in each of *columns* pixel columns.

    The worst, not the mean.  A reference several kilobases long is drawn a few
    bases to the pixel, and a position where every read disagrees would
    otherwise be averaged down to the rate of its neighbours -- which is the
    difference between a cloning error and a quiet stretch.  It is the same
    reasoning as the coverage profile below, which reports the thinnest column
    rather than the mean so that a narrow dropout survives.
    """
    depth, matches = group.depth, group.matches
    peaks = []
    for column in range(columns):
        start = group.ref_len * column // columns
        end = max(start + 1, group.ref_len * (column + 1) // columns)
        worst = 0.0
        for i in range(start, end):
            covered = depth[i]
            if covered:
                rate = (covered - matches[i]) / covered
                if rate > worst:
                    worst = rate
        peaks.append(worst)
    return peaks


def _mismatch_parts(group: GroupSummary, top: float) -> List[str]:
    """Per-position disagreement with the reference, drawn along the band.

    Disagreement is counted per position rather than per variant: a position
    where forty percent of reads disagree is drawn at forty percent whether they
    all read the same wrong base or three different ones.  Allele frequency is
    what the variant table reports, and it is blind to a position the reads
    disagree about in several directions at once, which is what a chimeric read
    and a bad basecall both look like.
    """
    if not group.depth or not group.matches:
        return []

    base = top + MISMATCH_H
    parts = [
        f'<rect class="sv-mm-bed" x="0" y="{top:.1f}" width="{WIDTH:.0f}" '
        f'height="{MISMATCH_H}" />'
    ]
    # Rules under the bars, labels over them: the one-percent rule sits where
    # the sequencing noise is densest, and a label drawn under that is unread.
    labels = []
    for mark in MISMATCH_MARKS:
        y = base - _mismatch_height(mark, MISMATCH_H)
        parts.append(
            f'<line class="sv-mm-rule" x1="0" y1="{y:.1f}" '
            f'x2="{WIDTH:.0f}" y2="{y:.1f}" />'
        )
        # Centred on its line, except within half a glyph of the top of the
        # track, where the label would hang over the row above.
        labels.append(
            f'<text class="sv-mm-tick" x="{TICK_X:.0f}" '
            f'y="{max(y, top + TICK_INSET):.1f}" '
            f'text-anchor="end" '
            f'dominant-baseline="middle">{mark:.0%}</text>'
        )

    labels.append(_row_label(ROW_LABELS[0], top, MISMATCH_H))
    columns = min(int(WIDTH), group.ref_len)
    step = WIDTH / columns
    quiet, loud = [], []
    for column, peak in enumerate(_mismatch_peaks(group, columns)):
        if peak <= 0.0:
            continue
        x = column * step + step / 2.0
        y = base - _mismatch_height(peak, MISMATCH_H)
        stroke = f"M{x:.1f},{base:.1f}V{y:.1f}"
        (loud if peak >= MISMATCH_ALERT else quiet).append(stroke)
    width = f'{max(1.0, step):.2f}'
    if quiet:
        parts.append(
            f'<path class="sv-mm-quiet" d="{"".join(quiet)}" '
            f'stroke-width="{width}" />'
        )
    if loud:
        parts.append(
            f'<path class="sv-mm" d="{"".join(loud)}" stroke-width="{width}" />'
        )
    return parts + labels


def _call_marks(group: GroupSummary, cell_w: float, top: float) -> List[str]:
    """Flat marks at the called positions, one size regardless of frequency.

    The frequency is in the table and the disagreement is in the bars.  Drawing
    it a third time as a height put the tallest mark on the page over a position
    that was not the worst one.
    """
    parts = []
    for variant in group.variants:
        x = _x(variant.pos, cell_w) + cell_w / 2.0
        tier = _TIER_OF.get(variant.consequence, "neutral")
        share = variant.count / variant.depth if variant.depth else 0.0
        parts.append(
            f'<g class="sv-mark sv-{tier}"><title>{_e(_mark_title(variant, share))}'
            f'</title><path d="M{x - 3.2:.1f},{top:.1f} h6.4 '
            f'l-3.2,{CALLS_H - 3:.1f} Z" /></g>'
        )
    return parts


def _mark_title(variant: Variant, share: float) -> str:
    """The hover text for a called position: where, what, and how many."""
    what = variant.alt if variant.kind == "snv" else f"{variant.kind} {variant.alt}"
    effect = f" \u00b7 {variant.effect}" if variant.effect else ""
    return (f"{variant.pos + 1}: {variant.ref}\u2192{what} \u00b7 "
            f"{share:.0%} of {variant.depth} reads{effect}")


def _upper_svg(group: GroupSummary, view: SummaryView, cell_w: float) -> str:
    """Disagreement, with a flat strip marking any called positions above it.

    A group with nothing called reserves no room for the strip, since keeping
    the height uniform spends the space on nothing and leaves the frame's
    boundary lines hanging in it.
    """
    marks_h = CALLS_H if group.variants else 0.0
    height = marks_h + MISMATCH_H

    parts = [
        f'<svg class="sv-band" viewBox="{-GUTTER:.0f} 0 {GUTTER + WIDTH:.0f} {height:.0f}" '
        f'width="{GUTTER + WIDTH:.0f}" height="{height:.0f}" role="img" '
        f'aria-label="{_e(group.name)} disagreement with the reference">'
    ]
    parts += _mismatch_parts(group, marks_h)
    if group.variants:
        parts += _call_marks(group, cell_w, 0.0)
    parts += _focus_parts(view.focus, cell_w, height)
    parts.append("</svg>")
    return "".join(parts)


def _lower_svg(group: GroupSummary, view: SummaryView, cell_w: float,
               ceiling: int) -> str:
    """The reads: coverage, and where this group falls short of the reference."""
    height = DEPTH_H
    parts = [
        f'<svg class="sv-band" viewBox="{-GUTTER:.0f} 0 {GUTTER + WIDTH:.0f} {height:.0f}" '
        f'width="{GUTTER + WIDTH:.0f}" height="{height:.0f}" role="img" '
        f'aria-label="{_e(group.name)} coverage">'
    ]
    peaks = _mismatch_peaks(group, min(int(WIDTH), group.ref_len)) \
        if group.depth else None
    parts += _depth_parts(group, cell_w, ceiling, 0.0, peaks)
    if group.ref_len < view.ref_len:
        # Stated against a shorter reference; say so rather than letting the
        # profile trail off as though coverage simply ran out.
        edge = _x(group.ref_len, cell_w)
        parts.append(
            f'<rect class="sv-ribbon-absent" x="{edge:.1f}" y="0" '
            f'width="{WIDTH - edge:.1f}" height="{height}" />'
        )
    parts += _focus_parts(view.focus, cell_w, height)
    parts.append("</svg>")
    return "".join(parts)


def flagged_runs(group: GroupSummary,
                 alert: float = MISMATCH_ALERT) -> List[Tuple[int, int, float]]:
    """Return the stretches disagreeing at or above *alert*.

    Each is ``(start, end, peak)`` over ``[start, end]`` inclusive.  Adjacent
    positions are joined: a chimeric block or a run of deleted bases covers
    hundreds of positions, and one stretch is a finding where several hundred
    rows are a list to scroll.
    """
    runs: List[Tuple[int, int, float]] = []
    start: Optional[int] = None
    peak = 0.0
    for i, (depth, matches) in enumerate(zip(group.depth, group.matches)):
        rate = (depth - matches) / depth if depth else 0.0
        if rate >= alert:
            if start is None:
                start, peak = i, rate
            else:
                peak = max(peak, rate)
        elif start is not None:
            runs.append((start, i - 1, peak))
            start = None
    if start is not None:
        runs.append((start, len(group.depth) - 1, peak))
    return runs


def _run_label(group: GroupSummary, run: Tuple[int, int, float]) -> str:
    """Name one flagged stretch: where it is, how wide, and how bad."""
    start, end, peak = run
    span = end - start + 1
    where = f"{start + 1:,}" if span == 1 else f"{start + 1:,}\u2013{end + 1:,}"
    width = "" if span == 1 else f" \u00b7 {span:,} bp"
    called = any(start <= v.pos <= end for v in group.variants)
    tail = " \u00b7 variant called" if called else ""
    return f"{where}{width} \u00b7 {peak:.0%} of reads{tail}"


def _call_count(group: GroupSummary, runs: Sequence) -> str:
    """The line on the collapsed detail, and the only statement of the count.

    A long list is not drawn in full — each window costs a hundred columns of
    reads — so the line says how many stretches are below it.
    """
    if not runs:
        return f"No position disagrees at {MISMATCH_ALERT:.0%} or more"
    n = len(runs)
    called = sum(1 for v in group.variants)
    plural = "" if n == 1 else "s"
    line = f"{n} position{plural} to check"
    if n > MAX_INSPECTORS:
        line += f", {MAX_INSPECTORS} shown"
    if called:
        line += f" \u00b7 {called} called variant{'' if called == 1 else 's'}"
    return line


# --------------------------------------------------------------------------
# The prose around each band
# --------------------------------------------------------------------------

def _chip(group: GroupSummary) -> str:
    verdict = group.verdict
    tier = "ok" if verdict in ("clean", "silent") else _TIER_OF.get(
        verdict, "neutral"
    )
    if verdict in ("variant", "noncoding"):
        tier = "neutral"
    text = _VERDICT_TEXT.get(verdict, verdict)
    return f'<span class="sv-chip sv-chip-{tier}">{_e(text)}</span>'


def _facts(group: GroupSummary, view: SummaryView) -> str:
    facts = []
    if group.n_reads and view.total_reads:
        facts.append(f"<b>{group.n_reads}</b> of <b>{view.total_reads}</b> reads")
    elif group.n_reads:
        facts.append(f"<b>{group.n_reads}</b> reads")
    if group.rows_drawn != group.n_reads:
        facts.append(f"<b>{group.rows_drawn}</b> drawn")
    if group.identity is not None:
        facts.append(f"<b>{group.identity:.1%}</b> identity")
    facts.append(f"<b>{group.mean_depth:.0f}&times;</b> mean depth")
    if group.covered < group.ref_len:
        facts.append(
            f"<b>{group.covered / group.ref_len:.0%}</b> of the reference covered"
        )
    return " &middot; ".join(facts)


def _reading_frame(view: SummaryView, start: int, end: int):
    """The frame to translate this window in.

    The focus region where one was named, and otherwise the longest coding
    feature covering the window: a reference that carries a CDS carries a
    translation, whether or not a caller singled a region out.
    """
    if view.focus is not None:
        return view.focus
    coding = [f for f in view.features
              if f.type.upper() in ("CDS", "GENE")
              and f.start < end and f.end > start
              and f.end - f.start >= 3]
    if not coding:
        return None
    widest = max(coding, key=lambda f: f.end - f.start)
    return (widest.start, widest.end)


def _inspectors(group: GroupSummary, view: SummaryView, index: int,
                runs: Sequence) -> str:
    """A base-resolution window per called variant, each collapsed until asked for.

    The map answers "where" and the table answers "what"; this answers "show me
    the actual bases", which is the one question neither can. Windows are
    ``<details>`` rather than JavaScript-driven panels: a summary should stay
    readable with nothing running, and a disclosure widget is the one interaction
    both engines implement identically.
    """
    source = view.source
    if source is None or not runs or index >= len(source.groups):
        return ""

    origin = source.groups[index]
    if not origin.rows:
        return ""

    fractions = [
        (d - m) / d if d else 0.0 for d, m in zip(group.depth, group.matches)
    ]
    # Worst first, so the windows drawn are the ones worth the space.
    ranked = sorted(runs, key=lambda r: (-r[2], r[0]))

    blocks = []
    for run in ranked[:MAX_INSPECTORS]:
        middle = (run[0] + run[1]) // 2
        start, end = window_bounds(middle, len(origin.ref_seq),
                                   INSPECTOR_COLUMNS)
        label = _run_label(group, run)
        window = window_svg(
            origin.ref_seq, origin.rows, start, end,
            frame=_reading_frame(view, start, end),
            max_read_rows=INSPECTOR_READS,
            label=label,
            features=[f for f in view.features
                      if f.start < end and f.end > start],
            mismatch=fractions, alert=MISMATCH_ALERT,
        )
        if not window.svg:
            continue
        hidden = ""
        if window.rows_hidden:
            hidden = (
                f'<p class="sv-none">{window.rows_shown} of '
                f"{window.rows_shown + window.rows_hidden} reads shown, those "
                f"disagreeing here first.</p>"
            )
        covering = [f.label or f.type for f in view.features
                    if f.start <= middle < f.end]
        place = ", ".join(dict.fromkeys(covering)) if covering else "intergenic"
        caption = (
            f'<span class="sv-zoom-at">{_e(label)}</span>'
            f'<span class="sv-zoom-in">{_e(place)}</span>'
        )
        blocks.append(
            '<div class="sv-zoom">'
            f'<div class="sv-zoom-head">{caption}</div>'
            f'<div class="sv-zoom-body">{window.svg}{hidden}</div>'
            "</div>"
        )

    if not blocks:
        return ""
    return f'<div class="sv-zooms">{"".join(blocks)}</div>'


# --------------------------------------------------------------------------
# The page shell
# --------------------------------------------------------------------------

#: The gutter as a share of the whole drawing, for the one track laid out in
#: CSS rather than in a viewBox.
GUTTER_PCT = 100.0 * GUTTER / (GUTTER + WIDTH)

#: Where the row's name ends, as a share of the whole drawing.  Percentage
#: padding resolves against the containing block rather than the element, so
#: this is stated over the drawing's full width and not over the gutter's.
LABEL_PCT = 100.0 * abs(LABEL_X) / (GUTTER + WIDTH)


def _shell(view: SummaryView, palette: dict, body: str, track_css: str,
           window_style: str) -> str:
    """Wrap *body* in a complete, self-contained document.

    THE SEAM.  Everything specific to being an HTML page rather than a drawing
    lives in this one function: the document skeleton, the reset, the palette
    emission, and the light/dark bridge that reads the host application's stored
    preference.  None of it is particular to a summary — the pileup page builds
    the same shell inside its own f-string — so when that markup is extracted
    into a shared asset this function is what gets replaced, and nothing above
    it has to change.
    """
    theme = view.theme
    prefix = theme.css_prefix

    palette_css = "\n".join(
        [":root {"]
        + [f"    --{prefix}-{k}: {v};" for k, v in palette["light"].items()]
        + ["}", '[data-theme="dark"] {']
        + [f"    --{prefix}-{k}: {v};" for k, v in palette["dark"].items()]
        + ["}"]
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{_e(view.title)}</title>
<style id="{theme.style_id}">
{palette_css}
:root {{
    --{prefix}-bg: #fafafa;
    --text: #1e293b;
    --muted: #94a3b8;
    --card-bg: #ffffff;
    --panel-line: #dfe3e8;
    --mono: 'SF Mono', SFMono-Regular, Menlo, Consolas, monospace;
}}
[data-theme="dark"] {{
    --{prefix}-bg: #1a1a2e;
    --text: #e0e0e0;
    --muted: #64748b;
    --card-bg: #16213e;
    --panel-line: #2c3a55;
}}
html, body {{
    background: var(--{prefix}-bg);
    color: var(--text);
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    margin: 0;
}}
/* The inset goes on body alone.  Setting it on both nests one inside the other
   and the page starts at twice the intended distance from every edge. */
body {{
    /* Room at the top for the fixed view toggle, and no more. */
    padding: 2.7rem 1.5rem 1.5rem;
}}
.sv-wrap {{ max-width: 760px; margin: 0 auto; }}
h1 {{ font-size: 1.3rem; margin: 0 0 0.2rem; letter-spacing: -0.01em; }}
.sv-sub {{ color: var(--muted); font-size: 0.85rem; margin-bottom: 1.1rem; }}
.sv-panel {{
    background: var(--card-bg);
    border: 1px solid var(--panel-line);
    border-radius: 3px;
    padding: 0.75rem 0.9rem;
    margin-bottom: 1.1rem;
}}
.sv-eyebrow {{
    font: 600 0.66rem/1 var(--mono);
    letter-spacing: 0.1em;
    text-transform: uppercase;
    color: var(--muted);
    margin-bottom: 0.45rem;
}}
.sv-group {{ margin: 0 0 1.6rem; }}
.sv-group-head {{
    display: flex; align-items: baseline; gap: 0.6rem;
    flex-wrap: wrap; margin-bottom: 0.15rem;
}}
.sv-name {{ font-size: 1.02rem; font-weight: 700; }}
.sv-star {{ color: var(--{prefix}-warn); cursor: help; }}
.sv-facts {{ color: var(--muted); font-size: 0.8rem; margin-bottom: 0.4rem; }}
.sv-chip {{
    font: 600 0.7rem/1 var(--mono);
    padding: 0.22rem 0.45rem;
    border-radius: 3px;
    white-space: nowrap;
}}
.sv-chip-ok {{ color: var(--{prefix}-ok); background: var(--{prefix}-ok-bg); }}
.sv-chip-warn {{ color: var(--{prefix}-warn); background: var(--{prefix}-warn-bg); }}
.sv-chip-bad {{ color: var(--{prefix}-bad); background: var(--{prefix}-bad-bg); }}
.sv-chip-neutral {{
    color: var(--{prefix}-neutral); background: var(--{prefix}-neutral-bg);
}}
/* An SVG scaled by its viewBox: laid out at a nominal width, drawn at whatever
   width the page gives it, crisp either way. */
svg.sv-band, svg.sv-map, svg.sv-annot {{
    display: block; width: 100%; height: auto; overflow: visible;
    /* A viewBox scales the type along with the drawing, so past a point the
       whole thing shrinks into illegibility. Below this the drawing keeps its
       size and its container scrolls instead — measured: at a 420px viewport an
       unfloored annotation track renders 6px tall. */
    min-width: {MIN_DRAW_WIDTH}px;
}}
.sv-scroll {{ overflow-x: auto; overflow-y: visible; }}
.sv-annot text {{ font: 10px var(--mono); }}
.sv-annot-out {{ fill: var(--muted); }}
.sv-annot path {{ stroke-width: 1; }}
.sv-axis, .sv-tick {{ stroke: var(--{prefix}-tick); stroke-width: 1; }}
.sv-tick-label {{ font: 9px var(--mono); fill: var(--muted); }}
.sv-tick-end {{ fill: var(--{prefix}-tick-label); }}
.sv-ribbon-absent {{ fill: var(--{prefix}-grid); opacity: 0.5; }}
.sv-mm-bed {{ fill: var(--{prefix}-mm-bed); }}
.sv-mm {{ stroke: var(--{prefix}-mm); fill: none; }}
.sv-mm-quiet {{ stroke: var(--{prefix}-mm-quiet); fill: none; }}
.sv-mark path {{ stroke: none; }}
.sv-mark.sv-ok path {{ fill: var(--{prefix}-ok); }}
.sv-mark.sv-warn path {{ fill: var(--{prefix}-warn); }}
.sv-mark.sv-bad path {{ fill: var(--{prefix}-bad); }}
.sv-mark.sv-neutral path {{ fill: var(--{prefix}-neutral); }}
.sv-mm-rule {{ stroke: var(--{prefix}-grid); stroke-width: 1;
    stroke-dasharray: 2 3; }}
.sv-row-label {{
    fill: var(--{prefix}-neutral); font: 9px var(--mono); opacity: 0.85;
}}
.sv-mm-tick {{
    fill: var(--{prefix}-neutral); font-size: 8px;
    paint-order: stroke; stroke: var(--{prefix}-mm-bed); stroke-width: 2.5px;
    stroke-linejoin: round;
}}
.sv-depth {{ fill: var(--{prefix}-depth); }}
.sv-depth-edge {{ fill: none; stroke: var(--{prefix}-depth-edge); stroke-width: 1; }}
.sv-depth-hot {{ stroke: var(--{prefix}-mm); fill: none; opacity: 0.85; }}
.sv-stem {{ stroke: var(--{prefix}-stem); stroke-width: 1.4; }}
.sv-head {{ stroke: var(--{prefix}-bg); stroke-width: 1; }}
.sv-head-letter {{
    font: 700 7px var(--mono); fill: #ffffff; pointer-events: none;
}}
.sv-focus-edge {{
    stroke: var(--{prefix}-boundary); stroke-width: 1; stroke-dasharray: 4 3;
    opacity: 0.75;
}}
.sv-ok {{ fill: var(--{prefix}-ok); }}
.sv-warn {{ fill: var(--{prefix}-warn); }}
.sv-bad {{ fill: var(--{prefix}-bad); }}
.sv-neutral {{ fill: var(--{prefix}-neutral); }}
.sv-dot {{
    display: inline-block; width: 7px; height: 7px; border-radius: 50%;
    margin-right: 0.4rem; vertical-align: middle;
}}
.sv-dot.sv-ok {{ background: var(--{prefix}-ok); }}
.sv-dot.sv-warn {{ background: var(--{prefix}-warn); }}
.sv-dot.sv-bad {{ background: var(--{prefix}-bad); }}
.sv-dot.sv-neutral {{ background: var(--{prefix}-neutral); }}
.sv-none {{ color: var(--muted); font-size: 0.8rem; font-style: italic; }}
.sv-sep {{ border: none; border-top: 1px solid var(--panel-line); margin: 0.9rem 0; }}
/* --- the three stacked tracks --- */
.sv-stack {{
    /* The tracks share one coordinate system, so they are flush: any gap
       between them reads as a gap in the sequence. */
    display: flex; flex-direction: column; align-items: center;
    margin: 0.35rem auto 0;
}}
.sv-stack > * {{ display: block; line-height: 0; }}
/* --- what this page is, and the page beside it --- */
/* --- The hover readout for a feature --- */
.sv-tip {{
    position: fixed; z-index: 30; pointer-events: none;
    display: none; max-width: 22rem;
    padding: 0.3rem 0.45rem; border-radius: 3px;
    border: 1px solid var(--panel-line); background: var(--card-bg);
    font: 0.72rem/1.35 var(--mono); color: var(--text);
    box-shadow: 0 2px 8px rgba(0, 0, 0, 0.12);
    white-space: pre-line;
}}
/* The label belongs to the glyph under it, so it must not take the hover
   itself: the pointer would then be over the text and not over the feature. */
.sv-annot text {{ pointer-events: none; }}
/* --- Views: fixed to the viewport, so the pair does not move between pages,
   and above the crossfade so it stays legible through it. --- */
.sv-views-fixed {{ position: fixed; top: 0.8rem; left: 1.5rem;
    z-index: 25; }}
/* The crossfade: the page's own background, drawn over the content.  It starts
   opaque and clears once the page has drawn, so a page is never seen filling in
   behind a fade that has already finished. */
.sv-fade {{
    position: fixed; inset: 0; z-index: 15; pointer-events: none;
    background: var(--{prefix}-bg); opacity: 1;
    transition: opacity 180ms ease-out;
}}
html.sv-ready .sv-fade {{ opacity: 0; }}
html.sv-leaving .sv-fade {{ opacity: 1; transition: opacity 90ms ease-in; }}
@media (prefers-reduced-motion: reduce) {{
    .sv-fade {{ transition: none; }}
}}
.sv-views {{
    display: inline-flex; border: 1px solid var(--panel-line);
    border-radius: 3px; overflow: hidden;
}}
.sv-view {{
    font: 600 0.68rem/1 var(--mono); letter-spacing: 0.04em;
    padding: 0.35rem 0.6rem; color: var(--muted); background: var(--card-bg);
    text-decoration: none; white-space: nowrap;
}}
.sv-view + .sv-view {{ border-left: 1px solid var(--panel-line); }}
.sv-view:hover {{ color: var(--text); background: var(--{prefix}-neutral-bg); }}
.sv-view[aria-current="page"] {{
    color: var(--text); background: var(--{prefix}-neutral-bg);
    cursor: default;
}}
.sv-stack > svg.sv-band:last-child {{ margin-top: 5px; }}
.sv-ref {{
    /* The other tracks carry the gutter inside their viewBox; this one holds
       two SVGs drawn over the sequence alone, so the gutter is a flex column.
       Both resolve to the same fraction of the container, which is what keeps
       a base at the same x in every track. */
    align-self: stretch; width: 100%; margin: 1px 0;
    display: flex; align-items: flex-end;
}}
.sv-ref-label {{
    flex: 0 0 {GUTTER_PCT:.4f}%; box-sizing: border-box;
    padding-right: {LABEL_PCT:.4f}%; padding-bottom: 4px; text-align: right;
    font: 9px var(--mono); color: var(--{prefix}-neutral); opacity: 0.85;
    line-height: 1;
}}
.sv-ref-body {{ flex: 1 1 auto; min-width: 0; }}
.sv-calls {{ margin-top: 0.5rem; font-size: 0.78rem; }}
.sv-calls > summary {{
    cursor: pointer; color: var(--muted); list-style: none;
    padding: 0.15rem 0;
}}
.sv-calls > summary::-webkit-details-marker {{ display: none; }}
.sv-calls > summary::before {{ content: "\\25B8 "; }}
.sv-calls[open] > summary::before {{ content: "\\25BE "; }}
.sv-calls-body {{ padding-top: 0.2rem; }}
/* --- what the page draws, at the foot of it --- */
.sv-note {{
    margin-top: 1.4rem; border-top: 1px solid var(--panel-line);
    padding-top: 0.5rem; font-size: 0.76rem; color: var(--muted);
}}
.sv-note > summary {{ cursor: pointer; list-style: none; }}
.sv-note > summary::-webkit-details-marker {{ display: none; }}
.sv-note > summary::before {{ content: "\\25B8 "; }}
.sv-note[open] > summary::before {{ content: "\\25BE "; }}
.sv-note-body {{ margin: 0.5rem 0 0; }}
.sv-note-row {{ display: flex; gap: 0.8rem; margin-bottom: 0.4rem; }}
.sv-note-row dt {{
    flex: 0 0 7.5rem; font-weight: 600; color: var(--text);
}}
.sv-note-row dd {{ margin: 0; }}
@media (max-width: 560px) {{
    .sv-note-row {{ flex-direction: column; gap: 0.1rem; }}
}}
/* --- base-resolution windows --- */
.sv-zooms {{ margin-top: 0.7rem; }}
/* One flat block per window.  A border and a disclosure of its own put two
   frames and two triangles around a drawing already inside a disclosure. */
.sv-zoom + .sv-zoom {{
    margin-top: 0.9rem;
    padding-top: 0.9rem;
    border-top: 1px solid var(--panel-line);
}}
.sv-zoom-head {{
    display: flex; align-items: baseline; gap: 0.5rem;
    margin-bottom: 0.3rem;
}}
.sv-zoom-at {{ font: 600 0.76rem/1.3 var(--mono); color: var(--text); }}
.sv-zoom-in {{ font-size: 0.72rem; color: var(--muted); }}
.sv-zoom-body {{ overflow-x: auto; }}
{window_style}
{track_css}
</style>
</head>
<body>
<div class="sv-fade"></div>
<div class="sv-tip" role="tooltip" aria-hidden="true"></div>
<noscript><style>.sv-fade {{ display: none; }}</style></noscript>
<div class="sv-wrap">
{body}
</div>
<script>
/* The crossfade to the pileup page.
 *
 * The drawing here is the markup, so the page is finished by the time this
 * runs and the overlay lifts on the next paint.  A modified click is left alone
 * so it can still open a tab, and a reader who asked for less motion navigates
 * with no fade at all. */
(function () {{
  var root = document.documentElement;
  requestAnimationFrame(function () {{ root.classList.add('sv-ready'); }});
  window.setTimeout(function () {{ root.classList.add('sv-ready'); }}, 700);

  /* Features name themselves at once, rather than after the delay a native
     tooltip takes.  The text is the glyph's own <title>, so the two cannot
     disagree, and it is set as textContent so a feature named after a tag
     cannot inject markup. */
  var tip = document.querySelector('.sv-tip');
  function moveTip(event) {{
    var pad = 14;
    var box = tip.getBoundingClientRect();
    var x = event.clientX + pad;
    var y = event.clientY + pad;
    if (x + box.width > window.innerWidth) {{ x = event.clientX - box.width - pad; }}
    if (y + box.height > window.innerHeight) {{ y = event.clientY - box.height - pad; }}
    tip.style.left = Math.max(4, x) + 'px';
    tip.style.top = Math.max(4, y) + 'px';
  }}
  document.addEventListener('mouseover', function (event) {{
    var host = event.target.closest
      && event.target.closest('.sv-annot g, .svz-feat');
    if (!host) {{ return; }}
    var label = host.querySelector('title');
    if (!label) {{ return; }}
    tip.textContent = label.textContent;
    tip.style.display = 'block';
    tip.setAttribute('aria-hidden', 'false');
    moveTip(event);
  }});
  document.addEventListener('mousemove', function (event) {{
    if (tip.style.display === 'block') {{ moveTip(event); }}
  }});
  document.addEventListener('mouseout', function (event) {{
    var host = event.target.closest
      && event.target.closest('.sv-annot g, .svz-feat');
    if (!host) {{ return; }}
    tip.style.display = 'none';
    tip.setAttribute('aria-hidden', 'true');
  }});

  var link = document.querySelector('.sv-views-fixed a.sv-view');
  if (!link) {{ return; }}
  var still = window.matchMedia
    && window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  link.addEventListener('click', function (event) {{
    if (event.button !== 0 || event.metaKey || event.ctrlKey
        || event.shiftKey || event.altKey) {{ return; }}
    if (still) {{ return; }}
    event.preventDefault();
    var href = link.getAttribute('href');
    root.classList.add('sv-leaving');
    window.setTimeout(function () {{ window.location.href = href; }}, 90);
  }});
}})();
</script>
<script id="{theme.script_id}">
(function () {{
  try {{
    var stored = localStorage.getItem('{theme.storage_key}');
    if (stored === 'dark') {{
      document.documentElement.setAttribute('data-theme', 'dark');
    }}
  }} catch (e) {{}}
}})();
</script>
</body>
</html>"""


def _one_lane(plan: TrackPlan) -> TrackPlan:
    """Flatten *plan* onto a single lane, labels on hover.

    Features are landmarks on this page rather than its subject, and a stack of
    lanes for them costs the height the tracks are read in.  Overlapping
    features are drawn over one another and each keeps its title, so a region
    covered twice reads as covered twice and hovering it names both.
    """
    if not plan.glyphs:
        return plan
    # A label that fits inside its own span is kept; one placed after the glyph
    # would run over whatever follows it now that there is one lane.
    flat = [replace(glyph, lane=0,
                    label_place=glyph.label_place if glyph.label_place == "in"
                    else "")
            for glyph in plan.glyphs]
    return replace(plan, glyphs=flat, lanes=1)


def _thresholds(view: SummaryView) -> str:
    """A collapsed note on what the page draws and what it called.

    The thresholds come from the view rather than from this module's defaults,
    so the note states the ones the page was actually built under.
    """
    rows = [
        ("Bars", "The share of covering reads that disagree with the reference "
                 "at each position. Counted per position, not per variant: a "
                 "position the reads disagree about in several directions is "
                 "reported at the total. A deletion counts as disagreement; an "
                 "uncovered position is not counted."),
        ("Bar colour", f"Muted below {MISMATCH_ALERT:.0%}, which is where every "
                       "position sits from sequencing error alone. Coloured at "
                       "or above it."),
        ("Bar height", "The worst position in each pixel column, not the mean. "
                       "Over a reference of several kilobases a column spans "
                       "several bases, and the mean of one disagreeing base "
                       "and its quiet neighbours is the rate of neither."),
        ("Positions to check",
         f"Every stretch disagreeing at {MISMATCH_ALERT:.0%} or more, worst "
         "first, whether or not a variant was called there. Adjacent positions "
         "are joined into one stretch."),
        ("Marks", "One flat mark per called position, sized the same "
                  "regardless of frequency. Hovering one reports the position, "
                  "the change, the reads supporting it, and its effect."),
        ("Variant calling",
         f"A variant needs {view.min_fraction:.0%} of covering reads, "
         f"{view.min_count} supporting read"
         f"{'' if view.min_count == 1 else 's'}, and a depth of "
         f"{view.min_depth}. Calling is per variant, so a position the reads "
         "disagree about in several directions can show a tall bar and call "
         "nothing. Chimeric reads and basecalling error both produce that "
         "pattern."),
        ("Reads", "Coverage, labelled with the depth it tops out at, on a log "
                  "scale. One region drawing many times the reads of the rest "
                  "-- primer dimer, a truncated product -- otherwise sets the "
                  "ceiling and presses the rest of the profile flat."),
    ]
    items = "".join(
        f'<div class="sv-note-row"><dt>{term}</dt><dd>{text}</dd></div>'
        for term, text in rows
    )
    return (
        '<details class="sv-note">'
        "<summary>How to read this page</summary>"
        f'<dl class="sv-note-body">{items}</dl>'
        "</details>"
    )


def _summary_title(title: str) -> str:
    """Name the page for what it is.

    The view is built from a pileup and carries its title, which the CLI
    composes as "Pileup: <reference>".  A caller's own ``--title`` is left as
    given.
    """
    if title.startswith("Pileup: "):
        return "Summary: " + title[len("Pileup: "):]
    return title


def _banner(counterpart: Optional[str]) -> str:
    """Which of the two pages this is, and a link to the other.

    Two pages are written for one alignment and they carry the same title, so
    the view is named rather than left to be inferred from the drawing.  Where
    no pileup was written the name stands alone.
    """
    current = '<span class="sv-view" aria-current="page">Summary</span>'
    other = ""
    if counterpart:
        other = (f'<a class="sv-view" href="{_e(counterpart)}" '
                 f'title="Every read, base by base, with the reads this page '
                 f'reduces.">Pileup</a>')
    return (f'<div class="sv-views-fixed"><div class="sv-views">'
            f"{current}{other}</div></div>")


def render_summary(view: SummaryView, max_lanes: int = 2,
                   pileup_href: Optional[str] = None) -> str:
    """Render *view* to a complete HTML document and return it as a string.

    Args:
        view: The reduction to draw.  Build one from an existing pileup with
            :meth:`~seqviewer.summary.SummaryView.from_view`.
        max_lanes: Feature lanes considered when placing the annotation track.
            Every glyph is then flattened onto one lane, so this only bounds
            what is kept: a feature that fits in no lane is named under the map
            rather than silently dropped.
        pileup_href: Relative link to the full pileup for the same alignment,
            when one was written beside this page.

    Returns:
        A self-contained page: no external stylesheets, scripts, or fonts.
    """
    cell_w = WIDTH / view.ref_len if view.ref_len else 1.0
    palette = _palette()

    plan = _one_lane(plan_track(view.features, view.ref_len, cell_w=cell_w,
                                max_lanes=max_lanes))

    # --- the map: one ruler and one annotation track for the whole page ---
    map_height = RULER_H
    map_parts = [
        f'<svg class="sv-map" viewBox="0 0 {WIDTH:.0f} {map_height}" '
        f'width="{WIDTH:.0f}" height="{map_height}" role="img" '
        f'aria-label="reference ruler">'
    ]
    map_parts += _ruler_parts(view.ref_len, cell_w, RULER_H)
    map_parts.append("</svg>")

    annotations = track_svg(plan)
    dropped = ""
    if plan.dropped:
        names = ", ".join(sorted({f.label or f.type for f in plan.dropped}))
        dropped = (
            f'<p class="sv-none">{len(plan.dropped)} feature(s) not drawn for '
            f"want of lanes: {_e(names)}.</p>"
        )

    highlighted = ""
    if view.highlight_ids:
        highlighted = (
            f'<div class="sv-sub">{_e(view.highlight_label)}: '
            f"{_e(', '.join(view.highlight_ids))}</div>"
        )

    # No heading: the group's own line carries its name, its reads and its
    # depth, and the document title names the reference for a bookmark or a tab.
    head = f"{_banner(pileup_href)}{highlighted}{dropped}"

    # The reference is drawn once and placed between each group's disagreement
    # track and its reads, so the coordinate both are read against sits between
    # them instead of at the top of the page.  Features carry their own titles,
    # so hovering one names it and no legend is spent on the vocabulary.
    reference_svg = "".join(map_parts) + annotations

    # One depth scale across every group, so their profiles are comparable.
    ceiling = max((g.max_depth for g in view.groups), default=0)

    sections = []
    for index, group in enumerate(view.groups):
        star = ' <span class="sv-star" title="Highlighted">&#9733;</span>' \
            if group.highlighted else ""
        status = (
            f'<span class="sv-facts">{_e(group.status)}</span>'
            if group.status else ""
        )
        runs = flagged_runs(group)
        sections.append(
            '<div class="sv-group">'
            f'<div class="sv-group-head"><span class="sv-name">'
            f"{_e(group.name)}{star}</span>{_chip(group)}{status}</div>"
            f'<div class="sv-facts">{_facts(group, view)}</div>'
            f'<div class="sv-stack">'
            f"{_upper_svg(group, view, cell_w)}"
            f'<div class="sv-ref"><span class="sv-ref-label">'
            f"{ROW_LABELS[1]}</span>"
            f'<div class="sv-ref-body">{reference_svg}</div></div>'
            f"{_lower_svg(group, view, cell_w, ceiling)}"
            f"</div>"
            f'<details class="sv-calls">'
            f'<summary>{_call_count(group, runs)}</summary>'
            f'<div class="sv-calls-body">'
            f"{_inspectors(group, view, index, runs)}</div>"
            f"</details>"
            "</div>"
        )

    body = (head + '<hr class="sv-sep">'.join(sections)
            + _thresholds(view))

    return _shell(view, palette, body, track_style(plan),
                  window_css(token_prefix=view.theme.css_prefix))
