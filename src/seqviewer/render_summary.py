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

import html as _html
from typing import List, Optional, Sequence, Tuple

from .annotate import plan_track, track_style, track_svg
from .summary import SEVERITY, GroupSummary, SummaryView, Variant
from .zoom import window_bounds, window_css, window_svg

# The pileup owns the package's colour table.  Imported rather than copied so
# the two views cannot drift apart; if it moves, this fails loudly at import
# time, which is the failure worth having.
from .render import _PALETTE as _PILEUP_PALETTE

__all__ = ["render_summary"]

#: Nominal drawing width.  Everything is laid out in this coordinate space and
#: the SVG is scaled to its container by its ``viewBox``, so the page is
#: responsive without measuring anything in the browser.
WIDTH = 1000.0

#: Narrowest the drawing is allowed to be scaled to before its container scrolls
#: instead.  A viewBox scales type along with geometry, so without a floor a
#: narrow window shrinks the whole band into illegibility.
MIN_DRAW_WIDTH = 680

RULER_H = 18
LOLLI_H = 64
RIBBON_H = 11
DEPTH_H = 30

#: Separation between the reference ribbon and the coverage profile, so a group
#: at full depth does not read as one thick bar.
RIBBON_GAP = 3

#: Lollipop head radius, and the stem lengths a 0% and a 100% allele get.
HEAD_R = 5.0
STEM_MIN = 14.0
STEM_SPAN = 24.0

#: Heads closer together than this are lifted onto another tier rather than
#: drawn on top of each other, up to this many tiers.
HEAD_SEPARATION = 11.0
TIER_RISE = 9.0
MAX_TIERS = 3

#: Rows of the variant table drawn before it is cut short.
MAX_TABLE_ROWS = 50

#: Base-resolution windows drawn per group, and their size.  Bounded because each
#: window is one SVG text element per base per row: a window is cheap, forty are
#: not, and a summary that carries forty has stopped summarising.
MAX_INSPECTORS = 6
INSPECTOR_COLUMNS = 60
INSPECTOR_READS = 6

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
    pos = step
    while pos < ref_len:
        x = _x(pos, cell_w)
        parts.append(
            f'<line class="sv-tick" x1="{x:.1f}" y1="{height - 4}" '
            f'x2="{x:.1f}" y2="{height}" />'
        )
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
# One group's band: lollipops, reference ribbon, coverage profile
# --------------------------------------------------------------------------

def _tiers(variants: Sequence[Variant], cell_w: float) -> List[int]:
    """Assign each variant a stagger tier so neighbouring heads do not overlap.

    Variants arrive sorted by position, so one pass keeping the last x used on
    each tier is enough.  Past the last tier a head is left to overlap rather
    than climbing out of the band — at that density the table is the readable
    account, and the lollipops are saying "a cluster, here".
    """
    last_x = [-1e9] * MAX_TIERS
    out = []
    for variant in variants:
        x = _x(variant.pos + 0.5, cell_w)
        for tier in range(MAX_TIERS):
            if x - last_x[tier] >= HEAD_SEPARATION:
                last_x[tier] = x
                out.append(tier)
                break
        else:
            out.append(MAX_TIERS - 1)
    return out


def _head(variant: Variant, x: float, y: float, tier_class: str) -> str:
    """The glyph at the top of a stem: shape by kind, colour by consequence."""
    if variant.kind == "del":
        half = HEAD_R + 0.5
        points = (f"{x - half:.1f},{y - half:.1f} {x + half:.1f},{y - half:.1f} "
                  f"{x:.1f},{y + half:.1f}")
        return f'<polygon class="sv-head {tier_class}" points="{points}" />'
    if variant.kind == "ins":
        half = HEAD_R + 0.5
        points = (f"{x - half:.1f},{y + half:.1f} {x:.1f},{y - half:.1f} "
                  f"{x + half:.1f},{y + half:.1f}")
        return f'<polygon class="sv-head {tier_class}" points="{points}" />'
    return (f'<circle class="sv-head {tier_class}" cx="{x:.1f}" cy="{y:.1f}" '
            f'r="{HEAD_R:.1f}" />')


def _lollipop_parts(group: GroupSummary, cell_w: float, baseline: float) -> List[str]:
    """One stem-and-head per called variant, rising from the ribbon at *baseline*."""
    parts: List[str] = []
    tiers = _tiers(group.variants, cell_w)
    for variant, tier in zip(group.variants, tiers):
        x = _x(variant.pos + 0.5, cell_w)
        stem = STEM_MIN + STEM_SPAN * min(1.0, variant.fraction) + TIER_RISE * tier
        top = baseline - stem
        tier_class = f"sv-{_TIER_OF.get(variant.consequence, 'neutral')}"
        title = (
            f"{variant.label} at {variant.pos + 1}\n"
            f"{variant.fraction:.0%} of {variant.depth} reads"
            + (f"\n{variant.effect}" if variant.effect else "")
        )
        parts.append("<g>")
        parts.append(f"<title>{_e(title)}</title>")
        parts.append(
            f'<line class="sv-stem" x1="{x:.1f}" y1="{baseline:.1f}" '
            f'x2="{x:.1f}" y2="{top:.1f}" />'
        )
        parts.append(_head(variant, x, top, tier_class))
        if variant.kind == "snv":
            parts.append(
                # dy in ems rather than dominant-baseline: the latter was
                # unsupported in WebKit for years, and the engines still
                # disagree about "central" versus "middle".
                f'<text class="sv-head-letter" x="{x:.1f}" y="{top:.1f}" '
                f'dy="0.35em" text-anchor="middle">'
                f"{_e(variant.alt)}</text>"
            )
        parts.append("</g>")
    return parts


def _ribbon_parts(group: GroupSummary, view: SummaryView,
                  cell_w: float, y: float) -> List[str]:
    """The reference bar, with the reading frame picked out along it."""
    width = _x(group.ref_len, cell_w)
    parts = [
        f'<rect class="sv-ribbon" x="0" y="{y}" width="{width:.1f}" '
        f'height="{RIBBON_H}" />'
    ]
    if view.focus is not None:
        start, end = view.focus
        left = _x(start, cell_w)
        parts.append(
            f'<rect class="sv-ribbon-focus" x="{left:.1f}" y="{y}" '
            f'width="{_x(end, cell_w) - left:.1f}" height="{RIBBON_H}" />'
        )
    if group.ref_len < view.ref_len:
        # This group is stated against a shorter reference; say so rather than
        # letting its band trail off as though coverage simply ran out.
        parts.append(
            f'<rect class="sv-ribbon-absent" x="{width:.1f}" y="{y}" '
            f'width="{WIDTH - width:.1f}" height="{RIBBON_H}" />'
        )
    return parts


def _depth_parts(group: GroupSummary, cell_w: float, ceiling: int,
                 top: float) -> List[str]:
    """A filled coverage profile under the ribbon.

    Each pixel column reports the *thinnest* coverage it spans, not the mean, so
    a dropout narrower than one pixel still shows as a notch.  On a page whose
    job is to be trusted at a glance, a coverage hole that averages away is the
    failure worth avoiding.
    """
    if not group.depth or not ceiling:
        return []
    base = top + DEPTH_H
    columns = min(int(WIDTH), group.ref_len)
    points = []
    for column in range(columns):
        start = group.ref_len * column // columns
        end = max(start + 1, group.ref_len * (column + 1) // columns)
        value = min(group.depth[start:end])
        x = WIDTH * column / columns
        y = base - DEPTH_H * min(1.0, value / ceiling)
        points.append(f"{x:.1f},{y:.1f}")
    # Close the polygon along the baseline so it fills.
    outline = " ".join(points)
    return [
        f'<polygon class="sv-depth" points="0,{base} {outline} '
        f'{WIDTH:.1f},{base} " />',
        f'<polyline class="sv-depth-edge" points="{outline}" />',
    ]


def _band_svg(group: GroupSummary, view: SummaryView, cell_w: float,
              ceiling: int) -> str:
    """One group's whole band as an ``<svg>`` element.

    A group with nothing called reserves no room for lollipops.  Keeping the
    height uniform would line the ribbons up, but it spends sixty pixels per
    clean group on empty space and leaves the frame's boundary lines hanging in
    it, which reads as a drawing that failed rather than a clone that passed.
    """
    lolli_h = LOLLI_H if group.variants else 0.0
    ribbon_y = lolli_h
    depth_y = ribbon_y + RIBBON_H + RIBBON_GAP
    height = depth_y + DEPTH_H

    parts = [
        f'<svg class="sv-band" viewBox="0 0 {WIDTH:.0f} {height:.0f}" '
        f'width="{WIDTH:.0f}" height="{height:.0f}" role="img" '
        f'aria-label="{_e(group.name)} summary">'
    ]
    parts += _depth_parts(group, cell_w, ceiling, depth_y)
    parts += _ribbon_parts(group, view, cell_w, ribbon_y)
    if group.variants:
        parts += _lollipop_parts(group, cell_w, lolli_h)
    parts += _focus_parts(view.focus, cell_w, height)
    parts.append("</svg>")
    return "".join(parts)


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


def _table(group: GroupSummary) -> str:
    """The called variants, worst first."""
    if not group.variants:
        return '<p class="sv-none">No variants cleared the calling thresholds.</p>'

    def rank(variant: Variant) -> tuple:
        consequence = variant.consequence
        order = (
            list(_TIER_OF).index(consequence)
            if consequence in _TIER_OF else len(_TIER_OF)
        )
        return (order, -variant.fraction, variant.pos)

    ordered = sorted(group.variants, key=rank)
    shown = ordered[:MAX_TABLE_ROWS]

    rows = []
    for variant in shown:
        tier = _TIER_OF.get(variant.consequence, "neutral")
        kind = {"snv": "SNV", "del": "Deletion", "ins": "Insertion"}[variant.kind]
        rows.append(
            "<tr>"
            f'<td><span class="sv-dot sv-{tier}"></span>{kind}</td>'
            f"<td class=\"sv-num\">{variant.pos + 1}</td>"
            f"<td class=\"sv-mono\">{_e(variant.label)}</td>"
            f'<td class="sv-num">{variant.fraction:.0%}</td>'
            f'<td class="sv-num">{variant.count}/{variant.depth}</td>'
            f"<td>{_e(variant.effect) or '&mdash;'}</td>"
            "</tr>"
        )

    more = ""
    if len(ordered) > len(shown):
        more = (
            f'<p class="sv-none">{len(ordered) - len(shown)} further '
            f"variants not listed.</p>"
        )

    return (
        '<table class="sv-table"><thead><tr>'
        "<th>Type</th><th>Position</th><th>Change</th><th>Frequency</th>"
        "<th>Reads</th><th>Effect</th>"
        "</tr></thead><tbody>" + "".join(rows) + "</tbody></table>" + more
    )


def _inspectors(group: GroupSummary, view: SummaryView, index: int) -> str:
    """A base-resolution window per called variant, each collapsed until asked for.

    The map answers "where" and the table answers "what"; this answers "show me
    the actual bases", which is the one question neither can. Windows are
    ``<details>`` rather than JavaScript-driven panels: a summary should stay
    readable with nothing running, and a disclosure widget is the one interaction
    both engines implement identically.
    """
    source = view.source
    if source is None or not group.variants or index >= len(source.groups):
        return ""

    origin = source.groups[index]
    if not origin.rows:
        return ""

    ranked = sorted(
        group.variants,
        key=lambda v: (SEVERITY.index(v.consequence)
                       if v.consequence in SEVERITY else len(SEVERITY),
                       -v.fraction),
    )

    blocks = []
    for variant in ranked[:MAX_INSPECTORS]:
        start, end = window_bounds(variant.pos, len(origin.ref_seq),
                                   INSPECTOR_COLUMNS)
        window = window_svg(
            origin.ref_seq, origin.rows, start, end,
            frame=view.focus, max_read_rows=INSPECTOR_READS,
            label=f"{variant.label} at {variant.pos + 1}",
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
        caption = (
            f"{variant.label} at {variant.pos + 1}"
            + (f" &middot; {_e(variant.effect)}" if variant.effect else "")
        )
        blocks.append(
            "<details class=\"sv-zoom\"><summary>"
            f"{caption} &middot; bases {start + 1}&ndash;{end}"
            "</summary>"
            f'<div class="sv-zoom-body">{window.svg}{hidden}</div>'
            "</details>"
        )

    if not blocks:
        return ""
    return f'<div class="sv-zooms">{"".join(blocks)}</div>'


# --------------------------------------------------------------------------
# The page shell
# --------------------------------------------------------------------------

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
    padding: 1.5rem;
}}
.sv-wrap {{ max-width: 1100px; margin: 0 auto; }}
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
.sv-ribbon {{ fill: var(--{prefix}-ribbon); }}
.sv-ribbon-focus {{ fill: var(--{prefix}-focus); }}
.sv-ribbon-absent {{ fill: var(--{prefix}-grid); opacity: 0.5; }}
.sv-depth {{ fill: var(--{prefix}-depth); }}
.sv-depth-edge {{ fill: none; stroke: var(--{prefix}-depth-edge); stroke-width: 1; }}
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
.sv-table {{
    border-collapse: collapse; font-size: 0.78rem; margin-top: 0.5rem;
    width: 100%; max-width: 640px;
}}
.sv-table th {{
    text-align: left; font-weight: 600; color: var(--muted);
    border-bottom: 1px solid var(--panel-line); padding: 0.25rem 0.6rem 0.25rem 0;
    font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.05em;
}}
.sv-table td {{
    padding: 0.25rem 0.6rem 0.25rem 0;
    border-bottom: 1px solid var(--panel-line);
}}
.sv-num, .sv-mono {{ font-family: var(--mono); }}
.sv-num {{ text-align: right; }}
.sv-dot {{
    display: inline-block; width: 7px; height: 7px; border-radius: 50%;
    margin-right: 0.4rem; vertical-align: middle;
}}
.sv-dot.sv-ok {{ background: var(--{prefix}-ok); }}
.sv-dot.sv-warn {{ background: var(--{prefix}-warn); }}
.sv-dot.sv-bad {{ background: var(--{prefix}-bad); }}
.sv-dot.sv-neutral {{ background: var(--{prefix}-neutral); }}
.sv-none {{ color: var(--muted); font-size: 0.8rem; font-style: italic; }}
.sv-key {{
    display: flex; flex-wrap: wrap; gap: 0.35rem 1rem;
    font-size: 0.75rem; color: var(--muted); margin-top: 0.3rem;
}}
.sv-kitem {{ display: flex; align-items: center; gap: 0.3rem; }}
.sv-sep {{ border: none; border-top: 1px solid var(--panel-line); margin: 1.3rem 0; }}
/* --- base-resolution windows --- */
.sv-zooms {{ margin-top: 0.6rem; }}
.sv-zoom {{
    border: 1px solid var(--panel-line);
    border-radius: 3px;
    margin-bottom: 0.3rem;
    background: var(--card-bg);
}}
.sv-zoom > summary {{
    cursor: pointer;
    padding: 0.3rem 0.5rem;
    font: 600 0.75rem/1.3 var(--mono);
    color: var(--text);
    /* Safari shows a default disclosure marker that ignores list-style; both
       engines honour ::-webkit-details-marker, so it is set for both. */
    list-style: none;
}}
.sv-zoom > summary::-webkit-details-marker {{ display: none; }}
.sv-zoom > summary::before {{
    content: "\\25B8 ";
    color: var(--muted);
}}
.sv-zoom[open] > summary::before {{ content: "\\25BE "; }}
.sv-zoom-body {{
    padding: 0.2rem 0.5rem 0.5rem;
    overflow-x: auto;
}}
{window_style}
{track_css}
</style>
</head>
<body>
<div class="sv-wrap">
{body}
</div>
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


def _key() -> str:
    """A legend for the glyph vocabulary, which is not self-evident."""
    shapes = [
        ("&#9679;", "substitution"),
        ("&#9660;", "deletion"),
        ("&#9650;", "insertion"),
    ]
    tiers = [
        ("bad", "frameshift or stop"),
        ("warn", "missense or in-frame indel"),
        ("ok", "silent"),
        ("neutral", "outside the reading frame"),
    ]
    items = [
        f'<span class="sv-kitem"><span class="sv-mono">{glyph}</span>{label}</span>'
        for glyph, label in shapes
    ]
    items += [
        f'<span class="sv-kitem"><span class="sv-dot sv-{tier}"></span>{label}</span>'
        for tier, label in tiers
    ]
    items.append(
        '<span class="sv-kitem">stem height &#8733; allele frequency</span>'
    )
    return f'<div class="sv-key">{"".join(items)}</div>'


def render_summary(view: SummaryView, max_lanes: int = 2) -> str:
    """Render *view* to a complete HTML document and return it as a string.

    Args:
        view: The reduction to draw.  Build one from an existing pileup with
            :meth:`~seqviewer.summary.SummaryView.from_view`.
        max_lanes: Feature lanes over the map.  Lower than the pileup page's
            default on purpose — a summary that grows a tall annotation stack
            has stopped being one.  Features that do not fit are named under the
            map rather than silently dropped.

    Returns:
        A self-contained page: no external stylesheets, scripts, or fonts.
    """
    cell_w = WIDTH / view.ref_len if view.ref_len else 1.0
    palette = _palette()

    plan = plan_track(view.features, view.ref_len, cell_w=cell_w,
                      max_lanes=max_lanes)

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
            f" &middot; {_e(view.highlight_label)}: "
            f"{_e(', '.join(view.highlight_ids))}"
        )

    head = (
        f"<h1>{_e(view.title)}</h1>"
        f'<div class="sv-sub">{view.total_reads} total reads &middot; '
        f"{len(view.groups)} group(s) &middot; {view.ref_len} bp"
        f"{highlighted}</div>"
        f'<div class="sv-panel">'
        f'<div class="sv-eyebrow">Reference</div>'
        f'<div class="sv-scroll">{"".join(map_parts)}{annotations}</div>'
        f"{dropped}"
        f"</div>"
    )

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
        sections.append(
            '<div class="sv-group">'
            f'<div class="sv-group-head"><span class="sv-name">'
            f"{_e(group.name)}{star}</span>{_chip(group)}{status}</div>"
            f'<div class="sv-facts">{_facts(group, view)}</div>'
            f'<div class="sv-scroll">'
            f"{_band_svg(group, view, cell_w, ceiling)}</div>"
            f"{_table(group)}"
            f"{_inspectors(group, view, index)}"
            "</div>"
        )

    body = head + _key() + '<hr class="sv-sep">' + \
        '<hr class="sv-sep">'.join(sections)

    return _shell(view, palette, body, track_style(plan),
                  window_css(token_prefix=view.theme.css_prefix))
