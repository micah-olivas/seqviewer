"""Lay annotations out over the reference bar and draw them as SVG.

The pileup itself is a pixel matrix — one rect per base per read — and canvas is
the right tool for that.  Feature glyphs are the opposite kind of drawing: a few
dozen outlined shapes with text in them, wanted crisp at any pixel ratio and
hoverable.  So this track is SVG, which buys real strokes, real text metrics,
native tooltips from ``<title>``, and per-theme fills without a JS colour table.

Everything here is pure geometry over a list of :class:`Feature`, which means the
packing and the label rules are testable without a browser.  The one number
shared with the canvas is :func:`cell_width`; both sides read it from here so a
glyph cannot land a pixel off the column it describes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .genbank import feature_spans
from .model import Feature

__all__ = [
    "cell_width", "plan_track", "track_svg", "TrackPlan", "Glyph",
    "LANE_HEIGHT", "FEATURE_PALETTE",
]

#: Height of one feature glyph, and the space between lanes.  Sized to hold a
#: 10px label with room to breathe, which is what makes it read as a plasmid
#: workspace rather than a heat-map row.
LANE_HEIGHT = 16
LANE_GAP = 4

#: Label type size, and the advance width of one character at that size.  A
#: monospace face is used deliberately: its advance is predictable, so labels can
#: be measured here instead of in the browser.
LABEL_SIZE = 10
_MONO_ADVANCE = 0.60

#: Padding inside a glyph before its label, and beside it for an outside label.
_LABEL_INSET = 4
_LABEL_OUTSET = 5

#: Longest arrowhead, in px.  Clamped to half the glyph so a short feature
#: degrades to a triangle rather than growing a head wider than its body.
_MAX_HEAD = 6

#: Fallback fills for files that carry no colour of their own — NCBI records
#: usually don't.  A file's own colour always wins over these.
FEATURE_PALETTE: Dict[str, Tuple[str, str]] = {
    "cds": ("#c9a227", "#d9b84a"),
    "promoter": ("#4a9d5f", "#5cba77"),
    "terminator": ("#b5452f", "#e0705a"),
    "origin": ("#7b5ea7", "#a68cd4"),
    "rbs": ("#2f7fa8", "#5aa9d6"),
    "bind": ("#c2708f", "#dd93ae"),
    "tag": ("#d2792b", "#e79a55"),
    "insert": ("#3f7f93", "#63a8bd"),
    "other": ("#7f8792", "#9aa3ad"),
}

#: GenBank feature types mapped onto those fills.
_TYPE_KEYS = {
    "cds": "cds", "orf": "cds", "gene": "cds", "mat_peptide": "cds",
    "sig_peptide": "cds",
    "promoter": "promoter", "regulatory": "promoter", "enhancer": "promoter",
    "terminator": "terminator", "polya_signal": "terminator",
    "rep_origin": "origin", "oric": "origin",
    "rbs": "rbs", "5'utr": "rbs", "3'utr": "rbs",
    "protein_bind": "bind", "primer_bind": "bind", "misc_binding": "bind",
    "insert": "insert",
}

#: Kept when lanes run out, most worth showing first.
_TYPE_PRIORITY = {
    "cds": 0, "insert": 0, "promoter": 1, "terminator": 1, "origin": 1,
    "rbs": 2, "tag": 2, "bind": 3, "other": 3,
}


def cell_width(n_cols: int) -> int:
    """Pixels per reference base.

    The canvas reads this too, rather than recomputing the rule in JavaScript,
    so the glyph track and the columns beneath it cannot disagree.
    """
    if n_cols < 200:
        return 4
    if n_cols < 500:
        return 3
    return 2


def _palette_key(feature_type: str) -> str:
    return _TYPE_KEYS.get(feature_type.strip().lower(), "other")


def _label_width(text: str) -> float:
    return len(text) * LABEL_SIZE * _MONO_ADVANCE


def _rgb(color: str) -> Tuple[int, int, int]:
    value = color.lstrip("#")
    if len(value) == 3:
        value = "".join(ch * 2 for ch in value)
    return int(value[0:2], 16), int(value[2:4], 16), int(value[4:6], 16)


def _hex(rgb: Sequence[float]) -> str:
    return "#" + "".join(f"{max(0, min(255, round(c))):02x}" for c in rgb)


def _luminance(color: str) -> float:
    """Perceived brightness on 0..1, the coefficients everyone uses for this."""
    r, g, b = _rgb(color)
    return (0.299 * r + 0.587 * g + 0.114 * b) / 255


def _mix(color: str, target: Sequence[int], amount: float) -> str:
    return _hex([c + (t - c) * amount for c, t in zip(_rgb(color), target)])


def _outline(color: str) -> str:
    """A darker edge for *color*, so a fill near the page ground still has shape."""
    return _mix(color, (0, 0, 0), 0.35)


def _for_dark(color: str) -> str:
    """Lift a fill that would disappear against the dark ground.

    A colour chosen in SnapGene was chosen against white.  Deep navy on a deep
    navy page is invisible, so anything below the floor is mixed toward white
    until it reads — while keeping its hue, which is the point of honouring the
    file at all.
    """
    if _luminance(color) >= 0.38:
        return color
    return _mix(color, (255, 255, 255), 0.45)


def _text_on(color: str) -> str:
    """Black or white, whichever the fill can carry."""
    return "#111111" if _luminance(color) >= 0.55 else "#ffffff"


@dataclass
class Glyph:
    """One drawn piece of one feature, already placed."""

    feature: Feature
    x: float
    width: float
    lane: int
    #: False on an edge that is a cut rather than the feature's real end — the
    #: second piece of an origin-crossing feature, or a clip to the reference.
    head_start: bool
    head_end: bool
    fill_light: str
    fill_dark: str
    label: str
    #: Where the label goes: "in", "after", or "" for hover-only.
    label_place: str
    title: str
    style_index: int = 0

    @property
    def strand(self) -> Optional[int]:
        return self.feature.strand


@dataclass
class TrackPlan:
    """A laid-out annotation track."""

    glyphs: List[Glyph]
    lanes: int
    dropped: List[Feature]
    width: float

    @property
    def height(self) -> int:
        if not self.lanes:
            return 0
        return self.lanes * LANE_HEIGHT + (self.lanes - 1) * LANE_GAP


def _describe(feature: Feature, ref_len: int, pieces: int) -> str:
    """Tooltip text: what it is, where it is, and which way it points."""
    arrow = {1: " →", -1: " ←"}.get(feature.strand, "")
    label = feature.label or feature.type
    if feature.wraps_origin or feature.start > feature.end:
        where = f"{feature.start + 1}..{ref_len}, 1..{feature.end}, crosses the origin"
    else:
        where = f"{feature.start + 1}..{feature.end}"
    bases = sum(e - s for s, e in feature_spans(feature, ref_len))
    detail = f"{feature.type} · {where} · {bases} bp"
    if pieces > 1:
        detail += " · drawn in 2 pieces"
    return f"{label}{arrow}\n{detail}"


def plan_track(
    features: Sequence[Feature],
    ref_len: int,
    cell_w: Optional[int] = None,
    max_lanes: int = 3,
) -> TrackPlan:
    """Place *features* into lanes over a reference of *ref_len* bases.

    Packing is greedy first-fit over features sorted by start, which is optimal
    for interval graphs, with ties broken by descending length so the longest
    feature of a cluster takes the top lane.  A feature that will not fit within
    *max_lanes* is dropped rather than silently growing the page, and comes back
    in :attr:`TrackPlan.dropped` so the caller can say so.
    """
    if cell_w is None:
        cell_w = cell_width(ref_len)
    width = ref_len * cell_w

    # Each feature becomes one or two spans: two when it crosses the origin.
    entries = []
    for feature in features:
        spans = feature_spans(feature, ref_len)
        if not spans:
            continue
        entries.append((feature, spans))

    entries.sort(key=lambda e: (
        min(s for s, _ in e[1]),
        -sum(end - s for s, end in e[1]),
    ))

    # A gap in pixels rather than bases: at 2px per base, a one-base gap is not
    # a visible seam between two glyphs.
    gap_px = 3.0
    lane_ends: List[float] = []
    glyphs: List[Glyph] = []
    dropped: List[Feature] = []
    style_index = 0

    for feature, spans in entries:
        pixel_spans = [(s * cell_w, e * cell_w) for s, e in spans]
        left = min(s for s, _ in pixel_spans)
        right = max(e for _, e in pixel_spans)

        lane = None
        for index, end in enumerate(lane_ends):
            if left >= end:
                lane = index
                break
        if lane is None:
            if len(lane_ends) < max_lanes:
                lane_ends.append(0.0)
                lane = len(lane_ends) - 1
            else:
                dropped.append(feature)
                continue
        lane_ends[lane] = right + gap_px

        key = _palette_key(feature.type)
        fallback_light, fallback_dark = FEATURE_PALETTE[key]
        if feature.color:
            fill_light = feature.color
            fill_dark = _for_dark(feature.color)
        else:
            fill_light, fill_dark = fallback_light, fallback_dark

        label = feature.label or feature.type
        title = _describe(feature, ref_len, len(pixel_spans))

        for piece, (start_px, end_px) in enumerate(sorted(pixel_spans)):
            wraps = len(pixel_spans) > 1
            # On a wrapped feature only the true ends get an arrowhead; the cuts
            # at base 1 and at the end of the reference are flat.
            first_piece = piece == 0
            glyphs.append(Glyph(
                feature=feature,
                x=start_px,
                width=end_px - start_px,
                lane=lane,
                head_start=not wraps or not first_piece,
                head_end=not wraps or first_piece,
                fill_light=fill_light,
                fill_dark=fill_dark,
                label=label,
                label_place="",
                title=title,
                style_index=style_index,
            ))
        style_index += 1

    _place_labels(glyphs, width)
    return TrackPlan(glyphs=glyphs, lanes=len(lane_ends), dropped=dropped,
                     width=width)


def _place_labels(glyphs: List[Glyph], width: float) -> None:
    """Decide, per glyph, whether its label goes inside, after, or nowhere.

    A label inside the glyph is the best outcome and the one a plasmid map leads
    with.  Failing that it sits just after the glyph, but only in the gap before
    the next thing in that lane, so labels never collide.  Anything narrower than
    that keeps its tooltip and no visible text — three characters and an ellipsis
    tells a reader less than the hover does.
    """
    by_lane: Dict[int, List[Glyph]] = {}
    for glyph in glyphs:
        by_lane.setdefault(glyph.lane, []).append(glyph)

    # One label per feature, on its widest piece.
    widest: Dict[int, Glyph] = {}
    for glyph in glyphs:
        current = widest.get(glyph.style_index)
        if current is None or glyph.width > current.width:
            widest[glyph.style_index] = glyph
    labelled = set(id(g) for g in widest.values())

    for lane_glyphs in by_lane.values():
        lane_glyphs.sort(key=lambda g: g.x)
        for index, glyph in enumerate(lane_glyphs):
            if id(glyph) not in labelled:
                continue
            text = _label_width(glyph.label)
            if text + 2 * _LABEL_INSET <= glyph.width:
                glyph.label_place = "in"
                continue
            following = lane_glyphs[index + 1].x if index + 1 < len(lane_glyphs) else width
            room = following - (glyph.x + glyph.width) - _LABEL_OUTSET
            if room >= text:
                glyph.label_place = "after"


def _path(glyph: Glyph) -> str:
    """The glyph outline: a block, with a point on whichever end is a real end."""
    height = LANE_HEIGHT
    x, w = glyph.x, glyph.width
    mid = height / 2
    strand = glyph.strand
    head = min(_MAX_HEAD, w / 2) if strand in (1, -1) else 0.0
    if strand == 1 and glyph.head_end and head > 0:
        return (f"M{x:.1f},0 H{x + w - head:.1f} L{x + w:.1f},{mid:.1f} "
                f"L{x + w - head:.1f},{height} H{x:.1f} Z")
    if strand == -1 and glyph.head_start and head > 0:
        return (f"M{x + w:.1f},0 H{x + head:.1f} L{x:.1f},{mid:.1f} "
                f"L{x + head:.1f},{height} H{x + w:.1f} Z")
    return f"M{x:.1f},0 H{x + w:.1f} V{height} H{x:.1f} Z"


def _escape(text: str) -> str:
    return (text.replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def track_style(plan: TrackPlan, prefix: str = "svf") -> str:
    """CSS for the track's fills, both themes.

    A feature's colour comes from the file, so it is a literal rather than a
    theme token and cannot be swapped by redefining a custom property.  Emitting
    one rule per feature per theme is what keeps a colour chosen against white
    from disappearing on the dark ground.
    """
    seen = {}
    for glyph in plan.glyphs:
        seen[glyph.style_index] = glyph
    lines = []
    for index, glyph in sorted(seen.items()):
        light, dark = glyph.fill_light, glyph.fill_dark
        lines.append(
            f".{prefix}{index}{{fill:{light};stroke:{_outline(light)}}}"
            f".{prefix}t{index}{{fill:{_text_on(light)}}}"
        )
        lines.append(
            f'[data-theme="dark"] .{prefix}{index}'
            f"{{fill:{dark};stroke:{_outline(dark)}}}"
            f'[data-theme="dark"] .{prefix}t{index}{{fill:{_text_on(dark)}}}'
        )
    return "\n".join(lines)


def track_svg(plan: TrackPlan, prefix: str = "svf") -> str:
    """The track as an ``<svg>`` element, or "" when there is nothing to draw."""
    if not plan.glyphs:
        return ""
    height = plan.height
    parts = [
        f'<svg class="sv-annot" width="{plan.width:.0f}" height="{height}" '
        f'viewBox="0 0 {plan.width:.0f} {height}" '
        f'role="img" aria-label="reference annotations">'
    ]
    for glyph in plan.glyphs:
        top = glyph.lane * (LANE_HEIGHT + LANE_GAP)
        parts.append(f'<g transform="translate(0,{top})">')
        parts.append(f'<title>{_escape(glyph.title)}</title>')
        parts.append(
            f'<path class="{prefix}{glyph.style_index}" d="{_path(glyph)}"/>'
        )
        if glyph.label_place == "in":
            parts.append(
                f'<text class="{prefix}t{glyph.style_index}" '
                f'x="{glyph.x + glyph.width / 2:.1f}" y="{LANE_HEIGHT / 2:.1f}" '
                f'text-anchor="middle" dominant-baseline="central">'
                f"{_escape(glyph.label)}</text>"
            )
        elif glyph.label_place == "after":
            parts.append(
                f'<text class="sv-annot-out" '
                f'x="{glyph.x + glyph.width + _LABEL_OUTSET:.1f}" '
                f'y="{LANE_HEIGHT / 2:.1f}" dominant-baseline="central">'
                f"{_escape(glyph.label)}</text>"
            )
        parts.append("</g>")
    parts.append("</svg>")
    return "".join(parts)
