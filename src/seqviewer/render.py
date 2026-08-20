"""Render a read pileup as a self-contained HTML page.

The page draws an HTML5 canvas matrix: one row per read, one cell per reference
position.  Matches are gray, mismatches take a per-base color, gaps are white.
Above the reads sit a ruler and a consensus row; below them, when the view marks
an insert, the reference and consensus translations of that insert.

Nothing here touches the filesystem, a subprocess, or an aligner.  The input is
a finished grid, so this module has no dependencies outside the standard library
and can be exercised with synthetic data — see :mod:`seqviewer.demo`.
"""

from __future__ import annotations

import html as _html
import json as _json
from importlib import resources as _resources

from .annotate import cell_width as _cell_width
from .annotate import plan_track as _plan_track
from .annotate import track_style as _track_style
from .annotate import region_band_svg as _region_band_svg
from .annotate import region_spans as _region_spans
from .annotate import track_svg as _track_svg
from .annotate import mismatch_track_svg as _mismatch_track_svg
from .annotate import MISMATCH_TRACK_HEIGHT as _MISMATCH_TRACK_HEIGHT
from .codon import translate as _translate
from .pileup import PileupView

__all__ = ["render"]


def _asset(name: str, substitutions=None) -> str:
    """Read a file from the package's ``assets`` directory.

    The stylesheet and the drawing code are real files rather than text inside
    this module's f-string.  An f-string demands every brace be doubled, which
    over 380 lines of JavaScript and 270 of CSS made each one a place to make a
    mistake, and no editor would lint either.

    Assets are inlined rather than linked, so a rendered page stays
    self-contained.  *substitutions* are literal token replacements: the files
    carry ``__PREFIX__`` and ``__STORAGE_KEY__`` rather than ``{}``
    placeholders, so their braces stay real CSS and JavaScript braces.  The
    trailing newline is dropped because the line the asset is spliced into
    supplies one.
    """
    text = _resources.files(__package__).joinpath("assets", name).read_text(
        encoding="utf-8")
    for token, value in (substitutions or {}).items():
        text = text.replace(token, value)
    return text.rstrip("\n")


#: Every colour the page draws, per theme.  The legend swatches read these as
#: CSS custom properties and the canvas reads the same values from an emitted
#: JS object, so a swatch can no longer disagree with the cell it describes.
_PALETTE = {
    "light": {
        "match": "#d4d8dc",
        "vector": "#e9ebee",
        "gap": "#ffffff",
        "ref": "#1e293b",
        "a": "#e03131",
        "t": "#1971c2",
        "c": "#e8590c",
        "g": "#e67700",
        "boundary": "#d97706",
        "region": "rgba(56,132,255,0.10)",
        "tick": "#94a3b8",
        "tick-label": "#1e293b",
        "flag": "#94a3b8",
        "aa-match": "#d1d5db",
        "aa-diff": "#e03131",
        "aa-diff-bg": "rgba(224,49,49,0.10)",
        "aa-bg": "#f8fafc",
        "aa-grid": "#e5e7eb",
    },
    "dark": {
        "match": "#4a5568",
        "vector": "#3a4455",
        "gap": "#ffffff",
        "ref": "#e0e0e0",
        "a": "#ff6b6b",
        "t": "#339af0",
        "c": "#ffa94d",
        "g": "#ffd43b",
        "boundary": "#f59e0b",
        "region": "rgba(120,170,255,0.16)",
        "tick": "#64748b",
        "tick-label": "#e0e0e0",
        "flag": "#64748b",
        "aa-match": "#4a5568",
        "aa-diff": "#ff6b6b",
        "aa-diff-bg": "rgba(255,107,107,0.18)",
        "aa-bg": "#1e293b",
        "aa-grid": "#334155",
    },
}

#: Substrings matched against a group's ``status``, in order, to pick the tier
#: it is drawn in.  "mismatch" precedes "match" so the longer word wins.
_TIERS = (
    ("perfect", "ok", "●"),
    ("silent", "warn", "◐"),
    ("missense", "bad", "▲"),
    ("nonsense", "bad", "▲"),
    ("frameshift", "bad", "▲"),
    ("stop", "bad", "▲"),
    ("mismatch", "bad", "▲"),
    ("match", "ok", "●"),
)


def _verdict(status: str, highlighted: bool):
    """Classify *status* into a ``(tier, glyph)`` pair, or None to show nothing.

    A blank status means the caller has no call to report, and renders no chip
    at all rather than an empty one.  A status the renderer has no opinion about
    gets the neutral tier: reporting something unrecognised should not paint a
    group as a failure, which is what the previous fall-through to the error
    colour did.
    """
    text = status.strip().lower()
    if not text:
        return None
    for needle, tier, glyph in _TIERS:
        if needle in text:
            return tier, glyph
    return ("ok", "●") if highlighted else ("neutral", "○")


def _chip(status: str, highlighted: bool) -> str:
    """Render *status* as a verdict chip, or "" when there is nothing to say."""
    verdict = _verdict(status, highlighted)
    if verdict is None:
        return ""
    tier, glyph = verdict
    return (
        f'<span class="sv-chip sv-chip-{tier}">'
        f'<span class="sv-glyph">{glyph}</span>'
        f'{_html.escape(status)}</span>'
    )


def render(view: PileupView) -> str:
    """Render *view* to a complete HTML document and return it as a string."""
    # The body below predates the typed model and reads its groups as plain
    # dicts.  Adapting here keeps the migrated renderer byte-for-byte identical
    # to the original while the public surface is dataclasses.
    groups = [
        {
            "ref_id": g.name,
            "n_reads": g.n_reads,
            "frac": g.fraction,
            "status": g.status,
            "is_recoverable": g.highlighted,
            "ref_seq": g.ref_seq,
            "pileup_rows": g.rows,
            "parent": g.parent,
        }
        for g in view.groups
    ]
    flank_lengths = view.flanks
    title = view.title

    _theme = view.theme
    _p = _theme.css_prefix
    _style_id = _theme.style_id
    _script_id = _theme.script_id
    _storage_key = _theme.storage_key

    pileup_css = _asset("pileup.css", {"__PREFIX__": _p})
    pileup_js = _asset("pileup.js")
    theme_js = _asset("theme.js", {"__STORAGE_KEY__": _storage_key})

    flanks_js = "null"
    if flank_lengths and (flank_lengths[0] or flank_lengths[1]):
        flanks_js = f"[{flank_lengths[0]},{flank_lengths[1]}]"

    single = len(groups) == 1
    translate = view.translate
    has_flanks = bool(flank_lengths and (flank_lengths[0] or flank_lengths[1]))
    any_flags = False
    any_starred = any(g["is_recoverable"] for g in groups)

    # The reference's geometry belongs to the run, not to a group, so with
    # several groups over one reference it is stated once in the masthead
    # instead of repeated in every group header.
    ref_lens = {len(g["ref_seq"]) for g in groups}
    hoist_geometry = not single and len(ref_lens) == 1

    # A fact states its number tersely and says what the number *is* on hover.
    # The tip defines the quantity and never restates the value: the line is
    # still the only place each fact is stated, and a tip cannot drift from a
    # number it does not contain.
    def _fact(text: str, tip: str) -> str:
        return f'<span class="sv-fact" title="{_html.escape(tip)}">{text}</span>'

    sections_html = []
    metas = []
    # One plan per group: groups may have different reference lengths, so each
    # gets its own layout and its own CSS prefix, keeping the per-feature rules
    # from colliding between groups.
    annot_plans = []
    for idx, g in enumerate(groups):
        star = ' <span class="sv-star" title="Highlighted">&#9733;</span>' \
            if g["is_recoverable"] else ""

        # Per-read identity across every called base in the group.
        identity = None
        if g["pileup_rows"]:
            total_bases = 0
            total_matches = 0
            for row in g["pileup_rows"]:
                aligned = [(b, m) for b, m in row if b != "-"]
                total_bases += len(aligned)
                total_matches += sum(1 for _, m in aligned if m)
            if total_bases > 0:
                identity = total_matches / total_bases

        ref_len = len(g["ref_seq"])
        drawn = len(g["pileup_rows"])

        # Each number appears once: the reads figure is a ratio rather than a
        # count and a separate percentage, and "drawn" shows up only when it
        # differs from the assigned count.
        geometry = [_fact(
            f"<b>{ref_len}</b> bp reference",
            "Length of the reference these reads are aligned to.",
        )]
        if has_flanks:
            geometry.append(_fact(
                f"insert <b>{flank_lengths[0] + 1}</b>&ndash;"
                f"<b>{ref_len - flank_lengths[1]}</b>",
                "Where the insert starts and ends in reference coordinates, "
                "counting from base 1 and including both ends. These are "
                "positions, not a count of bases. The flanks outside them are "
                "vector.",
            ))

        counts = []
        if g["n_reads"] and view.total_reads:
            counts.append(_fact(
                f"<b>{g['n_reads']}</b> of <b>{view.total_reads}</b> reads",
                "Reads assigned to this group, out of the page's total across "
                "all groups.",
            ))
        elif g["n_reads"]:
            counts.append(_fact(
                f"<b>{g['n_reads']}</b> reads",
                "Reads assigned to this group.",
            ))
        if drawn != g["n_reads"]:
            # Why this differs from the assigned count is the caller's business
            # — the grid arrives already built — so the tip names the quantity
            # without inventing a reason for the gap.
            counts.append(_fact(
                f"<b>{drawn}</b> rows drawn",
                "How many of this group's reads are plotted as rows below. "
                "Shown only when it differs from the number assigned.",
            ))
        if identity is not None:
            counts.append(_fact(
                f"<b>{identity:.1%}</b> identity",
                "Share of the bases the drawn reads actually called that match "
                "the reference. Gaps and uncovered positions count for neither "
                "side, so a read spanning half the reference with no mismatches "
                "is 100%.",
            ))

        facts = counts if hoist_geometry else geometry + counts

        metas.append({
            "name_html": f'{_html.escape(g["ref_id"])}{star}',
            "chip": _chip(g["status"], g["is_recoverable"]),
            "facts": " &middot; ".join(facts),
            "geometry": " &middot; ".join(geometry),
        })

        header = "" if single else (
            f'<div class="sv-group-head">'
            f'<span class="sv-group-name">{metas[-1]["name_html"]}</span>'
            f'{metas[-1]["chip"]}'
            f'<span class="sv-facts">{metas[-1]["facts"]}</span>'
            f'</div>'
        )

        # Encode pileup data compactly for JS:
        # '.' = match, base letter = mismatch, '-' = gap
        rows_encoded = []
        for row in g["pileup_rows"]:
            chars = []
            for base_char, is_match in row:
                if base_char == "-":
                    chars.append("-")
                elif is_match:
                    chars.append(".")
                else:
                    chars.append(base_char.upper())
            rows_encoded.append("".join(chars))

        # Build consensus from pileup: majority base at each position.
        # Also track positions where >10% of reads disagree with the
        # reference — these are flagged in the ruler as problem positions.
        from collections import Counter
        consensus_encoded = []
        flagged_cols = []
        mismatch_freqs = []
        ref_seq = g["ref_seq"]
        _MISMATCH_THRESHOLD = 0.10
        for col_idx in range(ref_len):
            counts = Counter()
            for row in g["pileup_rows"]:
                base, _ = row[col_idx]
                if base != "-":
                    counts[base.upper()] += 1
            total = sum(counts.values())
            freq = 0.0
            if counts:
                cons_base = counts.most_common(1)[0][0]
                if cons_base == ref_seq[col_idx].upper():
                    consensus_encoded.append(".")
                else:
                    consensus_encoded.append(cons_base)
                # Flag if >10% of reads differ from reference
                ref_base = ref_seq[col_idx].upper()
                ref_count = counts.get(ref_base, 0)
                freq = (total - ref_count) / total if total else 0.0
                if freq > _MISMATCH_THRESHOLD:
                    flagged_cols.append(col_idx)
            else:
                consensus_encoded.append("-")
            mismatch_freqs.append(freq)
        consensus_str = "".join(consensus_encoded)
        any_flags = any_flags or bool(flagged_cols)

        # Translate the focus region, when the view asks for a protein readout.
        # A focus region worth marking is not always a reading frame worth
        # translating, so this is gated on `translate` and not on `flanks`.
        cons_protein = ref_protein = ""
        if translate and has_flanks:
            consensus_dna = "".join(
                ref_seq[i] if c == "." else (c if c != "-" else "N")
                for i, c in enumerate(consensus_encoded)
            )
            _ins_start = flank_lengths[0]
            _ins_end = ref_len - flank_lengths[1]
            cons_protein = _translate(consensus_dna[_ins_start:_ins_end])
            ref_protein = _translate(ref_seq[_ins_start:_ins_end])

        # The parent row is encoded against the reference the same way the
        # reads are, so the designed change shows as the one column where the
        # parent disagrees with the variant built from it.
        parent_seq = g.get("parent") or ""
        parent_str = "".join(
            "." if parent_seq[i].upper() == ref_seq[i].upper() else parent_seq[i]
            for i in range(min(len(parent_seq), ref_len))
        ) if parent_seq else ""

        ref_seq_js = _json.dumps(ref_seq)
        rows_js = _json.dumps(rows_encoded)
        cons_js = _json.dumps(consensus_str)
        parent_js = _json.dumps(parent_str)
        flagged_js = _json.dumps(flagged_cols)
        n_rows = len(rows_encoded)

        # Translation data for canvas rendering
        ref_protein_js = _json.dumps(ref_protein) if ref_protein else "null"
        cons_protein_js = _json.dumps(cons_protein) if cons_protein else "null"

        # The annotation track shares the pileup's scroll container, so it
        # moves with the columns it describes.  cell_w comes from annotate
        # rather than being recomputed in the JS, so a glyph cannot land off
        # the column it names.
        cell_w = _cell_width(ref_len)
        annot_prefix = f"svf{idx}_"
        annot_plan = _plan_track(view.features, ref_len, cell_w=cell_w)
        annot_svg = _track_svg(annot_plan, prefix=annot_prefix)
        annot_plans.append((annot_plan, annot_prefix))
        annot_h = annot_plan.height + (5 if annot_plan.glyphs else 0)

        # The flank regions get a tint at the very back of the stack, so the
        # insert reads as a lit window between two shoulders before any text is
        # read.  It sits behind the glyphs rather than over them: it is context,
        # and tinting a feature's own colour would misreport it.
        region_tint = "".join(
            f'<div class="sv-region" style="left:{start * cell_w}px;'
            f'width:{(end - start) * cell_w}px"></div>'
            for start, end, _label, kind in _region_spans(flank_lengths, ref_len)
            if kind == "vec"
        )
        band_svg = _region_band_svg(flank_lengths, ref_len, cell_w=cell_w)
        band_h = 15 + 4 if band_svg else 0

        # Drawn only when asked for: another row of vertical space, spent only
        # when the page is worth it.  Sits right under the features track, so a
        # spike in disagreement can be read against what it falls in above it.
        mismatch_svg = (
            _mismatch_track_svg(mismatch_freqs, cell_w=cell_w)
            if view.mismatch_freq else ""
        )
        mismatch_h = _MISMATCH_TRACK_HEIGHT + 5 if mismatch_svg else 0

        if n_rows == 0:
            pileup_block = (
                f'<div class="pileup-empty">'
                f'No aligned reads available ({g["n_reads"]} reads unaligned)'
                f'</div>'
            )
        else:
            pileup_block = (
                f'<div class="pileup-container">'
                f'<div class="pileup-outer">'
                f'<div class="pileup-labels" id="labels-{idx}"></div>'
                f'<div class="pileup-scroll-wrap" id="wrap-{idx}">'
                f'<div class="pileup-scroll" id="scroll-{idx}">'
                f'{region_tint}'
                f'{band_svg}'
                f'{annot_svg}'
                f'{mismatch_svg}'
                f'<canvas id="ruler-{idx}" class="pileup-ruler"></canvas>'
                f'<canvas id="pileup-{idx}"></canvas>'
                f'</div>'
                f'</div>'
                f'</div>'
                f'</div>'
                f'<script>'
                f'(function(){{'
                f'var ref={ref_seq_js};'
                f'var cons={cons_js};'
                f'var rows={rows_js};'
                f'var flanks={flanks_js};'
                f'var refAA={ref_protein_js};'
                f'var consAA={cons_protein_js};'
                f'var flaggedCols={flagged_js};'
                f'var parent={parent_js};'
                f'drawPileup("pileup-{idx}","ruler-{idx}","labels-{idx}",ref,cons,rows,flanks,"scroll-{idx}","wrap-{idx}",refAA,consAA,flaggedCols,{cell_w},{annot_h + band_h},{mismatch_h},parent);'
                f'}})();'
                f'</script>'
            )

        sections_html.append(f'{header}\n{pileup_block}')

    body = '\n<hr class="group-sep">\n'.join(sections_html)

    # --- Masthead ---------------------------------------------------------
    # With one group the page and the group are the same object, so they share
    # one panel and the reference is named once.  With several, the panel
    # carries the run and each group keeps its own header below.
    if single:
        head_name = metas[0]["name_html"]
        head_chip = metas[0]["chip"]
        head_facts = metas[0]["facts"]
        # The title is usually just the group's name restated; show it only
        # when the caller put something else there.
        plain = view.title.strip().lower()
        name = groups[0]["ref_id"].strip().lower()
        eyebrow = "" if plain in (name, f"pileup: {name}") else (
            f'<div class="sv-eyebrow">{_html.escape(view.title)}</div>'
        )
    else:
        head_name = _html.escape(view.title)
        head_chip = ""
        eyebrow = ""
        tally = {}
        for g in groups:
            verdict = _verdict(g["status"], g["is_recoverable"])
            if verdict:
                tally[verdict] = tally.get(verdict, 0) + 1
        rollup = " ".join(
            f'<b>{count}</b>&#8202;<span class="sv-glyph sv-glyph-{tier}">{glyph}</span>'
            for (tier, glyph), count in tally.items()
        )
        head_facts = " &middot; ".join(
            bit for bit in (
                metas[0]["geometry"] if hoist_geometry else "",
                _fact(
                    f"<b>{view.total_reads}</b> reads",
                    "Reads across every group on this page.",
                ) if view.total_reads else "",
                _fact(
                    f"<b>{len(groups)}</b> groups",
                    "Subpopulations the reads were split into; each gets its "
                    "own header and pileup below.",
                ),
                rollup,
            ) if bit
        )

    # The star already marks the highlighted groups inline, so the masthead
    # defines what it means rather than listing the same names again.  When no
    # group carries the flag, the list is the only thing carrying the fact.
    highlight_line = ""
    if view.highlight_ids or any_starred:
        if any_starred:
            body_text = (
                f'<span class="sv-star">&#9733;</span> = '
                f'{_html.escape(view.highlight_label)}'
            )
        else:
            body_text = (
                f'{_html.escape(view.highlight_label)}: '
                f'{_html.escape(", ".join(view.highlight_ids))}'
            )
        highlight_line = f'<div class="sv-highlight">{body_text}</div>'

    # --- Key --------------------------------------------------------------
    # Grouped by what the mark means, and it decodes the things a reader
    # actually stalls on: region boundaries and disagreement flags are in, the
    # reference bar is out because the row label already names it.
    def _swatch(var: str, label: str, extra: str = "") -> str:
        return (
            f'<span class="sv-kitem">'
            f'<i class="sv-sw" style="background:var(--{_p}-{var});{extra}"></i>'
            f'{label}</span>'
        )

    clusters = [("bases", "".join(
        _swatch(base.lower(), base) for base in ("A", "T", "C", "G")
    ))]

    cells = [_swatch("match", "match")]
    if has_flanks:
        cells.append(_swatch("vector", "vector"))
    cells.append(_swatch("gap", "gap", f"border-color:var(--{_p}-tick)"))
    clusters.append(("cells", "".join(cells)))

    marks = []
    if any_flags:
        marks.append(
            '<span class="sv-kitem"><i class="sv-sw sv-sw-tri"></i>'
            '&gt;10% disagree</span>'
        )
    if has_flanks:
        marks.append(
            '<span class="sv-kitem"><i class="sv-sw sv-sw-dash"></i>'
            'boundary</span>'
        )
    if marks:
        clusters.append(("marks", "".join(marks)))

    key_html = "".join(
        f'<div class="sv-kgroup"><span class="sv-klabel">{label}</span>'
        f'<div class="sv-krow">{items}</div></div>'
        for label, items in clusters
    )

    # --- Palette emitted once, read by both the swatches and the canvas ----
    # Feature fills come from the file the reference was read from, so they are
    # literals rather than theme tokens and cannot be swapped by redefining a
    # custom property.  That makes them per-page, so they stay in the shell
    # rather than moving into the static stylesheet.
    annot_css = "\n".join(
        _track_style(plan, prefix=prefix) for plan, prefix in annot_plans
        if plan.glyphs
    )

    palette_css = "\n".join(
        [":root {"]
        + [f"    --{_p}-{k}: {v};" for k, v in _PALETTE["light"].items()]
        + ["}", '[data-theme="dark"] {']
        + [f"    --{_p}-{k}: {v};" for k, v in _PALETTE["dark"].items()]
        + ["}"]
    )
    palette_js = _json.dumps(_PALETTE)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{_html.escape(title)}</title>
<style id="{_style_id}">
{palette_css}
{pileup_css}
{annot_css}
</style>
<script>
/* The one palette, shared with the CSS custom properties above.  A swatch and
   the cell it describes cannot drift apart because both read this table. */
var SV_PALETTE = {palette_js};
{pileup_js}
</script>
</head>
<body>
<div class="sv-panel">
    <div class="sv-panel-id">
        {eyebrow}<div class="sv-idline"><span class="sv-name">{head_name}</span>{head_chip}</div>
        <div class="sv-facts">{head_facts}</div>
        {highlight_line}
    </div>
    <div class="sv-key">{key_html}</div>
</div>
{body}
<script id="{_script_id}">
{theme_js}
</script>
</body>
</html>"""
