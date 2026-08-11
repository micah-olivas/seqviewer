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

from .codon import translate as _translate
from .model import PileupView

__all__ = ["render"]


#: Every colour the page draws, per theme.  The legend swatches read these as
#: CSS custom properties and the canvas reads the same values from an emitted
#: JS object, so a swatch can no longer disagree with the cell it describes.
_PALETTE = {
    "light": {
        "match": "#c8ccd0",
        "vector": "#dfe2e6",
        "gap": "#ffffff",
        "ref": "#1e293b",
        "a": "#e03131",
        "t": "#1971c2",
        "c": "#e8590c",
        "g": "#e67700",
        "boundary": "#d97706",
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

    sections_html = []
    metas = []
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
        geometry = [f"<b>{ref_len}</b> bp"]
        if has_flanks:
            geometry.append(
                f"insert <b>{flank_lengths[0] + 1}</b>&ndash;"
                f"<b>{ref_len - flank_lengths[1]}</b>"
            )

        counts = []
        if g["n_reads"] and view.total_reads:
            counts.append(f"<b>{g['n_reads']}</b> of <b>{view.total_reads}</b> reads")
        elif g["n_reads"]:
            counts.append(f"<b>{g['n_reads']}</b> reads")
        if drawn != g["n_reads"]:
            counts.append(f"<b>{drawn}</b> drawn")
        if identity is not None:
            counts.append(f"<b>{identity:.1%}</b> identity")

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
        ref_seq = g["ref_seq"]
        _MISMATCH_THRESHOLD = 0.10
        for col_idx in range(ref_len):
            counts = Counter()
            for row in g["pileup_rows"]:
                base, _ = row[col_idx]
                if base != "-":
                    counts[base.upper()] += 1
            total = sum(counts.values())
            if counts:
                cons_base = counts.most_common(1)[0][0]
                if cons_base == ref_seq[col_idx].upper():
                    consensus_encoded.append(".")
                else:
                    consensus_encoded.append(cons_base)
                # Flag if >10% of reads differ from reference
                ref_base = ref_seq[col_idx].upper()
                ref_count = counts.get(ref_base, 0)
                if total and (total - ref_count) / total > _MISMATCH_THRESHOLD:
                    flagged_cols.append(col_idx)
            else:
                consensus_encoded.append("-")
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

        ref_seq_js = _json.dumps(ref_seq)
        rows_js = _json.dumps(rows_encoded)
        cons_js = _json.dumps(consensus_str)
        flagged_js = _json.dumps(flagged_cols)
        n_rows = len(rows_encoded)

        # Translation data for canvas rendering
        ref_protein_js = _json.dumps(ref_protein) if ref_protein else "null"
        cons_protein_js = _json.dumps(cons_protein) if cons_protein else "null"

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
                f'drawPileup("pileup-{idx}","ruler-{idx}","labels-{idx}",ref,cons,rows,flanks,"scroll-{idx}","wrap-{idx}",refAA,consAA,flaggedCols);'
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
                f"<b>{view.total_reads}</b> reads" if view.total_reads else "",
                f"<b>{len(groups)}</b> groups",
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
:root {{
    --{_p}-bg: #fafafa;
    --text: #1e293b;
    --muted: #94a3b8;
    --card-bg: #ffffff;
    --border: #e5e7eb;
    --panel-line: #dfe3e8;
    --ok: #0f7a52;
    --ok-bg: #e8f6ef;
    --warn: #a1650a;
    --warn-bg: #fdf3e2;
    --bad: #c22f2f;
    --bad-bg: #fdecec;
    --neutral: #5b6672;
    --neutral-bg: #eef1f4;
    /* One gutter width, used by the plot's row labels and by the indent on
       every text block, so text and plot share a single left edge. */
    --gutter: 2.5rem;
    --gutter-gap: 4px;
    --mono: 'SF Mono', SFMono-Regular, Menlo, Consolas, monospace;
}}
[data-theme="dark"] {{
    --{_p}-bg: #1a1a2e;
    --text: #e0e0e0;
    --muted: #64748b;
    --card-bg: #16213e;
    --border: #334155;
    --panel-line: #2c3a55;
    --ok: #35c493;
    --ok-bg: #14312a;
    --warn: #e0a355;
    --warn-bg: #33280f;
    --bad: #f0655f;
    --bad-bg: #35191c;
    --neutral: #93a1b0;
    --neutral-bg: #202b45;
}}
html, body {{
    background: var(--{_p}-bg);
    color: var(--text);
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    margin: 0;
    padding: 1.5rem;
}}
/* Every text block is indented by the plot's label gutter, so prose and
   pileup share one left edge instead of the two the old negative margin left. */
.sv-panel, .sv-group-head, .group-sep, .pileup-empty {{
    margin-left: calc(var(--gutter) + var(--gutter-gap));
}}
/* --- Masthead --- */
.sv-panel {{
    display: flex;
    flex-wrap: wrap;
    gap: 1rem 1.75rem;
    justify-content: space-between;
    align-items: flex-start;
    background: var(--card-bg);
    border: 1px solid var(--panel-line);
    border-radius: 3px;
    padding: 0.7rem 0.9rem;
    margin-bottom: 1.1rem;
}}
.sv-panel-id {{
    display: flex;
    flex-direction: column;
    gap: 0.3rem;
    min-width: 0;
}}
.sv-eyebrow {{
    font: 600 0.66rem/1 var(--mono);
    letter-spacing: 0.1em;
    text-transform: uppercase;
    color: var(--muted);
}}
.sv-idline {{
    display: flex;
    align-items: baseline;
    gap: 0.55rem;
    flex-wrap: wrap;
}}
.sv-name {{
    font-size: 1.15rem;
    font-weight: 700;
    letter-spacing: -0.01em;
}}
.sv-star {{
    color: var(--warn);
    cursor: help;
}}
.sv-facts {{
    font: 0.78rem/1.55 var(--mono);
    color: var(--muted);
    font-variant-numeric: tabular-nums;
}}
.sv-facts b {{
    color: var(--text);
    font-weight: 600;
}}
.sv-highlight {{
    font-size: 0.75rem;
    color: var(--muted);
}}
/* --- Verdict chip: hue, word, and glyph, so colour is never the only cue --- */
.sv-chip {{
    display: inline-flex;
    align-items: center;
    gap: 0.3rem;
    font-size: 0.72rem;
    font-weight: 600;
    padding: 0.28rem 0.45rem;
    border-radius: 2px;
    border: 1px solid currentColor;
    white-space: nowrap;
}}
.sv-glyph {{
    font-family: var(--mono);
    font-size: 0.85em;
}}
.sv-chip-ok {{ color: var(--ok); background: var(--ok-bg); }}
.sv-chip-warn {{ color: var(--warn); background: var(--warn-bg); }}
.sv-chip-bad {{ color: var(--bad); background: var(--bad-bg); }}
.sv-chip-neutral {{ color: var(--neutral); background: var(--neutral-bg); }}
.sv-glyph-ok {{ color: var(--ok); }}
.sv-glyph-warn {{ color: var(--warn); }}
.sv-glyph-bad {{ color: var(--bad); }}
.sv-glyph-neutral {{ color: var(--neutral); }}
/* --- Key: grouped by what the mark means --- */
.sv-key {{
    display: flex;
    gap: 1.1rem;
    flex-wrap: wrap;
    align-items: flex-start;
}}
.sv-kgroup {{
    display: flex;
    flex-direction: column;
    gap: 0.3rem;
}}
.sv-klabel {{
    font: 600 0.6rem/1 var(--mono);
    letter-spacing: 0.09em;
    text-transform: uppercase;
    color: var(--muted);
}}
.sv-krow {{
    display: flex;
    gap: 0.5rem;
    align-items: center;
    flex-wrap: wrap;
}}
.sv-kitem {{
    display: inline-flex;
    align-items: center;
    gap: 0.25rem;
    font-size: 0.7rem;
    color: var(--text);
    white-space: nowrap;
}}
.sv-sw {{
    width: 10px;
    height: 10px;
    border-radius: 1px;
    border: 1px solid transparent;
    display: inline-block;
    flex: none;
}}
.sv-sw-tri {{
    width: 0;
    height: 0;
    border-radius: 0;
    border: none;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 8px solid var(--{_p}-flag);
}}
.sv-sw-dash {{
    width: 0;
    height: 11px;
    border-radius: 0;
    border: none;
    border-left: 2px dashed var(--{_p}-boundary);
}}
/* --- Per-group header, shown only when there is more than one group --- */
.sv-group-head {{
    display: flex;
    align-items: baseline;
    gap: 0.7rem;
    flex-wrap: wrap;
    /* Longhand: a `margin` shorthand here would reset the gutter indent. */
    margin-top: 1.25rem;
    margin-bottom: 0.45rem;
}}
.sv-group-name {{
    font-weight: 700;
    font-size: 1rem;
}}
.pileup-container {{
    margin-bottom: 0.5rem;
}}
.pileup-outer {{
    display: flex;
    align-items: stretch;
}}
.pileup-scroll-wrap {{
    position: relative;
    flex: 1;
    min-width: 0;
}}
.pileup-scroll {{
    overflow-x: auto;
    scrollbar-width: none;
    background: transparent;
    border: none;
    border-radius: 0;
    padding: 0;
}}
.pileup-scroll::-webkit-scrollbar {{
    display: none;
}}
.pileup-mm-arrow {{
    position: absolute;
    top: 0;
    display: none;
    align-items: center;
    justify-content: center;
    width: 28px;
    font-size: 16px;
    pointer-events: none;
    z-index: 2;
    color: var(--text);
}}
.pileup-mm-arrow-l {{ left: 0; background: linear-gradient(to right, var(--{_p}-bg) 40%, transparent); }}
.pileup-mm-arrow-r {{ right: 0; background: linear-gradient(to left, var(--{_p}-bg) 40%, transparent); }}
.pileup-scroll canvas {{
    display: block;
}}
.pileup-labels {{
    display: flex;
    flex-direction: column;
    justify-content: flex-start;
    flex-shrink: 0;
    width: var(--gutter);
    padding-right: var(--gutter-gap);
    font: 9px/1 SF Mono, Menlo, Consolas, monospace;
    color: var(--muted);
    text-align: right;
    white-space: nowrap;
}}
.pileup-labels span {{
    display: flex;
    align-items: center;
    justify-content: flex-end;
}}
.pileup-empty {{
    font-size: 0.85rem;
    color: var(--muted);
    font-style: italic;
    padding: 1rem;
    background: var(--card-bg);
    border: 1px solid var(--border);
    border-radius: 6px;
}}
.group-sep {{
    border: none;
    border-top: 1px solid var(--border);
    margin-top: 1.5rem;
    margin-bottom: 1.5rem;
}}
</style>
<script>
/* The one palette, shared with the CSS custom properties above.  A swatch and
   the cell it describes cannot drift apart because both read this table. */
var SV_PALETTE = {palette_js};
function drawPileup(canvasId, rulerId, labelsId, refSeq, cons, rows, flanks, scrollId, wrapId, refAA, consAA, flaggedCols) {{
  var canvas = document.getElementById(canvasId);
  var rulerCanvas = document.getElementById(rulerId);
  var labelsEl = document.getElementById(labelsId);
  if (!canvas) return;
  var nCols = refSeq.length;
  var nRows = rows.length;
  var cellW = nCols < 200 ? 4 : nCols < 500 ? 3 : 2;
  var cellH = nRows < 100 ? 3 : nRows < 400 ? 2 : 1;
  var refH = Math.max(cellH, 6);
  var consH = refH;
  var gap = 4;
  var totalW = nCols * cellW;
  var dpr = window.devicePixelRatio || 1;

  // Translation rows below reads (aligned to insert region)
  var hasAA = refAA && consAA && flanks;
  var aaH = hasAA ? 14 : 0;       // height of each AA row
  var aaGap = hasAA ? 6 : 0;      // gap before AA section
  var aaCodonW = Math.max(3 * cellW, 8);  // min 8px so AA letters are legible

  // Canvas must be wide enough for both the nucleotide pileup and the AA section
  var aaEndPx = hasAA ? flanks[0] * cellW + refAA.length * aaCodonW : 0;
  var canvasW = Math.max(totalW, aaEndPx);

  var pileupH = refH + gap + consH + gap + nRows * cellH + aaGap + (hasAA ? aaH * 2 + 2 : 0);
  canvas.width = canvasW * dpr;
  canvas.height = pileupH * dpr;
  canvas.style.width = canvasW + 'px';
  canvas.style.height = pileupH + 'px';
  var ctx = canvas.getContext('2d');
  ctx.scale(dpr, dpr);
  var isDark = document.documentElement.getAttribute('data-theme') === 'dark';
  var P = SV_PALETTE[isDark ? 'dark' : 'light'];
  var matchColor = P.match;
  var vectorMatchColor = P.vector;
  var gapColor = P.gap;
  var refColor = P.ref;
  var consMatchColor = matchColor;
  var baseColors = {{'A': P.a, 'T': P.t, 'C': P.c, 'G': P.g}};
  // Flagged columns: positions where >10% of reads disagree with reference.
  // Falls back to consensus-derived mismatches when flaggedCols not provided.
  var mismatchCols = flaggedCols || [];
  if (!flaggedCols) {{
    for (var _mi = 0; _mi < cons.length; _mi++) {{
      var _ch = cons[_mi];
      if (_ch !== '.' && _ch !== '-') mismatchCols.push(_mi);
    }}
  }}
  var triRowH = mismatchCols.length > 0 ? 13 : 0;
  function isVector(col) {{
    return flanks && (col < flanks[0] || col >= nCols - flanks[1]);
  }}
  function pickMatch(col) {{
    return isVector(col) ? vectorMatchColor : matchColor;
  }}
  // --- Ruler ---
  var rulerH = (flanks ? 24 : 14) + triRowH;
  if (rulerCanvas) {{
    rulerCanvas.width = canvasW * dpr;
    rulerCanvas.height = rulerH * dpr;
    rulerCanvas.style.width = canvasW + 'px';
    rulerCanvas.style.height = rulerH + 'px';
    var rc = rulerCanvas.getContext('2d');
    rc.scale(dpr, dpr);
    var tickColor = P.tick;
    var labelColor = P['tick-label'];
    var boundaryColor = P.boundary;
    rc.clearRect(0, 0, canvasW, rulerH);
    var tickBottom = rulerH - triRowH;
    // Region labels on top row (if flanks present)
    var tickRowY = 0;
    if (flanks) {{
      tickRowY = 11;
      rc.fillStyle = boundaryColor;
      rc.font = '9px SF Mono,Menlo,Consolas,monospace';
      rc.textAlign = 'center';
      rc.textBaseline = 'top';
      var bLeft = flanks[0] * cellW;
      var bRight = (nCols - flanks[1]) * cellW;
      var minLabelPx = 40;
      if (flanks[0] > 0 && flanks[0] * cellW > minLabelPx) {{
        rc.fillText("5\u2032 vector", bLeft / 2, 0);
      }}
      var insertW = bRight - bLeft;
      if (insertW > minLabelPx) {{
        rc.fillText("insert", bLeft + insertW / 2, 0);
      }}
      if (flanks[1] > 0 && flanks[1] * cellW > minLabelPx) {{
        rc.fillText("3\u2032 vector", bRight + (totalW - bRight) / 2, 0);
      }}
      // Boundary dashed lines (start below both text rows, stop above triangle row)
      rc.setLineDash([3, 2]);
      rc.strokeStyle = boundaryColor;
      rc.lineWidth = 1;
      var dashY = tickRowY + 12;
      if (flanks[0] > 0 && dashY < tickBottom) {{
        rc.beginPath(); rc.moveTo(bLeft, dashY); rc.lineTo(bLeft, tickBottom); rc.stroke();
      }}
      if (flanks[1] > 0 && dashY < tickBottom) {{
        rc.beginPath(); rc.moveTo(bRight, dashY); rc.lineTo(bRight, tickBottom); rc.stroke();
      }}
      rc.setLineDash([]);
    }}
    // Tick labels + ticks (stop above triangle row)
    rc.strokeStyle = tickColor;
    rc.fillStyle = labelColor;
    rc.font = '10px SF Mono,Menlo,Consolas,monospace';
    rc.textBaseline = 'top';
    for (var i = 0; i < nCols; i++) {{
      var x = i * cellW + cellW / 2;
      if ((i + 1) % 100 === 0) {{
        rc.strokeStyle = tickColor;
        rc.beginPath(); rc.moveTo(x, tickRowY + 10); rc.lineTo(x, tickBottom); rc.stroke();
        rc.fillStyle = labelColor;
        rc.textAlign = 'center';
        rc.fillText(String(i + 1), x, tickRowY);
      }} else if ((i + 1) % 50 === 0) {{
        rc.strokeStyle = tickColor;
        rc.beginPath(); rc.moveTo(x, tickBottom - 3); rc.lineTo(x, tickBottom); rc.stroke();
      }}
    }}
    // --- Mismatch triangles (pointing down toward ref) ---
    if (mismatchCols.length > 0) {{
      var triH = 10, triW = Math.max(cellW * 2, 9);
      for (var _ti = 0; _ti < mismatchCols.length; _ti++) {{
        var mc = mismatchCols[_ti];
        rc.fillStyle = baseColors[cons[mc]] || P.flag;
        var cx = mc * cellW + cellW / 2;
        var ty = tickBottom + 1;
        rc.beginPath();
        rc.moveTo(cx - triW / 2, ty);
        rc.lineTo(cx + triW / 2, ty);
        rc.lineTo(cx, ty + triH);
        rc.closePath();
        rc.fill();
      }}
    }}
  }}
  // --- HTML row labels ---
  var consY = refH + gap;
  var readsY = consY + consH + gap;
  if (labelsEl) {{
    labelsEl.innerHTML = '';
    var rulerSpacer = document.createElement('span');
    rulerSpacer.style.height = rulerH + 'px';
    labelsEl.appendChild(rulerSpacer);
    var refLabel = document.createElement('span');
    refLabel.textContent = 'Ref';
    refLabel.style.height = refH + 'px';
    labelsEl.appendChild(refLabel);
    var gapSpacer1 = document.createElement('span');
    gapSpacer1.style.height = gap + 'px';
    labelsEl.appendChild(gapSpacer1);
    var consLabel = document.createElement('span');
    consLabel.textContent = 'Cons';
    consLabel.style.height = consH + 'px';
    labelsEl.appendChild(consLabel);
    var gapSpacer2 = document.createElement('span');
    gapSpacer2.style.height = gap + 'px';
    labelsEl.appendChild(gapSpacer2);
    if (nRows > 0) {{
      var readsLabel = document.createElement('span');
      readsLabel.textContent = 'Reads';
      readsLabel.style.height = (nRows * cellH) + 'px';
      labelsEl.appendChild(readsLabel);
    }}
    if (hasAA) {{
      var aaGapSpacer = document.createElement('span');
      aaGapSpacer.style.height = aaGap + 'px';
      labelsEl.appendChild(aaGapSpacer);
      var refAALabel = document.createElement('span');
      refAALabel.textContent = 'Ref AA';
      refAALabel.style.height = aaH + 'px';
      labelsEl.appendChild(refAALabel);
      var aaRowGap = document.createElement('span');
      aaRowGap.style.height = '2px';
      labelsEl.appendChild(aaRowGap);
      var consAALabel = document.createElement('span');
      consAALabel.textContent = 'Cons AA';
      consAALabel.style.height = aaH + 'px';
      labelsEl.appendChild(consAALabel);
    }}
  }}
  // --- Reference row ---
  ctx.fillStyle = refColor;
  for (var i = 0; i < nCols; i++) {{
    ctx.fillRect(i * cellW, 0, cellW, refH);
  }}
  // --- Consensus row ---
  for (var i = 0; i < cons.length; i++) {{
    var ch = cons[i];
    if (ch === '.') {{
      ctx.fillStyle = isVector(i) ? vectorMatchColor : consMatchColor;
    }} else if (ch === '-') {{
      ctx.fillStyle = gapColor;
    }} else {{
      ctx.fillStyle = baseColors[ch] || P.flag;
    }}
    ctx.fillRect(i * cellW, consY, cellW, consH);
  }}
  // --- Read rows ---
  for (var r = 0; r < nRows; r++) {{
    var row = rows[r];
    var y = readsY + r * cellH;
    for (var c = 0; c < row.length; c++) {{
      var ch = row[c];
      if (ch === '.') {{
        ctx.fillStyle = pickMatch(c);
      }} else if (ch === '-') {{
        ctx.fillStyle = gapColor;
      }} else {{
        ctx.fillStyle = baseColors[ch] || P.flag;
      }}
      ctx.fillRect(c * cellW, y, cellW, cellH);
    }}
  }}
  // --- Translation rows (aligned to insert region) ---
  var aaY = readsY + nRows * cellH + aaGap;
  if (hasAA) {{
    var insStart = flanks[0];
    var aaMatchColor = P['aa-match'];
    var aaDiffColor = P['aa-diff'];
    var aaDiffBg = P['aa-diff-bg'];
    var aaBg = P['aa-bg'];
    var aaFont = Math.min(aaH - 2, Math.max(7, aaCodonW - 2));
    ctx.font = aaFont + 'px SF Mono,Menlo,Consolas,monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

    // Draw background for the insert region
    var insX = insStart * cellW;
    var insW = refAA.length * aaCodonW;
    ctx.fillStyle = aaBg;
    ctx.fillRect(insX, aaY, insW, aaH * 2 + 2);

    for (var ai = 0; ai < refAA.length; ai++) {{
      var ax = insStart * cellW + ai * aaCodonW;
      var rAA = refAA[ai];
      var cAA = consAA[ai];
      var match = rAA === cAA;

      // Ref AA row
      if (!match) {{
        ctx.fillStyle = aaDiffBg;
        ctx.fillRect(ax, aaY, aaCodonW, aaH);
      }}
      ctx.fillStyle = match ? aaMatchColor : aaDiffColor;
      if (aaCodonW >= 7) {{
        ctx.fillText(rAA, ax + aaCodonW / 2, aaY + aaH / 2);
      }} else {{
        ctx.fillRect(ax + 1, aaY + 2, aaCodonW - 2, aaH - 4);
      }}

      // Cons AA row
      var caaY = aaY + aaH + 2;
      if (!match) {{
        ctx.fillStyle = aaDiffBg;
        ctx.fillRect(ax, caaY, aaCodonW, aaH);
      }}
      ctx.fillStyle = match ? aaMatchColor : aaDiffColor;
      if (aaCodonW >= 7) {{
        ctx.fillText(cAA, ax + aaCodonW / 2, caaY + aaH / 2);
      }} else {{
        ctx.fillRect(ax + 1, caaY + 2, aaCodonW - 2, aaH - 4);
      }}
    }}

    // Subtle codon grid lines
    ctx.strokeStyle = P['aa-grid'];
    ctx.lineWidth = 0.5;
    for (var ai = 1; ai < refAA.length; ai++) {{
      var lx = insStart * cellW + ai * aaCodonW;
      ctx.beginPath(); ctx.moveTo(lx, aaY); ctx.lineTo(lx, aaY + aaH * 2 + 2); ctx.stroke();
    }}
  }}
  // --- Region boundary dashed lines on pileup canvas ---
  if (flanks) {{
    ctx.save();
    ctx.setLineDash([4, 3]);
    ctx.strokeStyle = P.boundary;
    ctx.lineWidth = 1;
    var pH = pileupH;
    if (flanks[0] > 0) {{
      var bx = flanks[0] * cellW;
      ctx.beginPath(); ctx.moveTo(bx, 0); ctx.lineTo(bx, pH); ctx.stroke();
    }}
    if (flanks[1] > 0) {{
      var bx2 = (nCols - flanks[1]) * cellW;
      ctx.beginPath(); ctx.moveTo(bx2, 0); ctx.lineTo(bx2, pH); ctx.stroke();
    }}
    ctx.restore();
  }}
  // --- Mismatch overflow arrows ---
  if (mismatchCols.length > 0 && scrollId && wrapId) {{
    var scrollEl = document.getElementById(scrollId);
    var wrapEl = document.getElementById(wrapId);
    if (scrollEl && wrapEl) {{
      var leftArrow = document.createElement('div');
      leftArrow.className = 'pileup-mm-arrow pileup-mm-arrow-l';
      leftArrow.textContent = '\u25c4';
      leftArrow.style.height = rulerH + 'px';
      wrapEl.appendChild(leftArrow);
      var rightArrow = document.createElement('div');
      rightArrow.className = 'pileup-mm-arrow pileup-mm-arrow-r';
      rightArrow.textContent = '\u25ba';
      rightArrow.style.height = rulerH + 'px';
      wrapEl.appendChild(rightArrow);
      function updateMmArrows() {{
        var sl = scrollEl.scrollLeft;
        var vw = scrollEl.clientWidth;
        var hasL = false, hasR = false;
        for (var _ai = 0; _ai < mismatchCols.length; _ai++) {{
          var ax = mismatchCols[_ai] * cellW + cellW / 2;
          if (ax < sl + 4) hasL = true;
          if (ax > sl + vw - 4) hasR = true;
        }}
        leftArrow.style.display = hasL ? 'flex' : 'none';
        rightArrow.style.display = hasR ? 'flex' : 'none';
      }}
      scrollEl.addEventListener('scroll', updateMmArrows);
      updateMmArrows();
    }}
  }}
  // --- Tooltip ---
  function regionLabel(col) {{
    if (!flanks) return '';
    if (col < flanks[0]) return '[5\u2032 vector] ';
    if (col >= nCols - flanks[1]) return '[3\u2032 vector] ';
    return '[insert] ';
  }}
  var tooltip = document.createElement('div');
  tooltip.style.cssText = 'position:fixed;background:#1e293b;color:#fff;padding:4px 8px;'
    + 'border-radius:4px;font-size:11px;pointer-events:none;display:none;z-index:10;'
    + 'font-family:SF Mono,Menlo,Consolas,monospace;';
  document.body.appendChild(tooltip);
  canvas.addEventListener('mousemove', function(e) {{
    var rect = canvas.getBoundingClientRect();
    var x = e.clientX - rect.left;
    var yp = e.clientY - rect.top;
    var col = Math.floor(x / cellW);
    if (col < 0 || col >= nCols) {{ tooltip.style.display = 'none'; return; }}
    var rl = regionLabel(col);
    if (yp < refH) {{
      tooltip.textContent = rl + 'Ref pos ' + (col + 1) + ': ' + refSeq[col];
    }} else if (yp < consY + consH) {{
      var ch = cons[col];
      var base = ch === '.' ? refSeq[col] : ch;
      var note = ch === '.' ? ' (match)' : ch === '-' ? '' : ' (mismatch)';
      tooltip.textContent = rl + 'Consensus pos ' + (col + 1) + ': ' + base + note;
    }} else if (hasAA && yp >= aaY && yp < aaY + aaH * 2 + 2) {{
      var insStart = flanks[0];
      var aaIdx = Math.floor((col - insStart) / 3);
      if (aaIdx >= 0 && aaIdx < refAA.length) {{
        var isRefRow = yp < aaY + aaH;
        var which = isRefRow ? 'Ref' : 'Cons';
        var aa = isRefRow ? refAA[aaIdx] : consAA[aaIdx];
        var other = isRefRow ? consAA[aaIdx] : refAA[aaIdx];
        var note = aa === other ? ' (match)' : ' \u2260 ' + (isRefRow ? 'Cons' : 'Ref') + ': ' + other;
        tooltip.textContent = which + ' AA ' + (aaIdx + 1) + ': ' + aa + note;
      }} else {{
        tooltip.style.display = 'none'; return;
      }}
    }} else if (yp >= readsY && yp < readsY + nRows * cellH) {{
      var row_idx = Math.floor((yp - readsY) / cellH);
      if (row_idx >= 0 && row_idx < nRows) {{
        var ch = rows[row_idx][col];
        var label = ch === '.' ? refSeq[col] + ' (match)' : ch === '-' ? 'gap' : ch + ' (mismatch)';
        tooltip.textContent = rl + 'Read ' + (row_idx + 1) + ', pos ' + (col + 1) + ': ' + label;
      }} else {{
        tooltip.style.display = 'none'; return;
      }}
    }} else {{
      tooltip.style.display = 'none'; return;
    }}
    tooltip.style.display = 'block';
    tooltip.style.left = (e.clientX + 12) + 'px';
    tooltip.style.top = (e.clientY - 24) + 'px';
  }});
  canvas.addEventListener('mouseleave', function() {{
    tooltip.style.display = 'none';
  }});
}}
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
(function () {{
  try {{
    var stored = localStorage.getItem('{_storage_key}');
    if (stored === 'dark') {{
      document.documentElement.setAttribute('data-theme', 'dark');
    }}
  }} catch (e) {{}}
}})();
</script>
</body>
</html>"""
