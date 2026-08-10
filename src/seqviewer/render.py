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

    sections_html = []
    for idx, g in enumerate(groups):
        star = " &#9733;" if g["is_recoverable"] else ""
        if g["status"] == "Silent Mutation":
            status_class = "status-silent"
        elif g["is_recoverable"]:
            status_class = "status-correct"
        else:
            status_class = "status-other"

        # Compute per-read identity from pileup data
        identity_str = ""
        if g["pileup_rows"]:
            total_bases = 0
            total_matches = 0
            for row in g["pileup_rows"]:
                aligned = [(b, m) for b, m in row if b != "-"]
                total_bases += len(aligned)
                total_matches += sum(1 for _, m in aligned if m)
            if total_bases > 0:
                identity = total_matches / total_bases
                identity_str = f" &middot; Read identity: {identity:.1%}"

        ref_len = len(g["ref_seq"])

        header = (
            f'<div class="group-header">'
            f'<span class="ref-name">{_html.escape(g["ref_id"])}{star}</span>'
            f'<span class="group-meta">'
            f'{g["n_reads"]} reads ({g["frac"]:.0%}) &middot; '
            f'{ref_len} bp &middot; '
            f'Consensus: <span class="{status_class}">'
            f'{_html.escape(g["status"])}</span>'
            f'{identity_str}'
            f'</span></div>'
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

        # Reconstruct actual consensus DNA and translate the insert region
        consensus_dna = "".join(
            ref_seq[i] if c == "." else (c if c != "-" else "N")
            for i, c in enumerate(consensus_encoded)
        )
        _ins_start = flank_lengths[0] if flank_lengths else 0
        _ins_end = ref_len - (flank_lengths[1] if flank_lengths else 0)
        insert_dna = consensus_dna[_ins_start:_ins_end]
        ref_insert_dna = ref_seq[_ins_start:_ins_end]

        cons_protein = _translate(insert_dna)
        ref_protein = _translate(ref_insert_dna)

        ref_seq_js = _json.dumps(ref_seq)
        rows_js = _json.dumps(rows_encoded)
        cons_js = _json.dumps(consensus_str)
        flagged_js = _json.dumps(flagged_cols)
        n_rows = len(rows_encoded)
        n_cols = ref_len

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
                f'<div class="pileup-info">{n_rows} aligned reads &times; '
                f'{n_cols} bp</div>'
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

    recoverable_list = ", ".join(view.highlight_ids)
    recoverable_line = (
        f' &middot; {_html.escape(view.highlight_label)}: {_html.escape(recoverable_list)}'
        if recoverable_list else ""
    )

    vector_legend = ""
    if flank_lengths and (flank_lengths[0] or flank_lengths[1]):
        vector_legend = (
            '    <span class="legend-item">'
            '<span class="legend-swatch" style="background:#dfe2e6;"></span>'
            ' Vector Match</span>\n'
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{_html.escape(title)}</title>
<style id="{_style_id}">
:root {{
    --{_p}-bg: #fafafa;
    --text: #1e293b;
    --muted: #94a3b8;
    --card-bg: #ffffff;
    --border: #e5e7eb;
}}
[data-theme="dark"] {{
    --{_p}-bg: #1a1a2e;
    --text: #e0e0e0;
    --muted: #64748b;
    --card-bg: #16213e;
    --border: #334155;
}}
html, body {{
    background: var(--{_p}-bg);
    color: var(--text);
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
    margin: 0;
    padding: 1.5rem;
}}
h1 {{
    font-size: 1.4rem;
    margin: 0 0 0.25rem;
}}
.well-meta {{
    color: var(--muted);
    font-size: 0.9rem;
    margin-bottom: 1.5rem;
}}
.group-header {{
    display: flex;
    align-items: baseline;
    gap: 1rem;
    margin: 1rem 0 0.5rem;
}}
.ref-name {{
    font-weight: 700;
    font-size: 1.05rem;
}}
.group-meta {{
    color: var(--muted);
    font-size: 0.85rem;
}}
.status-correct {{
    color: #059669;
    font-weight: 600;
}}
.status-silent {{
    color: #d97706;
    font-weight: 600;
}}
.status-other {{
    color: #ef4444;
}}
.protein-seq {{
    margin-top: 0.5rem;
    font-family: 'Courier New', Courier, monospace;
    font-size: 10pt;
    white-space: nowrap;
    overflow-x: auto;
    color: var(--text);
    opacity: 0.85;
}}
.protein-label {{
    color: var(--muted);
    font-weight: 600;
    margin-right: 0.25rem;
    user-select: none;
}}
.pileup-container {{
    margin-bottom: 0.5rem;
    margin-left: -2.5rem;
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
    width: 2.5rem;
    padding-right: 4px;
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
.pileup-ruler {{
}}
.pileup-info {{
    font-size: 0.75rem;
    color: var(--muted);
    margin-top: 0.25rem;
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
.legend {{
    display: flex;
    gap: 1rem;
    align-items: center;
    font-size: 0.8rem;
    color: var(--muted);
    margin-bottom: 1rem;
    flex-wrap: wrap;
}}
.legend-item {{
    display: flex;
    align-items: center;
    gap: 0.3rem;
}}
.legend-swatch {{
    width: 12px;
    height: 12px;
    border-radius: 2px;
    border: 1px solid var(--border);
}}
.group-sep {{
    border: none;
    border-top: 1px solid var(--border);
    margin: 1.5rem 0;
}}
</style>
<script>
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
  var matchColor = isDark ? '#4a5568' : '#c8ccd0';
  var vectorMatchColor = isDark ? '#3a4455' : '#dfe2e6';
  var gapColor = isDark ? '#ffffff' : '#ffffff';
  var refColor = isDark ? '#e0e0e0' : '#1e293b';
  var consMatchColor = matchColor;
  var baseColors = isDark
    ? {{'A':'#ff6b6b','T':'#339af0','C':'#ffa94d','G':'#ffd43b'}}
    : {{'A':'#e03131','T':'#1971c2','C':'#e8590c','G':'#e67700'}};
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
    var tickColor = isDark ? '#64748b' : '#94a3b8';
    var labelColor = isDark ? '#e0e0e0' : '#1e293b';
    var boundaryColor = isDark ? '#f59e0b' : '#d97706';
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
        rc.fillStyle = baseColors[cons[mc]] || '#94a3b8';
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
      ctx.fillStyle = baseColors[ch] || '#94a3b8';
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
        ctx.fillStyle = baseColors[ch] || '#94a3b8';
      }}
      ctx.fillRect(c * cellW, y, cellW, cellH);
    }}
  }}
  // --- Translation rows (aligned to insert region) ---
  var aaY = readsY + nRows * cellH + aaGap;
  if (hasAA) {{
    var insStart = flanks[0];
    var aaMatchColor = isDark ? '#4a5568' : '#d1d5db';
    var aaDiffColor = isDark ? '#ff6b6b' : '#e03131';
    var aaDiffBg = isDark ? 'rgba(255,107,107,0.18)' : 'rgba(224,49,49,0.1)';
    var aaBg = isDark ? '#1e293b' : '#f8fafc';
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
    ctx.strokeStyle = isDark ? '#334155' : '#e5e7eb';
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
    ctx.strokeStyle = isDark ? '#f59e0b' : '#d97706';
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
<h1>{_html.escape(title)}</h1>
<div class="well-meta">
    {view.total_reads} total reads &middot;
    Top fraction: {view.top_fraction:.0%}{recoverable_line}
</div>
<div class="legend">
    <span style="font-weight:600;">Legend:</span>
    <span class="legend-item"><span class="legend-swatch" style="background:#c8ccd0;"></span> Match</span>
{vector_legend}    <span class="legend-item"><span class="legend-swatch" style="background:#e03131;"></span> A</span>
    <span class="legend-item"><span class="legend-swatch" style="background:#1971c2;"></span> T</span>
    <span class="legend-item"><span class="legend-swatch" style="background:#e8590c;"></span> C</span>
    <span class="legend-item"><span class="legend-swatch" style="background:#e67700;"></span> G</span>
    <span class="legend-item"><span class="legend-swatch" style="background:#ffffff;border:1px solid #d1d5db;"></span> Gap</span>
    <span class="legend-item"><span class="legend-swatch" style="background:#1e293b;"></span> Reference</span>

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
