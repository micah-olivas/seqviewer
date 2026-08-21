function drawPileup(canvasId, rulerId, labelsId, refSeq, cons, rows, flanks, scrollId, wrapId, refAA, consAA, flaggedCols, cellWIn, annotH, mismatchH, parent) {
  var canvas = document.getElementById(canvasId);
  var rulerCanvas = document.getElementById(rulerId);
  var labelsEl = document.getElementById(labelsId);
  if (!canvas) return;
  var nCols = refSeq.length;
  var nRows = rows.length;
  // Pixels per base is decided by annotate.cell_width and passed in, because
  // the SVG feature track measures itself with the same number. Recomputing the
  // rule here would let a glyph sit a pixel off the column it describes.
  var cellW = cellWIn;
  // A row is thinner the more of them there are, so a deep pileup still fits a
  // screen. The floor is 2px rather than 1: a one-pixel row is a hairline that
  // a fractional device pixel ratio can drop altogether.
  var cellH = nRows < 100 ? 4 : nRows < 400 ? 3 : 2;
  var refH = Math.max(cellH, 6);
  var consH = refH;
  // The parent sequence, when the caller supplied one: another row the
  // height of the reference, sitting directly under it so the designed
  // change reads as the column where the two disagree.
  var hasParent = !!(parent && parent.length);
  var parentH = hasParent ? refH : 0;
  var parentGap = hasParent ? 4 : 0;
  var gap = 4;
  var totalW = nCols * cellW;
  var dpr = window.devicePixelRatio || 1;

  // Translation rows below reads (aligned to insert region)
  var hasAA = refAA && consAA && flanks;
  var aaH = hasAA ? 14 : 0;       // height of each AA row
  var aaGap = hasAA ? 6 : 0;      // gap before AA section
  // A codon is exactly the three bases it translates. The old floor of 8px
  // made each glyph 2px wider than its codon whenever cellW was 2 -- every
  // reference >= 500 bp -- so the drift accumulated linearly and the last of
  // 266 residues on a 1000 bp insert sat 530px from the bases it came from.
  // A track that does not line up with the reference cannot be read against it.
  var aaCodonW = 3 * cellW;

  // Canvas must be wide enough for both the nucleotide pileup and the AA section
  var aaEndPx = hasAA ? flanks[0] * cellW + refAA.length * aaCodonW : 0;
  var canvasW = Math.max(totalW, aaEndPx);

  var pileupH = refH + parentGap + parentH + gap + consH + gap + nRows * cellH + aaGap + (hasAA ? aaH * 2 + 2 : 0);
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
  var baseColors = {'A': P.a, 'T': P.t, 'C': P.c, 'G': P.g};
  // Flagged columns: positions where >10% of reads disagree with reference.
  // Falls back to consensus-derived mismatches when flaggedCols not provided.
  var mismatchCols = flaggedCols || [];
  if (!flaggedCols) {
    for (var _mi = 0; _mi < cons.length; _mi++) {
      var _ch = cons[_mi];
      if (_ch !== '.' && _ch !== '-') mismatchCols.push(_mi);
    }
  }
  var triRowH = mismatchCols.length > 0 ? 13 : 0;
  function isVector(col) {
    return flanks && (col < flanks[0] || col >= nCols - flanks[1]);
  }
  function pickMatch(col) {
    return isVector(col) ? vectorMatchColor : matchColor;
  }
  // --- Ruler ---
  var rulerH = 14 + triRowH;
  // The flank tint stops above the reads. Run down the full height it sits
  // behind every row at a value close to the match colour, so a row has to be
  // read through it and the two greys compete. Over the tracks alone it still
  // says which span is which, which is all it was ever for; the reads get the
  // dashed boundaries instead, which cost no contrast.
  if (scrollId) {
    var scrollForRegions = document.getElementById(scrollId);
    if (scrollForRegions) {
      var regionEls = scrollForRegions.querySelectorAll('.sv-region');
      var headerH = (annotH || 0) + (mismatchH || 0) + rulerH;
      for (var ri = 0; ri < regionEls.length; ri++) {
        regionEls[ri].style.bottom = 'auto';
        regionEls[ri].style.height = headerH + 'px';
      }
    }
  }
  if (rulerCanvas) {
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
    // The regions are named by the band above the track now, not by text
    // floating here, so this only draws their boundaries.
    var tickRowY = 0;
    if (flanks) {
      rc.strokeStyle = boundaryColor;
      rc.lineWidth = 1;
      var bLeft = flanks[0] * cellW;
      var bRight = (nCols - flanks[1]) * cellW;
      if (flanks[0] > 0) {
        rc.beginPath(); rc.moveTo(bLeft, 0); rc.lineTo(bLeft, tickBottom); rc.stroke();
      }
      if (flanks[1] > 0) {
        rc.beginPath(); rc.moveTo(bRight, 0); rc.lineTo(bRight, tickBottom); rc.stroke();
      }
    }
    // Tick labels + ticks (stop above triangle row)
    rc.strokeStyle = tickColor;
    rc.fillStyle = labelColor;
    rc.font = '10px SF Mono,Menlo,Consolas,monospace';
    rc.textBaseline = 'top';
    for (var i = 0; i < nCols; i++) {
      var x = i * cellW + cellW / 2;
      if ((i + 1) % 100 === 0) {
        rc.strokeStyle = tickColor;
        rc.beginPath(); rc.moveTo(x, tickRowY + 10); rc.lineTo(x, tickBottom); rc.stroke();
        rc.fillStyle = labelColor;
        rc.textAlign = 'center';
        rc.fillText(String(i + 1), x, tickRowY);
      } else if ((i + 1) % 50 === 0) {
        rc.strokeStyle = tickColor;
        rc.beginPath(); rc.moveTo(x, tickBottom - 3); rc.lineTo(x, tickBottom); rc.stroke();
      }
    }
    // --- Mismatch triangles (pointing down toward ref) ---
    if (mismatchCols.length > 0) {
      var triH = 10, triW = Math.max(cellW * 2, 9);
      for (var _ti = 0; _ti < mismatchCols.length; _ti++) {
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
      }
    }
  }
  // --- HTML row labels ---
  var parentY = refH + parentGap;
  var consY = refH + parentGap + parentH + gap;
  var readsY = consY + consH + gap;
  if (labelsEl) {
    labelsEl.innerHTML = '';
    if (annotH) {
      var annotLabel = document.createElement('span');
      annotLabel.textContent = 'Features';
      annotLabel.style.height = annotH + 'px';
      labelsEl.appendChild(annotLabel);
    }
    if (mismatchH) {
      var mismatchLabel = document.createElement('span');
      mismatchLabel.textContent = 'Mismatches';
      mismatchLabel.title = 'Share of the reads called at this position that '
        + 'disagree with the reference. Gaps and uncalled reads count for '
        + 'neither side. Drawn on a 0-1% scale, so a position above 1% shows '
        + 'at full height.';
      mismatchLabel.style.height = mismatchH + 'px';
      labelsEl.appendChild(mismatchLabel);
    }
    var rulerSpacer = document.createElement('span');
    rulerSpacer.style.height = rulerH + 'px';
    labelsEl.appendChild(rulerSpacer);
    var refLabel = document.createElement('span');
    refLabel.textContent = 'Ref';
    refLabel.style.height = refH + 'px';
    labelsEl.appendChild(refLabel);
    if (hasParent) {
      var parentGapSpacer = document.createElement('span');
      parentGapSpacer.style.height = parentGap + 'px';
      labelsEl.appendChild(parentGapSpacer);
      var parentLabel = document.createElement('span');
      parentLabel.textContent = 'Parent';
      parentLabel.style.height = parentH + 'px';
      labelsEl.appendChild(parentLabel);
    }
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
    if (nRows > 0) {
      var readsLabel = document.createElement('span');
      readsLabel.textContent = 'Reads';
      readsLabel.style.height = (nRows * cellH) + 'px';
      labelsEl.appendChild(readsLabel);
    }
    if (hasAA) {
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
    }
  }
  // --- Reference row ---
  ctx.fillStyle = refColor;
  for (var i = 0; i < nCols; i++) {
    ctx.fillRect(i * cellW, 0, cellW, refH);
  }
  // --- Wild-type row ---
  if (hasParent) {
    for (var i = 0; i < parent.length; i++) {
      var wch = parent[i];
      if (wch === '.') {
        ctx.fillStyle = pickMatch(i);
      } else if (wch === '-') {
        ctx.fillStyle = gapColor;
      } else {
        ctx.fillStyle = baseColors[wch] || P.flag;
      }
      ctx.fillRect(i * cellW, parentY, cellW, parentH);
    }
  }
  // --- Consensus row ---
  for (var i = 0; i < cons.length; i++) {
    var ch = cons[i];
    if (ch === '.') {
      ctx.fillStyle = isVector(i) ? vectorMatchColor : consMatchColor;
    } else if (ch === '-') {
      ctx.fillStyle = gapColor;
    } else {
      ctx.fillStyle = baseColors[ch] || P.flag;
    }
    ctx.fillRect(i * cellW, consY, cellW, consH);
  }
  // --- Read rows ---
  for (var r = 0; r < nRows; r++) {
    var row = rows[r];
    var y = readsY + r * cellH;
    for (var c = 0; c < row.length; c++) {
      var ch = row[c];
      if (ch === '.') {
        ctx.fillStyle = pickMatch(c);
      } else if (ch === '-') {
        ctx.fillStyle = gapColor;
      } else {
        ctx.fillStyle = baseColors[ch] || P.flag;
      }
      ctx.fillRect(c * cellW, y, cellW, cellH);
    }
  }
  // --- Translation rows (aligned to insert region) ---
  var aaY = readsY + nRows * cellH + aaGap;
  if (hasAA) {
    var insStart = flanks[0];
    var aaMatchColor = P['aa-match'];
    var aaDiffColor = P['aa-diff'];
    var aaDiffBg = P['aa-diff-bg'];
    var aaBg = P['aa-bg'];
    var aaFont = Math.min(aaH - 2, aaCodonW - 1);
    var lettersFit = aaFont >= 7;
    ctx.font = aaFont + 'px SF Mono,Menlo,Consolas,monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

    // At 6px a codon cannot hold a legible letter, so weight carries the
    // meaning instead: a match is a quiet baseline tick, a mismatch is a full
    // block. This branch was unreachable while the 8px floor existed.
    function drawResidue(letter, x, y, isMatch) {
      if (lettersFit) {
        ctx.fillText(letter, x + aaCodonW / 2, y + aaH / 2);
      } else if (isMatch) {
        ctx.fillRect(x, y + aaH - 3, Math.max(1, aaCodonW - 1), 2);
      } else {
        ctx.fillRect(x, y + 1, Math.max(1, aaCodonW - 1), aaH - 2);
      }
    }

    // Draw background for the insert region
    var insX = insStart * cellW;
    var insW = refAA.length * aaCodonW;
    ctx.fillStyle = aaBg;
    ctx.fillRect(insX, aaY, insW, aaH * 2 + 2);

    for (var ai = 0; ai < refAA.length; ai++) {
      var ax = insStart * cellW + ai * aaCodonW;
      var rAA = refAA[ai];
      var cAA = consAA[ai];
      var match = rAA === cAA;

      // Ref AA row
      if (!match) {
        ctx.fillStyle = aaDiffBg;
        ctx.fillRect(ax, aaY, aaCodonW, aaH);
      }
      ctx.fillStyle = match ? aaMatchColor : aaDiffColor;
      drawResidue(rAA, ax, aaY, match);

      // Cons AA row
      var caaY = aaY + aaH + 2;
      if (!match) {
        ctx.fillStyle = aaDiffBg;
        ctx.fillRect(ax, caaY, aaCodonW, aaH);
      }
      ctx.fillStyle = match ? aaMatchColor : aaDiffColor;
      drawResidue(cAA, ax, caaY, match);
    }

    // Subtle codon grid lines, but only when the pitch is wide enough that a
    // line per codon reads as structure rather than as hatching.
    ctx.strokeStyle = P['aa-grid'];
    ctx.lineWidth = 0.5;
    for (var ai = 1; aaCodonW >= 9 && ai < refAA.length; ai++) {
      var lx = insStart * cellW + ai * aaCodonW;
      ctx.beginPath(); ctx.moveTo(lx, aaY); ctx.lineTo(lx, aaY + aaH * 2 + 2); ctx.stroke();
    }
  }
  // --- Region boundary dashed lines on pileup canvas ---
  if (flanks) {
    ctx.save();
    ctx.strokeStyle = P.boundary;
    ctx.lineWidth = 1;
    var pH = pileupH;
    if (flanks[0] > 0) {
      var bx = flanks[0] * cellW;
      ctx.beginPath(); ctx.moveTo(bx, 0); ctx.lineTo(bx, pH); ctx.stroke();
    }
    if (flanks[1] > 0) {
      var bx2 = (nCols - flanks[1]) * cellW;
      ctx.beginPath(); ctx.moveTo(bx2, 0); ctx.lineTo(bx2, pH); ctx.stroke();
    }
    ctx.restore();
  }
  // --- Mismatch overflow arrows ---
  if (mismatchCols.length > 0 && scrollId && wrapId) {
    var scrollEl = document.getElementById(scrollId);
    var wrapEl = document.getElementById(wrapId);
    if (scrollEl && wrapEl) {
      var arrowsTop = (annotH || 0) + (mismatchH || 0);
      var leftArrow = document.createElement('div');
      leftArrow.className = 'pileup-mm-arrow pileup-mm-arrow-l';
      leftArrow.textContent = '◄';
      leftArrow.style.top = arrowsTop + 'px';
      leftArrow.style.height = rulerH + 'px';
      wrapEl.appendChild(leftArrow);
      var rightArrow = document.createElement('div');
      rightArrow.className = 'pileup-mm-arrow pileup-mm-arrow-r';
      rightArrow.textContent = '►';
      rightArrow.style.top = arrowsTop + 'px';
      rightArrow.style.height = rulerH + 'px';
      wrapEl.appendChild(rightArrow);
      function updateMmArrows() {
        var sl = scrollEl.scrollLeft;
        var vw = scrollEl.clientWidth;
        var hasL = false, hasR = false;
        for (var _ai = 0; _ai < mismatchCols.length; _ai++) {
          var ax = mismatchCols[_ai] * cellW + cellW / 2;
          if (ax < sl + 4) hasL = true;
          if (ax > sl + vw - 4) hasR = true;
        }
        leftArrow.style.display = hasL ? 'flex' : 'none';
        rightArrow.style.display = hasR ? 'flex' : 'none';
      }
      scrollEl.addEventListener('scroll', updateMmArrows);
      updateMmArrows();
    }
  }
  // --- Tooltip ---
  function regionLabel(col) {
    if (!flanks) return '';
    if (col < flanks[0]) return '[5′ vector] ';
    if (col >= nCols - flanks[1]) return '[3′ vector] ';
    return '[insert] ';
  }
  var tooltip = document.createElement('div');
  tooltip.style.cssText = 'position:fixed;background:#1e293b;color:#fff;padding:4px 8px;'
    + 'border-radius:4px;font-size:11px;pointer-events:none;display:none;z-index:10;'
    + 'font-family:SF Mono,Menlo,Consolas,monospace;';
  document.body.appendChild(tooltip);
  canvas.addEventListener('mousemove', function(e) {
    var rect = canvas.getBoundingClientRect();
    var x = e.clientX - rect.left;
    var yp = e.clientY - rect.top;
    var col = Math.floor(x / cellW);
    if (col < 0 || col >= nCols) { tooltip.style.display = 'none'; return; }
    var rl = regionLabel(col);
    if (yp < refH) {
      tooltip.textContent = rl + 'Ref pos ' + (col + 1) + ': ' + refSeq[col];
    } else if (yp < consY + consH) {
      var ch = cons[col];
      var base = ch === '.' ? refSeq[col] : ch;
      var note = ch === '.' ? ' (match)' : ch === '-' ? '' : ' (mismatch)';
      tooltip.textContent = rl + 'Consensus pos ' + (col + 1) + ': ' + base + note;
    } else if (hasAA && yp >= aaY && yp < aaY + aaH * 2 + 2) {
      var insStart = flanks[0];
      // Derived from the drawn geometry, not from the nucleotide rate: reading
      // the residue index a different way than it was placed is what made the
      // tooltip name a third, differently wrong residue.
      var aaIdx = Math.floor((x - insStart * cellW) / aaCodonW);
      if (aaIdx >= 0 && aaIdx < refAA.length) {
        var isRefRow = yp < aaY + aaH;
        var which = isRefRow ? 'Ref' : 'Cons';
        var aa = isRefRow ? refAA[aaIdx] : consAA[aaIdx];
        var other = isRefRow ? consAA[aaIdx] : refAA[aaIdx];
        var note = aa === other ? ' (match)' : ' ≠ ' + (isRefRow ? 'Cons' : 'Ref') + ': ' + other;
        tooltip.textContent = which + ' AA ' + (aaIdx + 1) + ': ' + aa + note;
      } else {
        tooltip.style.display = 'none'; return;
      }
    } else if (yp >= readsY && yp < readsY + nRows * cellH) {
      var row_idx = Math.floor((yp - readsY) / cellH);
      if (row_idx >= 0 && row_idx < nRows) {
        var ch = rows[row_idx][col];
        var label = ch === '.' ? refSeq[col] + ' (match)' : ch === '-' ? 'gap' : ch + ' (mismatch)';
        tooltip.textContent = rl + 'Read ' + (row_idx + 1) + ', pos ' + (col + 1) + ': ' + label;
      } else {
        tooltip.style.display = 'none'; return;
      }
    } else {
      tooltip.style.display = 'none'; return;
    }
    tooltip.style.display = 'block';
    tooltip.style.left = (e.clientX + 12) + 'px';
    tooltip.style.top = (e.clientY - 24) + 'px';
  });
  canvas.addEventListener('mouseleave', function() {
    tooltip.style.display = 'none';
  });
}

/* Mirror horizontal scrolling across every pileup on the page.
 *
 * Each group scrolls in its own container, so without this a page of several
 * groups drifts out of register: the reader lines up a column in one pileup and
 * the group below it is showing somewhere else. Scrolling is what a reader does
 * to compare groups, so the groups have to move together.
 *
 * scrollLeft is mirrored as a raw pixel offset rather than as a fraction of the
 * scrollable width. Every group is drawn at the same pixels-per-base, so equal
 * offsets mean the same reference position even where the groups have different
 * reference lengths -- which is the register that matters. A fraction would put
 * them in different places on exactly those pages.
 *
 * Assigning scrollLeft here makes the browser fire scroll on the panes being
 * caught up, which re-enters this handler. Nothing guards against that beyond
 * the equality check, and nothing needs to: those panes are already at the
 * target offset, so the second pass assigns nothing and the echo stops. The
 * half-pixel tolerance is for fractional offsets at fractional zoom, where
 * assigning a value the pane cannot hold exactly would otherwise leave the two
 * nudging each other.
 */
function syncPileupScrolls() {
  var panes = [].slice.call(document.querySelectorAll('.pileup-scroll'));
  if (panes.length < 2) return;
  panes.forEach(function(pane) {
    pane.addEventListener('scroll', function() {
      var left = pane.scrollLeft;
      panes.forEach(function(other) {
        if (Math.abs(other.scrollLeft - left) > 0.5) other.scrollLeft = left;
      });
    });
  });
}
