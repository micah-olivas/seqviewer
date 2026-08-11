function drawPileup(canvasId, rulerId, labelsId, refSeq, cons, rows, flanks, scrollId, wrapId, refAA, consAA, flaggedCols, cellWIn, annotH) {
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
  var rulerH = (flanks ? 24 : 14) + triRowH;
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
    // Region labels on top row (if flanks present)
    var tickRowY = 0;
    if (flanks) {
      tickRowY = 11;
      rc.fillStyle = boundaryColor;
      rc.font = '9px SF Mono,Menlo,Consolas,monospace';
      rc.textAlign = 'center';
      rc.textBaseline = 'top';
      var bLeft = flanks[0] * cellW;
      var bRight = (nCols - flanks[1]) * cellW;
      var minLabelPx = 40;
      if (flanks[0] > 0 && flanks[0] * cellW > minLabelPx) {
        rc.fillText("5′ vector", bLeft / 2, 0);
      }
      var insertW = bRight - bLeft;
      if (insertW > minLabelPx) {
        rc.fillText("insert", bLeft + insertW / 2, 0);
      }
      if (flanks[1] > 0 && flanks[1] * cellW > minLabelPx) {
        rc.fillText("3′ vector", bRight + (totalW - bRight) / 2, 0);
      }
      // Boundary dashed lines (start below both text rows, stop above triangle row)
      rc.setLineDash([3, 2]);
      rc.strokeStyle = boundaryColor;
      rc.lineWidth = 1;
      var dashY = tickRowY + 12;
      if (flanks[0] > 0 && dashY < tickBottom) {
        rc.beginPath(); rc.moveTo(bLeft, dashY); rc.lineTo(bLeft, tickBottom); rc.stroke();
      }
      if (flanks[1] > 0 && dashY < tickBottom) {
        rc.beginPath(); rc.moveTo(bRight, dashY); rc.lineTo(bRight, tickBottom); rc.stroke();
      }
      rc.setLineDash([]);
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
  var consY = refH + gap;
  var readsY = consY + consH + gap;
  if (labelsEl) {
    labelsEl.innerHTML = '';
    if (annotH) {
      var annotLabel = document.createElement('span');
      annotLabel.textContent = 'Features';
      annotLabel.style.height = annotH + 'px';
      labelsEl.appendChild(annotLabel);
    }
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
    var aaFont = Math.min(aaH - 2, Math.max(7, aaCodonW - 2));
    ctx.font = aaFont + 'px SF Mono,Menlo,Consolas,monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';

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
      if (aaCodonW >= 7) {
        ctx.fillText(rAA, ax + aaCodonW / 2, aaY + aaH / 2);
      } else {
        ctx.fillRect(ax + 1, aaY + 2, aaCodonW - 2, aaH - 4);
      }

      // Cons AA row
      var caaY = aaY + aaH + 2;
      if (!match) {
        ctx.fillStyle = aaDiffBg;
        ctx.fillRect(ax, caaY, aaCodonW, aaH);
      }
      ctx.fillStyle = match ? aaMatchColor : aaDiffColor;
      if (aaCodonW >= 7) {
        ctx.fillText(cAA, ax + aaCodonW / 2, caaY + aaH / 2);
      } else {
        ctx.fillRect(ax + 1, caaY + 2, aaCodonW - 2, aaH - 4);
      }
    }

    // Subtle codon grid lines
    ctx.strokeStyle = P['aa-grid'];
    ctx.lineWidth = 0.5;
    for (var ai = 1; ai < refAA.length; ai++) {
      var lx = insStart * cellW + ai * aaCodonW;
      ctx.beginPath(); ctx.moveTo(lx, aaY); ctx.lineTo(lx, aaY + aaH * 2 + 2); ctx.stroke();
    }
  }
  // --- Region boundary dashed lines on pileup canvas ---
  if (flanks) {
    ctx.save();
    ctx.setLineDash([4, 3]);
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
      var leftArrow = document.createElement('div');
      leftArrow.className = 'pileup-mm-arrow pileup-mm-arrow-l';
      leftArrow.textContent = '◄';
      leftArrow.style.top = (annotH || 0) + 'px';
      leftArrow.style.height = rulerH + 'px';
      wrapEl.appendChild(leftArrow);
      var rightArrow = document.createElement('div');
      rightArrow.className = 'pileup-mm-arrow pileup-mm-arrow-r';
      rightArrow.textContent = '►';
      rightArrow.style.top = (annotH || 0) + 'px';
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
      var aaIdx = Math.floor((col - insStart) / 3);
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
