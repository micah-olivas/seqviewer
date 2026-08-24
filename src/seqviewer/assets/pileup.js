function drawPileup(canvasId, rulerId, labelsId, refSeq, cons, rows, flanks, scrollId, wrapId, refAA, consAA, flaggedCols, cellWIn, annotH, mismatchH, parent, mmRuns, parentAA) {
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
  // The parent is the library's baseline, so its translation is the row a
  // mutation is named against. Drawn last of the three because it is the one a
  // reader compares the other two back to.
  var hasParentAA = !!(hasAA && parentAA && parentAA.length);
  var aaRows = hasAA ? (hasParentAA ? 3 : 2) : 0;
  var aaH = hasAA ? 14 : 0;       // height of each AA row
  var aaBlockH = aaRows ? aaRows * aaH + (aaRows - 1) * 2 : 0;
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

  /* Mutation tabs: one per residue the consensus changes, named the way a
   * mutation is written -- ref residue, position, new residue, so V50I.
   *
   * Laid out before the canvas is sized, because how many rows of tabs are
   * needed decides how tall the canvas has to be. A tab is far wider than the
   * 6px codon it belongs to, so tabs are packed into rows by first fit and a
   * stem connects each back to its own codon. Widths are computed here rather
   * than measured: the face is monospace, so its advance is known, and the
   * layout has to exist before there is a context to measure with.
   */
  var TAB_FONT = 9, TAB_ADVANCE = 0.6, TAB_PAD = 4;
  var TAB_H = 12, TAB_ROW_GAP = 2, TAB_LEAD = 5, TAB_STEM = 4;
  var MAX_TAB_ROWS = 4;
  var mutTabs = [], mutRows = 0, mutDropped = 0;
  if (hasAA) {
    var rowEnds = [];
    /* What a tab names depends on what the page knows.
     *
     * With a parent, that is the library's baseline -- a WT -- and mutations are
     * named against it, which is the convention and the only numbering that
     * means anything to a reader. The reference is this well's assigned
     * identity, itself a variant, so a change named against it would be
     * numbered off an already-mutated sequence.
     *
     * Knowing both then answers the question the well exists to answer, and each
     * tab carries which of the three it is:
     *
     *   expected    the consensus moved off the parent, to what the assignment
     *               said it would be. The designed mutation, confirmed.
     *   unexpected  it moved off the parent to something else.
     *   missing     it did not move, but the assignment said it should have.
     *
     * With no parent there is only one baseline available, so a tab names the
     * reference-to-consensus change and has nothing to say about intent.
     */
    var namedAgainst = hasParentAA ? parentAA : refAA;
    var limit = Math.min(refAA.length, consAA.length,
                         hasParentAA ? parentAA.length : refAA.length);
    for (var mi = 0; mi < limit; mi++) {
      var base = namedAgainst[mi], obs = consAA[mi], want = refAA[mi];
      var kind, name;
      if (base !== obs) {
        name = base + (mi + 1) + obs;
        kind = !hasParentAA ? 'plain' : (obs === want ? 'expected' : 'unexpected');
      } else if (hasParentAA && want !== base) {
        // Expected but absent: the well reads as the parent where its assigned
        // identity says it should not.
        name = base + (mi + 1) + want;
        kind = 'missing';
      } else {
        continue;
      }
      var tabW = name.length * TAB_FONT * TAB_ADVANCE + TAB_PAD * 2;
      var stemX = flanks[0] * cellW + mi * aaCodonW + aaCodonW / 2;
      // Keep the tab on the canvas; the stem still points at the real codon.
      var left = Math.max(0, Math.min(stemX - tabW / 2, canvasW - tabW));
      var placed = -1;
      for (var r = 0; r < rowEnds.length; r++) {
        if (left >= rowEnds[r] + 2) { placed = r; break; }
      }
      if (placed < 0) {
        if (rowEnds.length >= MAX_TAB_ROWS) { mutDropped++; continue; }
        rowEnds.push(0);
        placed = rowEnds.length - 1;
      }
      rowEnds[placed] = left + tabW;
      mutTabs.push({name: name, stemX: stemX, left: left, w: tabW, row: placed,
                    kind: kind});
    }
    mutRows = rowEnds.length;
  }
  var mutH = mutRows ? TAB_LEAD + TAB_STEM + mutRows * TAB_H
                       + (mutRows - 1) * TAB_ROW_GAP : 0;

  var pileupH = refH + parentGap + parentH + gap + consH + gap + nRows * cellH + aaGap + aaBlockH + mutH;
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
  var consMatchColor = matchColor;
  var baseColors = {'A': P.a, 'T': P.t, 'C': P.c, 'G': P.g};
  // Flagged columns: positions past the flag threshold, from
  // summary.flagged_columns. Used only for the off-screen arrows now -- the
  // mismatch track above the reference draws the magnitude itself, so the
  // separate row of triangles that used to say "past the threshold here" is
  // gone. Falls back to consensus-derived mismatches when not provided.
  var mismatchCols = flaggedCols || [];
  if (!flaggedCols) {
    for (var _mi = 0; _mi < cons.length; _mi++) {
      var _ch = cons[_mi];
      if (_ch !== '.' && _ch !== '-') mismatchCols.push(_mi);
    }
  }
  function isVector(col) {
    return flanks && (col < flanks[0] || col >= nCols - flanks[1]);
  }
  function pickMatch(col) {
    return isVector(col) ? vectorMatchColor : matchColor;
  }
  // --- Ruler ---
  var rulerH = 14;
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
    var tickBottom = rulerH;
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
      if (hasParentAA) {
        var paaGap = document.createElement('span');
        paaGap.style.height = '2px';
        labelsEl.appendChild(paaGap);
        var parentAALabel = document.createElement('span');
        parentAALabel.textContent = 'Parent AA';
        parentAALabel.style.height = aaH + 'px';
        parentAALabel.title = 'The library baseline. Mutations are named '
          + 'against it, not against the well\u2019s assigned reference, which '
          + 'is itself a variant.';
        labelsEl.appendChild(parentAALabel);
      }
      if (mutRows) {
        var mutLabel = document.createElement('span');
        mutLabel.textContent = 'Changes';
        mutLabel.style.height = mutH + 'px';
        var named = mutTabs.length;
        // Say what was left out rather than let a page look complete when it
        // is not; the gutter labels already explain themselves on hover.
        mutLabel.title = named + (named === 1 ? ' residue changes'
                                             : ' residues change')
          + (mutDropped ? ', and ' + mutDropped + ' more not named for want '
                          + 'of room' : '');
        labelsEl.appendChild(mutLabel);
      }
    }
  }
  // --- Reference row ---
  // The match colour, not a dark bar of its own. Every cell agreeing with the
  // reference is drawn in it, so the reference being the same colour says what
  // the row is: the thing the greys below are agreeing with. A near-black bar
  // read as a separate kind of thing, and there is only one.
  //
  // Region-aware, like the read cells: a flank column takes the vector shade,
  // so a column is one colour from the reference down through the reads.
  for (var i = 0; i < nCols; i++) {
    ctx.fillStyle = pickMatch(i);
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
    ctx.fillRect(insX, aaY, insW, aaBlockH);

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

      // Parent AA row: the library baseline. Marked where the consensus has
      // moved away from it, which is what a mutation name describes.
      if (hasParentAA) {
        var pAA = parentAA[ai];
        var paaY = caaY + aaH + 2;
        var fromParent = pAA === cAA;
        if (!fromParent) {
          ctx.fillStyle = aaDiffBg;
          ctx.fillRect(ax, paaY, aaCodonW, aaH);
        }
        ctx.fillStyle = fromParent ? aaMatchColor : aaDiffColor;
        drawResidue(pAA === undefined ? '?' : pAA, ax, paaY, fromParent);
      }
    }

    // Subtle codon grid lines, but only when the pitch is wide enough that a
    // line per codon reads as structure rather than as hatching.
    ctx.strokeStyle = P['aa-grid'];
    ctx.lineWidth = 0.5;
    for (var ai = 1; aaCodonW >= 9 && ai < refAA.length; ai++) {
      var lx = insStart * cellW + ai * aaCodonW;
      ctx.beginPath(); ctx.moveTo(lx, aaY); ctx.lineTo(lx, aaY + aaBlockH); ctx.stroke();
    }

    // --- Mutation tabs ---
    // The AA rows say a residue changed; these say which change it was, in the
    // notation a person writes it in. A stem ties each tab to its own codon,
    // because a tab is five times the width of the codon it names and cannot
    // sit over it.
    if (mutRows) {
      var tabTop = aaY + aaBlockH + TAB_LEAD;
      ctx.font = TAB_FONT + 'px SF Mono,Menlo,Consolas,monospace';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      function tabTopOf(tab) {
        return tabTop + TAB_STEM + tab.row * (TAB_H + TAB_ROW_GAP);
      }
      // Stems first, boxes second. A tab on a lower row reaches past the rows
      // above it, and those tabs are wider than the codon pitch, so the stem
      // has to pass behind them -- drawing every box after every stem is what
      // makes it pass behind rather than through.
      function tabInk(kind) {
        if (kind === 'expected') return [P['tab-ok'], P['tab-ok-bg'], false];
        if (kind === 'missing') return [P['tab-warn'], P['tab-warn-bg'], true];
        return [aaDiffColor, aaDiffBg, false];
      }
      ctx.lineWidth = 1;
      for (var ti = 0; ti < mutTabs.length; ti++) {
        var tab = mutTabs[ti];
        var top = tabTopOf(tab);
        var centre = tab.left + tab.w / 2;
        ctx.strokeStyle = tabInk(tab.kind)[0];
        ctx.beginPath();
        ctx.moveTo(tab.stemX, aaY + aaBlockH);
        // Down to just above the tab, then across if the tab had to shift to
        // stay on the canvas. Stopping at the top edge rather than the middle
        // keeps the line out of the text.
        ctx.lineTo(tab.stemX, top - 2);
        ctx.lineTo(centre, top);
        ctx.stroke();
      }
      for (var ti = 0; ti < mutTabs.length; ti++) {
        var tab = mutTabs[ti];
        var top = tabTopOf(tab);
        // Opaque ground first, then the tint: aa-diff-bg is a translucent
        // rgba, so tinting alone would let the stems behind the tab show
        // through its text.
        var ink = tabInk(tab.kind);
        ctx.fillStyle = P['aa-bg'];
        ctx.fillRect(tab.left, top, tab.w, TAB_H);
        ctx.fillStyle = ink[1];
        ctx.fillRect(tab.left, top, tab.w, TAB_H);
        ctx.strokeStyle = ink[0];
        // A dashed edge for a change the assignment promised and the reads do
        // not show: the name is what was expected, not what is there.
        ctx.setLineDash(ink[2] ? [3, 2] : []);
        ctx.strokeRect(tab.left + 0.5, top + 0.5, tab.w - 1, TAB_H - 1);
        ctx.setLineDash([]);
        ctx.fillStyle = ink[0];
        ctx.fillText(tab.name, tab.left + tab.w / 2, top + TAB_H / 2);
      }
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
  // --- Opening position: the insert, not base 1 ---
  /* A page is opened to look at the insert; the vector around it is context. So
   * the first thing on screen is the insert rather than whatever happens to sit
   * at the left edge of the reference.
   *
   * Set before the overflow arrows attach, so their first update reads the
   * position the reader will actually see. Set on every group, not just one and
   * mirrored: the mirror only fires on a scroll event, and a group whose pane is
   * not scrollable never fires one.
   */
  if (flanks && scrollId) {
    var openEl = document.getElementById(scrollId);
    if (openEl) {
      var insLeft = flanks[0] * cellW;
      var insRight = (nCols - flanks[1]) * cellW;
      var viewW = openEl.clientWidth;
      var target;
      if (insRight - insLeft <= viewW) {
        // It fits: centre it, so both boundary lines are on screen at once.
        target = (insLeft + insRight) / 2 - viewW / 2;
      } else {
        // Wider than the pane: open at its 5' boundary, where reading starts,
        // rather than somewhere in the middle of it. The margin keeps the
        // boundary line itself visible instead of flush against the edge.
        target = insLeft - 12;
      }
      var maxLeft = openEl.scrollWidth - viewW;
      openEl.scrollLeft = Math.max(0, Math.min(target, maxLeft));
    }
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
    } else if (hasAA && yp >= aaY && yp < aaY + aaBlockH) {
      var insStart = flanks[0];
      // Derived from the drawn geometry, not from the nucleotide rate: reading
      // the residue index a different way than it was placed is what made the
      // tooltip name a third, differently wrong residue.
      var aaIdx = Math.floor((x - insStart * cellW) / aaCodonW);
      if (aaIdx >= 0 && aaIdx < refAA.length) {
        var band = Math.min(aaRows - 1,
                            Math.floor((yp - aaY) / (aaH + 2)));
        var names = ['Ref', 'Cons', 'Parent'];
        var seqs = [refAA, consAA, parentAA];
        var aa = seqs[band] ? seqs[band][aaIdx] : undefined;
        // Say the residue, then what the other rows read at the same codon, so
        // hovering any one row answers the comparison rather than half of it.
        var others = [];
        for (var bi = 0; bi < aaRows; bi++) {
          if (bi === band || !seqs[bi]) continue;
          others.push(names[bi] + ' ' + seqs[bi][aaIdx]);
        }
        tooltip.textContent = names[band] + ' AA ' + (aaIdx + 1) + ': ' + aa +
          (others.length ? '  (' + others.join(', ') + ')' : '');
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

  // --- Mismatch track readout ---
  /* The track shows how much a position disagrees; hovering says how much, in
   * reads rather than in pixels. A bar's height is a log-scaled fraction, which
   * is readable as "more" or "less" and not as a number, so the number has to be
   * available some other way.
   *
   * One handler over the whole track rather than an element per base. A
   * per-position rect with a title would be thousands of nodes on a plasmid, and
   * the track is already one path for exactly that reason.
   *
   * The counts arrive run-length encoded -- a mostly-clean reference is long runs
   * of the same pair -- and are expanded once here rather than searched per
   * mouse move.
   */
  if (mmRuns && mmRuns.length && scrollId) {
    var mmHost = document.getElementById(scrollId);
    var mmSvg = mmHost && mmHost.querySelector('svg.sv-mf');
    if (mmSvg) {
      var mmDis = new Array(nCols), mmCov = new Array(nCols), mmAt = 0;
      for (var mr = 0; mr < mmRuns.length; mr++) {
        for (var mk = 0; mk < mmRuns[mr][0] && mmAt < nCols; mk++) {
          mmDis[mmAt] = mmRuns[mr][1];
          mmCov[mmAt] = mmRuns[mr][2];
          mmAt++;
        }
      }
      mmSvg.addEventListener('mousemove', function(e) {
        var r = mmSvg.getBoundingClientRect();
        var col = Math.floor((e.clientX - r.left) / cellW);
        if (col < 0 || col >= nCols) { tooltip.style.display = 'none'; return; }
        var dis = mmDis[col], cov = mmCov[col];
        var text = regionLabel(col) + 'pos ' + (col + 1) + ': ';
        if (!cov) {
          text += 'no reads cover this position';
        } else {
          text += dis + ' of ' + cov + ' read' + (cov === 1 ? '' : 's') +
                  ' disagree' + (dis === 1 ? 's' : '') +
                  ' (' + (100 * dis / cov).toFixed(dis && dis / cov < 0.1 ? 2 : 1) + '%)';
        }
        tooltip.textContent = text;
        tooltip.style.display = 'block';
        tooltip.style.left = (e.clientX + 12) + 'px';
        tooltip.style.top = (e.clientY - 24) + 'px';
      });
      mmSvg.addEventListener('mouseleave', function() {
        tooltip.style.display = 'none';
      });
    }
  }
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
