"""Build docs/index.html from the package itself.

    python docs/build.py            # rebuild the page
    python docs/build.py --check    # fail if the committed page is stale

Anything the code already knows is read out of the code: the version, the
command list, and every flag of every command, taken from the argparse parsers
rather than transcribed.  A flag added to the CLI therefore appears on the page
at the next build, and a help string cannot drift out of step with the one the
terminal prints.

What cannot be derived is example output, since the runs need sequencing data
that is not in the repository.  Those are pasted runs held in ``examples/``.
Recapture them by running the command again and overwriting the file:

    seqview lengths reads.fastq --bins 14 > docs/examples/lengths.txt

``--check`` compares a fresh build against the committed page, ignoring the
date stamp, so CI reports a page that no longer matches the code.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

from seqviewer import __version__                              # noqa: E402
from seqviewer.cli import (                                     # noqa: E402
    build_lengths_parser, build_parser, build_quality_parser,
)

OUT = HERE / "index.html"
EXAMPLES = HERE / "examples"

#: The stamp line, matched so --check can ignore the date in it.
STAMP = re.compile(r'<p class="stamp">[^<]*</p>')


def escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def versioned(path: str) -> str:
    """*path* with a short hash of the file after it, so an edited page is
    fetched afresh rather than served from a browser's cache."""
    target = HERE / path
    if not target.exists():
        return path
    return f"{path}?v={hashlib.md5(target.read_bytes()).hexdigest()[:8]}"


def example(name: str) -> str:
    """One captured run, escaped for the page."""
    return escape(EXAMPLES.joinpath(f"{name}.txt").read_text().rstrip("\n"))


def arguments(parser):
    """Every argument *parser* takes, as ``(term, metavar, help)``.

    Read off the parser so the page lists what the command accepts rather than
    what someone last remembered it accepting.
    """
    out = []
    for action in parser._actions:
        if isinstance(action, argparse._HelpAction):
            continue
        if not action.option_strings:
            out.append((action.dest, "", action.help or ""))
            continue
        term = ", ".join(action.option_strings)
        if action.metavar:
            metavar = action.metavar
        elif action.nargs == 0:
            metavar = ""
        elif action.choices:
            metavar = "{" + ",".join(str(c) for c in action.choices) + "}"
        else:
            metavar = action.dest.upper()
        out.append((term, metavar, action.help or ""))
    return out


def flags(parser) -> str:
    """The argument list as a definition grid."""
    rows = ['<dl class="flags">']
    for term, metavar, text in arguments(parser):
        arg = f' <span class="arg">{escape(metavar)}</span>' if metavar else ""
        rows.append(f"<dt><code>{escape(term)}</code>{arg}</dt>")
        rows.append(f"<dd>{escape(text.strip())}</dd>")
    rows.append("</dl>")
    return "\n      ".join(rows)


TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<script>
try {
  var saved = localStorage.getItem("seqviewer-theme");
  if (saved === "dark" || saved === "light")
    document.documentElement.setAttribute("data-theme", saved);
} catch (e) {}
</script>
<title>seqviewer</title>
<style>
:root {
  color-scheme: light;
  --paper: #f8f6f2;
  --panel: #fffdfa;
  --ink: #221f1b;
  --muted: #58534b;
  --quiet: #766f65;
  --rule: #e5e0d8;
  --accent: var(--ink);    /* highlights are the text colour, not a hue */
  --flag: #c22f2f;         /* and the colour it flags a finding in */
  --term-bg: #efebe4;
  --edge: #cdc6bb;         /* border of code and output blocks */
  --term-ink: #2a2621;
  --rail: 7.5rem;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --paper: #2a2724; --panel: #35312c; --ink: #ece7df; --muted: #b5aea3;
    --quiet: #9d968a; --rule: #48433b; --edge: #615a50;
    --flag: #f0655f;
    --term-bg: #211f1c; --term-ink: #dcd6cb;
  }
}
:root[data-theme="light"] { color-scheme: light; }
:root[data-theme="dark"] {
  color-scheme: dark;
    --paper: #2a2724; --panel: #35312c; --ink: #ece7df; --muted: #b5aea3;
    --quiet: #9d968a; --rule: #48433b; --edge: #615a50;
    --flag: #f0655f;
    --term-bg: #211f1c; --term-ink: #dcd6cb;
}
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0; background: var(--paper); color: var(--ink);
  font: 16px/1.6 "Helvetica Neue", Helvetica, Arial, sans-serif;
}
.page { max-width: 54rem; margin: 0 auto; padding: 3rem 1.5rem 4rem; }

/* --- navigation: a row under the intro, and a ladder in the left margin
   where the window is wide enough to have one --- */
html { scroll-padding-top: 1.5rem; }
@media (prefers-reduced-motion: no-preference) {
  html { scroll-behavior: smooth; }
}
.ladder { display: flex; flex-wrap: wrap; gap: 0.3rem 1.3rem;
          padding: 0.85rem 0; border-bottom: 1px solid var(--rule); }
.ladder a { font-size: 0.9rem; line-height: 1.4; color: var(--muted);
            text-decoration: none; }
.ladder a:hover, .ladder a:focus-visible { color: var(--ink); }
.ladder a[aria-current="true"] { color: var(--ink); font-weight: 600; }
@media (min-width: 1180px) {
  .ladder { position: fixed; top: 3.1rem; width: 7rem;
            left: calc((100vw - 54rem) / 2 - 8.5rem);
            flex-direction: column; flex-wrap: nowrap; gap: 0;
            padding: 0; border-bottom: none;
            border-left: 2px solid var(--rule); }
  .ladder a { padding: 0.28rem 0 0.28rem 0.85rem; margin-left: -2px;
              border-left: 2px solid transparent; }
  .ladder a[aria-current="true"] { border-left-color: var(--accent); }
}

/* --- masthead --- */
.mast { position: relative; }
.theme { position: absolute; top: -0.2rem; right: 0; width: 2rem; height: 2rem;
         display: grid; place-items: center; padding: 0; cursor: pointer;
         color: var(--muted); background: none; border: 1px solid transparent;
         border-radius: 6px; }
.theme:hover, .theme:focus-visible { color: var(--ink);
                                     border-color: var(--edge); }
.theme svg { width: 1.1rem; height: 1.1rem; fill: none; stroke: currentColor;
             stroke-width: 1.6; stroke-linecap: round; stroke-linejoin: round; }
.mast { display: grid; grid-template-columns: var(--rail) 1fr; gap: 1.5rem;
        align-items: baseline; padding-bottom: 1.4rem;
        border-bottom: 2px solid var(--ink); }
.mast h1 { grid-column: 1; font-size: 1.5rem; line-height: 1.1;
           font-weight: 700; letter-spacing: -0.03em; margin: 0; }
/* The title sits in the rail and the intro in the content column beside it,
   on the same line as every tool's name and prose. */
.mast p { grid-column: 2; margin: 0; color: var(--muted);
          max-width: 38rem; min-width: 0; }

/* --- the rail: one row per tool --- */
.tool { display: grid; grid-template-columns: var(--rail) 1fr; gap: 1.5rem;
        padding: 2.1rem 0; border-bottom: 1px solid var(--rule); }
.tool > .label { position: sticky; top: 1.5rem; align-self: start; }
.tool > .label .name {
  display: block; font-size: 0.95rem; line-height: 1.3; font-weight: 600;
  color: var(--ink);
}
/* A grid item defaults to min-width:auto, which is a refusal to shrink below
   its widest content.  A wide transcript inside one therefore widens the
   column, the grid, and the page, and the overflow-x below never engages.
   Zero here is what lets the scroll container be the transcript itself. */
.tool > .label, .tool > .body { min-width: 0; }
.tool > .body > :first-child { margin-top: 0; }
.tool > .body > :last-child { margin-bottom: 0; }

p { margin: 0 0 0.95rem; max-width: 38rem; }

/* --- a command and what it printed --- */
.run { margin: 1.3rem 0 1.4rem; min-width: 0; overflow: hidden;
        border: 1px solid var(--edge); border-radius: 10px; }
/* Output never wraps.  It is fixed-pitch, and a line folded at the window edge
   reads as two: the histogram stops being a histogram.  It is bounded by its
   column and scrolls sideways within it, with the scrollbar hidden.  A command
   does wrap, since it has a copy button and nothing in it depends on its
   columns. */
.run pre {
  font: 0.78rem/1.55 ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas,
        monospace;
  margin: 0; padding: 0.9rem 1.15rem;
  white-space: pre; overflow-x: auto; overscroll-behavior-x: contain;
  scrollbar-width: none;
}
.run pre::-webkit-scrollbar { display: none; }
.run .cmd { background: var(--panel); color: var(--ink);
            white-space: pre-wrap; overflow-wrap: anywhere; }
.run { position: relative; }
.copy { position: absolute; top: 0.5rem; right: 0.5rem; z-index: 1;
        font: 0.72rem/1 "Helvetica Neue", Helvetica, Arial, sans-serif;
        color: var(--muted); background: var(--paper);
        border: 1px solid var(--edge); border-radius: 6px;
        width: 1.9rem; height: 1.9rem; padding: 0; cursor: pointer;
        display: grid; place-items: center; }
.copy svg { width: 1rem; height: 1rem; fill: none; stroke: currentColor;
            stroke-width: 1.6; stroke-linecap: round; stroke-linejoin: round; }
.copy { opacity: 0; transition: opacity 0.12s; }
.run .cmd:hover ~ .copy, .copy:hover, .copy:focus-visible,
.copy[data-state="done"] { opacity: 1; }
@media (hover: none) { .copy { opacity: 1; } }
.copy:hover, .copy:focus-visible { color: var(--ink); border-color: var(--quiet); }
.copy[data-state="done"] { color: var(--accent); border-color: var(--accent); }
.run .cmd { padding-right: 3.2rem; }
.run .cmd b { font-weight: 600; color: var(--accent); }
.run .out { background: var(--term-bg); color: var(--term-ink);
            border-top: 1px solid var(--rule); }

figure { margin: 1.2rem 0; }
figure iframe { display: block; width: 100%;
                border: 1px solid var(--edge); border-radius: 10px;
                background: #fff; }
figcaption { margin-top: 0.5rem; font-size: 0.86rem; color: var(--muted);
             max-width: 38rem; }

/* --- what one run turned out to be, as against how the tool works --- */
.reading { margin: 1.2rem 0; padding-left: 0.9rem;
           border-left: 2px solid var(--flag);
           color: var(--muted); font-size: 0.95rem; max-width: 38rem; }
.reading p:last-child { margin-bottom: 0; }

/* --- every argument the command takes, read off its parser --- */
.ref { margin-top: 1.6rem; min-width: 0; }
.ref h3 { font: 600 0.78rem/1.3 ui-monospace, SFMono-Regular, "SF Mono",
          Menlo, Consolas, monospace; color: var(--quiet);
          margin: 0 0 0.6rem; }
/* The term column takes what the longest flag needs, but never more than a
   third of the row: --min-overlap-pos MIN_OVERLAP_POS would otherwise set the
   column width and squeeze the descriptions into a ribbon. */
dl.flags { display: grid;
           grid-template-columns: minmax(0, min(max-content, 33%)) 1fr;
           gap: 0 1.2rem; margin: 0; font-size: 0.9rem; }
dl.flags dt, dl.flags dd { padding: 0.34rem 0; min-width: 0;
                           overflow-wrap: anywhere;
                           border-bottom: 1px solid var(--rule); }
dl.flags dd { margin: 0; color: var(--muted); }
dl.flags dt:last-of-type, dl.flags dd:last-of-type { border-bottom: none; }
dl.flags dt code { background: none; padding: 0; font-weight: 600;
                   color: var(--ink); }
dl.flags .arg { font: 0.78rem/1 ui-monospace, SFMono-Regular, "SF Mono", Menlo,
                Consolas, monospace; color: var(--quiet); }

code { font: 0.87em/1 ui-monospace, SFMono-Regular, "SF Mono", Menlo,
       Consolas, monospace;
       background: var(--term-bg); padding: 0.08em 0.3em; }
a { color: var(--accent); text-decoration: underline;
    text-underline-offset: 0.15em; text-decoration-thickness: 0.06em; }
a:focus-visible { outline: 2px solid var(--accent); outline-offset: 3px; }

.stamp { margin: 1.6rem 0 0; font-size: 0.8rem; color: var(--quiet); }

@media (max-width: 720px) {
  .mast, .tool { grid-template-columns: 1fr; gap: 0.5rem; }
  .mast p { grid-column: 1; }
  .mast h1 { font-size: 1.6rem; margin-bottom: 0.2rem; }
  .tool > .label { position: static; }
  .page { padding: 2.2rem 1.1rem 3rem; }
  dl.flags { grid-template-columns: 1fr; }
  dl.flags dt { border-bottom: none; padding-bottom: 0; }
  dl.flags dd { padding-top: 0.1rem; padding-bottom: 0.7rem; }
}
</style>
</head>
<body>
<main class="page">

<header class="mast">
  <h1>seqviewer</h1>
  <p>Lightweight tools for interpreting <code>fastqs</code></p>
</header>

<nav class="ladder" aria-label="Sections">
  <a href="#install">install</a>
  <a href="#lengths">lengths</a>
  <a href="#qualities">qualities</a>
  <a href="#pileup">pileup</a>
</nav>

<section class="tool" id="install">
  <div class="label">
    <span class="name">install</span>
  </div>
  <div class="body">
    <p>Aligning reads also needs <code>minimap2</code> and
    <code>samtools</code> on your <code>PATH</code>.</p>
    <div class="run">
      <pre class="cmd">$ <b>uv tool install</b> --python 3.13 'seqviewer[all] @ git+https://github.com/micah-olivas/seqviewer'</pre>
    </div>
    <p>Run <code>seqview</code> with no arguments to list the commands.</p>
    <div class="run">
      <pre class="cmd">$ <b>seqview</b></pre>
      <pre class="out">__SEQVIEW__</pre>
    </div>
  </div>
</section>

<section class="tool" id="lengths">
  <div class="label">
    <span class="name">lengths</span>
  </div>
  <div class="body">
    <p>This command needs no reference or aligner. It tallies lengths as it
    reads, so memory use depends on the range of lengths, not the number of
    reads, and a multi-gigabyte FASTQ takes one pass.</p>

    <div class="run">
      <pre class="cmd">$ <b>seqview lengths</b> sample_reads.fastq --bins 14</pre>
      <pre class="out">__LENGTHS__</pre>
    </div>

    <div class="reading">
      <p>This run is two populations: 4,324 reads of 100&ndash;245&nbsp;bp, and
      260 that span the full 3.5&nbsp;kb construct. The short ones are 81% of
      the reads and 29% of the bases, which is why the mean, 453, sits so far
      from the N50, 1,304.</p>
    </div>

    <p>The axis covers the central 99% of reads, not the full range, because a few
    concatemers would otherwise stretch it on their own. Reads outside the axis
    get a row of their own, labelled with the extreme they reach. The figures
    under the histogram cover the whole run. The longest read here is
    4,700&nbsp;bp, which is past the last bin.</p>

    <div class="ref">
      <h3>seqview lengths [options] reads</h3>
      __LENGTHS_FLAGS__
    </div>
  </div>
</section>

<section class="tool" id="qualities">
  <div class="label">
    <span class="name">qualities</span>
  </div>
  <div class="body">
    <p>This command plots how many bases the sequencer scored at each Phred
    value. Scores are per base, not per read, so a read that is good for most of
    its length and poor at the end still shows its poor bases. Like
    <code>lengths</code> it keeps only a tally, one counter per score, so memory
    use is the same for any file size and a multi-gigabyte FASTQ takes one
    pass.</p>

    <div class="run">
      <pre class="cmd">$ <b>seqview qualities</b> sample_reads.fastq --bins 14</pre>
      <pre class="out">__QUALITY__</pre>
    </div>

    <div class="reading">
      <p>Most of this run sits at Q30&ndash;38, with 83.6% of bases at Q30 or
      above. The tail toward low scores holds the poor reads and the ends of
      reads, where quality falls.</p>
    </div>

    <p>By default every base is one count. <code>--by read</code> counts each
    read once instead, at the mean of its bases rounded to a whole score, so the
    plot says how many reads are good rather than how many bases. The mean is
    taken over error rates, as NanoPlot does, since a score is the logarithm of
    an error rate. <code>--mean phred</code> averages the scores instead.</p>

    <div class="run">
      <pre class="cmd">$ <b>seqview qualities</b> sample_reads.fastq --bins 14 --by read</pre>
      <pre class="out">__QUALITY_BY_READ__</pre>
    </div>

    <div class="reading">
      <p>By read, only 47.8% are at Q30 or above, against 83.6% of bases. Quality
      falls along each read, and the poor end pulls the read's mean down.
      <code>--mean phred</code> gives 91.0% here.</p>
    </div>

    <p>Scores are read as Phred+33, which Illumina, nanopore and most other FASTQ
    use. Pass <code>--offset 64</code> for the older Illumina 1.3&ndash;1.7
    files. If any base scores below the offset, the summary says so, which
    usually means the wrong offset.</p>

    <div class="ref">
      <h3>seqview qualities [options] reads</h3>
      __QUALITY_FLAGS__
    </div>
  </div>
</section>

<section class="tool" id="pileup">
  <div class="label">
    <span class="name">pileup</span>
  </div>
  <div class="body">
    <p>This command aligns a run with <code>minimap2</code> and writes two
    pages. The pileup shows every read base by base. The summary condenses the
    same alignment to one screen and tells you whether the clone has an error,
    where the pileup only shows what each read says. The pages link to each
    other.</p>

    <div class="run">
      <pre class="cmd">$ <b>seqview pileup</b> reads.fastq vector.dna out.html --max 300</pre>
      <pre class="out">__PILEUP__</pre>
    </div>

    <p>Each page is a single self-contained HTML file, with no server or
    external assets, so it opens in any browser and can be emailed or kept with
    the run. The command also writes a log beside them. The summary is shown
    below as it opens, drawn from synthetic reads.</p>

    <figure>
      <iframe title="Summary page, rendered from synthetic data"
              src="__SUMMARY_PAGE__" style="height: 20rem" loading="lazy"></iframe>
      <figcaption>The summary page for one clone, live, from synthetic reads
      with a missense change and an in-frame deletion planted in it. Three bands, top to bottom: how far the reads
      disagree with the reference at each position, the annotated reference,
      and read depth. Bars are muted below 10%, the rate sequencing error alone
      produces, and coloured at or above it.
      <a href="demo/summary.html">Open it full size</a>, or see
      <a href="demo/pileup.html">the pileup</a> from the same kind of
      data.</figcaption>
    </figure>

    <p>A reference of several kilobases is drawn at a few bases per pixel, so each
    column shows the worst position it covers, not the mean. Averaging would
    hide a single position where every read disagrees.</p>

    <p>The pileup page below is drawn from synthetic reads against a 1&nbsp;kb
    reference. All of the clone's reads disagree with the reference at one
    position.</p>

    <figure>
      <iframe title="Pileup page, rendered from synthetic data"
              src="__PILEUP_PAGE__" style="height: 46rem" loading="lazy"></iframe>
      <figcaption><a href="demo/pileup.html">Open it full size</a>.</figcaption>
    </figure>

    <div class="ref">
      <h3>seqview pileup [options] reads reference out</h3>
      __PILEUP_FLAGS__
    </div>
  </div>
</section>

<p class="stamp">__STAMP__</p>

</main>
<script>
/* Light and dark.  Follows the system until the button is used, then keeps the
   choice.  The embedded pages read the same key, and are told directly so they
   change without a reload. */
(function () {
  var root = document.documentElement;
  var MOON = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M13.5 9.6A5.6 5.6 0 0 1 6.4 2.5a5.6 5.6 0 1 0 7.1 7.1z"/></svg>';
  var SUN = '<svg viewBox="0 0 16 16" aria-hidden="true"><circle cx="8" cy="8" r="3"/><path d="M8 1.5v1.6M8 12.9v1.6M1.5 8h1.6M12.9 8h1.6M3.4 3.4l1.1 1.1M11.5 11.5l1.1 1.1M3.4 12.6l1.1-1.1M11.5 4.5l1.1-1.1"/></svg>';
  var button = document.createElement("button");
  button.type = "button"; button.className = "theme";

  function current() {
    return root.getAttribute("data-theme") ||
      (window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches
        ? "dark" : "light");
  }
  function tell(frame) {
    var doc = frame.contentDocument;
    if (!doc || !doc.documentElement) return;
    if (current() === "dark") doc.documentElement.setAttribute("data-theme", "dark");
    else doc.documentElement.removeAttribute("data-theme");
  }
  function show() {
    var dark = current() === "dark";
    button.innerHTML = dark ? SUN : MOON;
    button.setAttribute("aria-label", dark ? "Switch to light mode" : "Switch to dark mode");
    document.querySelectorAll("figure iframe").forEach(tell);
  }
  /* Colours cross-fade for a moment while the theme changes.  The rule exists
     only for that moment, so nothing else on the page is slowed by it. */
  var FADE = "html.theme-fade, html.theme-fade *, html.theme-fade *::before, " +
    "html.theme-fade *::after { transition: background-color .3s ease, " +
    "color .3s ease, border-color .3s ease, fill .3s ease, stroke .3s ease " +
    "!important; }";
  function fade(doc) {
    try {
      if (!doc || !doc.documentElement) return;
      if (window.matchMedia &&
          matchMedia("(prefers-reduced-motion: reduce)").matches) return;
      if (!doc.getElementById("theme-fade-style")) {
        var style = doc.createElement("style");
        style.id = "theme-fade-style"; style.textContent = FADE;
        (doc.head || doc.documentElement).appendChild(style);
      }
      var html = doc.documentElement;
      html.classList.add("theme-fade");
      clearTimeout(html.fadeTimer);
      html.fadeTimer = setTimeout(function () {
        html.classList.remove("theme-fade");
      }, 450);
    } catch (e) {}
  }
  button.addEventListener("click", function () {
    fade(document);
    document.querySelectorAll("figure iframe").forEach(function (frame) {
      fade(frame.contentDocument);
    });
    var next = current() === "dark" ? "light" : "dark";
    root.setAttribute("data-theme", next);
    try { localStorage.setItem("seqviewer-theme", next); } catch (e) {}
    show();
  });
  document.querySelectorAll("figure iframe").forEach(function (frame) {
    frame.addEventListener("load", function () { tell(frame); });
  });
  var mast = document.querySelector(".mast");
  if (mast) mast.appendChild(button);
  show();
})();
/* Mark the section in view in the navigation. */
(function () {
  var links = Array.prototype.slice.call(
      document.querySelectorAll(".ladder a"));
  var sections = links.map(function (a) {
    return document.querySelector(a.getAttribute("href"));
  });
  function mark() {
    var line = window.innerHeight * 0.3, here = 0;
    sections.forEach(function (sec, i) {
      if (sec && sec.getBoundingClientRect().top <= line) here = i;
    });
    links.forEach(function (a, i) {
      if (i === here) a.setAttribute("aria-current", "true");
      else a.removeAttribute("aria-current");
    });
  }
  window.addEventListener("scroll", mark, { passive: true });
  window.addEventListener("resize", mark);
  mark();
})();
/* Size each embedded page to its content, and again when it changes size, as
   when a section of it is opened. */
document.querySelectorAll("figure iframe").forEach(function (frame) {
  function fit() {
    var doc = frame.contentDocument;
    if (!doc || !doc.body) return;
    /* The body's own height, since the document's cannot be less than the
       frame's and so would never let it shrink. */
    var style = frame.contentWindow.getComputedStyle(doc.body);
    var margins = parseFloat(style.marginTop) + parseFloat(style.marginBottom);
    var height = Math.max(doc.body.getBoundingClientRect().height,
                          doc.body.scrollHeight);
    frame.style.height = Math.ceil(height + margins) + 2 + "px";
    if (window.ResizeObserver && !frame.watched) {
      frame.watched = true;
      new ResizeObserver(fit).observe(doc.body);
    }
  }
  frame.addEventListener("load", fit);
  if (frame.contentDocument && frame.contentDocument.readyState === "complete") fit();
});
/* Copy button on each command.  Takes the text without the "$ " prompt. */
document.querySelectorAll(".run").forEach(function (run) {
  var cmd = run.querySelector(".cmd");
  if (!cmd || !navigator.clipboard) return;
  var btn = document.createElement("button");
  var ICON = '<svg viewBox="0 0 16 16" aria-hidden="true"><rect x="5.5" y="5.5" width="8" height="8" rx="1.5"/><path d="M10.5 3.5v-.5a1.5 1.5 0 0 0-1.5-1.5H3.5A1.5 1.5 0 0 0 2 3v5.5A1.5 1.5 0 0 0 3.5 10H4"/></svg>';
  var DONE = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M3 8.5l3.2 3.2L13 4.8"/></svg>';
  btn.type = "button"; btn.className = "copy"; btn.innerHTML = ICON;
  btn.setAttribute("aria-label", "Copy command");
  btn.addEventListener("click", function () {
    var text = cmd.textContent.replace(/^\\$ /, "").trim();
    navigator.clipboard.writeText(text).then(function () {
      btn.innerHTML = DONE; btn.dataset.state = "done";
      btn.setAttribute("aria-label", "Copied");
      setTimeout(function () {
        btn.innerHTML = ICON; delete btn.dataset.state;
        btn.setAttribute("aria-label", "Copy command");
      }, 1600);
    });
  });
  run.appendChild(btn);
});
</script>
</body>
</html>
"""


def build(today=None) -> str:
    stamp = (f"seqviewer {__version__}. Built from the package by "
             f"docs/build.py on "
             f"{today or datetime.date.today().isoformat()}.")
    return (TEMPLATE
            .replace("__SEQVIEW__", example("seqview"))
            .replace("__LENGTHS__", example("lengths"))
            .replace("__QUALITY__", example("quality"))
            .replace("__QUALITY_BY_READ__", example("quality_by_read"))
            .replace("__SUMMARY_PAGE__", versioned("demo/summary.html"))
            .replace("__PILEUP_PAGE__", versioned("demo/pileup.html"))
            .replace("__PILEUP__", example("pileup"))
            .replace("__LENGTHS_FLAGS__", flags(build_lengths_parser(
                "seqview lengths")))
            .replace("__QUALITY_FLAGS__", flags(build_quality_parser(
                "seqview qualities")))
            .replace("__PILEUP_FLAGS__", flags(build_parser("seqview pileup")))
            .replace("__STAMP__", escape(stamp)))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="docs/build.py", description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true",
                        help="report whether the committed page still matches "
                             "the code, without writing it. The date stamp is "
                             "ignored, so a rebuild on a later day is not a "
                             "difference")
    args = parser.parse_args(argv)
    fresh = build()

    if args.check:
        if not OUT.exists():
            print(f"{OUT} does not exist; run: python docs/build.py",
                  file=sys.stderr)
            return 1
        if STAMP.sub("", OUT.read_text()) != STAMP.sub("", fresh):
            print(f"{OUT} is out of date; run: python docs/build.py",
                  file=sys.stderr)
            return 1
        print(f"{OUT.name} is up to date")
        return 0

    OUT.write_text(fresh)
    print(f"wrote {OUT}  ({len(fresh) / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
