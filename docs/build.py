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
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

from seqviewer import __version__                              # noqa: E402
from seqviewer.cli import build_lengths_parser, build_parser   # noqa: E402

OUT = HERE / "index.html"
EXAMPLES = HERE / "examples"

#: The stamp line, matched so --check can ignore the date in it.
STAMP = re.compile(r'<p class="stamp">[^<]*</p>')


def escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


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
<title>seqviewer</title>
<style>
:root {
  color-scheme: light dark;
  --paper: #f7f8f9;
  --panel: #ffffff;
  --ink: #16212b;
  --muted: #4a5560;
  --quiet: #667380;
  --rule: #dde1e5;
  --accent: #1f77b4;       /* the colour the package plots its figures in */
  --flag: #c22f2f;         /* and the colour it flags a finding in */
  --term-bg: #eef1f3;
  --term-ink: #1d262e;
  --rail: 7.5rem;
}
@media (prefers-color-scheme: dark) {
  :root {
    --paper: #10151a; --panel: #161d24; --ink: #e4e9ed; --muted: #a9b4be;
    --quiet: #8795a1; --rule: #26313b; --accent: #6cb6e8; --flag: #f0655f;
    --term-bg: #0c1116; --term-ink: #d3dae0;
  }
}
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0; background: var(--paper); color: var(--ink);
  font: 15.5px/1.62 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
        Helvetica, Arial, sans-serif;
}
.page { max-width: 54rem; margin: 0 auto; padding: 3rem 1.5rem 4rem; }

/* --- masthead --- */
.mast { display: grid; grid-template-columns: var(--rail) 1fr; gap: 1.5rem;
        align-items: baseline; padding-bottom: 1.4rem;
        border-bottom: 2px solid var(--ink); }
.mast h1 { grid-column: 1 / -1; font: 600 1.6rem/1.1 ui-monospace,
           SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
           letter-spacing: -0.02em; margin: 0 0 0.7rem; }
/* The intro starts in the content column, on the same line as every tool's
   prose, rather than in the rail. */
.mast p { grid-column: 2 / -1; margin: 0; color: var(--muted);
          max-width: 38rem; min-width: 0; }

/* --- the rail: one row per tool, and room for the next --- */
.tool { display: grid; grid-template-columns: var(--rail) 1fr; gap: 1.5rem;
        padding: 2.1rem 0; border-bottom: 1px solid var(--rule); }
.tool > .label { position: sticky; top: 1.5rem; align-self: start; }
.tool > .label .name {
  display: block; font: 600 0.87rem/1.3 ui-monospace, SFMono-Regular,
  "SF Mono", Menlo, Consolas, monospace; color: var(--ink);
}
.tool > .label .what { display: block; margin-top: 0.2rem;
                       font-size: 0.76rem; color: var(--quiet); }
/* A grid item defaults to min-width:auto, which is a refusal to shrink below
   its widest content.  A wide transcript inside one therefore widens the
   column, the grid, and the page, and the overflow-x below never engages.
   Zero here is what lets the scroll container be the transcript itself. */
.tool > .label, .tool > .body { min-width: 0; }
.tool > .body > :first-child { margin-top: 0; }
.tool > .body > :last-child { margin-bottom: 0; }
.tool.next .label .name { color: var(--quiet); }
.tool.next .body { color: var(--quiet); font-size: 0.95rem; }

p { margin: 0 0 0.95rem; max-width: 38rem; }

/* --- a command and what it printed --- */
.run { margin: 1.1rem 0 1.2rem; min-width: 0; }
/* Nothing here wraps.  A command and its output are both fixed-pitch, and a
   line folded at the window edge reads as two: the histogram stops being a
   histogram, and a wrapped command stops being one you can copy by eye.  Each
   block is bounded by its column and scrolls sideways within it. */
.run pre {
  font: 0.78rem/1.55 ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas,
        monospace;
  margin: 0; padding: 0.7rem 0.9rem;
  white-space: pre; overflow-x: auto; overscroll-behavior-x: contain;
}
.run .cmd { background: var(--panel); color: var(--ink);
            border: 1px solid var(--rule); }
.run .cmd b { font-weight: 600; color: var(--accent); }
.run .out { background: var(--term-bg); color: var(--term-ink);
            border: 1px solid var(--rule); border-top: none; }

figure { margin: 1.2rem 0; }
figure img { display: block; width: 100%; height: auto;
             border: 1px solid var(--rule); }
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

ul.plain { margin: 0.4rem 0 0.95rem; padding-left: 1.05rem; max-width: 38rem; }
ul.plain li { margin-bottom: 0.22rem; }

.stamp { margin: 1.6rem 0 0; font-size: 0.8rem; color: var(--quiet); }

@media (max-width: 720px) {
  .mast, .tool { grid-template-columns: 1fr; gap: 0.5rem; }
  .mast p { grid-column: 1; }
  .tool > .label { position: static; }
  .tool > .label .what { display: inline; margin-left: 0.5rem; }
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
  <p>Lightweight tools for working with sequencing files. Each one does a
  single job on the formats a run already arrives in, and writes something you
  can open or paste. They share an entry point, <code>seqview</code>, and
  nothing else.</p>
</header>

<section class="tool">
  <div class="label">
    <span class="name">install</span>
    <span class="what">and the entry point</span>
  </div>
  <div class="body">
    <p>The core needs nothing beyond the standard library. Aligning reads needs
    pysam and biopython, with <code>minimap2</code> and <code>samtools</code> on
    PATH. The <code>fast</code> extra adds numpy, which about halves the time to
    scan a large FASTQ, and <code>plot</code> adds matplotlib for figure
    output.</p>
    <div class="run">
      <pre class="cmd">$ <b>uv tool install</b> --python 3.13 'seqviewer[cli,fast,plot] @ git+https://github.com/micah-olivas/seqviewer'</pre>
    </div>
    <p>One entry point carries the commands.</p>
    <div class="run">
      <pre class="cmd">$ <b>seqview</b></pre>
      <pre class="out">__SEQVIEW__</pre>
    </div>
  </div>
</section>

<section class="tool">
  <div class="label">
    <span class="name">lengths</span>
    <span class="what">read-length distribution</span>
  </div>
  <div class="body">
    <p>Needs no reference and no aligner, so it runs on a file straight off the
    sequencer. Lengths are tallied as the file is read, which keeps memory
    bounded by the range of lengths rather than by the number of reads: a
    multi-gigabyte FASTQ is scanned in one pass.</p>

    <div class="run">
      <pre class="cmd">$ <b>seqview lengths</b> CRNGS7_1_sample_1.fastq --bins 14</pre>
      <pre class="out">__LENGTHS__</pre>
    </div>

    <div class="reading">
      <p>This run is two populations: 4,386 reads of 128&ndash;372&nbsp;bp, and
      259 that span the full 3.5&nbsp;kb construct. The short ones are 80% of
      the file and 2% of its bases, which is the gap between the mean, 408, and
      the N50, 1,338.</p>
    </div>

    <p>The axis covers the central 99% of reads rather than the full range,
    because a few concatemers otherwise reach the top of it on their own. Reads
    outside the axis are counted in a row of their own, labelled with the
    extreme they reach, and every figure under the histogram covers the whole
    run: the longest read here is 4,700&nbsp;bp, past the last bin.</p>

    <div class="ref">
      <h3>seqview lengths [options] reads</h3>
      __LENGTHS_FLAGS__
    </div>
  </div>
</section>

<section class="tool">
  <div class="label">
    <span class="name">pileup</span>
    <span class="what">reads against a reference</span>
  </div>
  <div class="body">
    <p>Aligns a run with minimap2 and writes two pages. The pileup carries every
    read base by base. The summary reduces the same alignment to one screen, and
    reports whether the clone carries an error where the pileup reports what
    each read says. Each page links to the other.</p>

    <div class="run">
      <pre class="cmd">$ <b>seqview pileup</b> reads.fastq vector.dna out.html --max 300</pre>
      <pre class="out">__PILEUP__</pre>
    </div>

    <figure>
      <img alt="The summary page for a 3,549 bp vector: a per-position
                disagreement track, the annotated reference, and the coverage
                profile."
           src="img/summary.png">
      <figcaption>Three bands, top to bottom: how far the reads disagree with
      the reference at each position, the annotated reference, and read depth.
      Bars are muted below 10%, the rate sequencing error alone produces, and
      coloured at or above it.</figcaption>
    </figure>

    <p>A reference of several kilobases is drawn a few bases to the pixel, so
    each column reports the worst position it spans rather than the mean. One
    position where every read disagrees would otherwise be averaged into its
    quiet neighbours.</p>

    <div class="ref">
      <h3>seqview pileup [options] reads reference out</h3>
      __PILEUP_FLAGS__
    </div>
  </div>
</section>

<section class="tool next">
  <div class="label">
    <span class="name">next</span>
    <span class="what">not written yet</span>
  </div>
  <div class="body">
    <p>A tool belongs here when it does one job over a file a run already
    produces, and finishes without a project directory or a config.
    Candidates:</p>
    <ul class="plain">
      <li>per-barcode yield across a demultiplexed plate</li>
      <li>quality profile along the read</li>
      <li>adapter and primer occupancy</li>
    </ul>
    <p>Sections still to write: reading the summary page, reading the pileup,
    calling thresholds and what they miss, reference and feature formats,
    library use, theming, tests.</p>
  </div>
</section>

<p class="stamp">__STAMP__</p>

</main>
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
            .replace("__PILEUP__", example("pileup"))
            .replace("__LENGTHS_FLAGS__", flags(build_lengths_parser(
                "seqview lengths")))
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
