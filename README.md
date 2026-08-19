# seqviewer

Viewers for sequencing alignments.

## Install

The core is dependency-free. Reading alignments needs pysam; aligning reads also
needs `minimap2` and `samtools` on PATH; reading an annotated reference needs
biopython.

```bash
pip install 'seqviewer @ git+https://github.com/micah-olivas/seqviewer'
pip install 'seqviewer[cli] @ git+https://github.com/micah-olivas/seqviewer'
```

Not on PyPI. From a checkout, `pip install -e '.[cli]'` installs the same extras
editable. With uv, no install step is needed: `uv run` builds the project into
its own environment first, so the commands below work in a fresh clone.

## Input levels

Input comes in at three levels. A caller pays only for the level it uses.

| Level | Entry point | Requires |
| --- | --- | --- |
| A grid you already have | `render` | nothing |
| A SAM or BAM | `grid_from_alignment` | pysam |
| Raw reads | `grid_from_reads` | pysam, minimap2, samtools |

The alignment helpers are imported lazily, so importing `seqviewer` does
not require pysam.

## Use

A grid is a list of rows, one per read, each as wide as the reference. Every cell
is `(base, is_match)`, where a base of `"-"` marks a gap or an uncovered
position.

```python
from pathlib import Path
from seqviewer import PileupGroup, PileupView, render

group = PileupGroup(
    name="pUC19-WT",
    ref_seq=ref_seq,
    rows=rows,
    n_reads=65,
    fraction=0.65,
    status="Perfect Match",
    highlighted=True,
)
view = PileupView(title="Well A1", groups=[group], total_reads=100,
                  flanks=(100, 100))
Path("pileup.html").write_text(render(view))
```

Starting from an alignment instead:

```python
from seqviewer import grid_from_alignment

rows = grid_from_alignment("well_A1.bam", ref_seq, read_names={"read_1", "read_2"})
```

`grid_from_alignment` drops reads that do not cross the reference midpoint.
Concatemer fragments cover one flank and stop, so this filter removes them while
keeping full-length reads. Pass `min_overlap_pos=0` to disable it, or an explicit
position to move it.

To see both pages without reads or an aligner:

```bash
python -m seqviewer.demo demo.html
python -m seqviewer.demo demo.html --summary
```

## Summarized view

`SummaryView.from_view` reduces a `PileupView` to per-position depth, agreement,
and deletions, plus a list of called variants classified against the reading
frame. `render_summary` draws that.

```python
from pathlib import Path
from seqviewer import SummaryView, render_summary

Path("summary.html").write_text(render_summary(SummaryView.from_view(view)))
```

Glyph shape carries the kind of change — substitution, deletion, insertion — and
color carries the consequence: frameshift or premature stop, missense or in-frame
indel, silent, or outside the reading frame. Stem height is the allele fraction.
Each called variant also gets a base-resolution window, drawn as one letter per
base with each codon bracketed under the three bases it is translated from.

A variant is called when at least 25% of the reads covering its position carry it
and at least two of them do. `min_fraction` and `min_count` move both floors. The
read floor matters at the depths these pages are made for: with ten reads, one
read is 10%, so a fraction alone admits the per-base error rate as an allele.

Substitutions and deletions are recovered from the grid. Insertions are not. An
inserted base has no reference position, so `align` drops it and a row stays
exactly as wide as the reference, leaving nowhere to hold one. `Variant` models
`kind="ins"` and `summarize_group` accepts insertion evidence passed alongside
the grid, so a caller holding it can supply it; a grid alone reports none.
`tests/test_cli.py` pins that against minimap2 output.

The two pages scale differently. A pileup holds one row per read, so its size
follows the run's depth; a summary holds one band per group and one window per
called variant, so its size follows what was called rather than how deep the run
was.

`summary.flagged_columns` is the per-column disagreement statistic as a pure
function of a grid and a reference, for callers that want the flagged positions
without a page.

## Command line

`seqviewer-pileup` aligns reads to a reference and writes a page. It takes a
directory of FASTQs, which is how a sequencing run arrives:

```bash
uv run --extra cli seqviewer-pileup reads/ reference.dna pileup.html
```

The directory's files are pooled into one pileup; `--per-file` draws a group per
file instead, so samples can be compared on one page. A single FASTQ still works
in place of the directory, and `.gz` is read directly. The reference is FASTA,
GenBank, ApE, or SnapGene.

`--extra cli` covers what the command needs beyond the core: pysam to read the
alignment and biopython to read an annotated reference. `minimap2` and
`samtools` still have to be on PATH. Without uv, run the same thing as
`python -m seqviewer.cli`.

That form has to be run from a checkout. To get the command on PATH and run it
from wherever the reads are:

```bash
uv tool install --python 3.13 'seqviewer[cli] @ git+https://github.com/micah-olivas/seqviewer'
seqviewer-pileup reads/ reference.dna pileup.html      # from anywhere
```

Pass `--editable '/path/to/seqviewer[cli]'` instead to track a local checkout, so
a `git pull` is the whole upgrade. The Python has to be named: `uv tool install`
otherwise takes the default interpreter, and a conda base of 3.8 does not satisfy
this package's `requires-python`. `uv tool uninstall seqviewer` undoes it.

A page draws 500 reads by default, sampled uniformly from across the whole of
every file. A row costs roughly 2 KB of HTML, so a 35,000-read pool drawn whole
is a 70 MB page — too large to open, and no more readable for holding every
read. `--max` moves the number and `--max 0` draws all of them. The sample is
seeded, so the same directory gives the same page twice.

`--summary` writes the summarized page beside the pileup, named from its stem:
`pileup.html` and `pileup.summary.html`. `--variant-freq` and `--variant-reads`
move the two calling floors. The log written beside the page records which
thresholds produced it and what was called.

Worth knowing: `--insert LABEL` marks a feature as the focus region, which is
what draws the boundary lines, the translation rows, and the frame that variants
are classified against. `--order` sets the row order; `cluster` is hierarchical,
average linkage over what each read disagrees about, which groups a
subpopulation that the cheaper `mismatch` ordering splits when a read carries an
unrelated error further left. `--help` lists the rest.

A track above the reference shows, per position, the share of covering reads
that disagree with it, on a log scale marked at 1% and 10%. A deletion counts as
disagreement, so a column half the reads have deleted reads as half disagreeing
rather than as clean. `seqviewer.summary.mismatch_fractions` is the one
definition of that number, and both pages read it.

## Reference and Feature

`Reference` and `Feature` describe a construct and are the substrate the
package's viewers share. `Feature` records a half-open `[start, end)` span, a
strand, and whether the feature wraps the origin of a circular construct — a
wrapped feature's real extent is not `[start, end)`, so callers refuse it rather
than read it as one stretch.

The pileup's vector/insert split is exposed as `PileupView.flanks`, a
`(5' length, 3' length)` pair. `Reference.flank_lengths()` derives that pair from
an annotated insert, and `PileupView.from_reference()` applies it:

```python
from seqviewer import Feature, PileupView, Reference

reference = Reference(seq=seq, name="pUC19", topology="circular",
                      features=[Feature("insert", 100, 900, strand=1)])
view = PileupView.from_reference("Well A1", reference, groups)
assert view.flanks == (100, 100)
```

Feature glyphs are laid out by `seqviewer.annotate` and drawn as inline SVG.
Overlapping features pack into lanes; a feature that will not fit within the lane
limit is dropped and returned, so a page can say what it left out.

## Theming

A rendered page carries its own light and dark palettes and reads a stored
preference on load. An application that already stores one passes its own names
so the pages follow the same setting as the rest of its output:

```python
from seqviewer import Theme

theme = Theme(storage_key="app-theme", css_prefix="app",
              style_id="app-bridge", script_id="app-sync")
```

## Translation

`seqviewer.codon` implements the standard genetic code (NCBI table 1)
without Biopython. Ambiguous codons resolve to a single residue when every base
they could stand for gives the same one, so `ACN` is threonine; gaps are written
as `N` before translation, which makes such codons routine. `tests/test_codon.py`
checks all 4,096 codons over the full IUPAC alphabet against
`Bio.Seq.Seq.translate()`.

## Out of scope

These are excluded so the package does not become a plasmid editor:

- sequence editing
- restriction-enzyme and REBASE site calculation
- primer design
- ORF finding
- implementing aligners — `grid_from_reads` invokes minimap2 and reads the result
- chromatogram viewing
- plate and well-layout maps

## Tests

```bash
pip install -e '.[dev,sam]'
pytest
```

Alignment tests build a SAM by hand and run without minimap2 or samtools; the
tests that need both are skipped when they are absent.

## License

MIT. See `LICENSE`.
