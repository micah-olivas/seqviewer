"""Read-length distributions rendered as a terminal histogram.

Reads a FASTQ for the length of each record and renders the distribution as one
row per length bin, with bars drawn in eighth-block characters so bar length
resolves to an eighth of a character cell.

Lengths are read as a stream and sequences are discarded as they are counted, so
memory does not scale with file size.  :func:`seqviewer.cli.read_fastq` is not
reused here because it retains a record per read.

The axis spans a central share of the reads rather than the full range.  Reads
outside that range are counted in a row of their own, and are included in every
figure :func:`summarise` reports, so the reported minimum, maximum and N50
describe the whole file rather than the part in view.

Bars run horizontally because the labelled axis is the one that needs width: a
length range and a read count are each several characters.
"""

from __future__ import annotations

import gzip
import math
from dataclasses import dataclass, field
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple

__all__ = [
    "Bin",
    "Binning",
    "DEFAULT_BINS",
    "DEFAULT_BULK",
    "Summary",
    "bin_lengths",
    "distribution",
    "histogram",
    "read_lengths",
    "summarise",
    "summary_lines",
]

#: Left-to-right eighth blocks, so a bar can end part way through a cell.
#: Without them bar length resolves only to a whole cell.
_EIGHTHS = " ▏▎▍▌▋▊▉█"

#: Bins used when the caller does not specify a count.  Fits a histogram, a
#: summary and a prompt in a 30-row terminal.
DEFAULT_BINS = 24

#: Percent of reads the axis covers, centred.  Ultra-long reads such as
#: concatemers otherwise set the upper limit, which compresses the rest of the
#: distribution into the lowest bins.
DEFAULT_BULK = 99.0


@dataclass(frozen=True)
class Bin:
    """One bar: the half-open length range ``[low, high)`` and its read count.

    The last bin of the axis is closed at both ends, so that the longest read in
    view falls inside a bin.
    """

    low: int
    high: int
    count: int


@dataclass(frozen=True)
class Binning:
    """The axis and the reads outside it.

    ``below`` and ``above`` count reads shorter than ``low`` and longer than
    ``high``.  ``shortest`` and ``longest`` are the extremes over all reads, not
    over the axis.
    """

    bins: List[Bin] = field(default_factory=list)
    below: int = 0
    above: int = 0
    low: int = 0
    high: int = 0
    longest: int = 0
    shortest: int = 0

    @property
    def clipped(self) -> bool:
        """Whether any read falls outside the axis."""
        return bool(self.below or self.above)


@dataclass(frozen=True)
class Summary:
    """Read count, total bases, and length statistics over every read."""

    reads: int
    bases: int
    shortest: int
    longest: int
    median: int
    mean: float
    n50: int

    @property
    def empty(self) -> bool:
        return self.reads == 0


def _open(path):
    """Open a FASTQ, plain or gzipped, as text."""
    return (gzip.open if str(path).endswith(".gz") else open)(path, "rt")


def read_lengths(paths: Iterable) -> Iterator[int]:
    """Yield the sequence length of every record in *paths*.

    Reads four-line FASTQ records and retains only the length.  A record whose
    fourth line is absent is dropped rather than counted, which is the state a
    file interrupted mid-write is left in.
    """
    for path in paths:
        with _open(path) as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                seq = handle.readline()
                handle.readline()              # the '+' line
                qual = handle.readline()
                if not qual:
                    break
                if header[0] == "@":
                    yield len(seq.strip())


def _percentile(ordered: Sequence[int], pct: float) -> int:
    """Return the *pct* percentile of a sorted sequence, linearly interpolated."""
    if not ordered:
        return 0
    if pct <= 0:
        return ordered[0]
    if pct >= 100:
        return ordered[-1]
    pos = (len(ordered) - 1) * pct / 100.0
    low, high = math.floor(pos), math.ceil(pos)
    if low == high:
        return ordered[low]
    return int(round(ordered[low] + (ordered[high] - ordered[low]) * (pos - low)))


def summarise(lengths: Sequence[int], presorted: bool = False) -> Summary:
    """Return read count, total bases and length statistics for *lengths*.

    N50 is the length at which reads of that length or longer account for half of
    all bases.  It differs from the mean where a run mixes many short reads with
    a few long ones: the mean sits near the short reads and the N50 near the
    long.  Both are reported so the difference is visible.

    Pass ``presorted=True`` when *lengths* is already ascending, to skip the
    sort.
    """
    if not lengths:
        return Summary(0, 0, 0, 0, 0, 0.0, 0)

    ordered = list(lengths) if presorted else sorted(lengths)
    n = len(ordered)
    total = sum(ordered)

    mid = n // 2
    median = ordered[mid] if n % 2 else (ordered[mid - 1] + ordered[mid]) // 2

    half, run, n50 = total / 2, 0, ordered[-1]
    for length in reversed(ordered):
        run += length
        if run >= half:
            n50 = length
            break

    return Summary(reads=n, bases=total, shortest=ordered[0],
                   longest=ordered[-1], median=median, mean=total / n, n50=n50)


def bin_lengths(
    lengths: Sequence[int],
    count: int = DEFAULT_BINS,
    bulk: float = DEFAULT_BULK,
    presorted: bool = False,
) -> Binning:
    """Bin *lengths* over the central *bulk* percent of reads.

    Bins are equal width, which keeps the length axis linear and bar heights
    comparable between bins.  Equal-population bins would give every bar the same
    height.

    *bulk* is split between the two ends: at 99, the axis runs from the 0.5th to
    the 99.5th percentile.  Reads outside it are counted in ``Binning.below`` and
    ``Binning.above`` and are not binned.  ``bulk=100`` spans the full range.

    Lengths that are all equal give a single bin.  Where clipping would leave no
    reads inside the axis, the full range is used instead.
    """
    if not lengths:
        return Binning()

    ordered = list(lengths) if presorted else sorted(lengths)
    tail = max(0.0, 100.0 - bulk) / 2.0
    low = _percentile(ordered, tail)
    high = _percentile(ordered, 100.0 - tail)

    below = sum(1 for v in ordered if v < low)
    above = sum(1 for v in ordered if v > high)
    inside = [v for v in ordered if low <= v <= high]
    if not inside:
        inside, below, above = ordered, 0, 0
        low, high = ordered[0], ordered[-1]

    if low == high:
        bins = [Bin(low, low + 1, len(inside))]
    else:
        count = max(1, count)
        width = max(1, math.ceil((high - low + 1) / count))
        edges = list(range(low, high + 1, width))
        counts = [0] * len(edges)
        for value in inside:
            counts[min((value - low) // width, len(edges) - 1)] += 1
        bins = [Bin(e, min(e + width, high + 1), counts[i])
                for i, e in enumerate(edges)]

    return Binning(bins=bins, below=below, above=above, low=low, high=high,
                   longest=ordered[-1], shortest=ordered[0])


def distribution(
    lengths: Sequence[int],
    count: int = DEFAULT_BINS,
    bulk: float = DEFAULT_BULK,
) -> Tuple[Summary, Binning]:
    """Return the summary and the binning for *lengths*, sorting once."""
    ordered = sorted(lengths)
    return (summarise(ordered, presorted=True),
            bin_lengths(ordered, count, bulk, presorted=True))


def _bar(eighths: int) -> str:
    """Return a bar of *eighths* eighth-cells."""
    if eighths <= 0:
        return ""
    full, part = divmod(eighths, 8)
    return "█" * full + (_EIGHTHS[part] if part else "")


def histogram(binning: Binning, width: int = 80, log: bool = False) -> List[str]:
    """Render *binning* as lines of text, none wider than *width*.

    Clipped tails are drawn as rows above and below the axis, labelled with
    their read count and the extreme length they reach.  Tail counts do not
    scale the bars; the tallest binned count does.

    *log* scales bar length by ``log(1 + count)``, which keeps the smaller bins
    of a peaked distribution distinguishable.  Axis labels stay linear under
    either scale.
    """
    if not binning.bins:
        return ["no reads"]

    rows = []
    if binning.below:
        rows.append((f"<{binning.low:,}", binning.below,
                     f"shorter, down to {binning.shortest:,}"))
    for b in binning.bins:
        rows.append((f"{b.low:,}–{b.high - 1:,}", b.count, ""))
    if binning.above:
        rows.append((f">{binning.high:,}", binning.above,
                     f"longer, up to {binning.longest:,}"))

    counts = [f"{r[1]:,}" for r in rows]
    notes = [r[2] for r in rows]
    label_w = max(len(r[0]) for r in rows + [("bp", 0, "")])
    count_w = max(len(x) for x in counts + ["reads"])
    note_w = max(len(x) for x in notes)

    bar_w = max(8, width - label_w - count_w - note_w - (6 if note_w else 4))

    def height(count: int) -> float:
        return math.log1p(count) if log else float(count)

    ceiling = max((height(b.count) for b in binning.bins), default=1.0) or 1.0

    lines = [f"{'bp':>{label_w}}  {'':{bar_w}}  {'reads':>{count_w}}"]
    for (label, count, note), shown in zip(rows, counts):
        eighths = int(round(bar_w * 8 * min(1.0, height(count) / ceiling)))
        if count and eighths < 1:
            eighths = 1                        # a row with reads is never blank
        line = f"{label:>{label_w}}  {_bar(eighths):{bar_w}}  {shown:>{count_w}}"
        if note:
            line += f"  {note}"
        lines.append(line)
    return lines


def summary_lines(summary: Summary,
                  binning: Optional[Binning] = None) -> List[str]:
    """Return the distribution's figures as lines of text.

    The figures cover every read.  When *binning* is clipped, a third line gives
    the axis range and how many reads it holds, so a maximum past the last bin
    reads as the file's own rather than as an error.
    """
    if summary.empty:
        return ["no reads"]
    lines = [
        f"{summary.reads:,} reads · {summary.bases:,} bases",
        f"min {summary.shortest:,} · median {summary.median:,} · "
        f"mean {summary.mean:,.0f} · max {summary.longest:,} · "
        f"N50 {summary.n50:,}",
    ]
    if binning is not None and binning.clipped:
        inside = summary.reads - binning.below - binning.above
        lines.append(
            f"axis {binning.low:,}–{binning.high:,} bp, holding {inside:,} of "
            f"{summary.reads:,} reads ({100 * inside / summary.reads:.1f}%). "
            f"Pass --bulk 100 for the full range."
        )
    return lines
