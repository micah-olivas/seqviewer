"""Read-length distributions rendered as a terminal histogram.

Reads a FASTQ for the length of each record and renders the distribution as one
row per length bin, with bars drawn in eighth-block characters so bar length
resolves to an eighth of a character cell.

Lengths are tallied as they are read, into one count per distinct length rather
than one value per read.  A run's lengths are integers bounded by the longest
read, so the tally is bounded by the length range and not by the number of
reads, and a multi-gigabyte file is scanned in one pass holding nothing but the
tally.  Every figure reported here is counted off the tally: quantiles, the
median and N50 are exact rather than estimated, and no read is sampled away.

Records are found by their line boundaries and their contents are not decoded.
The length of a sequence line is its count of bytes, which is its count of bases
for the ASCII a FASTQ holds.  Two scanners find those boundaries, described at
:func:`_tally_python` and :func:`_tally_numpy`; the second is used where numpy is
installed and is faster on a large file.  They return the same tally, which the
tests check.

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
import os
from bisect import bisect_right
from collections import Counter
from dataclasses import dataclass, field
from typing import (Callable, Dict, Iterable, Iterator, List, Optional,
                    Sequence, Tuple)

try:                                    # an accelerator, not a requirement
    import numpy as _np
except ImportError:                     # pragma: no cover - numpy is usual
    _np = None

#: Whether the array scanner is available.  It bears on how long a scan takes
#: and on nothing it reports.
HAVE_NUMPY = _np is not None

__all__ = [
    "Bin",
    "Binning",
    "DEFAULT_BINS",
    "DEFAULT_BULK",
    "HAVE_NUMPY",
    "LengthCounts",
    "PALETTE",
    "Palette",
    "Summary",
    "bin_counts",
    "bin_lengths",
    "count_lengths",
    "distribution",
    "histogram",
    "read_lengths",
    "summarise",
    "summarise_counts",
    "summary_lines",
]

#: Left-to-right eighth blocks, so a bar can end part way through a cell.
#: Without them bar length resolves only to a whole cell.
_EIGHTHS = " ▏▎▍▌▋▊▉█"

#: Bytes read per block.  Large enough that the per-block work in Python is
#: spread over thousands of records, small enough to stay in cache.
_BLOCK = 1 << 22

#: Lengths the array tally starts with room for.  It grows to fit a longer read.
_TALLY_START = 4096

#: Bins used when the caller does not specify a count.  Fits a histogram, a
#: summary and a prompt in a 30-row terminal.
DEFAULT_BINS = 24

#: Percent of reads the axis covers, centred.  Ultra-long reads such as
#: concatemers otherwise set the upper limit, which compresses the rest of the
#: distribution into the lowest bins.
DEFAULT_BULK = 99.0


@dataclass(frozen=True)
class Palette:
    """ANSI codes for the histogram.  The default renders as plain text."""

    bar: str = ""
    peak: str = ""
    tail: str = ""
    head: str = ""
    reset: str = ""


#: Colours for a terminal that reports one.  The peak bin is brightened and the
#: clipped tails are dimmed, so the bulk of the distribution is what the eye
#: reaches first.
PALETTE = Palette(
    bar="\x1b[38;5;37m",
    peak="\x1b[1;38;5;44m",
    tail="\x1b[2m",
    head="\x1b[2m",
    reset="\x1b[0m",
)


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


class LengthCounts:
    """Reads tallied by length.

    Holds one count per distinct length, so its size follows the run's length
    range rather than its read count.  Statistics are counted off it and are
    exact.  The sorted lengths and their running read totals are built on first
    use and reused until the tally changes.
    """

    __slots__ = ("_counts", "_keys", "_cum", "_bases")

    def __init__(self, counts: Optional[Dict[int, int]] = None):
        self._counts: Counter = Counter()
        if counts:
            self._counts.update(counts)
        self._invalidate()

    @classmethod
    def from_lengths(cls, lengths: Iterable[int]) -> "LengthCounts":
        """Tally an iterable of lengths."""
        counts = cls()
        counts.update(lengths)
        return counts

    def _invalidate(self) -> None:
        self._keys = None
        self._cum = None
        self._bases = None

    def add(self, length: int, times: int = 1) -> None:
        """Tally *times* reads of *length*."""
        self._counts[length] += times
        self._invalidate()

    def update(self, lengths: Iterable[int]) -> None:
        """Tally every length in *lengths*."""
        self._counts.update(lengths)
        self._invalidate()

    def merge(self, tally) -> None:
        """Add a tally from a scanner: a mapping, or an array indexed by length."""
        if _np is not None and isinstance(tally, _np.ndarray):
            found = _np.flatnonzero(tally)
            tally = dict(zip(found.tolist(), tally[found].tolist()))
        self._counts.update(tally)
        self._invalidate()

    def _prepare(self) -> Tuple[List[int], List[int]]:
        """Return the sorted lengths and the running read total at each."""
        if self._keys is None:
            keys = sorted(self._counts)
            cum, run = [], 0
            for key in keys:
                run += self._counts[key]
                cum.append(run)
            self._keys, self._cum = keys, cum
        return self._keys, self._cum

    def items(self) -> List[Tuple[int, int]]:
        """Return ``(length, count)`` pairs in ascending length."""
        keys, _ = self._prepare()
        return [(k, self._counts[k]) for k in keys]

    @property
    def empty(self) -> bool:
        return not self._counts

    @property
    def total(self) -> int:
        """Reads tallied."""
        _, cum = self._prepare()
        return cum[-1] if cum else 0

    @property
    def bases(self) -> int:
        """Bases tallied."""
        if self._bases is None:
            self._bases = sum(k * c for k, c in self._counts.items())
        return self._bases

    @property
    def shortest(self) -> int:
        keys, _ = self._prepare()
        return keys[0] if keys else 0

    @property
    def longest(self) -> int:
        keys, _ = self._prepare()
        return keys[-1] if keys else 0

    def value_at(self, index: int) -> int:
        """Return the length of the *index*-th shortest read, counting from 0."""
        keys, cum = self._prepare()
        if not keys:
            return 0
        index = min(max(index, 0), cum[-1] - 1)
        return keys[bisect_right(cum, index)]

    def quantile(self, pct: float) -> int:
        """Return the *pct* percentile, interpolated between adjacent reads."""
        total = self.total
        if not total:
            return 0
        if pct <= 0:
            return self.shortest
        if pct >= 100:
            return self.longest
        pos = (total - 1) * pct / 100.0
        low, high = math.floor(pos), math.ceil(pos)
        below = self.value_at(low)
        if low == high:
            return below
        above = self.value_at(high)
        return int(round(below + (above - below) * (pos - low)))

    @property
    def median(self) -> int:
        total = self.total
        if not total:
            return 0
        mid = total // 2
        if total % 2:
            return self.value_at(mid)
        return (self.value_at(mid - 1) + self.value_at(mid)) // 2

    @property
    def mean(self) -> float:
        total = self.total
        return self.bases / total if total else 0.0

    @property
    def n50(self) -> int:
        """The length at which reads that long or longer hold half the bases."""
        keys, _ = self._prepare()
        if not keys:
            return 0
        half, run = self.bases / 2, 0
        for key in reversed(keys):
            run += key * self._counts[key]
            if run >= half:
                return key
        return keys[0]

    def between(self, low: int, high: int) -> Tuple[int, int, int]:
        """Return the reads below *low*, within ``[low, high]``, and above it."""
        below = above = inside = 0
        for key, count in self._counts.items():
            if key < low:
                below += count
            elif key > high:
                above += count
            else:
                inside += count
        return below, inside, above


def _open_pair(path):
    """Return the file handle and the byte stream to read records from.

    The handle is returned alongside the stream so a caller tracking progress
    can read a position from the file even when the stream decompresses.
    """
    handle = open(path, "rb")
    if str(path).endswith(".gz"):
        return handle, gzip.GzipFile(fileobj=handle, mode="rb")
    return handle, handle


def _length_blocks(stream) -> Iterator[List[int]]:
    """Yield the record lengths in *stream*, one list per block read.

    Each block is split on line boundaries in one call and every fourth line
    taken, so the work per record is done in C rather than in a Python loop.  A
    block is trimmed back to a whole number of records and the remainder carried
    forward, so a record split across two blocks is counted once, when its
    fourth line arrives.  A final record whose fourth line is absent is dropped,
    which is the state a file interrupted mid-write is left in.

    One list is yielded per block, empty where a block completed no record, so a
    caller reporting progress hears from every block.
    """
    buf = b""
    sep = b"\n"
    first = True
    while True:
        data = stream.read(_BLOCK)
        if not data:
            break
        if first:
            first = False
            if b"\r\n" in data:
                sep = b"\r\n"           # split off the CR rather than count it
        buf += data
        lines = buf.split(sep)
        partial = lines.pop()           # the line the block ended part way into
        whole = len(lines) - len(lines) % 4
        leftover = lines[whole:]        # complete lines of an unfinished record
        buf = sep.join(leftover) + sep + partial if leftover else partial
        yield list(map(len, lines[1:whole:4])) if whole else []

    if buf:
        lines = buf.split(sep)
        if not lines[-1]:
            lines.pop()                 # a trailing separator, not a line
        whole = len(lines) - len(lines) % 4
        if whole:
            yield list(map(len, lines[1:whole:4]))


def _tally_python(
    stream,
    on_block: Optional[Callable[[int], None]] = None,
) -> Tuple[Counter, int]:
    """Return the record lengths in *stream* as a counter, and the read count."""
    tally: Counter = Counter()
    reads = 0
    for block in _length_blocks(stream):
        if block:
            tally.update(block)         # counted in C, one call per block
            reads += len(block)
        if on_block is not None:
            on_block(reads)
    return tally, reads


def _grow(tally, top: int):
    """Return *tally* with room for a read of length *top*."""
    if top < tally.size:
        return tally
    grown = _np.zeros(top + 1, dtype=_np.int64)
    grown[:tally.size] = tally
    return grown


def _tally_numpy(
    stream,
    on_block: Optional[Callable[[int], None]] = None,
):
    """Return the record lengths in *stream* as an array indexed by length.

    The newlines in a block are located in one pass and the gaps between them
    are the line lengths, which keeps the work per record inside the array layer.
    Blocks are trimmed and carried the same way as in :func:`_length_blocks`, and
    a final record whose fourth line is absent is dropped the same way.
    """
    tally = _np.zeros(_TALLY_START, dtype=_np.int64)
    carry = b""
    crlf = False
    first = True
    reads = 0
    drop = 1                            # the newline itself

    while True:
        data = stream.read(_BLOCK)
        if not data:
            break
        if first:
            first = False
            crlf = b"\r\n" in data
            drop = 2 if crlf else 1     # the CR as well
        if carry:
            data = carry + data
        marks = _np.flatnonzero(_np.frombuffer(data, dtype=_np.uint8) == 10)
        whole = len(marks) - len(marks) % 4
        if whole < 4:
            carry = data                # a record longer than one block
        else:
            spans = _np.diff(marks[:whole], prepend=-1) - drop
            seqs = _np.maximum(spans[1::4], 0)
            if seqs.size:
                tally = _grow(tally, int(seqs.max()))
                tally += _np.bincount(seqs, minlength=tally.size)
                reads += int(seqs.size)
            carry = data[marks[whole - 1] + 1:]
        if on_block is not None:
            on_block(reads)

    if carry:
        marks = _np.flatnonzero(_np.frombuffer(carry, dtype=_np.uint8) == 10)
        # Three newlines and bytes after the last: a final record whose fourth
        # line arrived without one.  A trailing newline instead means the file
        # stops mid-record, and the record is dropped.
        if len(marks) == 3 and marks[-1] != len(carry) - 1:
            spans = _np.diff(marks, prepend=-1) - drop
            length = max(int(spans[1]), 0)
            tally = _grow(tally, length)
            tally[length] += 1
            reads += 1

    return tally, reads


def read_lengths(paths: Iterable) -> Iterator[int]:
    """Yield the sequence length of every complete record in *paths*.

    A record is four lines and its header is not inspected, so a file whose
    records are not in fours is miscounted rather than resynchronised.
    """
    for path in paths:
        handle, stream = _open_pair(path)
        try:
            for block in _length_blocks(stream):
                yield from block
        finally:
            if stream is not handle:
                stream.close()
            handle.close()


def count_lengths(
    paths: Iterable,
    progress: Optional[Callable[[int, int, int], None]] = None,
    fast: Optional[bool] = None,
) -> LengthCounts:
    """Tally the record lengths in *paths* in one pass.

    Nothing is held but the tally, so the cost in memory follows the run's
    length range rather than its size on disk.

    *progress* is called with the bytes read, the bytes to read and the reads
    counted so far, after each block and once at the end of each file.  For a
    gzipped file the figures are compressed bytes, which is what the file's size
    reports.

    *fast* chooses the scanner: the array one where None and numpy is installed,
    and the pure-Python one where False.
    """
    if fast is None:
        fast = HAVE_NUMPY
    elif fast and not HAVE_NUMPY:
        raise RuntimeError("the array scanner needs numpy installed")
    scan = _tally_numpy if fast else _tally_python

    paths = list(paths)
    sizes = []
    for path in paths:
        try:
            sizes.append(os.path.getsize(path))
        except OSError:
            sizes.append(0)
    total = sum(sizes)

    counts = LengthCounts()
    done_bytes = 0
    done_reads = 0

    for path, size in zip(paths, sizes):
        handle, stream = _open_pair(path)
        report = None
        if progress is not None:
            def report(reads, _handle=handle, _bytes=done_bytes,
                       _reads=done_reads):
                progress(_bytes + _handle.tell(), total, _reads + reads)
        try:
            tally, reads = scan(stream, report)
        finally:
            if stream is not handle:
                stream.close()
            handle.close()
        counts.merge(tally)
        done_bytes += size
        done_reads += reads
        if progress is not None:
            progress(done_bytes, total, done_reads)

    return counts


def summarise_counts(counts: LengthCounts) -> Summary:
    """Return read count, total bases and length statistics for a tally.

    N50 is the length at which reads of that length or longer account for half of
    all bases.  It differs from the mean where a run mixes many short reads with
    a few long ones: the mean sits near the short reads and the N50 near the
    long.  Both are reported so the difference is visible.
    """
    if counts.empty:
        return Summary(0, 0, 0, 0, 0, 0.0, 0)
    return Summary(reads=counts.total, bases=counts.bases,
                   shortest=counts.shortest, longest=counts.longest,
                   median=counts.median, mean=counts.mean, n50=counts.n50)


def bin_counts(
    counts: LengthCounts,
    count: int = DEFAULT_BINS,
    bulk: float = DEFAULT_BULK,
) -> Binning:
    """Bin a tally over the central *bulk* percent of reads.

    Bins are equal width, which keeps the length axis linear and bar heights
    comparable between bins.  Equal-population bins would give every bar the same
    height.

    *bulk* is split between the two ends: at 99, the axis runs from the 0.5th to
    the 99.5th percentile.  Reads outside it are counted in ``Binning.below`` and
    ``Binning.above`` and are not binned.  ``bulk=100`` spans the full range.

    Lengths that are all equal give a single bin.  Where clipping would leave no
    reads inside the axis, the full range is used instead.
    """
    if counts.empty:
        return Binning()

    tail = max(0.0, 100.0 - bulk) / 2.0
    low = counts.quantile(tail)
    high = counts.quantile(100.0 - tail)
    below, inside, above = counts.between(low, high)
    if not inside:
        low, high = counts.shortest, counts.longest
        below, inside, above = 0, counts.total, 0

    if low == high:
        bins = [Bin(low, low + 1, inside)]
    else:
        count = max(1, count)
        width = max(1, math.ceil((high - low + 1) / count))
        edges = list(range(low, high + 1, width))
        tallies = [0] * len(edges)
        last = len(edges) - 1
        for length, n in counts.items():
            if low <= length <= high:
                tallies[min((length - low) // width, last)] += n
        bins = [Bin(e, min(e + width, high + 1), tallies[i])
                for i, e in enumerate(edges)]

    return Binning(bins=bins, below=below, above=above, low=low, high=high,
                   longest=counts.longest, shortest=counts.shortest)


def summarise(lengths: Sequence[int]) -> Summary:
    """Return read count, total bases and length statistics for *lengths*.

    See :func:`summarise_counts`, which this tallies for.
    """
    return summarise_counts(LengthCounts.from_lengths(lengths))


def bin_lengths(
    lengths: Sequence[int],
    count: int = DEFAULT_BINS,
    bulk: float = DEFAULT_BULK,
) -> Binning:
    """Bin *lengths* over the central *bulk* percent of reads.

    See :func:`bin_counts`, which this tallies for.
    """
    return bin_counts(LengthCounts.from_lengths(lengths), count, bulk)


def distribution(
    lengths: Sequence[int],
    count: int = DEFAULT_BINS,
    bulk: float = DEFAULT_BULK,
) -> Tuple[Summary, Binning]:
    """Return the summary and the binning for *lengths*, tallying once."""
    counts = LengthCounts.from_lengths(lengths)
    return summarise_counts(counts), bin_counts(counts, count, bulk)


def _bar(eighths: int) -> str:
    """Return a bar of *eighths* eighth-cells."""
    if eighths <= 0:
        return ""
    full, part = divmod(eighths, 8)
    return "█" * full + (_EIGHTHS[part] if part else "")


def histogram(
    binning: Binning,
    width: int = 80,
    log: bool = False,
    palette: Optional[Palette] = None,
) -> List[str]:
    """Render *binning* as lines of text, none wider than *width*.

    Clipped tails are drawn as rows above and below the axis, labelled with
    their read count and the extreme length they reach.  Tail counts do not
    scale the bars; the tallest binned count does.

    *log* scales bar length by ``log(1 + count)``, which keeps the smaller bins
    of a peaked distribution distinguishable.  Axis labels stay linear under
    either scale.

    *palette* colours the output.  Its codes are added after the columns are laid
    out, so they neither shift the alignment nor count toward *width*.
    """
    if not binning.bins:
        return ["no reads"]

    pal = palette or Palette()
    rows: List[Tuple[str, int, str, bool]] = []
    if binning.below:
        rows.append((f"<{binning.low:,}", binning.below,
                     f"shorter, down to {binning.shortest:,}", True))
    for b in binning.bins:
        rows.append((f"{b.low:,}–{b.high - 1:,}", b.count, "", False))
    if binning.above:
        rows.append((f">{binning.high:,}", binning.above,
                     f"longer, up to {binning.longest:,}", True))

    counts = [f"{r[1]:,}" for r in rows]
    label_w = max(len(r[0]) for r in rows + [("bp", 0, "", False)])
    count_w = max(len(x) for x in counts + ["reads"])
    note_w = max(len(r[2]) for r in rows)
    bar_w = max(8, width - label_w - count_w - note_w - (6 if note_w else 4))

    def height(count: int) -> float:
        return math.log1p(count) if log else float(count)

    ceiling = max((height(b.count) for b in binning.bins), default=1.0) or 1.0
    peak = max(b.count for b in binning.bins)

    head = f"{'bp':>{label_w}}  {'':{bar_w}}  {'reads':>{count_w}}"
    lines = [f"{pal.head}{head}{pal.reset}" if pal.head else head]

    for (label, count, note, tail), shown in zip(rows, counts):
        eighths = int(round(bar_w * 8 * min(1.0, height(count) / ceiling)))
        if count and eighths < 1:
            eighths = 1                        # a row with reads is never blank
        bar = f"{_bar(eighths):{bar_w}}"
        if tail:
            body = f"{label:>{label_w}}  {bar}  {shown:>{count_w}}  {note}"
            lines.append(f"{pal.tail}{body}{pal.reset}" if pal.tail else body)
            continue
        colour = pal.peak if count == peak else pal.bar
        if colour:
            bar = f"{colour}{bar}{pal.reset}"
        lines.append(f"{label:>{label_w}}  {bar}  {shown:>{count_w}}")
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
