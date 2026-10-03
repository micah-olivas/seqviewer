"""Base-quality distributions rendered as a terminal histogram.

Reads a FASTQ for Phred scores and renders the distribution as one row per score
bin, with bars drawn in eighth-block characters, the same way
:mod:`seqviewer.lengths` draws read lengths.  Two things can be scored:

``base``
    Every base is one count, at its own score.  A read that is good for most of
    its length and poor at the end still shows its poor bases.

``read``
    Every read is one count, at a mean of its bases' scores rounded to the
    nearest whole score.  This is the view that says how many reads are usable,
    and it hides where in a read the poor bases are.  The mean is taken over
    error probabilities by default, as nanopore tools such as NanoPlot take it:
    a score is ``-10 log10`` of an error rate, so a read's mean score is that of
    its mean error rate, and a few bad bases pull it down much further than a
    plain average of the scores would.  The plain average is available too.

Either way scores are tallied as they are read, into one count per distinct
score.  A Phred score is a single byte, so the tally holds at most a few dozen
counters whatever the size of the run, and a multi-gigabyte file is scanned in
one pass holding nothing but the tally.  Every figure reported here is counted
off the tally, so the median and the shares above a threshold are exact and
nothing is sampled away.

Records are found by their line boundaries, as in :mod:`seqviewer.lengths`, and
only the fourth line of each is looked at.  Where numpy is installed the lines
are found and summed as arrays, which is faster on a large file; the pure-Python
scanners return the same tallies, which the tests check.
"""

from __future__ import annotations

import math
import os
from collections import Counter
from dataclasses import dataclass
from typing import Callable, Iterable, Iterator, List, Optional, Tuple

from .lengths import HAVE_NUMPY, PALETTE, Palette, _BLOCK, _bar, _open_pair

try:                                    # an accelerator, not a requirement
    import numpy as _np
except ImportError:                     # pragma: no cover - numpy is usual
    _np = None


__all__ = [
    "BY",
    "Bin",
    "Binning",
    "DEFAULT_BINS",
    "DEFAULT_OFFSET",
    "MEANS",
    "PALETTE",
    "QualityCounts",
    "Summary",
    "bin_counts",
    "count_qualities",
    "histogram",
    "summarise",
    "summary_lines",
]

#: What a count can be of: every base, or every read at its mean score.
BY = ("base", "read")

#: How a read's mean is taken: over error probabilities, or over the scores.
MEANS = ("error", "phred")

#: Bins used when the caller does not specify a count.  Illumina scores run to
#: about 41, so this gives bins two scores wide.
DEFAULT_BINS = 21

#: The ASCII value of score 0.  Sanger and current Illumina FASTQ use 33; the
#: Illumina 1.3-1.7 files use 64.
DEFAULT_OFFSET = 33

#: The thresholds the summary reports the share at or above.
THRESHOLDS = (20, 30)


def _unit(by: str) -> str:
    return "bases" if by == "base" else "reads"


@dataclass(frozen=True)
class Bin:
    """A half-open range of scores, ``low <= q < high``, and its count."""

    low: int
    high: int
    count: int


@dataclass(frozen=True)
class Binning:
    """The tally cut into equal-width bins from score 0 to the highest seen."""

    bins: List[Bin]
    width: int
    unit: str = "bases"

    @property
    def empty(self) -> bool:
        return not self.bins


@dataclass(frozen=True)
class Summary:
    """Figures over every base or read counted, as the run was scored."""

    by: str
    reads: int
    bases: int
    total: int
    lowest: int
    highest: int
    median: int
    mean: float
    shares: Tuple[Tuple[int, float], ...]
    skipped: int
    read_mean: Optional[str] = None

    @property
    def unit(self) -> str:
        return _unit(self.by)

    @property
    def empty(self) -> bool:
        return not self.total


class QualityCounts:
    """How many bases, or reads, were seen at each Phred score.

    Built from a scanner's tally.  For ``by="base"`` the tally is indexed by byte
    value and the offset is taken off each; bytes below it cannot be scores, and
    are counted in :attr:`skipped` rather than dropped without a word, since they
    usually mean the wrong offset was given.  For ``by="read"`` the scanner has
    already reduced each read to a score, and :attr:`skipped` counts the reads
    left out for holding such a byte.
    """

    def __init__(self, counts: Optional[dict] = None, reads: int = 0,
                 skipped: int = 0, by: str = "base",
                 bases: Optional[int] = None,
                 read_mean: Optional[str] = None):
        if by not in BY:
            raise ValueError(f"by must be one of {', '.join(BY)}")
        self._counts: Counter = Counter()
        if counts:
            self._counts.update(counts)
        self.by = by
        self.reads = reads
        self.skipped = skipped
        self._bases = bases
        self.read_mean = read_mean

    @classmethod
    def from_tally(cls, tally, reads: int, offset: int = DEFAULT_OFFSET):
        """Build a per-base count from a mapping or array indexed by byte."""
        if _np is not None and isinstance(tally, _np.ndarray):
            found = _np.flatnonzero(tally)
            tally = dict(zip(found.tolist(), tally[found].tolist()))
        counts, skipped = {}, 0
        for byte, count in tally.items():
            if byte < offset:
                skipped += count
            else:
                counts[byte - offset] = counts.get(byte - offset, 0) + count
        return cls(counts, reads, skipped)

    @classmethod
    def from_read_tally(cls, tally, reads: int, bases: int, skipped: int,
                        read_mean: str = "error"):
        """Build a per-read count from a mapping or array indexed by score."""
        if _np is not None and isinstance(tally, _np.ndarray):
            found = _np.flatnonzero(tally)
            tally = dict(zip(found.tolist(), tally[found].tolist()))
        return cls(dict(tally), reads, skipped, by="read", bases=bases,
                   read_mean=read_mean)

    def merge(self, other: "QualityCounts") -> None:
        """Add another run's counts to this one."""
        if other.by != self.by or other.read_mean != self.read_mean:
            raise ValueError("counts scored differently cannot be merged")
        self._counts.update(other._counts)
        self.reads += other.reads
        self.skipped += other.skipped
        if self.by == "read":
            self._bases = (self._bases or 0) + (other._bases or 0)

    def items(self) -> List[Tuple[int, int]]:
        """Return ``(score, count)`` pairs in ascending score."""
        return sorted(self._counts.items())

    @property
    def unit(self) -> str:
        return _unit(self.by)

    @property
    def empty(self) -> bool:
        return not self._counts

    @property
    def total(self) -> int:
        """Units counted: bases, or reads."""
        return sum(self._counts.values())

    @property
    def bases(self) -> int:
        """Bases in the reads counted."""
        return self.total if self._bases is None else self._bases

    @property
    def lowest(self) -> int:
        return min(self._counts) if self._counts else 0

    @property
    def highest(self) -> int:
        return max(self._counts) if self._counts else 0

    @property
    def mean(self) -> float:
        total = self.total
        if not total:
            return 0.0
        return sum(q * n for q, n in self._counts.items()) / total

    def median(self) -> int:
        """The score half the counted units fall at or below."""
        total = self.total
        run = 0
        for score, count in self.items():
            run += count
            if run * 2 >= total:
                return score
        return 0

    def share_at_least(self, score: int) -> float:
        """The fraction of units scoring *score* or more."""
        total = self.total
        if not total:
            return 0.0
        return sum(n for q, n in self._counts.items() if q >= score) / total


def _quality_blocks(stream) -> Iterator[List[bytes]]:
    """Yield the quality lines in *stream*, one list per block read.

    A block is split on line boundaries and trimmed back to a whole number of
    records, the remainder carried forward, so a record split across two blocks
    is counted once, when its fourth line arrives.  A final record whose fourth
    line is absent is dropped, which is the state a file interrupted mid-write
    is left in.  One list is yielded per block, empty where a block completed no
    record, so a caller reporting progress hears from every block.
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
        yield lines[3:whole:4] if whole else []

    if buf:
        lines = buf.split(sep)
        if not lines[-1]:
            lines.pop()                 # a trailing separator, not a line
        whole = len(lines) - len(lines) % 4
        if whole:
            yield lines[3:whole:4]


def _tally_python(
    stream,
    on_block: Optional[Callable[[int], None]] = None,
) -> Tuple[Counter, int]:
    """Return the byte values on the quality lines in *stream*, and the reads.

    The lines of a block are joined and counted in one call, so the work per
    base is done in C rather than in a Python loop.  *on_block* is called after
    each block with the reads counted so far.
    """
    tally: Counter = Counter()
    reads = 0
    for block in _quality_blocks(stream):
        if block:
            tally.update(b"".join(block))
            reads += len(block)
        if on_block is not None:
            on_block(reads)
    return tally, reads


def _mean_score(total: int, length: int, offset: int) -> int:
    """The mean score of a read, rounded half up, in exact integer arithmetic.

    *total* is the sum of the read's bytes.  Integers rather than floats, so the
    two scanners cannot disagree over a mean that lands on a half.
    """
    return (2 * (total - offset * length) + length) // (2 * length)


def _error_table(offset: int) -> List[float]:
    """The error probability of each byte value, ``10 ** (-score / 10)``.

    Bytes below *offset* are not scores and get 0; a read holding one is set
    aside before its sum is used.
    """
    return [10.0 ** (-(b - offset) / 10.0) if b >= offset else 0.0
            for b in range(256)]


def _error_score(error_sum: float, length: int) -> int:
    """The score of a read's mean error rate, rounded half up and kept in range."""
    return min(255, max(0, int(math.floor(
        -10.0 * math.log10(error_sum / length) + 0.5))))


def _reads_python(
    stream,
    offset: int,
    mean: str = "error",
    on_block: Optional[Callable[[int], None]] = None,
) -> Tuple[Counter, int, int, int]:
    """Return each read's mean score tallied, with reads, bases and skipped.

    *mean* is ``"error"`` to average error probabilities, or ``"phred"`` to
    average the scores themselves.

    A read with no bases has no mean and is not counted.  A read holding a byte
    below *offset* has a mean that means nothing and is counted in *skipped*
    instead.  Each line is summed in C, so the loop is per read, not per base.
    """
    tally: Counter = Counter()
    reads = bases = skipped = 0
    table = _error_table(offset)
    for block in _quality_blocks(stream):
        for line in block:
            length = len(line)
            if not length:
                continue
            if min(line) < offset:
                skipped += 1
                continue
            if mean == "error":
                score = _error_score(sum(map(table.__getitem__, line)), length)
            else:
                score = _mean_score(sum(line), length, offset)
            tally[score] += 1
            reads += 1
            bases += length
        if on_block is not None:
            on_block(reads + skipped)
    return tally, reads, bases, skipped


def _numpy_spans(stream):
    """Yield ``(data, starts, ends)`` for the quality lines in each block.

    The newlines in a block are located in one pass, and the quality line of each
    record lies between its third newline and its fourth.  Blocks are trimmed
    and carried the same way as in :func:`_quality_blocks`, and a final record
    whose fourth line is absent is dropped the same way.  A block that completes
    no record yields empty arrays, so a caller reporting progress hears from
    every block.  *ends* leaves out the CR of a CRLF file.
    """
    empty = _np.zeros(0, dtype=_np.int64)
    carry = b""
    first = True
    drop = 0                            # the CR before each newline, if any

    while True:
        data = stream.read(_BLOCK)
        if not data:
            break
        if first:
            first = False
            drop = 1 if b"\r\n" in data else 0
        if carry:
            data = carry + data
        marks = _np.flatnonzero(_np.frombuffer(data, dtype=_np.uint8) == 10)
        whole = len(marks) - len(marks) % 4
        if whole < 4:
            carry = data                # a record longer than one block
            yield data, empty, empty
        else:
            carry = data[marks[whole - 1] + 1:]
            yield data, marks[2:whole:4] + 1, marks[3:whole:4] - drop

    if carry:
        marks = _np.flatnonzero(_np.frombuffer(carry, dtype=_np.uint8) == 10)
        # Three newlines and bytes after the last: a final record whose fourth
        # line arrived without one.  A trailing newline instead means the file
        # stops mid-record, and the record is dropped.
        if len(marks) == 3 and marks[-1] != len(carry) - 1:
            yield carry, marks[2:3] + 1, _np.array([len(carry)])


def _tally_numpy(
    stream,
    on_block: Optional[Callable[[int], None]] = None,
):
    """Return the byte values on the quality lines in *stream*, as an array.

    A mask over the block selects the bytes inside the quality lines, so the work
    per base stays inside the array layer.
    """
    tally = _np.zeros(256, dtype=_np.int64)
    reads = 0
    for data, starts, ends in _numpy_spans(stream):
        if starts.size:
            edges = _np.zeros(len(data) + 1, dtype=_np.int8)
            edges[starts] += 1
            edges[ends] -= 1
            inside = _np.cumsum(edges[:-1], dtype=_np.int8) > 0
            tally += _np.bincount(
                _np.frombuffer(data, dtype=_np.uint8)[inside], minlength=256)
            reads += int(starts.size)
        if on_block is not None:
            on_block(reads)
    return tally, reads


def _reads_numpy(
    stream,
    offset: int,
    mean: str = "error",
    on_block: Optional[Callable[[int], None]] = None,
):
    """Return each read's mean score tallied as an array, as :func:`_reads_python`.

    The bytes of a block are summed once, as a running total, so the sum over any
    line is a difference of two entries; the same is done for the bytes below
    *offset*, which is how a read holding one is found.  Error probabilities are
    summed per read from a table lookup of each byte, which keeps every sum local
    to its own read.
    """
    tally = _np.zeros(256, dtype=_np.int64)
    reads = bases = skipped = 0
    table = _np.array(_error_table(offset))
    for data, starts, ends in _numpy_spans(stream):
        if starts.size:
            raw = _np.frombuffer(data, dtype=_np.uint8)
            total = _np.concatenate(([0], _np.cumsum(raw, dtype=_np.int64)))
            low = _np.concatenate(
                ([0], _np.cumsum(raw < offset, dtype=_np.int64)))
            length = ends - starts
            sums = total[ends] - total[starts]
            has_low = (low[ends] - low[starts]) > 0
            scored = (length > 0) & ~has_low
            skipped += int(((length > 0) & has_low).sum())
            n = length[scored]
            if mean == "error":
                edges = _np.zeros(len(data) + 1, dtype=_np.int8)
                edges[starts] += 1
                edges[ends] -= 1
                inside = _np.cumsum(edges[:-1], dtype=_np.int8) > 0
                owner = _np.repeat(_np.arange(starts.size), length)
                errors = _np.bincount(owner, weights=table[raw[inside]],
                                      minlength=starts.size)
                means = _np.clip(_np.floor(
                    -10.0 * _np.log10(errors[scored] / n) + 0.5), 0, 255
                ).astype(_np.int64)
            else:
                means = (2 * (sums[scored] - offset * n) + n) // (2 * n)
            tally += _np.bincount(means, minlength=256)
            reads += int(scored.sum())
            bases += int(n.sum())
        if on_block is not None:
            on_block(reads + skipped)
    return tally, reads, bases, skipped


def count_qualities(
    paths: Iterable,
    offset: int = DEFAULT_OFFSET,
    progress: Optional[Callable[[int, int, int], None]] = None,
    fast: Optional[bool] = None,
    by: str = "base",
    mean: str = "error",
) -> QualityCounts:
    """Tally the base qualities in *paths* in one pass.

    *by* is ``"base"`` to count every base at its score, or ``"read"`` to count
    every read at the mean of its scores.  *mean* says how that mean is taken, and
    is used only for reads: ``"error"`` averages error probabilities and reports
    the score of the result, and ``"phred"`` averages the scores.  Nothing is held
    but the tally, so the cost in memory is the same for a run of any size.

    *progress* is called after each block and once at the end of each file, with
    the bytes read, the bytes to read and the reads counted so far.  For a
    gzipped file the byte figures are compressed bytes, which is what the file's
    size reports.

    *fast* chooses the scanner: the array one where None and numpy is installed,
    and the pure-Python one where False.
    """
    if by not in BY:
        raise ValueError(f"by must be one of {', '.join(BY)}")
    if mean not in MEANS:
        raise ValueError(f"mean must be one of {', '.join(MEANS)}")
    if fast is None:
        fast = HAVE_NUMPY
    elif fast and not HAVE_NUMPY:
        raise RuntimeError("the array scanner needs numpy installed")

    paths = list(paths)
    sizes = []
    for path in paths:
        try:
            sizes.append(os.path.getsize(path))
        except OSError:
            sizes.append(0)
    total = sum(sizes)

    counts = QualityCounts(by=by, read_mean=mean if by == "read" else None)
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
            if by == "base":
                scan = _tally_numpy if fast else _tally_python
                tally, reads = scan(stream, report)
                found = QualityCounts.from_tally(tally, reads, offset)
            else:
                scan = _reads_numpy if fast else _reads_python
                tally, reads, bases, skipped = scan(stream, offset, mean, report)
                found = QualityCounts.from_read_tally(tally, reads, bases,
                                                      skipped, mean)
                reads += skipped
        finally:
            if stream is not handle:
                stream.close()
            handle.close()
        counts.merge(found)
        done_bytes += size
        done_reads += reads
        if progress is not None:
            progress(done_bytes, total, done_reads)
    return counts


def bin_counts(counts: QualityCounts, bins: int = DEFAULT_BINS) -> Binning:
    """Cut *counts* into equal-width bins, from the lowest score seen upward.

    The width is the smallest whole number of scores that fits the range from
    zero to the highest seen in *bins* bins, so every bin holds the same number
    of scores and its label is exact.  Bins are aligned to multiples of the
    width, and the axis starts at the one holding the lowest score, so a run with
    no poor scores does not open on rows of nothing.
    """
    if counts.empty:
        return Binning([], 1, counts.unit)
    top = counts.highest
    width = max(1, math.ceil((top + 1) / max(1, bins)))
    tally = dict(counts.items())
    out = []
    for low in range(counts.lowest // width * width, top + 1, width):
        high = low + width
        out.append(Bin(low, high, sum(tally.get(q, 0)
                                      for q in range(low, high))))
    return Binning(out, width, counts.unit)


def summarise(counts: QualityCounts) -> Summary:
    """Return the figures over every base or read counted."""
    return Summary(
        by=counts.by,
        reads=counts.reads,
        bases=counts.bases,
        total=counts.total,
        lowest=counts.lowest,
        highest=counts.highest,
        median=counts.median(),
        mean=counts.mean,
        shares=tuple((q, counts.share_at_least(q)) for q in THRESHOLDS),
        skipped=counts.skipped,
        read_mean=counts.read_mean,
    )


def histogram(
    binning: Binning,
    width: int = 80,
    log: bool = False,
    palette: Optional[Palette] = None,
) -> List[str]:
    """Render *binning* as lines of text, none wider than *width*.

    *log* scales bar length by ``log(1 + count)``, which keeps the smaller bins
    of a peaked distribution distinguishable.  Axis labels stay linear under
    either scale.  *palette* colours the output; its codes are added after the
    columns are laid out, so they neither shift the alignment nor count toward
    *width*.
    """
    if not binning.bins:
        return [f"no {binning.unit}"]

    pal = palette or Palette()
    rows = []
    for b in binning.bins:
        top = b.high - 1
        label = f"Q{b.low}" if binning.width == 1 else f"Q{b.low}–{top}"
        rows.append((label, b.count))

    counts = [f"{n:,}" for _, n in rows]
    label_w = max(len(r[0]) for r in rows + [("Q", 0)])
    count_w = max(len(x) for x in counts + [binning.unit])
    bar_w = max(8, width - label_w - count_w - 4)

    def height(count: int) -> float:
        return math.log1p(count) if log else float(count)

    ceiling = max((height(n) for _, n in rows), default=1.0) or 1.0
    peak = max(n for _, n in rows)

    head = f"{'Q':>{label_w}}  {'':{bar_w}}  {binning.unit:>{count_w}}"
    lines = [f"{pal.head}{head}{pal.reset}" if pal.head else head]

    for (label, count), shown in zip(rows, counts):
        eighths = int(round(bar_w * 8 * min(1.0, height(count) / ceiling)))
        if count and eighths < 1:
            eighths = 1                        # a row with counts is never blank
        bar = f"{_bar(eighths):{bar_w}}"
        colour = pal.peak if count == peak else pal.bar
        if colour:
            bar = f"{colour}{bar}{pal.reset}"
        lines.append(f"{label:>{label_w}}  {bar}  {shown:>{count_w}}")
    return lines


def _how(summary: Summary) -> str:
    """How the reads were scored, for the line of totals; nothing for bases."""
    if summary.by != "read":
        return ""
    if summary.read_mean == "phred":
        return " · scored by each read's mean score"
    return " · scored by each read's mean error rate"


def summary_lines(summary: Summary) -> List[str]:
    """Return the distribution's figures as lines of text."""
    if summary.empty:
        return [f"no {summary.unit}"]
    unit = summary.unit
    lines = [
        f"{summary.reads:,} reads · {summary.bases:,} bases"
        + _how(summary),
        f"min Q{summary.lowest} · median Q{summary.median} · "
        f"mean Q{summary.mean:.1f} · max Q{summary.highest}",
        " · ".join(f"{100 * share:.1f}% of {unit} at Q{q} or above"
                   for q, share in summary.shares),
    ]
    if summary.skipped:
        if summary.by == "base":
            lines.append(f"{summary.skipped:,} bases scored below the offset "
                         f"and were left out. Check --offset.")
        else:
            lines.append(f"{summary.skipped:,} reads held a base scored below "
                         f"the offset and were left out. Check --offset.")
    return lines
