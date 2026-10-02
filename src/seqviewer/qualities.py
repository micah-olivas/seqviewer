"""Read quality as a histogram: one bar per whole Q, one count per read.

Each read is reduced to one number, its mean quality, and the run is tallied by
that number the way :mod:`seqviewer.lengths` tallies it by length.  The tally,
the binning, the bars and the streaming are that module's; what is here is how
a quality line becomes a number, and the words the histogram uses.

A read's mean quality is the mean of its per-base *error rates*, expressed as a
Phred score: ``-10 log10(mean(10 ** (-Q / 10)))``.  Averaging the scores
themselves gives a larger number, because Q is a logarithm and the mean of
logarithms is at least the logarithm of the mean -- by more, the more a read's
scores vary.  ``docs/index.html`` works an example.

Scores are read as Phred+33, which Illumina 1.8 and later, nanopore and PacBio
write.  A record whose quality line is empty has nothing to average and is
counted apart rather than binned.

The tally is kept in tenths of Q, so the figures are exact to a tenth, and the
bars are drawn by whole Q.
"""
from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from typing import Callable, Iterable, List, Optional, Sequence, Tuple

from . import lengths
from .lengths import (Axis, Bin, Binning, HAVE_NUMPY, LengthCounts, Palette,
                      histogram, record_lines, scan_files)

try:                                    # the array scanner is optional
    import numpy as _np
except ImportError:                     # pragma: no cover - exercised without it
    _np = None

__all__ = [
    "DEFAULT_BINS",
    "DEFAULT_THRESHOLDS",
    "QUALITY_AXIS",
    "QualitySummary",
    "bin_qualities",
    "count_qualities",
    "mean_quality",
    "render",
    "summarise_qualities",
    "summary_lines",
]

#: The offset scores are stored at: Phred+33.
PHRED_OFFSET = 33

#: The highest score a FASTQ can hold: ``~`` is 126, less the offset.
MAX_Q = 93

#: The tally's resolution: one key per tenth of Q.
SCALE = 10

#: Bins across the axis.  The run in ``docs/examples`` spans thirty-one Q, so
#: this is one bar per whole Q there; a run wider than it has neighbouring Q
#: merged into a bar.
DEFAULT_BINS = 32

#: Scores the summary reports the share of reads at or above.  The Phred
#: decades: one error in ten bases, a hundred, a thousand.
DEFAULT_THRESHOLDS = (10, 20, 30)

#: Error rate per score byte.  Bytes below the offset are not a score and are
#: read as Q0, the worst there is, rather than raising mid-file; bytes above
#: ``~`` are read as Q93.
_ERR = [10.0 ** (-max(0, min(b - PHRED_OFFSET, MAX_Q)) / 10.0)
        for b in range(256)]
_ERR_NP = _np.array(_ERR, dtype=_np.float64) if _np is not None else None

#: Tally length for the array scanner: every tenth of Q up to Q93.
_KEYS = MAX_Q * SCALE + 1


def mean_quality(quality: bytes) -> Optional[float]:
    """Return the mean quality of one quality line, or None if it is empty.

    The mean is of error rates, not of scores; see the module docstring.
    """
    if not quality:
        return None
    total = sum(n * _ERR[b] for b, n in Counter(quality).items())
    return -10.0 * math.log10(total / len(quality))


#: Added before flooring a key.  A read whose bases share one score lands
#: exactly on a tenth, and the two scanners reach it by different float paths
#: that can fall either side of it by 1e-12; this keeps both on the same key.
_EPS = 1e-6


def _key(q: float) -> int:
    """The tally key for a mean quality: tenths of Q, clamped to the scale.

    Floored, not rounded, so that a key is the tenth the read is *in*: bar 19
    holds Q19.0 up to Q20.0, and a read counted at or above Q20 is one.
    Rounding put a read at Q19.96 in bar 20 and over the Q20 threshold.
    """
    return min(max(int(math.floor(q * SCALE + _EPS)), 0), MAX_Q * SCALE)


# --------------------------------------------------------------------------
# Scanners.  Both return ``(tally, scored, unscored)``; the pass that runs them
# wants ``(tally, reads)`` and is handed a closure that keeps the third.
# --------------------------------------------------------------------------

def _tally_python(
    stream,
    on_block: Optional[Callable[[int, Counter], None]] = None,
    stop: Optional[Callable[[], bool]] = None,
) -> Tuple[Counter, int, int]:
    """Tally the reads in *stream* by mean quality, in tenths of Q.

    Each quality line is counted by byte value in one C call, so the work in
    Python per read is one step per distinct score rather than per base.
    """
    tally: Counter = Counter()
    scored = unscored = 0
    for block in record_lines(stream, 3):
        for quality in block:
            q = mean_quality(quality)
            if q is None:
                unscored += 1
                continue
            tally[_key(q)] += 1
            scored += 1
        if on_block is not None:
            on_block(scored, tally)
        if stop is not None and stop():
            break
    return tally, scored, unscored


def _tally_numpy(
    stream,
    on_block: Optional[Callable[[int, "object"], None]] = None,
    stop: Optional[Callable[[], bool]] = None,
):
    """Tally the reads in *stream* by mean quality, with the work in arrays.

    Every byte of a block is mapped to its error rate in one lookup and the
    rates are summed once, cumulatively; a read's total is then the difference
    of two entries of that sum at the ends of its quality line.  Blocks are
    trimmed and carried as :func:`seqviewer.lengths._tally_numpy` trims them,
    and an unterminated final record is read the same way.
    """
    tally = _np.zeros(_KEYS, dtype=_np.int64)
    carry = b""
    crlf = False
    first = True
    scored = unscored = 0

    def take(data, starts, ends):
        nonlocal tally, scored, unscored
        rates = _ERR_NP[_np.frombuffer(data, dtype=_np.uint8)]
        # Summed into a buffer one longer than the block, with a zero in front,
        # so a read's total is running[end] - running[start] with no special
        # case at the first byte -- and without concatenating a second copy.
        running = _np.empty(rates.size + 1, dtype=_np.float64)
        running[0] = 0.0
        _np.cumsum(rates, out=running[1:])
        spans = ends - starts
        has = spans > 0
        unscored += int((~has).sum())
        if not has.any():
            return
        sums = running[ends[has]] - running[starts[has]]
        means = sums / spans[has]
        keys = _np.clip(_np.floor(-10.0 * SCALE * _np.log10(means) + _EPS),
                        0, MAX_Q * SCALE).astype(_np.int64)
        tally += _np.bincount(keys, minlength=_KEYS)
        scored += int(keys.size)

    while True:
        # Read at call time, not import time, so the block size the length
        # scanner uses -- and a test that shrinks it -- governs this one too.
        data = stream.read(lengths._BLOCK)
        if not data:
            break
        if first:
            first = False
            crlf = b"\r\n" in data
        if carry:
            data = carry + data
        marks = _np.flatnonzero(_np.frombuffer(data, dtype=_np.uint8) == 10)
        whole = len(marks) - len(marks) % 4
        if whole < 4:
            carry = data                # a record longer than one block
        else:
            starts = marks[2:whole:4] + 1
            ends = marks[3:whole:4] - (1 if crlf else 0)
            take(data, starts, ends)
            carry = data[marks[whole - 1] + 1:]
        if on_block is not None:
            on_block(scored, tally)
        if stop is not None and stop():
            return tally, scored, unscored

    if carry:
        marks = _np.flatnonzero(_np.frombuffer(carry, dtype=_np.uint8) == 10)
        # A final record whose fourth line arrived without a newline: the same
        # test the length scanner applies.
        if len(marks) == 3 and marks[-1] != len(carry) - 1:
            end = len(carry)
            if crlf and carry.endswith(b"\r"):
                end -= 1
            take(carry, _np.array([marks[2] + 1]), _np.array([end]))

    return tally, scored, unscored


@dataclass
class QualityTally:
    """Reads tallied by mean quality, and the records that had none.

    ``counts`` is keyed in tenths of Q; divide a key by :data:`SCALE` for the
    score.
    """

    counts: LengthCounts = field(default_factory=LengthCounts)
    unscored: int = 0

    @property
    def empty(self) -> bool:
        return self.counts.empty


def count_qualities(
    paths: Iterable,
    progress: Optional[Callable[[int, int, int, Callable], None]] = None,
    fast: Optional[bool] = None,
    stop: Optional[Callable[[], bool]] = None,
) -> QualityTally:
    """Tally the reads in *paths* by mean quality in one pass.

    *progress*, *fast* and *stop* mean what they mean to
    :func:`seqviewer.lengths.count_lengths`, with the read count in *progress*
    being reads scored.
    """
    if fast is None:
        fast = HAVE_NUMPY
    elif fast and not HAVE_NUMPY:
        raise RuntimeError("the array scanner needs numpy installed")
    scanner = _tally_numpy if fast else _tally_python
    unscored = 0

    def scan(stream, on_block, should_stop):
        nonlocal unscored
        tally, scored, missing = scanner(stream, on_block, should_stop)
        unscored += missing
        return tally, scored

    counts = scan_files(paths, scan, progress, stop)
    return QualityTally(counts, unscored)


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class QualitySummary:
    """Per-read mean quality over every scored read.

    ``at_least`` pairs each threshold with the reads scoring at or above it.
    """

    reads: int
    unscored: int
    lowest: float
    median: float
    highest: float
    at_least: Tuple[Tuple[int, int], ...] = ()

    @property
    def empty(self) -> bool:
        return self.reads == 0


def summarise_qualities(
    tally: QualityTally,
    thresholds: Sequence[int] = DEFAULT_THRESHOLDS,
) -> QualitySummary:
    """Return the figures for *tally*: extremes, median, and shares over Q."""
    counts = tally.counts
    if counts.empty:
        return QualitySummary(0, tally.unscored, 0.0, 0.0, 0.0)
    total = counts.total
    at_least = []
    for q in thresholds:
        below, _, _ = counts.between(q * SCALE, MAX_Q * SCALE)
        at_least.append((q, total - below))
    return QualitySummary(
        reads=total,
        unscored=tally.unscored,
        lowest=counts.shortest / SCALE,
        median=counts.median / SCALE,
        highest=counts.longest / SCALE,
        at_least=tuple(at_least),
    )


def summary_lines(summary: QualitySummary) -> List[str]:
    """Return the figures as lines of text."""
    if summary.empty:
        if summary.unscored:
            return [f"no reads with qualities; {summary.unscored:,} records "
                    "had an empty quality line"]
        return ["no reads"]
    lines = [f"{summary.reads:,} reads",
             f"min Q{summary.lowest:.1f} · median Q{summary.median:.1f} · "
             f"max Q{summary.highest:.1f}"]
    if summary.at_least:
        lines.append(" · ".join(
            f"≥Q{q} {n:,} ({100 * n / summary.reads:.1f}%)"
            for q, n in summary.at_least))
    if summary.unscored:
        lines.append(f"{summary.unscored:,} records with an empty quality line "
                     "are left out")
    return lines


# --------------------------------------------------------------------------
# Drawing
# --------------------------------------------------------------------------

class QualityAxis(Axis):
    """The histogram's words for quality: whole Q per bar, from tenths."""

    unit = "Q"

    def bin_label(self, b: Bin) -> str:
        low, high = b.low // SCALE, (b.high - 1) // SCALE
        return f"{low}" if low == high else f"{low}–{high}"

    def below(self, binning: Binning) -> Tuple[str, str]:
        return (f"<{binning.low // SCALE}",
                f"lower, down to Q{binning.shortest / SCALE:.1f}")

    def above(self, binning: Binning) -> Tuple[str, str]:
        return (f">{binning.high // SCALE}",
                f"higher, up to Q{binning.longest / SCALE:.1f}")


#: The axis the quality histogram is drawn on.
QUALITY_AXIS = QualityAxis()


def bin_qualities(tally: QualityTally, count: int = DEFAULT_BINS) -> Binning:
    """Bin *tally* by whole Q, over every read.

    Unlike lengths, nothing is clipped: scores are bounded, so there are no
    concatemer-like outliers to stretch the axis.
    """
    return lengths.bin_counts(tally.counts, count, bulk=100.0, step=SCALE)


def render(
    tally: QualityTally,
    bins: int = DEFAULT_BINS,
    width: int = 80,
    log: bool = False,
    palette: Optional[Palette] = None,
    thresholds: Sequence[int] = DEFAULT_THRESHOLDS,
    stats: bool = True,
) -> Tuple[List[str], List[str]]:
    """Return the histogram's lines and, with *stats*, the figures under it."""
    hist = histogram(bin_qualities(tally, bins), width, log, palette,
                     axis=QUALITY_AXIS)
    texts = (summary_lines(summarise_qualities(tally, thresholds))
             if stats else [])
    return hist, texts
