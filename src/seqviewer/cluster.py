"""Order pileup rows by hierarchical clustering of what they disagree about.

The pileup's other orderings sort rows by a key.  ``--order mismatch`` sorts the
per-read mismatch pattern as a string, which puts identical patterns together but
is dominated by the leftmost difference: one wild-type read carrying an early
sequencing error sorts between two halves of a real subpopulation and splits the
block.  This module clusters instead, so reads group by how similar they are
overall rather than by where they first differ.

Three things make it affordable in pure Python, with no new dependency:

**Reads collapse onto signatures.**
    Two reads disagreeing with the reference at exactly the same positions are one
    row of the distance matrix, not two.  Clustering runs over the distinct
    signatures and the reads are expanded back afterwards.

**A signature is an integer.**
    One bit per column that varies anywhere in the group.  Distance is then
    ``popcount(a ^ b)`` — the Hamming distance over those columns — which is a
    single machine operation rather than a loop over the reference.

**Linkage is nearest-neighbour chain.**
    Average linkage is a reducible linkage, so Müllner's NN-chain algorithm
    (2011) computes the exact same dendrogram as the naive method in O(n²) rather
    than O(n³).  Measured on this implementation, the naive method needs about
    150s for 2000 distinct signatures and this needs a few seconds.

For very long runs the linkage reports progress, but only once it has been going
long enough to be worth saying so — see :data:`PROGRESS_AFTER`.
"""

from __future__ import annotations

import sys
import time
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from .grid import Row

__all__ = [
    "PROGRESS_AFTER",
    "cluster_rows",
    "leaf_order",
    "nn_chain",
    "signatures",
    "variable_columns",
]

#: Seconds a linkage runs before it starts reporting progress.  Silence under
#: this: a progress line for an operation that finishes in half a second is
#: noise, and the common case here is well under that.
PROGRESS_AFTER = 1.5

#: Minimum gap between progress reports once they start.
PROGRESS_EVERY = 0.2

_bit_count = getattr(int, "bit_count", None)
if _bit_count is None:                          # Python 3.9 and 3.10
    def _popcount(value: int) -> int:
        return bin(value).count("1")
else:
    def _popcount(value: int) -> int:
        return value.bit_count()


def variable_columns(rows: Sequence[Row], ref_seq: str) -> List[int]:
    """Reference positions where at least one read disagrees.

    The clustering space.  A column every read agrees on carries no information
    about which reads belong together, and dropping it is what keeps a signature
    to a handful of bits on a mostly-clean reference.
    """
    ref = ref_seq.upper()
    out = []
    for col in range(len(ref)):
        for row in rows:
            base = row[col][0]
            if base != "-" and base.upper() != ref[col]:
                out.append(col)
                break
    return out


def signatures(
    rows: Sequence[Row], ref_seq: str, columns: Sequence[int],
) -> List[int]:
    """One integer per read: bit *k* set where the read disagrees at ``columns[k]``.

    A deletion is not a disagreement here even though
    :func:`seqviewer.summary.mismatch_fractions` counts it as one.  The question
    is different: that one asks how much a column disagrees, this one asks which
    reads carry the same substitutions.  A read that stops short would otherwise
    look like a read that deleted everything it never covered.
    """
    ref = ref_seq.upper()
    out = []
    for row in rows:
        mask = 0
        for bit, col in enumerate(columns):
            base = row[col][0]
            if base != "-" and base.upper() != ref[col]:
                mask |= 1 << bit
        out.append(mask)
    return out


def _distance_matrix(masks: Sequence[int]) -> List[List[float]]:
    """Full symmetric Hamming distances between signatures."""
    n = len(masks)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        mi = masks[i]
        Di = D[i]
        for j in range(i + 1, n):
            d = float(_popcount(mi ^ masks[j]))
            Di[j] = d
            D[j][i] = d
    return D


def nn_chain(
    D: List[List[float]],
    progress: Optional[Callable[[int, int], None]] = None,
) -> List[Tuple[int, int, float]]:
    """Average-linkage merges of *D*, by nearest-neighbour chain.

    *D* is a full symmetric distance matrix and is consumed: distances are
    updated in place by the Lance-Williams rule as clusters merge, and the
    surviving cluster of each merge reuses the lower of the two slots, so the
    matrix never grows.

    Returns one ``(survivor, absorbed, height)`` per merge, in merge order.
    *progress* is called with ``(merges done, merges total)`` when supplied.
    """
    n = len(D)
    if n < 2:
        return []

    active = [True] * n
    size = [1] * n
    chain: List[int] = []
    merges: List[Tuple[int, int, float]] = []

    for step in range(n - 1):
        if not chain:
            chain = [next(i for i in range(n) if active[i])]

        # Walk the chain until its last two entries are each other's nearest
        # neighbour.  A nearest-neighbour chain always reaches such a pair, and
        # because average linkage is reducible the rest of the chain stays valid
        # across the merge — which is what makes this quadratic and not cubic.
        while True:
            x = chain[-1]
            Dx = D[x]
            best, y = float("inf"), -1
            for i in range(n):
                if not active[i] or i == x:
                    continue
                if Dx[i] < best:
                    best, y = Dx[i], i
            if len(chain) > 1 and y == chain[-2]:
                break
            chain.append(y)

        x = chain.pop()
        y = chain.pop()
        survivor, absorbed = (x, y) if x < y else (y, x)
        merges.append((survivor, absorbed, best))

        ns, na = size[survivor], size[absorbed]
        total = ns + na
        Ds, Da = D[survivor], D[absorbed]
        for i in range(n):
            if not active[i] or i == survivor or i == absorbed:
                continue
            merged = (Ds[i] * ns + Da[i] * na) / total
            Ds[i] = merged
            D[i][survivor] = merged
        size[survivor] = total
        active[absorbed] = False

        if progress is not None:
            progress(step + 1, n - 1)

    return merges


def leaf_order(merges: Sequence[Tuple[int, int, float]], n: int) -> List[int]:
    """Leaves of the dendrogram, left to right.

    Each merge puts the absorbed cluster's leaves immediately after the
    survivor's, so the order a reader sees is the dendrogram's own.
    """
    blocks: Dict[int, List[int]] = {i: [i] for i in range(n)}
    for survivor, absorbed, _ in merges:
        blocks[survivor] = blocks[survivor] + blocks.pop(absorbed)
    # One block survives unless the matrix was empty.
    return next(iter(blocks.values())) if blocks else []


def _reporter(label: str):
    """Progress to stderr, but only for a run long enough to warrant it.

    Silent on a pipe or a file: a redraw sequence in a log is noise, and the log
    written beside a page is read afterwards rather than watched.
    """
    if not sys.stderr.isatty():
        return None

    state = {"started": 0.0, "last": 0.0}
    begin = time.perf_counter()

    def report(done: int, total: int) -> None:
        now = time.perf_counter()
        if not state["started"]:
            if now - begin < PROGRESS_AFTER:
                return
            state["started"] = now
        if done < total and now - state["last"] < PROGRESS_EVERY:
            return
        state["last"] = now
        width = 24
        filled = int(width * done / total) if total else width
        bar = "#" * filled + "-" * (width - filled)
        sys.stderr.write(f"\r{label} [{bar}] {done}/{total}")
        if done >= total:
            sys.stderr.write("\r" + " " * (width + len(label) + 24) + "\r")
        sys.stderr.flush()

    return report


def cluster_rows(
    rows: Sequence[Row],
    ref_seq: str,
    columns: Optional[Sequence[int]] = None,
    progress: bool = True,
) -> List[Row]:
    """Reorder *rows* so that similar reads sit together.

    Args:
        rows: The grid, one row per read, each as wide as *ref_seq*.
        ref_seq: The reference the rows are stated against.
        columns: Positions to cluster on.  Defaults to every column some read
            disagrees at.  Pass a narrower set — the flagged columns, say — to
            cluster on the positions a subpopulation is defined by and ignore
            each read's own errors.
        progress: Whether a long linkage may report progress to a terminal.

    Returns:
        The same rows, reordered.  Reads sharing a signature stay in the order
        they arrived, so the result is deterministic.
    """
    if len(rows) < 3:
        return list(rows)

    if columns is None:
        columns = variable_columns(rows, ref_seq)
    if not columns:
        return list(rows)                      # nothing disagrees; nothing to order

    masks = signatures(rows, ref_seq, columns)

    # Distinct signatures, first-seen order, with the reads that carry each.
    groups: Dict[int, List[int]] = {}
    for index, mask in enumerate(masks):
        groups.setdefault(mask, []).append(index)
    distinct = list(groups)
    if len(distinct) < 2:
        return list(rows)

    D = _distance_matrix(distinct)
    reporter = _reporter("clustering") if progress else None
    merges = nn_chain(D, reporter)
    order = leaf_order(merges, len(distinct))

    out: List[Row] = []
    for slot in order:
        for index in groups[distinct[slot]]:
            out.append(rows[index])
    return out
