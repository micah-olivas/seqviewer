"""Reduce a pileup to the few numbers a one-glance map draws.

The pileup viewer answers "what does every read say".  This module answers the
other question — "what does this alignment say, in total" — by collapsing a grid
of reads into three per-position arrays and a list of called variants.  Nothing
here draws anything; :mod:`seqviewer.render_summary` consumes what this produces.

The reduction is deliberately the whole of the statistics.  A renderer that
recomputed disagreement per column would be a second implementation of the same
number, which is how a legend comes to disagree with the cells it describes.

Two properties of the grid shape what can be recovered from it:

**Deletions and uncovered positions are the same cell.**
    :func:`seqviewer.align._row_for` starts every row as ``("-", True)`` and
    writes ``("-", True)`` for a deletion, so the two are byte-identical.  They
    are separated here by position rather than by content: a run of ``"-"``
    strictly inside a read's first and last called base is a deletion, and one
    running off either end is no coverage.  That rule is load-bearing on reads
    being contiguous in reference coordinates, which holds because ``_row_for``
    only ever writes where ``get_aligned_pairs`` reports.  Soft-clip or
    supplementary-alignment handling in ``align`` would break it.

**Insertions are not in the grid at all.**
    An inserted base has no reference position, and ``_row_for`` drops it.  A
    row is exactly ``len(ref_seq)`` wide by construction, so there is nowhere to
    put one.  :class:`Variant` therefore models ``kind="ins"`` and
    :func:`summarize_group` accepts an *insertions* sidecar, but a plain grid
    supplies none and a summary of one reports no insertions — which is the
    honest answer, not the absence of any.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field, replace
from typing import List, Mapping, Optional, Sequence, Tuple

from .codon import translate_codon
from .construct import Feature
from .grid import Row
from .pileup import PileupGroup, PileupView
from .theme import Theme

__all__ = [
    "DEFAULT_FLAG_THRESHOLD",
    "DEFAULT_MIN_COUNT",
    "DEFAULT_MIN_DEPTH",
    "DEFAULT_MIN_FRACTION",
    "GroupSummary",
    "SummaryView",
    "Variant",
    "flagged_columns",
    "mismatch_fractions",
    "summarize_group",
]

#: Fraction of covering reads an allele needs before it is called.
#:
#: Higher than the 10% the pileup page flags columns at, and deliberately so:
#: flagging marks a column worth a human's eye, while calling asserts an allele
#: is really there.  At the shallow depths these pages are made for, 10% is one
#: read — on a 2% per-base error rate over 10 reads that fires on hundreds of
#: columns, which buries the one mutation that matters.
DEFAULT_MIN_FRACTION = 0.25

#: Reads that must support an allele before it is called, whatever the fraction.
#: A single read is never a variant; it is the error rate.
DEFAULT_MIN_COUNT = 2

#: Reads that must cover a position before it is called at all.
DEFAULT_MIN_DEPTH = 3

#: Share of covering reads that must disagree with the reference before a column
#: is worth a reader's eye.  Lower than the calling threshold on purpose: marking
#: a column and asserting an allele are different claims.
DEFAULT_FLAG_THRESHOLD = 0.10

#: Consequences, worst first.  A group's verdict is the worst one it carries,
#: and the renderer styles on these rather than on any free text.
SEVERITY: Tuple[str, ...] = (
    "frameshift",
    "nonsense",
    "inframe_indel",
    "missense",
    "silent",
    "noncoding",
)


@dataclass(frozen=True)
class Variant:
    """One called difference from the reference.

    ``count`` is the reads supporting it and ``depth`` the reads covering the
    position, both counted over the group this variant was called in.

    ``consequence`` is the machine-readable classification — one of
    :data:`SEVERITY`, or ``""`` when no reading frame was known — and ``effect``
    is its human form, ``"T40N"`` or ``"Δ1 aa at 151"``.  Styling reads the
    former; only display reads the latter.
    """

    pos: int
    kind: str
    ref: str
    alt: str
    count: int
    depth: int
    consequence: str = ""
    effect: str = ""

    @property
    def fraction(self) -> float:
        """Share of covering reads carrying this allele."""
        return self.count / self.depth if self.depth else 0.0

    @property
    def length(self) -> int:
        """Reference bases affected: a deletion's run, otherwise one."""
        return len(self.ref) if self.kind == "del" else 1

    @property
    def label(self) -> str:
        """Short form for a table cell or a glyph label."""
        if self.kind == "del":
            return f"Δ{len(self.ref)} bp"
        if self.kind == "ins":
            return f"+{len(self.alt)} bp"
        return f"{self.ref}→{self.alt}"


@dataclass
class GroupSummary:
    """One subpopulation, reduced.

    The three arrays are all ``ref_len`` long and are the whole of the
    per-position statistics: ``depth`` counts reads covering a position whether
    they called a base or deleted it, ``matches`` counts those agreeing with the
    reference, and ``deletions`` counts those deleting it.  Bases actually called
    at a position are ``depth - deletions``.

    No reference sequence is held here.  It belongs to the view, which draws one
    coordinate system; keeping a copy per group is how the pileup page comes to
    serialise the same 1000 bases once for every subpopulation it shows.
    """

    name: str
    ref_len: int
    depth: List[int] = field(default_factory=list)
    matches: List[int] = field(default_factory=list)
    deletions: List[int] = field(default_factory=list)
    variants: List[Variant] = field(default_factory=list)
    n_reads: int = 0
    fraction: float = 0.0
    rows_drawn: int = 0
    highlighted: bool = False
    #: Whatever the caller's pipeline called this group.  Display only — the
    #: verdict below is computed, so a page never styles on free text.
    status: str = ""

    @property
    def verdict(self) -> str:
        """The worst consequence carried by any called variant.

        ``"clean"`` when nothing was called.  ``"variant"`` when something was
        called but no reading frame was known to classify it, which is the
        honest answer for a construct with no annotated insert.
        """
        if not self.variants:
            return "clean"
        ranked = [v.consequence for v in self.variants if v.consequence]
        if not ranked:
            return "variant"
        return min(ranked, key=SEVERITY.index)

    @property
    def covered(self) -> int:
        """Reference positions any read reached."""
        return sum(1 for d in self.depth if d)

    @property
    def mean_depth(self) -> float:
        """Mean coverage over covered positions, not over the whole reference."""
        covered = [d for d in self.depth if d]
        return sum(covered) / len(covered) if covered else 0.0

    @property
    def max_depth(self) -> int:
        return max(self.depth) if self.depth else 0

    @property
    def identity(self) -> Optional[float]:
        """Share of called bases agreeing with the reference, or None if none were."""
        called = sum(self.depth) - sum(self.deletions)
        return sum(self.matches) / called if called else None


def _covered_span(row: Row) -> Tuple[Optional[int], Optional[int]]:
    """First and last position this read called a base at, or ``(None, None)``."""
    first: Optional[int] = None
    last: Optional[int] = None
    for i, (base, _) in enumerate(row):
        if base != "-":
            if first is None:
                first = i
            last = i
    return first, last


def _deletion_runs(row: Row, first: int, last: int) -> List[Tuple[int, int]]:
    """Half-open runs of ``"-"`` lying strictly inside the read's covered span.

    Anything outside that span is absence of evidence rather than evidence of a
    deletion, which is the only thing separating the two in a grid.
    """
    runs: List[Tuple[int, int]] = []
    start: Optional[int] = None
    for i in range(first, last + 1):
        if row[i][0] == "-":
            if start is None:
                start = i
        elif start is not None:
            runs.append((start, i))
            start = None
    if start is not None:  # unreachable: row[last] is a called base
        runs.append((start, last + 1))
    return runs


def _deletion_variants(
    passing: Sequence[bool],
    runs_by_row: Sequence[Sequence[Tuple[int, int]]],
    depth: Sequence[int],
    deletions: Sequence[int],
    ref_seq: str,
    min_count: int,
) -> List[Variant]:
    """Merge adjacent deleted positions into one variant per run.

    Support is the reads whose own deletion covers the whole merged run, which
    is the exact count for the allele being reported.  Reads deleting staggered,
    partly-overlapping spans support no single run; when that leaves a run with
    no whole-run support, the thinnest per-position count stands in, so a real
    deletion is still reported rather than vanishing between two tallies.
    """
    out: List[Variant] = []
    i, n = 0, len(ref_seq)
    while i < n:
        if not passing[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and passing[j + 1]:
            j += 1
        start, end = i, j + 1
        whole = sum(
            1 for runs in runs_by_row
            if any(s <= start and e >= end for s, e in runs)
        )
        count = whole or min(deletions[start:end])
        if count >= min_count:
            out.append(Variant(
                pos=start,
                kind="del",
                ref=ref_seq[start:end],
                alt="",
                count=count,
                depth=depth[start],
            ))
        i = end
    return out


def summarize_group(
    group: PileupGroup,
    *,
    min_fraction: float = DEFAULT_MIN_FRACTION,
    min_count: int = DEFAULT_MIN_COUNT,
    min_depth: int = DEFAULT_MIN_DEPTH,
    insertions: Optional[Mapping[int, Mapping[str, int]]] = None,
) -> GroupSummary:
    """Collapse *group*'s reads into per-position arrays and called variants.

    Args:
        group: The reads to reduce, as the pileup viewer takes them.
        min_fraction: Share of covering reads an allele needs to be called.
        min_count: Reads that must support an allele, whatever the fraction.
        min_depth: Reads that must cover a position before it is called.
        insertions: Optional sidecar evidence, ``position -> {sequence: count}``,
            for insertions the grid cannot carry.  A grid alone supplies none.

    Returns:
        The reduction.  Variants come back unclassified — assigning a
        consequence needs a reading frame, which belongs to the view, so
        :meth:`SummaryView.from_view` is what fills ``consequence`` and
        ``effect`` in.
    """
    ref_seq = group.ref_seq.upper()
    n = len(ref_seq)

    calls: List[Counter] = [Counter() for _ in range(n)]
    deletions = [0] * n
    runs_by_row: List[Sequence[Tuple[int, int]]] = []

    for row in group.rows:
        first, last = _covered_span(row)
        if first is None or last is None:
            continue                        # a read that called nothing
        runs = _deletion_runs(row, first, last)
        runs_by_row.append(runs)
        for start, end in runs:
            for i in range(start, end):
                deletions[i] += 1
        for i in range(first, last + 1):
            base = row[i][0]
            if base != "-":
                calls[i][base.upper()] += 1

    depth = [sum(calls[i].values()) + deletions[i] for i in range(n)]
    matches = [calls[i].get(ref_seq[i], 0) for i in range(n)]

    variants: List[Variant] = []
    deleted_enough = [False] * n
    for i in range(n):
        if depth[i] < min_depth:
            continue
        for base, count in calls[i].items():
            if base == ref_seq[i] or count < min_count:
                continue
            if count / depth[i] >= min_fraction:
                variants.append(Variant(
                    pos=i, kind="snv", ref=ref_seq[i], alt=base,
                    count=count, depth=depth[i],
                ))
        if deletions[i] / depth[i] >= min_fraction:
            deleted_enough[i] = True

    variants.extend(_deletion_variants(
        deleted_enough, runs_by_row, depth, deletions, ref_seq, min_count,
    ))

    for pos, observed in (insertions or {}).items():
        if not 0 <= pos < n:
            continue
        covering = depth[pos]
        for seq, count in observed.items():
            if covering < min_depth or count < min_count:
                continue
            if count / covering >= min_fraction:
                variants.append(Variant(
                    pos=pos, kind="ins", ref="", alt=seq.upper(),
                    count=count, depth=covering,
                ))

    variants.sort(key=lambda v: (v.pos, v.kind, v.alt))

    return GroupSummary(
        name=group.name,
        ref_len=n,
        depth=depth,
        matches=matches,
        deletions=deletions,
        variants=variants,
        n_reads=group.n_reads,
        fraction=group.fraction,
        rows_drawn=len(group.rows),
        highlighted=group.highlighted,
        status=group.status,
    )


def mismatch_fractions(rows: Sequence[Row], ref_seq: str) -> List[float]:
    """Per position, the share of covering reads that disagree with the reference.

    The one definition of disagreement in the package.  Both the pileup's track
    and this module's calls read it, so they cannot report different numbers for
    the same column — which they did while each computed its own.

    Two choices in it, both load-bearing:

    * **A deletion is disagreement.** Counting only called bases makes a column
      where half the reads deleted the base read as perfectly clean, because the
      deleted reads leave both the numerator and the denominator.
    * **The denominator is reads that reached the position**, not every read in
      the group, so a position at the edge of a short read's span is not diluted
      by reads that never covered it.

    Absence of coverage is neither: a ``"-"`` outside a read's own covered span
    contributes to nothing.  Separating that from a deletion is what
    :func:`_covered_span` is for.
    """
    ref = ref_seq.upper()
    n = len(ref)
    agree = [0] * n
    covering = [0] * n

    for row in rows:
        first, last = _covered_span(row)
        if first is None or last is None:
            continue
        for i in range(first, last + 1):
            covering[i] += 1
            base = row[i][0]
            if base != "-" and base.upper() == ref[i]:
                agree[i] += 1

    return [
        (covering[i] - agree[i]) / covering[i] if covering[i] else 0.0
        for i in range(n)
    ]


def flagged_columns(
    rows: Sequence[Row],
    ref_seq: str,
    threshold: float = DEFAULT_FLAG_THRESHOLD,
) -> List[int]:
    """Reference positions where more than *threshold* of covering reads disagree.

    The "worth a look" positions, thresholded out of
    :func:`mismatch_fractions` so the boolean and the magnitude cannot disagree.
    """
    fractions = mismatch_fractions(rows, ref_seq)
    return [i for i, f in enumerate(fractions) if f > threshold]


def _classify(
    variant: Variant, ref_seq: str, focus: Optional[Tuple[int, int]],
) -> Variant:
    """Assign *variant* a consequence and a human effect string.

    *focus* is the ``[start, end)`` reading frame.  A variant outside it, or a
    view with no frame at all, is ``"noncoding"`` with no effect text: the
    position is real but there is nothing to say about a protein.
    """
    if focus is None:
        return variant
    start, end = focus
    if not start <= variant.pos < end:
        return replace(variant, consequence="noncoding")

    if variant.kind in ("del", "ins"):
        size = len(variant.ref) if variant.kind == "del" else len(variant.alt)
        residue = (variant.pos - start) // 3 + 1
        if size % 3:
            return replace(
                variant, consequence="frameshift",
                effect=f"frameshift at {residue}",
            )
        verb = "Δ" if variant.kind == "del" else "+"
        return replace(
            variant, consequence="inframe_indel",
            effect=f"{verb}{size // 3} aa at {residue}",
        )

    offset = variant.pos - start
    codon_start = start + (offset // 3) * 3
    if codon_start + 3 > end:
        return replace(variant, consequence="noncoding")

    ref_codon = ref_seq[codon_start:codon_start + 3].upper()
    within = offset % 3
    alt_codon = ref_codon[:within] + variant.alt + ref_codon[within + 1:]
    was, now = translate_codon(ref_codon), translate_codon(alt_codon)
    residue = offset // 3 + 1

    if was == now:
        return replace(variant, consequence="silent", effect=f"silent ({was}{residue})")
    if now == "*":
        return replace(variant, consequence="nonsense", effect=f"{was}{residue}*")
    return replace(variant, consequence="missense", effect=f"{was}{residue}{now}")


@dataclass
class SummaryView:
    """Everything one summarized page shows.

    The reference is held once, here, rather than once per group: a summary
    draws one coordinate system, and every group's arrays are stated against it.

    ``focus`` is the reading frame variants are classified against, as a
    ``[start, end)`` pair on the reference.  It is derived from the pileup's
    ``flanks`` and is None when no region was singled out or when the view asked
    for no translation, in which case every variant is reported as noncoding.
    """

    title: str
    ref_seq: str
    groups: List[GroupSummary] = field(default_factory=list)
    total_reads: int = 0
    highlight_ids: List[str] = field(default_factory=list)
    highlight_label: str = "Highlighted"
    focus: Optional[Tuple[int, int]] = None
    features: List[Feature] = field(default_factory=list)
    theme: Theme = field(default_factory=Theme)
    #: The view this was reduced from, when there was one.  A reference to it, not
    #: a copy — kept because base-resolution detail needs the reads themselves,
    #: which a reduction deliberately does not carry.  None when a summary was
    #: assembled directly, and then a page draws no such detail.
    source: Optional[PileupView] = None

    @property
    def ref_len(self) -> int:
        return len(self.ref_seq)

    @classmethod
    def from_view(
        cls,
        view: PileupView,
        *,
        min_fraction: float = DEFAULT_MIN_FRACTION,
        min_count: int = DEFAULT_MIN_COUNT,
        min_depth: int = DEFAULT_MIN_DEPTH,
        insertions: Optional[Mapping[str, Mapping[int, Mapping[str, int]]]] = None,
    ) -> "SummaryView":
        """Reduce an existing :class:`~seqviewer.pileup.PileupView`.

        This is the whole of the plumbing: anything that can draw a pileup can
        draw its summary, with the same features, the same focus region, and the
        same theme.

        *insertions* is keyed by group name, so a caller holding sidecar
        evidence can pass it for the groups it has it for and omit the rest.

        Groups may be stated against references of different lengths — the
        pileup view allows it — so the longest is taken as the page's coordinate
        system and shorter groups simply stop early.
        """
        if not view.groups:
            raise ValueError("a summary needs at least one group")

        widest = max(view.groups, key=lambda g: len(g.ref_seq))
        ref_seq = widest.ref_seq
        ref_len = len(ref_seq)

        focus: Optional[Tuple[int, int]] = None
        if view.flanks and view.translate:
            five, three = view.flanks
            start, end = five, ref_len - three
            # A frame must hold at least one whole codon to say anything.
            if end - start >= 3:
                focus = (start, end)

        by_group = insertions or {}
        groups = []
        for group in view.groups:
            summary = summarize_group(
                group,
                min_fraction=min_fraction,
                min_count=min_count,
                min_depth=min_depth,
                insertions=by_group.get(group.name),
            )
            summary.variants = [
                _classify(v, group.ref_seq, focus) for v in summary.variants
            ]
            groups.append(summary)

        return cls(
            title=view.title,
            ref_seq=ref_seq,
            groups=groups,
            total_reads=view.total_reads,
            highlight_ids=list(view.highlight_ids),
            highlight_label=view.highlight_label,
            focus=focus,
            features=list(view.features),
            theme=view.theme,
            source=view,
        )
