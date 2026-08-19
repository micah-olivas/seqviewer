"""Render a pileup from synthetic data.

Useful for looking at the viewer without reads, an aligner, or a reference on
disk, and as a smoke test that the renderer still produces a page.

    python -m seqviewer.demo [output.html]
    python -m seqviewer.demo [output.html] --summary

``--summary`` writes the summarized view instead, from a second construct built
to exercise it: several annotated features, and subpopulations carrying a
missense substitution, an in-frame deletion, a frameshift, and a coverage
dropout, which the pileup demo has none of.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .codon import translate_codon
from .construct import Feature, Reference
from .grid import Cell, Row
from .pileup import PileupGroup, PileupView
from .render import render

BASES = "ACGT"


def _random_seq(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(BASES) for _ in range(length))


def _other_base(exclude: str, rng: random.Random) -> str:
    return rng.choice([b for b in BASES if b != exclude])


def synthetic_rows(
    ref_seq: str,
    n_rows: int,
    rng: random.Random,
    mismatch_rate: float = 0.02,
    gap_rate: float = 0.01,
    mutations: Optional[Dict[int, str]] = None,
) -> List[Row]:
    """Synthesize *n_rows* reads against *ref_seq*.

    *mutations* maps a reference position to a base forced into every read
    there, which is how a real subpopulation shows up: one column that
    disagrees in every read rather than scattered sequencing noise.
    """
    mutations = mutations or {}
    rows: List[Row] = []
    for _ in range(n_rows):
        row: List[Cell] = []
        for i, ref_base in enumerate(ref_seq):
            ref_upper = ref_base.upper()
            roll = rng.random()
            if i in mutations:
                row.append((mutations[i], False))
            elif roll < gap_rate:
                row.append(("-", False))
            elif roll < gap_rate + mismatch_rate:
                row.append((_other_base(ref_upper, rng), False))
            else:
                row.append((ref_upper, True))
        rows.append(row)
    return rows


def build_view(seed: int = 42) -> PileupView:
    """Build a two-group demo view: a wild-type majority and a mutant minority."""
    rng = random.Random(seed)

    flank_5p, insert_len, flank_3p = 100, 800, 100
    ref_seq = (
        _random_seq(flank_5p, rng)
        + _random_seq(insert_len, rng)
        + _random_seq(flank_3p, rng)
    )
    reference = Reference(
        seq=ref_seq,
        name="pUC19",
        topology="circular",
        features=[
            Feature("insert", flank_5p, flank_5p + insert_len, strand=1, label="CDS"),
        ],
    )

    wt = synthetic_rows(ref_seq, 10, rng)

    # One column every mutant read disagrees on, inside the insert.
    mutant_pos = flank_5p + 20
    mutant = synthetic_rows(
        ref_seq, 8, rng,
        mismatch_rate=0.03,
        mutations={mutant_pos: _other_base(ref_seq[mutant_pos].upper(), rng)},
    )

    groups = [
        PileupGroup("pUC19-WT", ref_seq, wt, n_reads=65, fraction=0.65,
                    status="Perfect Match", highlighted=True),
        PileupGroup("pUC19-K44A", ref_seq, mutant, n_reads=35, fraction=0.35,
                    status="Mismatch", highlighted=False),
    ]

    return PileupView.from_reference(
        title="Pileup: pUC19 demo",
        reference=reference,
        groups=groups,
        total_reads=100,
        highlight_ids=["pUC19-WT"],
        highlight_label="Recoverable",
    )


def _planted_rows(
    ref_seq: str,
    n_rows: int,
    rng: random.Random,
    snvs: Optional[Dict[int, str]] = None,
    deletions: Sequence[Tuple[int, int]] = (),
    carriers: float = 1.0,
    covers: Optional[Tuple[int, int]] = None,
    error_rate: float = 0.004,
) -> List[Row]:
    """Reads carrying planted variants over a low background error rate.

    *carriers* is the share of reads carrying the planted alleles, which is what
    makes a subpopulation rather than a fixed difference.  *covers* restricts
    every read to a window, so a group can show a coverage dropout — the one
    thing a summary should never average away.
    """
    snvs = snvs or {}
    first, last = covers or (0, len(ref_seq))
    rows: List[Row] = []
    for index in range(n_rows):
        carrier = index < round(n_rows * carriers)
        deleted = set()
        if carrier:
            for start, length in deletions:
                deleted.update(range(start, start + length))
        row: List[Cell] = []
        for i, ref_base in enumerate(ref_seq):
            ref_upper = ref_base.upper()
            if not first <= i < last or i in deleted:
                row.append(("-", True))
            elif carrier and i in snvs:
                row.append((snvs[i], False))
            elif rng.random() < error_rate:
                row.append((_other_base(ref_upper, rng), False))
            else:
                row.append((ref_upper, True))
        rows.append(row)
    return rows


def _missense_site(ref_seq: str, frame_start: int, near: int) -> Tuple[int, str]:
    """The first position at or after *near* where one base change is a missense.

    Searched rather than assumed: a demo of a summarized view should show a
    missense call, and not every site can produce one.  Third codon positions in
    particular are often wholly degenerate — every base at the end of ``GT_``
    reads valine — so picking a position blind gives a silent change instead.
    """
    for pos in range(near, len(ref_seq) - 3):
        offset = pos - frame_start
        codon_start = frame_start + (offset // 3) * 3
        codon = ref_seq[codon_start:codon_start + 3].upper()
        within = offset % 3
        was = translate_codon(codon)
        for base in BASES:
            if base == codon[within]:
                continue
            now = translate_codon(codon[:within] + base + codon[within + 1:])
            if now != was and now != "*":
                return pos, base
    raise ValueError("no missense site found in the reading frame")


def build_summary_view(seed: int = 7) -> PileupView:
    """A construct annotated enough, and broken enough, to exercise the summary.

    Three clones over one plasmid: one clean, one carrying a missense change and
    an in-frame deletion, and one carrying a frameshift and a premature stop
    over a read set that does not reach the 3' end.
    """
    rng = random.Random(seed)

    flank_5p, insert_len, flank_3p = 420, 900, 480
    ref_seq = (
        _random_seq(flank_5p, rng)
        + _random_seq(insert_len, rng)
        + _random_seq(flank_3p, rng)
    )
    insert_start = flank_5p
    reference = Reference(
        seq=ref_seq,
        name="pSV-demo",
        topology="circular",
        features=[
            Feature("promoter", 250, 400, strand=1, label="T7"),
            Feature("insert", insert_start, insert_start + insert_len,
                    strand=1, label="target CDS"),
            Feature("CDS", 1400, 1700, strand=-1, label="AmpR"),
            Feature("rep_origin", 1750, 1795, label="ori"),
        ],
    )

    # A site inside the insert where one substitution really is a missense.
    missense_at, missense_to = _missense_site(ref_seq, insert_start,
                                              insert_start + 147)
    clean = _planted_rows(ref_seq, 24, rng)
    inframe = _planted_rows(
        ref_seq, 18, rng,
        snvs={missense_at: missense_to},
        deletions=[(insert_start + 300, 3)],
    )
    broken = _planted_rows(
        ref_seq, 14, rng,
        deletions=[(insert_start + 501, 1)],
        carriers=0.8,
        covers=(0, 1500),
    )

    groups = [
        PileupGroup("clone-A", ref_seq, clean, n_reads=24, fraction=0.43,
                    status="Perfect Match", highlighted=True),
        PileupGroup("clone-B", ref_seq, inframe, n_reads=18, fraction=0.32),
        PileupGroup("clone-C", ref_seq, broken, n_reads=14, fraction=0.25),
    ]

    return PileupView.from_reference(
        title="Summary: pSV-demo",
        reference=reference,
        groups=groups,
        total_reads=56,
        highlight_ids=["clone-A"],
        highlight_label="Recoverable",
    )


def main(argv: Optional[List[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    summary = "--summary" in argv
    argv = [a for a in argv if a != "--summary"]

    if summary:
        from .render_summary import render_summary
        from .summary import SummaryView

        out = Path(argv[0]) if argv else Path("demo_summary.html")
        out.write_text(render_summary(SummaryView.from_view(build_summary_view())))
    else:
        out = Path(argv[0]) if argv else Path("demo_pileup.html")
        out.write_text(render(build_view()))
    print(f"Wrote {out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
