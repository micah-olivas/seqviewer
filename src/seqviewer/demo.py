"""Render a pileup from synthetic data.

Useful for looking at the viewer without reads, an aligner, or a reference on
disk, and as a smoke test that the renderer still produces a page.

    python -m seqviewer.demo [output.html]
"""

from __future__ import annotations

import random
import sys
from pathlib import Path
from typing import Dict, List, Optional

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
        top_fraction=0.65,
        highlight_ids=["pUC19-WT"],
        highlight_label="Recoverable",
    )


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    out = Path(argv[0]) if argv else Path("demo_pileup.html")
    out.write_text(render(build_view()))
    print(f"Wrote {out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
