"""Turning reads into the grid a pileup is drawn from.

Two entry points, with different requirements:

``grid_from_alignment``
    Reads an existing SAM or BAM.  Needs :mod:`pysam` and nothing else.

``grid_from_reads``
    Aligns reads first, then reads the result.  Needs :mod:`pysam` plus
    ``minimap2`` and ``samtools`` on PATH.

Callers that already have a grid never import this module, and so never pay for
either.  Nothing here implements alignment; ``grid_from_reads`` invokes
minimap2 and reads what comes back.
"""

from __future__ import annotations

import logging
import os
import subprocess
import tempfile
from typing import Iterable, List, NamedTuple, Optional, Sequence, Set

from .model import Row

logger = logging.getLogger(__name__)

__all__ = ["Read", "grid_from_alignment", "grid_from_reads", "reads_from_alignment"]


class Read(NamedTuple):
    """One read, in the neutral form both entry points accept."""

    name: str
    seq: str
    qual: Optional[str] = None


def _require_pysam():
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - exercised only without pysam
        raise ImportError(
            "reading alignments needs pysam: pip install 'seqviewer[sam]'"
        ) from exc
    return pysam


def _cluster(rows: List[Row]) -> List[Row]:
    """Order rows so reads sharing a mismatch pattern sit together.

    The key writes ``.`` for a match or a gap and the base itself for a
    mismatch, so identical patterns sort adjacent and the cleanest reads — all
    dots — come first.  Subpopulations then read as blocks rather than noise.
    """
    return sorted(
        rows,
        key=lambda row: "".join(
            "." if is_match or base == "-" else base.upper() for base, is_match in row
        ),
    )


def _spans(read, midpoint: int) -> bool:
    """Whether *read* crosses *midpoint* from both sides.

    Concatemer fragments cover one flank and stop; only a full-length read
    crosses the middle, so this is what keeps partial reads out of the grid.
    Pass ``midpoint=0`` to accept everything.
    """
    if not midpoint:
        return True
    start, end = read.reference_start, read.reference_end
    return start is not None and end is not None and end > midpoint > start


def _row_for(read, ref_seq: str, ref_len: int) -> Row:
    """Project one aligned read onto reference coordinates."""
    row: List = [("-", True)] * ref_len
    for qpos, rpos in read.get_aligned_pairs():
        if rpos is None or rpos >= ref_len:
            continue
        if qpos is None:
            row[rpos] = ("-", True)          # deletion
        else:
            base = read.query_sequence[qpos]
            row[rpos] = (base, base.upper() == ref_seq[rpos].upper())
    return row


def _parse(path: str, ref_seq: str, read_names: Optional[Set[str]],
           min_overlap_pos: int) -> List[Row]:
    pysam = _require_pysam()
    ref_len = len(ref_seq)
    midpoint = ref_len // 2 if min_overlap_pos < 0 else min_overlap_pos

    rows: List[Row] = []
    try:
        with pysam.AlignmentFile(path, check_sq=False) as handle:
            for read in handle.fetch(until_eof=True):
                if read_names is not None and read.query_name not in read_names:
                    continue
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if not _spans(read, midpoint):
                    continue
                rows.append(_row_for(read, ref_seq, ref_len))
    except Exception as exc:
        logger.warning("could not read alignment %s: %s", path, exc)
        return []

    return _cluster(rows)


def grid_from_alignment(
    path: str,
    ref_seq: str,
    read_names: Optional[Iterable[str]] = None,
    min_overlap_pos: int = -1,
) -> List[Row]:
    """Build a pileup grid from an existing SAM or BAM.

    Args:
        path: The SAM or BAM to read.  Reads are taken as aligned, so whatever
            reference they were aligned against is the one *ref_seq* must be.
        ref_seq: The reference sequence, which sets the grid's width.
        read_names: Restrict to these read names.  None takes every read.
        min_overlap_pos: Drop reads that do not cross this reference position.
            Negative means the reference midpoint; 0 disables the filter.

    Returns:
        One row per surviving read, clustered by mismatch pattern.
    """
    names = set(read_names) if read_names is not None else None
    return _parse(path, ref_seq, names, min_overlap_pos)


def reads_from_alignment(path: str, read_names: Optional[Iterable[str]] = None) -> List[Read]:
    """Pull reads back out of a SAM or BAM as sequences.

    Use this to re-align a subset against a different reference: an alignment
    carries the reference it was made against, so reads assigned to one variant
    but stored against another show systematic mismatches until they are moved.
    """
    pysam = _require_pysam()
    names = set(read_names) if read_names is not None else None

    out: List[Read] = []
    try:
        with pysam.AlignmentFile(path, check_sq=False) as handle:
            for read in handle.fetch(until_eof=True):
                if names is not None and read.query_name not in names:
                    continue
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                seq = read.query_sequence or ""
                if not seq:
                    continue
                if read.query_qualities is not None:
                    qual = "".join(chr(q + 33) for q in read.query_qualities)
                else:
                    qual = "I" * len(seq)
                out.append(Read(read.query_name, seq, qual))
    except Exception as exc:
        logger.warning("could not read alignment %s: %s", path, exc)
    return out


def grid_from_reads(
    reads: Sequence[Read],
    ref_fasta: str,
    ref_seq: str,
    minimap2: str = "minimap2",
    samtools: str = "samtools",
    ref_index: Optional[str] = None,
    min_overlap_pos: int = -1,
) -> List[Row]:
    """Align *reads* to a reference and build a pileup grid from the result.

    Args:
        reads: The reads to align.
        ref_fasta: FASTA to align against.
        ref_seq: That FASTA's sequence, which sets the grid's width.
        minimap2: minimap2 executable.
        samtools: samtools executable.
        ref_index: A prebuilt ``.mmi`` to use in place of *ref_fasta*.
        min_overlap_pos: As in :func:`grid_from_alignment`.

    Returns:
        One row per aligned read, clustered by mismatch pattern.  Empty if
        nothing aligned or an external tool failed.
    """
    if not reads:
        return []

    with tempfile.TemporaryDirectory() as tmp:
        fastq = os.path.join(tmp, "reads.fastq")
        bam = os.path.join(tmp, "aligned.bam")

        with open(fastq, "w") as handle:
            for read in reads:
                qual = read.qual or "I" * len(read.seq)
                handle.write(f"@{read.name}\n{read.seq}\n+\n{qual}\n")

        try:
            aligner = subprocess.Popen(
                [minimap2, "-a", "--MD", "--secondary=no", ref_index or ref_fasta, fastq],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            )
            subprocess.run([samtools, "sort", "-o", bam], stdin=aligner.stdout,
                           stderr=subprocess.DEVNULL, check=False)
            aligner.wait()
            subprocess.run([samtools, "index", bam],
                           stderr=subprocess.DEVNULL, check=False)
        except Exception as exc:
            logger.warning("alignment failed: %s", exc)
            return []

        rows = _parse(bam, ref_seq, None, min_overlap_pos)

    if not rows:
        logger.warning("alignment of %d reads against %s produced no rows",
                       len(reads), ref_fasta)
    return rows
