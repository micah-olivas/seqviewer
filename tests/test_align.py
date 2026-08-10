"""Alignment-layer tests.

These build a SAM by hand rather than shelling out, so they exercise
``grid_from_alignment`` without minimap2 or samtools.  ``grid_from_reads``
needs both binaries and is skipped when they are absent.
"""

import shutil

import pytest

pytest.importorskip("pysam")

from seqviewer import grid_from_alignment, reads_from_alignment  # noqa: E402
from seqviewer.align import Read, grid_from_reads  # noqa: E402

REF = "ACGTACGTACGTACGTACGT"  # 20 bp, midpoint 10


def _sam(tmp_path, reads):
    """Write a SAM. Each read is (name, pos_1based, cigar, seq, flag)."""
    lines = ["@HD\tVN:1.6\tSO:unsorted", f"@SQ\tSN:ref\tLN:{len(REF)}"]
    for name, pos, cigar, seq, flag in reads:
        qual = "I" * len(seq)
        lines.append(f"{name}\t{flag}\tref\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}")
    path = tmp_path / "reads.sam"
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_perfect_read_is_all_matches(tmp_path):
    sam = _sam(tmp_path, [("r1", 1, "20M", REF, 0)])
    rows = grid_from_alignment(sam, REF)
    assert len(rows) == 1
    assert all(is_match for _, is_match in rows[0])
    assert "".join(base for base, _ in rows[0]) == REF


def test_mismatch_is_flagged_at_the_right_position(tmp_path):
    mutated = REF[:5] + "T" + REF[6:]      # ref has C at index 5
    sam = _sam(tmp_path, [("r1", 1, "20M", mutated, 0)])
    row = grid_from_alignment(sam, REF)[0]
    assert row[5] == ("T", False)
    assert all(is_match for i, (_, is_match) in enumerate(row) if i != 5)


def test_deletion_becomes_a_gap(tmp_path):
    seq = REF[:5] + REF[6:]                # 19 bp, one base deleted
    sam = _sam(tmp_path, [("r1", 1, "5M1D14M", seq, 0)])
    row = grid_from_alignment(sam, REF)[0]
    assert row[5] == ("-", True)


def test_rows_are_the_reference_width_even_for_partial_coverage(tmp_path):
    """A read covering the middle still occupies the full row; the rest is gaps."""
    sam = _sam(tmp_path, [("r1", 6, "10M", REF[5:15], 0)])
    row = grid_from_alignment(sam, REF)[0]
    assert len(row) == len(REF)
    assert row[0] == ("-", True)
    assert row[5][1] is True


def test_reads_not_spanning_the_midpoint_are_dropped(tmp_path):
    """Concatemer fragments cover one flank and stop; only full reads cross."""
    sam = _sam(tmp_path, [
        ("full", 1, "20M", REF, 0),
        ("five_prime", 1, "8M", REF[:8], 0),
        ("three_prime", 13, "8M", REF[12:], 0),
    ])
    assert len(grid_from_alignment(sam, REF)) == 1


def test_the_midpoint_filter_can_be_disabled(tmp_path):
    sam = _sam(tmp_path, [
        ("full", 1, "20M", REF, 0),
        ("short", 1, "8M", REF[:8], 0),
    ])
    assert len(grid_from_alignment(sam, REF, min_overlap_pos=0)) == 2


def test_unmapped_and_secondary_reads_are_skipped(tmp_path):
    sam = _sam(tmp_path, [
        ("mapped", 1, "20M", REF, 0),
        ("unmapped", 0, "*", REF, 4),
        ("secondary", 1, "20M", REF, 256),
    ])
    assert len(grid_from_alignment(sam, REF)) == 1


def test_read_names_restrict_the_grid(tmp_path):
    sam = _sam(tmp_path, [("r1", 1, "20M", REF, 0), ("r2", 1, "20M", REF, 0)])
    assert len(grid_from_alignment(sam, REF)) == 2
    assert len(grid_from_alignment(sam, REF, read_names={"r1"})) == 1
    assert len(grid_from_alignment(sam, REF, read_names=set())) == 0


def test_rows_cluster_by_mismatch_pattern(tmp_path):
    """Cleanest reads first, and identical patterns adjacent."""
    mutated = REF[:5] + "T" + REF[6:]
    sam = _sam(tmp_path, [
        ("mut1", 1, "20M", mutated, 0),
        ("clean", 1, "20M", REF, 0),
        ("mut2", 1, "20M", mutated, 0),
    ])
    rows = grid_from_alignment(sam, REF)
    mismatch_counts = [sum(1 for _, m in row if not m) for row in rows]
    assert mismatch_counts == [0, 1, 1]


def test_missing_file_returns_no_rows_rather_than_raising(tmp_path):
    assert grid_from_alignment(str(tmp_path / "nope.sam"), REF) == []


def test_reads_can_be_recovered_for_realignment(tmp_path):
    sam = _sam(tmp_path, [("r1", 1, "20M", REF, 0), ("r2", 1, "20M", REF, 0)])
    reads = reads_from_alignment(sam)
    assert [r.name for r in reads] == ["r1", "r2"]
    assert reads[0].seq == REF
    assert len(reads[0].qual) == len(REF)
    assert [r.name for r in reads_from_alignment(sam, {"r2"})] == ["r2"]


def test_empty_read_list_short_circuits():
    assert grid_from_reads([], "missing.fasta", REF) == []


@pytest.mark.skipif(
    not (shutil.which("minimap2") and shutil.which("samtools")),
    reason="needs minimap2 and samtools on PATH",
)
def test_grid_from_reads_round_trip(tmp_path):
    long_ref = REF * 5   # minimap2 needs something longer than 20 bp to seed on
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(f">ref\n{long_ref}\n")
    rows = grid_from_reads([Read("r1", long_ref, None)], str(fasta), long_ref)
    assert len(rows) == 1
    assert all(is_match for _, is_match in rows[0])
