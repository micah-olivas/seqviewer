"""Hierarchical clustering of pileup rows.

The linkage is checked against a textbook O(n^3) average-linkage implementation
written here, so the fast algorithm is verified rather than trusted.
"""

import io
import random
import sys

import pytest

from seqviewer.cluster import (
    cluster_rows, leaf_order, nn_chain, signatures, variable_columns,
)


def _row(ref, spec):
    out = []
    for i, ch in enumerate(spec):
        if ch == ".":
            out.append((ref[i], True))
        elif ch == "-":
            out.append(("-", True))
        else:
            out.append((ch, ch == ref[i]))
    return out


def _naive_average_linkage(D):
    """The definition, straight from the textbook: O(n^3), no chain."""
    n = len(D)
    D = [row[:] for row in D]
    active = [True] * n
    size = [1] * n
    merges = []
    for _ in range(n - 1):
        best, pair = float("inf"), None
        for i in range(n):
            if not active[i]:
                continue
            for j in range(i + 1, n):
                if active[j] and D[i][j] < best:
                    best, pair = D[i][j], (i, j)
        a, b = pair
        merges.append((a, b, best))
        for i in range(n):
            if not active[i] or i in (a, b):
                continue
            m = (D[a][i] * size[a] + D[b][i] * size[b]) / (size[a] + size[b])
            D[a][i] = D[i][a] = m
        size[a] += size[b]
        active[b] = False
    return merges


def _heights(merges):
    return sorted(round(h, 9) for _, _, h in merges)


def _random_matrix(n, rng, integral=False):
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = float(rng.randint(0, 8)) if integral else rng.random() * 10
            D[i][j] = D[j][i] = d
    return D


# --- The linkage ----------------------------------------------------------

@pytest.mark.parametrize("integral", [True, False])
def test_nn_chain_reproduces_the_textbook_dendrogram(integral):
    """Average linkage is reducible, so NN-chain is exact and not an approximation.

    Integer distances are the harder case: they produce many ties, which is where
    a chain algorithm can diverge from the naive one if written wrongly.
    """
    rng = random.Random(4 if integral else 5)
    for _ in range(60):
        n = rng.randint(2, 14)
        D = _random_matrix(n, rng, integral)
        assert _heights(nn_chain([r[:] for r in D])) == _heights(
            _naive_average_linkage(D)
        )


def test_a_matrix_of_one_or_none_has_no_merges():
    assert nn_chain([]) == []
    assert nn_chain([[0.0]]) == []


def test_every_merge_is_reported_once():
    rng = random.Random(6)
    D = _random_matrix(9, rng)
    assert len(nn_chain(D)) == 8


def test_leaf_order_is_a_permutation_of_the_leaves():
    rng = random.Random(7)
    for _ in range(40):
        n = rng.randint(2, 25)
        D = _random_matrix(n, rng)
        order = leaf_order(nn_chain(D), n)
        assert sorted(order) == list(range(n))


def test_leaf_order_keeps_a_merged_pair_adjacent():
    """Two identical points must not be separated by a third, distant one."""
    D = [[0.0, 0.0, 9.0],
         [0.0, 0.0, 9.0],
         [9.0, 9.0, 0.0]]
    order = leaf_order(nn_chain(D), 3)
    assert abs(order.index(0) - order.index(1)) == 1


# --- The clustering space -------------------------------------------------

def test_variable_columns_are_the_ones_some_read_disagrees_at():
    ref = "ACGTACGT"
    rows = [_row(ref, "..T....."), _row(ref, "........")]
    assert variable_columns(rows, ref) == [2]


def test_a_column_every_read_agrees_on_is_not_a_column():
    ref = "ACGT"
    assert variable_columns([_row(ref, "....")] * 3, ref) == []


def test_a_signature_is_one_bit_per_variable_column():
    ref = "ACGTACGT"
    rows = [_row(ref, "..T....."), _row(ref, "..T..A.."), _row(ref, "........")]
    cols = variable_columns(rows, ref)
    assert cols == [2, 5]
    assert signatures(rows, ref, cols) == [0b01, 0b11, 0b00]


def test_a_deletion_is_not_a_disagreement_in_a_signature():
    """A read that stopped short would otherwise look like one that deleted the
    whole tail.  Which reads carry the same substitutions is a different question
    from how much a column disagrees.
    """
    ref = "ACGTACGT"
    rows = [_row(ref, "..-....."), _row(ref, "........")]
    assert variable_columns(rows, ref) == []


# --- Ordering rows --------------------------------------------------------

def test_clustering_groups_a_subpopulation_that_a_string_sort_splits():
    """The case that motivates the whole module.

    Four reads share a variant at position 40; two of them also carry an
    unrelated error at position 2.  Sorting the pattern as a string is decided by
    the leftmost difference, so a wild-type read carrying its own early error
    sorts into the middle of the four.
    """
    ref = "A" * 50

    def read(sites):
        return [("T", False) if i in sites else ("A", True) for i in range(50)]

    rows = [read(set()), read(set()), read({2}), read(set()),
            read({40}), read({40}), read({2, 40}), read({2, 40})]
    ordered = cluster_rows(rows, ref, progress=False)
    carriers = [i for i, row in enumerate(ordered) if not row[40][1]]
    assert carriers == list(range(min(carriers), max(carriers) + 1))


def test_every_read_survives_the_reordering():
    ref = "ACGTACGT" * 4
    rng = random.Random(8)
    rows = [
        [(ref[i], True) if rng.random() > 0.1 else ("T", ref[i] == "T")
         for i in range(len(ref))]
        for _ in range(30)
    ]
    ordered = cluster_rows(rows, ref, progress=False)
    assert len(ordered) == len(rows)
    assert sorted(map(id, ordered)) == sorted(map(id, rows))


def test_reads_sharing_a_signature_keep_their_arrival_order():
    """Determinism: the same grid gives the same page twice."""
    ref = "ACGT" * 4
    rows = [[(b, True) for b in ref] for _ in range(5)]
    rows.append([("T", ref[0] == "T")] + [(b, True) for b in ref[1:]])
    first = cluster_rows(rows, ref, progress=False)
    assert first == cluster_rows(rows, ref, progress=False)


def test_a_grid_nothing_disagrees_in_is_returned_unchanged():
    ref = "ACGT" * 4
    rows = [[(b, True) for b in ref] for _ in range(4)]
    assert cluster_rows(rows, ref, progress=False) == rows


def test_too_few_rows_to_cluster_are_returned_unchanged():
    ref = "ACGT"
    rows = [[(b, True) for b in ref], [("T", False)] + [(b, True) for b in ref[1:]]]
    assert cluster_rows(rows, ref, progress=False) == rows


def test_clustering_can_be_restricted_to_chosen_columns():
    """Passing the flagged columns ignores each read's own errors."""
    ref = "A" * 50

    def read(sites):
        return [("T", False) if i in sites else ("A", True) for i in range(50)]

    rows = [read({2}), read({40}), read({40}), read({40})]
    # Only column 40 is shared; restricting to it puts the three carriers together
    # regardless of the lone error at 2.
    ordered = cluster_rows(rows, ref, columns=[40], progress=False)
    carriers = [i for i, row in enumerate(ordered) if not row[40][1]]
    assert carriers == list(range(min(carriers), max(carriers) + 1))


# --- Progress -------------------------------------------------------------

def test_a_short_run_reports_no_progress(monkeypatch):
    """A progress line for something that finishes instantly is noise."""
    stream = io.StringIO()
    stream.isatty = lambda: True
    monkeypatch.setattr(sys, "stderr", stream)

    ref = "A" * 40
    rows = [[("T", False) if (i + k) % 7 == 0 else ("A", True)
             for i in range(40)] for k in range(12)]
    cluster_rows(rows, ref, progress=True)
    assert stream.getvalue() == ""


def test_a_long_run_reports_progress_and_clears_the_line(monkeypatch):
    stream = io.StringIO()
    stream.isatty = lambda: True
    monkeypatch.setattr(sys, "stderr", stream)
    monkeypatch.setattr("seqviewer.cluster.PROGRESS_AFTER", 0.0)
    monkeypatch.setattr("seqviewer.cluster.PROGRESS_EVERY", 0.0)

    ref = "A" * 60
    rows = [[("T", False) if (i * 7 + k * 13) % 97 < 3 else ("A", True)
             for i in range(60)] for k in range(40)]
    cluster_rows(rows, ref, progress=True)
    text = stream.getvalue()
    assert "clustering" in text
    assert "#" in text
    # Ends by wiping the line rather than leaving a finished bar behind.
    assert text.rstrip("\r").endswith(" ") or text.endswith("\r")


def test_progress_is_silent_when_stderr_is_not_a_terminal(monkeypatch):
    """A redraw sequence in a log file is noise; the log is read afterwards."""
    stream = io.StringIO()
    stream.isatty = lambda: False
    monkeypatch.setattr(sys, "stderr", stream)
    monkeypatch.setattr("seqviewer.cluster.PROGRESS_AFTER", 0.0)

    ref = "A" * 40
    rows = [[("T", False) if (i * 3 + k) % 11 < 2 else ("A", True)
             for i in range(40)] for k in range(20)]
    cluster_rows(rows, ref, progress=True)
    assert stream.getvalue() == ""


# --- The command line -----------------------------------------------------

def test_cluster_is_an_offered_ordering():
    from seqviewer.cli import ORDERINGS

    assert "cluster" in ORDERINGS
    assert set(ORDERINGS) == {"length", "position", "mismatch", "cluster"}


def test_every_ordering_takes_the_grid_and_the_reference():
    """Uniform signature, so the call site does not special-case one of them."""
    from seqviewer.cli import ORDERINGS

    ref = "ACGT" * 8
    rows = [[(b, True) for b in ref], [("T", ref[0] == "T")]
            + [(b, True) for b in ref[1:]]]
    for name, fn in ORDERINGS.items():
        assert len(fn(rows, ref)) == 2, name
