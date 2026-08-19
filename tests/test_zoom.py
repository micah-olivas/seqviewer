"""Base-resolution windows: register, bounds, and what a window admits it hid.

The register assertions are the point of the module — a codon must sit over the
three bases it is translated from — so they are checked against the emitted
geometry rather than trusted.
"""

import re
import xml.etree.ElementTree as ET

import pytest

from seqviewer.zoom import (
    CELL_W, DEFAULT_COLUMNS, Window, window_bounds, window_css, window_svg,
)

# frame 3..18: ATG GTG AAA CCC TTT = M V K P F
REF = "GGG" + "ATGGTGAAACCCTTT" + "AAA"
FRAME = (3, 18)


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


def _rows(*specs):
    return [_row(REF, s) for s in specs]


def _brackets(svg):
    """Left and right x of every codon bracket, from its path."""
    lefts = [float(m) - 1 for m in re.findall(r'd="M(-?[\d.]+),', svg)]
    rights = [float(m) + 1 for m in re.findall(r"H(-?[\d.]+)", svg)]
    return list(zip(lefts, rights))


def _residues(svg):
    return re.findall(r'class="svz-aa-(?:same|diff)"[^>]*>([^<]+)<', svg)


# --- Bounds ---------------------------------------------------------------

def test_a_window_is_centred_on_its_position():
    assert window_bounds(50, 1000, 20) == (40, 60)


def test_a_window_near_the_start_shifts_rather_than_shrinks():
    """A variant at base 2 still deserves a full window of context."""
    assert window_bounds(2, 1000, 20) == (0, 20)


def test_a_window_near_the_end_shifts_rather_than_shrinks():
    assert window_bounds(998, 1000, 20) == (980, 1000)


def test_a_window_wider_than_the_reference_is_the_reference():
    assert window_bounds(5, 10, 40) == (0, 10)


def test_the_default_window_is_the_documented_size():
    start, end = window_bounds(500, 1000)
    assert end - start == DEFAULT_COLUMNS


# --- Register: the whole point -------------------------------------------

@pytest.mark.parametrize("start", [0, 3, 4, 5, 6, 7, 12])
def test_a_codon_is_exactly_three_bases_wide_at_every_window_offset(start):
    end = min(len(REF), start + 12)
    svg = window_svg(REF, _rows("." * len(REF)), start, end, frame=FRAME).svg
    pairs = _brackets(svg)
    assert pairs
    for left, right in pairs:
        assert right - left == pytest.approx(3 * CELL_W)


@pytest.mark.parametrize("start", [0, 3, 4, 5, 6, 7, 12])
def test_every_codon_starts_on_the_reading_frame(start):
    """A window opened mid-codon must not re-phase the frame to its own edge."""
    end = min(len(REF), start + 12)
    svg = window_svg(REF, _rows("." * len(REF)), start, end, frame=FRAME).svg
    for left, _ in _brackets(svg):
        base = left / CELL_W + start
        assert base == pytest.approx(round(base))
        assert (round(base) - FRAME[0]) % 3 == 0


def test_a_codon_sits_over_its_own_three_bases():
    """The bracket for codon k covers exactly the bases the residue translates."""
    window = window_svg(REF, _rows("." * len(REF)), 0, 18, frame=FRAME)
    first_left = _brackets(window.svg)[0][0]
    # Codon 0 of the frame starts at reference base 3, and the window starts at 0.
    assert first_left == pytest.approx(3 * CELL_W)


def test_residues_are_the_translation_of_the_reference_when_reads_agree():
    window = window_svg(REF, _rows("." * len(REF)), 0, len(REF), frame=FRAME)
    assert _residues(window.svg) == list("MVKPF")


def test_a_codon_straddling_the_edge_keeps_its_bracket_but_loses_its_letter():
    """The cut bracket says "continues past here"; an invisible letter says nothing."""
    window = window_svg(REF, _rows("." * len(REF)), 5, 17, frame=FRAME)
    assert len(_brackets(window.svg)) > len(_residues(window.svg))


def test_no_frame_means_no_codons_at_all():
    window = window_svg(REF, _rows("." * len(REF)), 0, 18, frame=None)
    assert _brackets(window.svg) == []
    assert window.residues == 0


def test_a_frame_too_short_for_a_codon_draws_none():
    window = window_svg(REF, _rows("." * len(REF)), 0, 18, frame=(3, 5))
    assert window.residues == 0


def test_a_changed_codon_is_marked_as_differing():
    # ATG -> GTG at the first base of codon 0.
    spec = "..." + "G" + "." * 14 + "..."
    window = window_svg(REF, _rows(spec, spec), 0, 18, frame=FRAME)
    assert "svz-aa-diff" in window.svg
    assert _residues(window.svg)[0] == "V"


# --- Deletions versus no coverage ----------------------------------------

def test_a_deletion_in_the_consensus_is_drawn_as_a_dash_not_a_letter():
    spec = "..." + "---" + "." * 12 + "..."
    window = window_svg(REF, _rows(spec, spec), 0, 18, frame=FRAME)
    assert 'class="svz-del"' in window.svg


def test_a_deletion_inside_a_read_is_drawn_but_missing_coverage_is_not():
    """The same distinction the reducer makes, made visible."""
    deleted = "..." + "-" + "." * 14 + "..."      # gap flanked by called bases
    truncated = "-" * 6 + "." * (len(REF) - 6)     # gap running off the 5' end
    with_del = window_svg(REF, _rows(deleted), 0, 18, frame=FRAME).svg
    with_gap = window_svg(REF, _rows(truncated), 0, 18, frame=FRAME).svg
    assert 'class="svz-rdel"' in with_del
    assert 'class="svz-rdel"' not in with_gap


# --- Reads ---------------------------------------------------------------

def test_reads_disagreeing_in_the_window_are_drawn_first():
    """An inspector is opened because something is wrong here."""
    clean = "." * len(REF)
    dirty = "..." + "G" + "." * 14 + "..."
    window = window_svg(REF, _rows(clean, clean, dirty), 0, 18, frame=FRAME,
                        max_read_rows=1)
    # Only one read is drawn, and it must be the informative one.
    assert window.rows_shown == 1
    assert "svz-b" in window.svg


def test_a_window_says_how_many_reads_it_left_out():
    rows = _rows(*["." * len(REF)] * 20)
    window = window_svg(REF, rows, 0, 18, frame=FRAME, max_read_rows=6)
    assert window.rows_shown == 6
    assert window.rows_hidden == 14


def test_a_window_with_no_reads_still_draws_the_reference():
    window = window_svg(REF, [], 0, 18, frame=FRAME)
    assert window.rows_shown == 0
    assert "svz-ref" in window.svg


# --- The element ---------------------------------------------------------

def test_the_window_is_well_formed_xml():
    window = window_svg(REF, _rows("." * len(REF)), 0, 18, frame=FRAME)
    ET.fromstring(window.svg)


def test_the_window_avoids_dominant_baseline():
    """Long broken in WebKit; a dy in ems is exact in both engines."""
    window = window_svg(REF, _rows("." * len(REF)), 0, 18, frame=FRAME)
    assert "dominant-baseline" not in window.svg
    assert 'dy="0.35em"' in window.svg


def test_the_window_is_exactly_as_wide_as_its_columns():
    window = window_svg(REF, _rows("." * len(REF)), 2, 14)
    assert window.columns == 12
    assert window.width == pytest.approx(12 * CELL_W)


def test_an_empty_span_draws_nothing():
    window = window_svg(REF, _rows("." * len(REF)), 5, 5)
    assert window.svg == ""
    assert window.columns == 0


def test_a_span_beyond_the_reference_is_clipped():
    window = window_svg(REF, _rows("." * len(REF)), 15, 999)
    assert window.end == len(REF)


def test_labels_are_escaped():
    window = window_svg(REF, _rows("." * len(REF)), 0, 12,
                        label='<script>alert(1)</script>')
    assert "<script>" not in window.svg
    ET.fromstring(window.svg)


def test_the_prefix_keeps_two_windows_from_colliding():
    a = window_svg(REF, _rows("." * len(REF)), 0, 12, prefix="one")
    b = window_svg(REF, _rows("." * len(REF)), 0, 12, prefix="two")
    assert "one-ref" in a.svg and "two-ref" not in a.svg
    assert "two-ref" in b.svg


def test_css_follows_the_hosts_token_prefix():
    assert "var(--app-a)" in window_css(token_prefix="app")
    assert "var(--cv-a)" in window_css()


def test_a_window_reports_its_own_geometry():
    window = window_svg(REF, _rows("." * len(REF), "." * len(REF)), 0, 12,
                        frame=FRAME)
    assert isinstance(window, Window)
    assert window.height > 0
    assert window.width > 0
