"""The plate map: reading a well, choosing a format, and drawing it."""
import xml.etree.ElementTree as ET

import pytest

from seqviewer.cli import build_parser
from seqviewer.pileup import PileupGroup, PileupView
from seqviewer.plate import FORMATS, Well, parse_well, plate_svg
from seqviewer.render import render
from seqviewer.render_summary import render_summary
from seqviewer.summary import SummaryView


# --- reading a well --------------------------------------------------------

@pytest.mark.parametrize("text, row, col", [
    ("A1", 0, 0), ("a1", 0, 0), ("A01", 0, 0), (" H12 ", 7, 11),
    ("P24", 15, 23), ("b7", 1, 6),
])
def test_wells_are_read_in_any_case_with_or_without_padding(text, row, col):
    well = parse_well(text)
    assert (well.row, well.col) == (row, col)


def test_the_format_is_the_smallest_the_well_fits():
    """Anything inside A1-H12 is read as 96; past either edge is 384."""
    assert parse_well("H12").plate == 96
    assert parse_well("I1").plate == 384       # row past H
    assert parse_well("A13").plate == 384      # column past 12
    assert parse_well("P24").plate == 384


def test_an_explicit_format_is_kept_even_when_a_smaller_one_fits():
    """A1 on a 384 plate is still A1 on a 384 plate."""
    assert parse_well("A1", plate=384).plate == 384


def test_the_label_round_trips():
    for text in ("A1", "H12", "I13", "P24"):
        assert parse_well(text).label == text


@pytest.mark.parametrize("text", ["", "1A", "AA1", "A0", "A", "12", "Q1", "A25"])
def test_things_that_are_not_wells_are_refused(text):
    with pytest.raises(ValueError):
        parse_well(text)


def test_a_well_off_the_named_plate_is_refused_with_the_plate_edge():
    with pytest.raises(ValueError, match="A1 to H12"):
        parse_well("I1", plate=96)


def test_an_unknown_format_is_refused():
    with pytest.raises(ValueError, match="no 1536-well plate"):
        parse_well("A1", plate=1536)


# --- drawing it -----------------------------------------------------------

def _wells_drawn(svg: str) -> int:
    """Wells in the empty-well path plus the one filled circle.

    The SVG is inline, so it declares no namespace and its tags parse bare.
    Each empty well starts with one absolute ``M`` move.
    """
    root = ET.fromstring(svg)
    empty = root.find(".//*[@class='sv-plate-well']")
    return empty.get("d").count("M") + len(root.findall(".//circle"))


@pytest.mark.parametrize("fmt", sorted(FORMATS))
def test_every_well_of_the_plate_is_drawn(fmt):
    """A 384 map missing its interior would read as a 96."""
    rows, cols = FORMATS[fmt]
    svg = plate_svg(Well(0, 0, fmt))
    assert _wells_drawn(svg) == rows * cols


def test_the_map_is_well_formed_svg_and_names_the_well():
    svg = plate_svg(parse_well("B7"))
    root = ET.fromstring(svg)
    assert root.tag.endswith("svg")
    assert "Well B7 of 96" in svg
    assert ">B7<" in svg


def test_the_filled_well_is_where_the_grid_puts_it():
    """The marked circle is one cell right and one cell down of A1's."""
    a1 = ET.fromstring(plate_svg(Well(0, 0, 96))).find(".//circle")
    b2 = ET.fromstring(plate_svg(Well(1, 1, 96))).find(".//circle")
    dx = float(b2.get("cx")) - float(a1.get("cx"))
    dy = float(b2.get("cy")) - float(a1.get("cy"))
    assert dx == pytest.approx(dy)
    assert dx > 0


def test_both_formats_come_out_about_the_same_size():
    """384 is packed tighter so it does not take twice the corner."""
    def size(fmt):
        root = ET.fromstring(plate_svg(Well(0, 0, fmt)))
        return float(root.get("width")), float(root.get("height"))
    w96, h96 = size(96)
    w384, h384 = size(384)
    assert 0.8 < w384 / w96 < 1.25
    assert 0.8 < h384 / h96 < 1.25


# --- on the pages -----------------------------------------------------------

def _view(**kwargs):
    ref = "ACGT" * 25
    rows = [[(b, True) for b in ref], [("A", b == "A") for b in ref]]
    kwargs.setdefault("title", "well page")
    kwargs.setdefault("groups", [PileupGroup("g", ref, rows, n_reads=2,
                                             fraction=1.0)])
    return PileupView(**kwargs)


def test_no_well_no_map():
    """The stylesheet always carries the rules; the markup appears only when
    there is a well to place."""
    html = render(_view())
    assert 'class="sv-plate-fixed"' not in html
    assert 'class="sv-plated"' not in html
    summary = render_summary(SummaryView.from_view(_view(), min_depth=1))
    assert 'class="sv-plate-fixed"' not in summary
    assert 'class="sv-plated"' not in summary


def test_the_pileup_page_carries_the_map_and_makes_room_for_it():
    html = render(_view(well=parse_well("C4")))
    assert 'class="sv-plate-fixed"' in html
    assert "Well C4 of 96" in html
    assert '<body class="sv-plated">' in html


def test_the_summary_carries_the_same_well_as_its_pileup():
    """from_view forwards it, so the two pages cannot disagree."""
    view = _view(well=parse_well("O23"))
    summary = SummaryView.from_view(view, min_depth=1)
    assert summary.well == view.well
    html = render_summary(summary)
    assert "Well O23 of 384" in html
    assert '<body class="sv-plated">' in html


def test_the_map_takes_the_corner_opposite_the_toggle():
    """Both are fixed; one is left and one is right, on both pages."""
    view = _view(well=parse_well("A1"))
    for html in (render(view, summary_href="s.html"),
                 render_summary(SummaryView.from_view(view, min_depth=1),
                                pileup_href="p.html")):
        assert "sv-views-fixed" in html and "sv-plate-fixed" in html
        assert "right: 1.5rem" in html


# --- the flag -----------------------------------------------------------

def test_the_cli_accepts_a_well_and_a_format():
    args = build_parser().parse_args(
        ["reads", "ref", "out", "--well", "P24", "--plate", "384"])
    assert (args.well, args.plate) == ("P24", 384)


def test_the_cli_refuses_a_format_it_cannot_draw():
    with pytest.raises(SystemExit):
        build_parser().parse_args(["reads", "ref", "out", "--well", "A1",
                                   "--plate", "1536"])
