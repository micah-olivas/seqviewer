"""Feature layout and SVG emission.

Pure geometry over Feature records, so none of this needs a browser or Biopython.
"""

import re

from seqviewer.annotate import (
    LANE_GAP, LANE_HEIGHT, cell_width, plan_track, track_style, track_svg,
)
from seqviewer.model import Feature


def _f(start, end, **kw):
    kw.setdefault("type", "CDS")
    kw.setdefault("label", "f")
    return Feature(start=start, end=end, **kw)


def _lanes(plan):
    return {g.feature.label: g.lane for g in plan.glyphs}


# --- The shared column width ---------------------------------------------

def test_cell_width_follows_the_column_count():
    """The canvas reads this same function instead of recomputing the rule, so a
    glyph cannot land a pixel off the column it describes.
    """
    assert cell_width(100) == 4
    assert cell_width(199) == 4
    assert cell_width(200) == 3
    assert cell_width(499) == 3
    assert cell_width(500) == 2
    assert cell_width(5000) == 2


# --- Packing --------------------------------------------------------------

def test_features_that_do_not_overlap_share_a_lane():
    plan = plan_track([_f(0, 100, label="a"), _f(200, 300, label="b")], 1000)
    assert _lanes(plan) == {"a": 0, "b": 0}
    assert plan.lanes == 1


def test_overlapping_features_take_separate_lanes():
    plan = plan_track([_f(0, 500, label="a"), _f(400, 900, label="b")], 1000)
    assert _lanes(plan) == {"a": 0, "b": 1}
    assert plan.lanes == 2


def test_abutting_features_are_separated_by_a_pixel_gap_not_a_base_gap():
    """At 2px per base an adjacent glyph would share an edge and read as one
    feature, so the gap is measured in pixels.
    """
    plan = plan_track([_f(0, 500, label="a"), _f(500, 1000, label="b")], 1000)
    assert plan.lanes == 2


def test_the_longest_feature_of_a_cluster_takes_the_top_lane():
    plan = plan_track([_f(10, 40, label="short"), _f(0, 900, label="long")], 1000)
    assert _lanes(plan)["long"] == 0


def test_features_beyond_the_lane_limit_are_dropped_and_reported():
    """Silently growing the page is worse than saying what was left out."""
    crowd = [_f(0, 900, label=f"f{i}") for i in range(5)]
    plan = plan_track(crowd, 1000, max_lanes=2)
    assert plan.lanes == 2
    assert len(plan.glyphs) == 2
    assert len(plan.dropped) == 3


def test_a_feature_outside_the_reference_is_not_placed_at_all():
    """Groups may have shorter references than the view's features describe."""
    plan = plan_track([_f(2000, 2400, label="far")], 1000)
    assert plan.glyphs == []
    assert plan.lanes == 0


def test_track_height_accounts_for_lanes_and_their_gaps():
    plan = plan_track([_f(0, 500, label="a"), _f(400, 900, label="b")], 1000)
    assert plan.height == 2 * LANE_HEIGHT + LANE_GAP


def test_an_empty_track_has_no_height_and_no_svg():
    plan = plan_track([], 1000)
    assert plan.height == 0
    assert track_svg(plan) == ""


# --- Origin wrapping ------------------------------------------------------

def test_a_wrapping_feature_is_drawn_in_two_pieces_on_one_lane():
    wrapper = Feature("misc_feature", start=900, end=100, label="span",
                      wraps_origin=True)
    plan = plan_track([wrapper], 1000)
    assert len(plan.glyphs) == 2
    assert {g.lane for g in plan.glyphs} == {0}
    assert plan.lanes == 1


def test_a_wrapping_feature_carries_one_label_not_two():
    wrapper = Feature("misc_feature", start=900, end=100, label="span",
                      wraps_origin=True)
    plan = plan_track([wrapper], 1000)
    assert sum(1 for g in plan.glyphs if g.label_place) == 1


def test_the_cut_edge_of_a_wrapping_feature_gets_no_arrowhead():
    """A + strand feature crossing the origin runs off the end of the reference
    and finishes just after base 1.  So the piece drawn at the start of the
    reference owns the real 3' end and its arrowhead, and the piece at the end of
    the reference is cut flat — a point there would claim it stops.
    """
    wrapper = Feature("CDS", start=900, end=100, label="span", strand=1,
                      wraps_origin=True)
    plan = plan_track([wrapper], 1000)
    after_origin, before_origin = sorted(plan.glyphs, key=lambda g: g.x)
    assert after_origin.x == 0
    assert after_origin.head_end is True
    assert before_origin.head_end is False


def test_a_wrapping_feature_says_so_in_its_tooltip():
    wrapper = Feature("misc_feature", start=900, end=100, label="span",
                      wraps_origin=True)
    title = plan_track([wrapper], 1000).glyphs[0].title
    assert "crosses the origin" in title
    assert "200 bp" in title


# --- Labels ---------------------------------------------------------------

def test_a_wide_feature_carries_its_label_inside():
    plan = plan_track([_f(0, 500, label="AmpR")], 1000)
    assert plan.glyphs[0].label_place == "in"


def test_a_narrow_feature_puts_its_label_after_itself_when_there_is_room():
    plan = plan_track([_f(0, 6, label="a long name")], 1000)
    assert plan.glyphs[0].label_place == "after"


def test_a_narrow_feature_with_no_room_keeps_only_its_tooltip():
    plan = plan_track(
        [_f(0, 6, label="a long name"), _f(8, 400, label="neighbour")], 1000)
    crowded = next(g for g in plan.glyphs if g.label == "a long name")
    assert crowded.label_place == ""
    assert "a long name" in crowded.title


def test_an_outside_label_never_runs_past_the_reference():
    plan = plan_track([_f(996, 1000, label="a very long trailing name")], 1000)
    assert plan.glyphs[0].label_place == ""


# --- Colour ---------------------------------------------------------------

def test_a_file_colour_beats_the_type_palette():
    plan = plan_track([_f(0, 400, color="#ffd281")], 1000)
    assert plan.glyphs[0].fill_light == "#ffd281"


def test_a_feature_with_no_colour_falls_back_by_type():
    cds = plan_track([_f(0, 400, type="CDS")], 1000).glyphs[0]
    ori = plan_track([_f(0, 400, type="rep_origin")], 1000).glyphs[0]
    assert cds.fill_light != ori.fill_light


def test_an_unknown_type_still_gets_a_colour():
    plan = plan_track([_f(0, 400, type="mobile_element")], 1000)
    assert plan.glyphs[0].fill_light.startswith("#")


def test_a_dark_file_colour_is_lifted_for_the_dark_theme():
    """A colour chosen in SnapGene was chosen against white; deep navy on the
    dark ground would vanish.
    """
    glyph = plan_track([_f(0, 400, color="#101830")], 1000).glyphs[0]
    assert glyph.fill_dark != glyph.fill_light
    assert glyph.fill_dark != "#101830"


def test_a_bright_file_colour_is_left_alone_in_dark():
    glyph = plan_track([_f(0, 400, color="#ffd281")], 1000).glyphs[0]
    assert glyph.fill_dark == "#ffd281"


def test_label_text_takes_its_colour_from_the_fill_it_sits_on():
    light_fill = track_style(plan_track([_f(0, 400, color="#ffe000")], 1000))
    dark_fill = track_style(plan_track([_f(0, 400, color="#203040")], 1000))
    assert "#111111" in light_fill
    assert "#ffffff" in dark_fill


def test_every_feature_gets_a_rule_in_both_themes():
    css = track_style(plan_track(
        [_f(0, 400, label="a"), _f(500, 900, label="b")], 1000))
    assert css.count('[data-theme="dark"] .svf') == 4  # 2 features x fill + text


# --- SVG ------------------------------------------------------------------

def test_a_stranded_feature_is_drawn_with_a_point():
    forward = track_svg(plan_track([_f(0, 400, strand=1)], 1000))
    plain = track_svg(plan_track([_f(0, 400, strand=None)], 1000))
    # A pointed glyph needs diagonal segments; a block is only H and V moves.
    assert " L" in forward
    assert " L" not in plain


def test_a_reverse_feature_points_the_other_way():
    forward = track_svg(plan_track([_f(100, 400, strand=1)], 1000))
    reverse = track_svg(plan_track([_f(100, 400, strand=-1)], 1000))
    assert forward != reverse


def test_a_short_feature_does_not_grow_a_head_wider_than_itself():
    svg = track_svg(plan_track([_f(0, 2, strand=1)], 1000))
    assert svg  # drawn at all
    coordinates = [float(v) for v in re.findall(r"[HL](-?\d+\.?\d*)", svg)]
    assert max(coordinates) <= 8.0


def test_each_feature_carries_a_hover_title():
    svg = track_svg(plan_track([_f(0, 400, label="AmpR", strand=1)], 1000))
    assert "<title>" in svg
    assert "AmpR" in svg
    assert "→" in svg


def test_a_label_from_a_file_cannot_inject_markup():
    """Labels are free text out of a user's file, so they reach the page escaped."""
    nasty = '</title><script>alert(1)</script>'
    svg = track_svg(plan_track([_f(0, 800, label=nasty)], 1000))
    assert "<script>" not in svg
    assert "&lt;script&gt;" in svg


def test_the_svg_is_sized_to_the_reference():
    plan = plan_track([_f(0, 400)], 1000)
    svg = track_svg(plan)
    assert plan.width == 1000 * cell_width(1000)
    assert f'width="{plan.width:.0f}"' in svg


def test_lanes_are_stacked_by_transform():
    plan = plan_track([_f(0, 500, label="a"), _f(400, 900, label="b")], 1000)
    svg = track_svg(plan)
    assert "translate(0,0)" in svg
    assert f"translate(0,{LANE_HEIGHT + LANE_GAP})" in svg
