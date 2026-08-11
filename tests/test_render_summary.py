"""The summarized page: a complete, self-contained document of well-formed SVG."""

import re
import xml.etree.ElementTree as ET

from seqviewer.construct import Feature, Reference
from seqviewer.pileup import PileupGroup, PileupView
from seqviewer.render_summary import MAX_TABLE_ROWS, render_summary
from seqviewer.summary import SummaryView
from seqviewer.theme import Theme
from seqviewer.demo import build_summary_view

REF = "GGG" + "ATGGTGAAA" + "CCC"


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


def _view(*specs, features=(), flanks=(3, 3), ref=REF, **kwargs):
    group = PileupGroup(
        kwargs.pop("name", "g"), ref, [_row(ref, s) for s in specs],
        n_reads=len(specs), fraction=1.0, **kwargs
    )
    view = PileupView(title="test view", groups=[group], flanks=flanks,
                      features=list(features), ref_len=len(ref))
    return SummaryView.from_view(view, min_depth=1, min_count=2,
                                 min_fraction=0.25)


def _render(*specs, **kwargs):
    return render_summary(_view(*specs, **kwargs))


def _svgs(html):
    return re.findall(r"<svg\b.*?</svg>", html, re.S)


# --- The document ---------------------------------------------------------

def test_renders_a_complete_document():
    html = _render("." * len(REF))
    assert html.startswith("<!DOCTYPE html>")
    assert html.rstrip().endswith("</html>")
    assert "<title>test view</title>" in html


def test_page_is_self_contained():
    """A strict reading: no network requests of any kind."""
    html = _render("." * len(REF))
    assert "http://" not in html
    assert "https://" not in html
    assert not re.search(r"<script[^>]+src=", html)
    assert not re.search(r"<link[^>]+stylesheet", html)


def test_every_drawing_is_well_formed_xml():
    """SVG is generated as text, so nothing else checks that it parses."""
    html = render_summary(SummaryView.from_view(build_summary_view()))
    blocks = _svgs(html)
    assert blocks
    for block in blocks:
        ET.fromstring(block)


def test_the_page_carries_no_javascript_data_payload():
    """The only script is the theme bridge.

    Drawing in SVG rather than canvas is what keeps free-text feature labels out
    of a <script> block, where json.dumps would not have escaped ``</script>``.
    """
    scripts = re.findall(r"<script\b[^>]*>(.*?)</script>", _render("." * len(REF)), re.S)
    assert len(scripts) == 1
    assert "localStorage" in scripts[0]


def test_deterministic():
    assert _render("." * len(REF)) == _render("." * len(REF))


# --- Escaping -------------------------------------------------------------

def test_title_is_escaped():
    view = _view("." * len(REF))
    view.title = "<script>alert(1)</script>"
    html = render_summary(view)
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;" in html


def test_group_names_are_escaped():
    assert "<b>x</b>" not in _render("." * len(REF), name="<b>x</b>")


def test_feature_labels_from_a_file_are_escaped():
    """A label is free text out of a user's GenBank or SnapGene file."""
    nasty = Feature("CDS", 0, 15, label='</title><script>alert(1)</script>')
    html = _render("." * len(REF), features=[nasty])
    assert "<script>alert(1)</script>" not in html
    for block in _svgs(html):
        ET.fromstring(block)


def test_a_status_string_is_escaped():
    html = _render("." * len(REF), status="<img src=x>")
    assert "<img src=x>" not in html


# --- What the drawing contains -------------------------------------------

def test_a_clean_group_says_so_and_draws_no_lollipops():
    html = _render("." * len(REF), "." * len(REF))
    assert "No variants cleared" in html
    assert "No variants called" in html
    assert 'class="sv-head' not in html


def test_a_called_variant_becomes_a_lollipop_and_a_table_row():
    spec = "..." + "G........" + "..."
    html = _render(spec, spec)
    assert html.count('class="sv-head ') == 1
    assert "M1V" in html
    assert "Missense" in html


def test_a_substitution_carries_its_alternate_base_as_a_letter():
    spec = "..." + "G........" + "..."
    assert 'class="sv-head-letter"' in _render(spec, spec)


def test_a_deletion_is_drawn_as_a_different_shape_from_a_substitution():
    deletion = "..." + "...---..." + "..."
    substitution = "..." + "G........" + "..."
    assert "<polygon class=\"sv-head" in _render(deletion, deletion)
    assert "<circle class=\"sv-head" in _render(substitution, substitution)


def test_severity_reaches_the_glyph_as_a_class():
    frameshift = "..." + "...-....." + "..."
    html = _render(frameshift, frameshift)
    assert "sv-bad" in html
    assert "Frameshift" in html


def test_the_reading_frame_is_marked_on_the_band():
    with_frame = _render("." * len(REF), flanks=(3, 3))
    without = _render("." * len(REF), flanks=None)
    # The rule is always in the stylesheet; only the element is conditional.
    assert '<line class="sv-focus-edge"' in with_frame
    assert '<line class="sv-focus-edge"' not in without
    assert '<rect class="sv-ribbon-focus"' in with_frame


def test_features_are_drawn_by_the_shared_annotation_module():
    """The track is annotate.py's, not a second implementation."""
    html = _render("." * len(REF), features=[Feature("CDS", 0, 15, label="orf")])
    assert 'class="sv-annot"' in html
    assert re.search(r'<path class="svf\d+"', html)
    assert "orf" in html


def test_a_page_with_no_features_draws_no_track():
    assert 'class="sv-annot"' not in _render("." * len(REF))


def test_features_that_do_not_fit_are_named_rather_than_dropped_silently():
    crowd = [Feature("protein_bind", 0, 15, label=f"site{i}") for i in range(4)]
    html = _render("." * len(REF), features=crowd)
    assert "not drawn for want of lanes" in html


def test_a_coverage_profile_is_drawn():
    assert 'class="sv-depth"' in _render("." * len(REF))


def test_the_legend_explains_the_glyph_vocabulary():
    html = _render("." * len(REF))
    assert "substitution" in html
    assert "deletion" in html
    assert "allele frequency" in html


# --- Long output ----------------------------------------------------------

def test_a_very_long_variant_list_is_cut_short_and_says_so():
    ref = "A" * 400
    specs = ["T" * 400, "T" * 400]
    view = PileupView(title="t", groups=[
        PileupGroup("g", ref, [_row(ref, s) for s in specs], n_reads=2)
    ])
    html = render_summary(
        SummaryView.from_view(view, min_depth=1, min_count=2, min_fraction=0.25)
    )
    assert "further variants not listed" in html
    assert html.count("<tr>") <= MAX_TABLE_ROWS + 1


# --- Theme bridge ---------------------------------------------------------

def test_theme_names_are_namespaced_by_default():
    html = _render("." * len(REF))
    assert 'id="cv-theme-bridge"' in html
    assert "--cv-bg" in html
    assert "seqviewer-theme" in html


def test_host_application_can_supply_its_own_theme_names():
    view = _view("." * len(REF))
    view.theme = Theme(storage_key="app-theme", css_prefix="app",
                       style_id="app-bridge", script_id="app-sync")
    html = render_summary(view)
    assert 'id="app-bridge"' in html
    assert "--app-bg" in html
    assert "localStorage.getItem('app-theme')" in html
    assert "--cv-bg" not in html


def test_the_palette_is_shared_with_the_pileup_rather_than_copied():
    """Both views must reach the same hex for the same thing."""
    from seqviewer.render import _PALETTE

    html = _render("." * len(REF))
    for key in ("a", "t", "c", "g", "boundary"):
        assert f"--cv-{key}: {_PALETTE['light'][key]};" in html


# --- The demo -------------------------------------------------------------

def test_the_summary_demo_renders_and_shows_each_kind_of_finding():
    html = render_summary(SummaryView.from_view(build_summary_view()))
    assert html.startswith("<!DOCTYPE html>")
    for name in ("clone-A", "clone-B", "clone-C"):
        assert name in html
    assert "No variants called" in html      # the clean clone
    assert "In-frame indel" in html
    assert "Frameshift" in html


def test_the_summary_demo_is_deterministic():
    from seqviewer.demo import build_summary_view as build

    first = render_summary(SummaryView.from_view(build(seed=3)))
    assert first == render_summary(SummaryView.from_view(build(seed=3)))
    assert first != render_summary(SummaryView.from_view(build(seed=4)))


def test_the_pileup_demo_is_unchanged_by_this_work():
    """The seeded pileup page is a regression check other work relies on."""
    from seqviewer.demo import build_view
    from seqviewer.render import render

    assert render(build_view(seed=1)) == render(build_view(seed=1))
