import json
import re

import pytest

from seqviewer import PileupGroup, PileupView, Theme, render
from seqviewer.demo import build_view


def _view(**kwargs):
    ref = "ACGT" * 25
    rows = [[(b, True) for b in ref], [("A", b == "A") for b in ref]]
    kwargs.setdefault("title", "test view")
    kwargs.setdefault("groups", [PileupGroup("ref-1", ref, rows, n_reads=2, fraction=1.0)])
    return PileupView(**kwargs)


def test_renders_a_complete_document():
    html = render(_view())
    assert html.startswith("<!DOCTYPE html>")
    assert html.rstrip().endswith("</html>")
    assert "<title>test view</title>" in html


def test_page_is_self_contained():
    """A strict reading: no network requests of any kind."""
    html = render(_view())
    assert "http://" not in html
    assert "https://" not in html
    assert not re.search(r"<script[^>]+src=", html)
    assert not re.search(r"<link[^>]+stylesheet", html)


def test_title_is_escaped():
    html = render(_view(title="<script>alert(1)</script>"))
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;" in html


def test_group_names_are_escaped():
    ref = "ACGT"
    view = _view(groups=[PileupGroup("<b>x</b>", ref, [], n_reads=0, fraction=0.0)])
    assert "<b>x</b>" not in render(view)


def test_theme_names_are_namespaced_by_default():
    html = render(_view())
    assert 'id="cv-theme-bridge"' in html
    assert "--cv-bg" in html
    assert "seqviewer-theme" in html
    assert "usortm" not in html


def test_host_application_can_supply_its_own_theme_names():
    theme = Theme(storage_key="app-theme", css_prefix="app",
                  style_id="app-bridge", script_id="app-sync")
    html = render(_view(theme=theme))
    assert 'id="app-bridge"' in html
    assert "--app-bg" in html
    assert "localStorage.getItem('app-theme')" in html
    assert "--cv-bg" not in html


def test_highlight_label_is_configurable():
    view = _view(highlight_ids=["ref-1"], highlight_label="Recoverable")
    assert "Recoverable: ref-1" in render(view)


def test_no_highlight_line_when_nothing_is_highlighted():
    assert "Highlighted:" not in render(_view())


def test_highlighted_groups_get_a_star():
    ref = "ACGT"
    plain = render(_view(groups=[PileupGroup("g", ref, [], highlighted=False)]))
    starred = render(_view(groups=[PileupGroup("g", ref, [], highlighted=True)]))
    assert "&#9733;" not in plain
    assert "&#9733;" in starred


def test_empty_group_renders_a_notice_rather_than_failing():
    html = render(_view(groups=[PileupGroup("g", "ACGT", [], n_reads=7)]))
    assert "No aligned reads available" in html


def test_flanks_add_the_vector_key_and_region_data():
    without = render(_view())
    with_flanks = render(_view(flanks=(10, 10)))
    assert ">vector</span>" not in without
    assert ">vector</span>" in with_flanks
    assert "var flanks=[10,10];" in with_flanks


def test_flanks_name_the_insert_span_once():
    html = render(_view(flanks=(10, 10)))
    assert "insert <b>11</b>&ndash;<b>90</b>" in html


def test_every_fact_says_what_its_number_is():
    """The line is deliberately terse, which left "10 rows drawn" sitting next
    to "65 of 100 reads" as two read counts in similar words with nothing to
    tell them apart.  Each fact carries the distinction in a tooltip instead of
    spending line width on it.
    """
    ref = "ACGT" * 25
    rows = [[(b, True) for b in ref]]
    html = render(_view(flanks=(10, 10), total_reads=9,
                        groups=[PileupGroup("g", ref, rows, n_reads=4)]))

    facts = re.findall(r'<span class="sv-fact" title="([^"]*)">(.*?)</span>', html)
    assert facts, "no fact carries a tooltip"

    labelled = {re.sub(r"<[^>]+>", "", text) for _, text in facts}
    assert {"4 of 9 reads", "1 rows drawn", "100.0% identity"} <= labelled

    for tip, text in facts:
        assert tip.strip(), f"{text} has an empty tooltip"


# --- The masthead ---------------------------------------------------------

def test_swatches_and_canvas_read_one_palette():
    """The legend used to hardcode light-mode hexes while the canvas swapped
    palettes in JS, so in dark mode a swatch could contradict the cell it
    described.

    Asserted as the invariant — every colour the canvas fills with equals the
    custom property the matching swatch reads, in both themes — rather than on
    the shape of the emission, so reshaping the JS cannot make this pass while
    the two silently drift apart.
    """
    html = render(_view(flanks=(10, 10)))

    payload = re.search(r"var SV_PALETTE = (\{.*?\});", html, re.S)
    assert payload, "the page emits no palette for the canvas"
    canvas = json.loads(payload.group(1))

    def custom_properties(selector):
        block = html.split(selector, 1)[1].split("}", 1)[0]
        return dict(re.findall(r"--cv-([a-z-]+):\s*([^;]+);", block))

    for theme, selector in (("light", ":root {"), ("dark", '[data-theme="dark"] {')):
        declared = custom_properties(selector)
        for key, value in canvas[theme].items():
            assert declared.get(key) == value, (
                f"{theme} {key}: swatch reads {declared.get(key)!r}, "
                f"canvas fills {value!r}"
            )

    # And the swatches read those properties rather than carrying literals.
    swatches = re.findall(r'class="sv-sw"[^>]*style="background:([^;"]+)', html)
    assert swatches, "no swatches rendered"
    assert all(value.startswith("var(--cv-") for value in swatches), swatches


def test_text_and_plot_share_one_left_edge():
    assert "margin-left: -2.5rem" not in render(_view())


def test_facts_are_not_restated_under_the_plot():
    html = render(_view())
    assert "aligned reads &times;" not in html
    assert "Top fraction" not in html


def test_a_single_group_gets_no_separate_header():
    html = render(_view(groups=[PileupGroup("only-ref", "ACGT", [])]))
    assert '<div class="sv-group-head">' not in html
    assert '<span class="sv-name">only-ref' in html


def test_several_groups_each_get_their_own_header():
    ref = "ACGT"
    html = render(_view(groups=[PileupGroup("a", ref, []), PileupGroup("b", ref, [])]))
    assert html.count('<div class="sv-group-head">') == 2


def test_a_title_that_only_restates_the_group_name_is_not_repeated():
    html = render(_view(groups=[PileupGroup("puc19", "ACGT", [])],
                        title="Pileup: puc19"))
    assert 'class="sv-eyebrow"' not in html


def test_shared_reference_geometry_is_stated_once_across_groups():
    """bp and insert span describe the reference, not a group, so several groups
    over one reference state them in the masthead rather than in each header.
    """
    ref = "ACGT" * 25
    html = render(_view(groups=[PileupGroup("a", ref, []), PileupGroup("b", ref, [])],
                        flanks=(10, 10)))
    assert html.count("insert <b>11</b>&ndash;<b>90</b>") == 1
    assert html.count("</b> bp") == 1


def test_groups_over_different_references_keep_their_own_geometry():
    html = render(_view(groups=[PileupGroup("a", "ACGT" * 25, []),
                                PileupGroup("b", "ACGT" * 30, [])]))
    assert html.count("</b> bp") == 2


def test_the_star_is_defined_rather_than_restated():
    """The star already marks the group inline; repeating the name adds nothing,
    but the label explaining what the star means was previously missing.
    """
    html = render(_view(groups=[PileupGroup("g", "ACGT", [], highlighted=True)],
                        highlight_ids=["g"], highlight_label="Recoverable"))
    assert "&#9733;</span> = Recoverable" in html
    assert "Recoverable: g" not in html


def test_a_title_that_adds_information_is_kept():
    html = render(_view(groups=[PileupGroup("puc19", "ACGT", [])],
                        title="Plate 3 well A1"))
    assert 'class="sv-eyebrow"' in html
    assert "Plate 3 well A1" in html


# --- Verdict chips --------------------------------------------------------

def _chip_html(status, **kwargs):
    return render(_view(groups=[PileupGroup("g", "ACGT", [], status=status, **kwargs)]))


def test_a_blank_status_renders_no_chip():
    assert '<span class="sv-chip' not in _chip_html("")


def test_an_unrecognised_status_is_not_drawn_as_a_failure():
    """A caller reporting something the renderer has no opinion about used to
    fall through to the error colour, so every group from a driver that passed
    a constant like "Aligned" printed in alarm red.
    """
    html = _chip_html("Aligned")
    assert 'class="sv-chip sv-chip-neutral"' in html
    assert "sv-chip-bad" not in html.split("</style>")[1]


def test_a_mismatch_status_is_drawn_as_a_failure():
    assert 'class="sv-chip sv-chip-bad"' in _chip_html("Mismatch")


def test_perfect_match_outranks_the_mismatch_needle():
    """"Perfect Match" contains "match"; "Mismatch" must not win on it."""
    assert 'class="sv-chip sv-chip-ok"' in _chip_html("Perfect Match")


def test_a_silent_mutation_is_drawn_as_a_warning():
    assert 'class="sv-chip sv-chip-warn"' in _chip_html("Silent Mutation")


def test_the_verdict_carries_a_glyph_as_well_as_a_colour():
    """Colour alone fails for colourblind readers and in print."""
    assert "&#9733;" not in _chip_html("Mismatch")
    assert "▲" in _chip_html("Mismatch")
    assert "●" in _chip_html("Perfect Match")


def test_translation_track_is_drawn_only_when_an_insert_is_known():
    """The payload always carries a translation; the page draws it only with flanks.

    Placing the amino-acid row needs the insert's start, so the canvas gates on
    ``flanks`` being present rather than on the translation being non-empty.
    """
    assert "var flanks=null;" in render(_view())
    with_flanks = render(_view(flanks=(10, 10)))
    assert "var flanks=[10,10];" in with_flanks
    assert "var refAA=null;" not in with_flanks


def test_translation_is_suppressed_when_the_view_asks_for_none():
    """A focus region worth marking is not always a reading frame worth
    translating — a vector's payload can be one without being the other — so the
    protein readout is its own flag rather than a consequence of `flanks`.
    """
    html = render(_view(flanks=(10, 10), translate=False))
    assert "var refAA=null;" in html
    assert "var consAA=null;" in html
    # The focus region still gets its boundary lines.
    assert "var flanks=[10,10];" in html


def test_translation_is_drawn_by_default_when_a_focus_region_exists():
    html = render(_view(flanks=(10, 10)))
    assert "var refAA=null;" not in html


def test_asking_for_a_translation_without_a_focus_region_draws_none():
    """There is nowhere to place amino-acid rows without a region to translate."""
    html = render(_view(translate=True))
    assert "var refAA=null;" in html


def test_the_translation_flag_does_not_disturb_the_reads():
    """Guards the refactor: the consensus reconstruction moved inside the flag."""
    ref = "ACGT" * 25
    rows = [[(b, True) for b in ref]]
    on = render(_view(groups=[PileupGroup("g", ref, rows)], flanks=(10, 10)))
    off = render(_view(groups=[PileupGroup("g", ref, rows)], flanks=(10, 10),
                       translate=False))
    assert '"' + "." * 100 + '"' in on
    assert '"' + "." * 100 + '"' in off


def test_reads_are_encoded_compactly():
    """Matches collapse to '.', so the payload stays small on wide references."""
    ref = "ACGT" * 25
    rows = [[(b, True) for b in ref]]
    html = render(_view(groups=[PileupGroup("g", ref, rows)]))
    assert '"' + "." * 100 + '"' in html


# --- Mismatch frequency track ----------------------------------------------

def test_the_mismatch_track_is_off_by_default():
    assert 'class="sv-mf"' not in render(_view())


def test_asking_for_it_draws_the_track():
    assert 'class="sv-mf"' in render(_view(mismatch_freq=True))


def test_the_track_is_positioned_after_the_features_track():
    """"Below the features" means in the markup between the annotation SVG and
    the ruler canvas, since the page stacks them top to bottom in document flow.
    """
    from seqviewer.construct import Feature

    ref = "ACGT" * 25
    html = render(_view(mismatch_freq=True, features=[Feature("CDS", 0, 10, label="f")],
                        ref_len=len(ref)))
    annot_pos = html.index('class="sv-annot"')
    mf_pos = html.index('class="sv-mf"')
    ruler_pos = html.index('class="pileup-ruler"')
    assert annot_pos < mf_pos < ruler_pos


def test_a_group_with_no_rows_draws_no_pileup_at_all():
    """The empty-group notice replaces the whole pileup, mismatch track included."""
    html = render(_view(mismatch_freq=True,
                        groups=[PileupGroup("g", "ACGT", [], n_reads=7)]))
    assert 'class="sv-mf"' not in html
    assert "No aligned reads available" in html


def test_demo_view_renders():
    html = render(build_view())
    assert html.startswith("<!DOCTYPE html>")
    assert "pUC19-WT" in html
    assert "pUC19-K44A" in html


def test_demo_is_deterministic():
    assert render(build_view(seed=1)) == render(build_view(seed=1))
    assert render(build_view(seed=1)) != render(build_view(seed=2))


def test_demo_derives_flanks_from_its_reference():
    assert build_view().flanks == (100, 100)


def test_renderer_needs_no_third_party_imports():
    """The core is dependency-free; this fails loudly if that stops being true."""
    import ast
    import importlib
    import pathlib

    # seqviewer.render is the re-exported function, so ask for the module.
    render_mod = importlib.import_module("seqviewer.render")

    tree = ast.parse(pathlib.Path(render_mod.__file__).read_text())
    roots = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            roots.add(node.module.split(".")[0])
    # importlib is stdlib: the renderer reads its CSS and JS out of the
    # package's assets directory with importlib.resources.
    assert roots <= {"html", "json", "collections", "importlib", "__future__"}, roots


def test_scroll_sync_runs_after_the_groups():
    """The mirror has to find every pileup, so it cannot run before them."""
    html = render(_view())
    assert "syncPileupScrolls();" in html
    assert html.index("syncPileupScrolls();") > html.rindex("drawPileup(")


def test_scroll_sync_mirrors_a_raw_offset():
    """Equal pixel offsets are the same base; equal fractions are not."""
    html = render(_view())
    call = html[html.index("function syncPileupScrolls"):]
    call = call[:call.index("\n}")]
    assert "other.scrollLeft = left" in call
    assert "scrollWidth" not in call


def test_scroll_sync_leaves_a_lone_group_alone():
    """One group has nothing to stay in register with."""
    html = render(_view())
    body = html[html.index("function syncPileupScrolls"):]
    assert "panes.length < 2" in body[:body.index("\n}")]


# --- Opening position ------------------------------------------------------

def _opening_block(html):
    """The initial-scroll block out of the inlined drawPileup source."""
    start = html.index("// --- Opening position")
    return html[start:html.index("// --- Mismatch overflow arrows ---", start)]


def test_the_pileup_opens_on_the_insert_when_there_is_one():
    """A page is opened to look at the insert; the vector is context."""
    block = _opening_block(render(_view(flanks=(10, 10))))
    assert "if (flanks && scrollId)" in block
    assert "openEl.scrollLeft" in block


def test_the_opening_position_is_set_before_the_overflow_arrows_attach():
    """Their first update has to read the position the reader will see."""
    html = render(_view(flanks=(10, 10)))
    assert html.index("// --- Opening position") < html.index(
        "// --- Mismatch overflow arrows ---"
    )


def test_a_centred_insert_is_measured_against_the_pane_not_the_page():
    block = _opening_block(render(_view(flanks=(10, 10))))
    assert "openEl.clientWidth" in block
    assert "(insLeft + insRight) / 2 - viewW / 2" in block


def test_an_insert_wider_than_the_pane_opens_at_its_start():
    """Reading starts at the 5' boundary, not in the middle of a long insert."""
    block = _opening_block(render(_view(flanks=(10, 10))))
    assert "insRight - insLeft <= viewW" in block
    assert "target = insLeft - 12" in block


def test_the_opening_position_is_clamped_to_what_can_be_scrolled():
    block = _opening_block(render(_view(flanks=(10, 10))))
    assert "openEl.scrollWidth - viewW" in block
    assert "Math.max(0, Math.min(target, maxLeft))" in block


def test_a_page_with_no_insert_opens_where_it_always_did():
    """flanks is null, so the block is inert rather than absent."""
    assert "var flanks=null;" in render(_view())


# --- The help cursor promises a tooltip -----------------------------------

def test_the_help_cursor_is_gated_on_carrying_a_title():
    """A help cursor is a promise of a tooltip.

    The legend's own star is the case that got this wrong: it is the definition
    of the mark, so it has nothing to explain on hover, and it drew a question
    mark that did nothing.
    """
    css = render(_view())
    assert ".sv-fact[title], .sv-star[title] {" in css
    assert re.search(r"\.sv-star \{\s*color: var\(--warn\);\s*\}", css)


def test_the_legends_star_carries_no_title_and_so_claims_no_cursor():
    html = render(_view(highlight_ids=["ref-1"], highlight_label="Recoverable",
                        groups=[PileupGroup("ref-1", "ACGT" * 25,
                                            [[(b, True) for b in "ACGT" * 25]],
                                            n_reads=1, highlighted=True)]))
    # The body element, not the stylesheet rule of the same name.
    legend = html[html.index('class="sv-highlight"'):]
    legend = legend[:legend.index("</div>")]
    assert '<span class="sv-star">' in legend
    assert "title=" not in legend
