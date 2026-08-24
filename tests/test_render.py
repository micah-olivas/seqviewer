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


# --- Mismatch track --------------------------------------------------------

def test_the_mismatch_track_is_always_drawn():
    """It replaced the ruler's flag triangles, which were unconditional.

    Gating it would leave a page with no account of disagreement by default.
    """
    assert 'class="sv-mf"' in render(_view())


def test_the_ruler_no_longer_draws_a_separate_row_of_flag_triangles():
    """One tall row, not several. The magnitude is in the track."""
    from seqviewer.summary import DEFAULT_FLAG_THRESHOLD

    ref = "ACGT" * 25
    rows = [[(b, True) for b in ref], [("A", b == "A") for b in ref]]
    html = render(_view(groups=[PileupGroup("g", ref, rows, n_reads=2)]))
    assert "triRowH" not in html
    assert DEFAULT_FLAG_THRESHOLD == 0.10


def test_the_page_and_a_summary_agree_about_disagreement():
    """Both read summary.mismatch_fractions, so they cannot report different
    numbers for the same column. They did while each computed its own: the old
    loop counted only called bases, so a column where half the reads had deleted
    the base read as perfectly clean.
    """
    from seqviewer.summary import mismatch_fractions

    ref = "ACGTACGT"
    deleted = [(ref[i], True) if i != 2 else ("-", True) for i in range(8)]
    clean = [(b, True) for b in ref]
    rows = [deleted, deleted, clean, clean]
    assert mismatch_fractions(rows, ref)[2] == 0.5


def test_the_track_is_positioned_after_the_features_track():
    """"Below the features" means in the markup between the annotation SVG and
    the ruler canvas, since the page stacks them top to bottom in document flow.
    """
    from seqviewer.construct import Feature

    ref = "ACGT" * 25
    html = render(_view(features=[Feature("CDS", 0, 10, label="f")],
                        ref_len=len(ref)))
    annot_pos = html.index('class="sv-annot"')
    mf_pos = html.index('class="sv-mf"')
    ruler_pos = html.index('class="pileup-ruler"')
    assert annot_pos < mf_pos < ruler_pos


def test_a_group_with_no_rows_draws_no_pileup_at_all():
    """The empty-group notice replaces the whole pileup, mismatch track included."""
    html = render(_view(groups=[PileupGroup("g", "ACGT", [], n_reads=7)]))
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


# --- Mutation tabs ---------------------------------------------------------

def _tab_block(html):
    start = html.index("// --- Mutation tabs ---")
    return html[start:html.index("// --- Region boundary", start)]


def test_a_changed_residue_is_named_the_way_a_mutation_is_written():
    """Baseline residue, position, new residue: V50I."""
    html = render(_view(flanks=(10, 10)))
    assert "name = base + (mi + 1) + obs" in html


def test_a_mutation_is_named_against_the_parent_when_there_is_one():
    """The parent is the library baseline; the reference is this well's assigned
    identity and is itself a variant, so numbering against it would count off an
    already-mutated sequence.
    """
    html = render(_view(flanks=(10, 10)))
    assert "var namedAgainst = hasParentAA ? parentAA : refAA;" in html


def test_the_tab_band_adds_canvas_height_only_when_something_changed():
    html = render(_view(flanks=(10, 10)))
    assert "var mutH = mutRows ?" in html
    assert "+ mutH;" in html


def test_tab_rows_are_packed_by_first_fit_and_capped():
    """A tab is far wider than its 6px codon, so they cannot all share a row."""
    html = render(_view(flanks=(10, 10)))
    assert "MAX_TAB_ROWS" in html
    assert "mutDropped++" in html


def test_a_tab_that_will_not_fit_is_counted_rather_than_dropped_silently():
    html = render(_view(flanks=(10, 10)))
    assert "not named for want" in html


def test_tab_widths_are_computed_rather_than_measured():
    """The layout must exist before there is a context to measure with, and the
    face is monospace, so the advance is known.
    """
    html = render(_view(flanks=(10, 10)))
    assert "TAB_ADVANCE" in html
    assert "name.length * TAB_FONT * TAB_ADVANCE" in html


def test_stems_are_drawn_before_the_boxes():
    """A lower row's stem passes the rows above it, so boxes must occlude."""
    block = _tab_block(render(_view(flanks=(10, 10))))
    assert block.index("ctx.moveTo(tab.stemX") < block.index("ctx.fillRect(tab.left")


def test_a_tab_face_is_opaque_before_it_is_tinted():
    """The tint is translucent, so tint alone would show the stems through."""
    block = _tab_block(render(_view(flanks=(10, 10))))
    opaque = block.index("ctx.fillStyle = P['aa-bg'];")
    assert opaque < block.index("ctx.fillStyle = ink[1];")


def test_a_stem_stops_at_the_top_of_its_tab_not_its_middle():
    """Through the middle would strike the text out."""
    block = _tab_block(render(_view(flanks=(10, 10))))
    assert "ctx.lineTo(tab.stemX, top - 2)" in block
    assert "ctx.lineTo(centre, top)" in block


def test_the_tab_band_gets_a_gutter_label_carrying_the_count():
    html = render(_view(flanks=(10, 10)))
    assert "mutLabel.textContent = 'Changes'" in html
    assert "mutLabel.title" in html


def test_no_tabs_without_the_translation_rows():
    """There is no residue to name if no reading frame was translated."""
    assert "var flanks=null;" in render(_view())


# --- Mismatch track readout ------------------------------------------------

def test_the_track_carries_its_counts_for_the_readout():
    """A bar's height is a log-scaled fraction, readable as more or less but not
    as a number, so the number travels alongside it.
    """
    html = render(_view())
    assert "var mmRuns=" in html
    assert "mmDis[col]" in html and "mmCov[col]" in html


def test_the_counts_are_run_length_encoded():
    """A mostly-clean reference is long runs of the same pair."""
    from seqviewer.render import _run_lengths

    assert _run_lengths([(0, 3), (0, 3), (1, 3)]) == [[2, 0, 3], [1, 1, 3]]
    assert _run_lengths([]) == []


def test_the_encoded_counts_are_shorter_than_the_positions_they_cover():
    ref = "ACGT" * 50
    rows = [[(b, True) for b in ref] for _ in range(4)]
    html = render(_view(groups=[PileupGroup("g", ref, rows, n_reads=4)]))
    runs = json.loads(re.search(r"var mmRuns=(\[.*?\]);", html).group(1))
    assert sum(r[0] for r in runs) == len(ref)
    assert len(runs) < len(ref)


def test_the_readout_reports_reads_and_not_only_a_percentage():
    """1 of 3 and 100 of 300 are both 33%, and only one of them is noise."""
    html = render(_view())
    assert "' of ' + cov + ' read'" in html


def test_the_readout_says_so_when_nothing_covers_a_position():
    assert "no reads cover this position" in render(_view())


def test_one_handler_covers_the_track_rather_than_an_element_per_base():
    """Thousands of titled rects would be thousands of nodes on a plasmid."""
    html = render(_view())
    assert "mmSvg.addEventListener('mousemove'" in html
    assert html.count("mmSvg.addEventListener") == 2      # mousemove, mouseleave


def test_the_track_shows_a_crosshair_not_a_help_cursor():
    """A help cursor on this page means a title is present; here there is none."""
    css = render(_view())
    assert "cursor: crosshair;" in css


# --- Parent AA row ---------------------------------------------------------

def _parent_view(**kwargs):
    """A well whose reference is the parent carrying one designed change."""
    parent = "ATG" * 20                      # 60 bases
    # Flanks of 3 leave an 18-codon frame; the designed change is its first
    # codon, ATG -> GTG, so M1V against the parent.
    ref = parent[:3] + "G" + parent[4:]
    rows = [[(b, b == ref[i]) for i, b in enumerate(parent)]] * 3
    kwargs.setdefault("title", "parent view")
    kwargs.setdefault("flanks", (3, 3))
    kwargs.setdefault("groups", [PileupGroup("well", ref, rows, n_reads=3,
                                             fraction=1.0, parent=parent)])
    return PileupView(**kwargs)


def test_the_parent_is_translated_and_reaches_the_page():
    html = render(_parent_view())
    assert "var parentAA=" in html
    assert "var parentAA=null;" not in html


def test_a_page_without_a_parent_carries_no_parent_translation():
    assert "var parentAA=null;" in render(_view(flanks=(10, 10)))


def test_the_amino_acid_block_grows_to_three_rows_for_a_parent():
    html = render(_parent_view())
    assert "var aaRows = hasAA ? (hasParentAA ? 3 : 2) : 0;" in html
    assert "aaRows * aaH + (aaRows - 1) * 2" in html


def test_the_parent_row_gets_its_own_gutter_label():
    html = render(_parent_view())
    assert "parentAALabel.textContent = 'Parent AA'" in html
    assert "parentAALabel.title" in html


def test_a_change_is_classified_against_both_baselines():
    """Expected, unexpected, or missing -- which is the question a well asks."""
    html = render(_parent_view())
    for kind in ("'expected'", "'unexpected'", "'missing'"):
        assert kind in html


def test_a_change_matching_the_assignment_is_expected():
    html = render(_parent_view())
    assert "obs === want ? 'expected' : 'unexpected'" in html


def test_a_change_the_assignment_promised_but_the_reads_lack_is_missing():
    html = render(_parent_view())
    assert "hasParentAA && want !== base" in html
    assert "name = base + (mi + 1) + want" in html


def test_the_three_kinds_take_three_inks():
    html = render(_parent_view())
    assert "P['tab-ok']" in html
    assert "P['tab-warn']" in html
    for theme in ("light", "dark"):
        from seqviewer.render import _PALETTE
        assert "tab-ok" in _PALETTE[theme]
        assert "tab-warn" in _PALETTE[theme]


def test_a_missing_change_is_drawn_with_a_dashed_edge():
    """The name is what was expected, not what is there."""
    html = render(_parent_view())
    assert "ctx.setLineDash(ink[2] ? [3, 2] : []);" in html


def test_hovering_one_amino_acid_row_reports_the_others():
    html = render(_parent_view())
    assert "var names = ['Ref', 'Cons', 'Parent'];" in html
    assert "others.join(', ')" in html


# --- The pinned scale's gutter ---------------------------------------------

def test_the_gutter_is_sized_to_the_track_it_labels():
    """One panel the height of the track, not a box per label: the two marks sit
    about nine pixels apart and nine-pixel boxes would overlap.
    """
    from seqviewer.annotate import MISMATCH_TRACK_HEIGHT

    html = render(_view())
    assert f"height:{MISMATCH_TRACK_HEIGHT}px" in html


def test_the_gutter_fades_to_transparent_over_a_fixed_length():
    """A percentage stop would shrink the runway for a shorter label."""
    html = render(_view())
    assert "var(--cv-bg) 26px," in html
    assert "transparent 40px" in html


def test_the_labels_carry_no_background_of_their_own():
    """The panel is the ground; a per-label box is what overlapped."""
    css = render(_view())
    block = css[css.index(".pileup-mf-marks span {"):]
    assert "background" not in block[:block.index("}")]


def test_the_focus_region_is_not_filled_the_same_as_the_flanks():
    """Both are outlined the same way; the fill is what separates them, and the
    insert reads warm rather than as a second blue block.
    """
    from seqviewer.render import _PALETTE

    for theme in ("light", "dark"):
        assert "region-focus" in _PALETTE[theme]
        assert _PALETTE[theme]["region-focus"] != _PALETTE[theme]["region"]
    css = render(_view())
    assert "fill: var(--cv-region-focus);" in css
