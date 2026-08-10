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


def test_flanks_add_the_vector_legend_and_region_data():
    without = render(_view())
    with_flanks = render(_view(flanks=(10, 10)))
    assert "Vector Match" not in without
    assert "Vector Match" in with_flanks
    assert "var flanks=[10,10];" in with_flanks


def test_translation_track_is_drawn_only_when_an_insert_is_known():
    """The payload always carries a translation; the page draws it only with flanks.

    Placing the amino-acid row needs the insert's start, so the canvas gates on
    ``flanks`` being present rather than on the translation being non-empty.
    """
    assert "var flanks=null;" in render(_view())
    with_flanks = render(_view(flanks=(10, 10)))
    assert "var flanks=[10,10];" in with_flanks
    assert "var refAA=null;" not in with_flanks


def test_reads_are_encoded_compactly():
    """Matches collapse to '.', so the payload stays small on wide references."""
    ref = "ACGT" * 25
    rows = [[(b, True) for b in ref]]
    html = render(_view(groups=[PileupGroup("g", ref, rows)]))
    assert '"' + "." * 100 + '"' in html


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
    assert roots <= {"html", "json", "collections", "__future__"}, roots
