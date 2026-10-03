"""Write the two example pages the docs embed.

    python docs/demo/make_demo.py

Both are the package's synthetic demos cut down to a single clone, so each page
shows one clone's findings and nothing else.  In the summary it is the clone
carrying a missense change and an in-frame deletion; in the pileup, the mutant
group, whose reads all disagree with the reference at one position.
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "src"))

from seqviewer.demo import build_summary_view, build_view      # noqa: E402
from seqviewer.render import render                             # noqa: E402
from seqviewer.render_summary import render_summary             # noqa: E402
from seqviewer.summary import SummaryView                       # noqa: E402


def single_clone_view():
    view = build_summary_view()
    clone = next(g for g in view.groups if g.name == "clone-B")
    clone = dataclasses.replace(
        clone, name="clone-1", n_reads=clone.n_reads, fraction=1.0,
        status="", highlighted=False,
    )
    return dataclasses.replace(
        view, groups=[clone], total_reads=clone.n_reads, highlight_ids=[],
    )


#: The drawing is laid out in a 778-unit space and scaled to its container, so
#: its type shrinks with a narrow column.  Embedded in a page, the labels read
#: better a few units larger and the wrapper has no need of its own margin.
LARGER_TYPE = """
.sv-wrap { padding-bottom: 0.5rem; }
.sv-annot text { font-size: 13px; }
.sv-tick-label, .sv-row-label { font-size: 11.5px; }
.sv-mm-tick { font-size: 9px; }
"""


#: Embedded, a page scrolls with the one around it, and a scrollbar of its own
#: is only in the way.  Trackpad and wheel scrolling are unaffected.
#: An embedded page has no page transition to fade in from, and the theme is
#: chosen before the first paint (below) so it never shows one and then the
#: other.  The overlay the pages fade in from is removed.
NO_FADE = """
.sv-fade { display: none; }
"""

#: Run first in the page's head: take the theme of the page around it, the saved
#: choice or else the system's, so an embed is drawn in it from the start.
THEME_FIRST = """<script>
(function () {
  try {
    var up = window.parent.document.documentElement.getAttribute("data-theme");
    var saved = localStorage.getItem("seqviewer-theme");
    var choice = up || (saved === "dark" || saved === "light" ? saved : "");
    var dark = choice ? choice === "dark" : window.matchMedia(
      "(prefers-color-scheme: dark)").matches;
    if (dark) document.documentElement.setAttribute("data-theme", "dark");
    else document.documentElement.removeAttribute("data-theme");
  } catch (e) {}
})();
</script>"""

NO_SCROLLBARS = """
* { scrollbar-width: none; }
*::-webkit-scrollbar { display: none; }
"""


#: The pages' own dark palette is a deep navy, which sits below the docs page
#: rather than on it.  Embedded, the dark theme is raised to the docs page's
#: warm grey, a good deal lighter, with the lines and the quiet text lifted to
#: match.
LIGHTER_DARK = """
[data-theme="dark"] {
    --cv-bg: #302d29;
    --card-bg: #3b3732;
    --panel-line: #554f47;
    --muted: #b0a99e;
    --cv-tick: #928b80;
    --cv-tick-label: #ece7df;
    --cv-flag: #928b80;
    --cv-match: #675f55;
    --cv-vector: #565046;
    --cv-aa-match: #675f55;
    --cv-aa-bg: #34302b;
    --cv-aa-grid: #554f47;
    --cv-depth: #5a5247;
    --cv-depth-edge: #9a9286;
    --cv-stem: #9a9286;
    --cv-ribbon: #675f55;
    --cv-focus: #a39b8f;
    --cv-grid: #554f47;
    --cv-mm-bed: #45403a;
    --cv-mm-quiet: #7d756a;
    --cv-neutral-bg: #5c554b;
}
"""


def embed(page: str, extra: str = "") -> str:
    page = page.replace(
        "</style>", extra + LIGHTER_DARK + NO_FADE + NO_SCROLLBARS + "</style>", 1)
    return page.replace("<head>", "<head>" + THEME_FIRST, 1)


def single_clone_pileup():
    """The pileup demo cut to its mutant group, which has a finding to show."""
    view = build_view()
    clone = next(g for g in view.groups if g.name == "pUC19-K44A")
    clone = dataclasses.replace(
        clone, name="clone-1", fraction=1.0, highlighted=False)
    return dataclasses.replace(
        view, groups=[clone], total_reads=clone.n_reads, highlight_ids=[])


def main() -> None:
    (HERE / "pileup.html").write_text(embed(render(single_clone_pileup())))
    page = render_summary(SummaryView.from_view(single_clone_view()))
    (HERE / "summary.html").write_text(embed(page, LARGER_TYPE))


if __name__ == "__main__":
    main()
