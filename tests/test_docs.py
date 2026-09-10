"""The documentation page is generated, so it can be checked against the code.

Its flag tables are read off the argparse parsers, which is what makes a
committed page that no longer matches the CLI a test failure rather than
something a reader discovers.
"""
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
BUILD = ROOT / "docs" / "build.py"


def load_build():
    spec = importlib.util.spec_from_file_location("docs_build", BUILD)
    module = importlib.util.module_from_spec(spec)
    sys.modules["docs_build"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def docs():
    return load_build()


def test_committed_page_matches_the_code(docs):
    """A rebuild reproduces the committed page, the date stamp aside."""
    committed = docs.OUT.read_text()
    fresh = docs.build()
    assert docs.STAMP.sub("", committed) == docs.STAMP.sub("", fresh), (
        "docs/index.html is out of date; run: python docs/build.py")


def test_every_flag_of_every_command_is_documented(docs):
    """No argument the CLI accepts is missing from the page."""
    page = docs.build()
    parsers = [docs.build_parser("seqview pileup"),
               docs.build_lengths_parser("seqview lengths")]
    for parser in parsers:
        for term, _, _ in docs.arguments(parser):
            assert f"<code>{docs.escape(term)}</code>" in page, term


def test_every_flag_carries_help_text(docs):
    """An undocumented flag renders as a blank row, so require the text."""
    parsers = [docs.build_parser("seqview pileup"),
               docs.build_lengths_parser("seqview lengths")]
    missing = [term for parser in parsers
               for term, _, text in docs.arguments(parser) if not text.strip()]
    assert not missing, f"no help text for: {', '.join(missing)}"


def test_captured_output_is_shown_verbatim(docs):
    """The transcripts reach the page unedited, so they say what the tool says."""
    page = docs.build()
    for name in ("seqview", "lengths", "pileup"):
        for line in docs.example(name).splitlines():
            assert line in page, line


def test_check_mode_passes_for_the_committed_page(docs):
    assert docs.main(["--check"]) == 0
