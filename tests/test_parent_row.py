"""Tests for the parent row and the extent of the region tint."""

import pytest

from seqviewer import PileupGroup, PileupView, render

REF = "ACGT" * 30                    # 120 bp
WT = REF[:50] + ("A" if REF[50] != "A" else "C") + REF[51:]


def _rows(n=3):
    return [[("A", True)] * len(REF) for _ in range(n)]


def _page(**kw):
    group = PileupGroup(name="v1", ref_seq=REF, rows=_rows(), n_reads=3,
                        **kw)
    return render(PileupView(title="Pileup: Plate 1 Well A1",
                             groups=[group], total_reads=3,
                             flanks=(20, 30)))


class TestWildTypeRow:

    def test_absent_unless_given(self):
        assert 'var parent="";' in _page().replace("var parent='';", 'var parent="";')

    def test_encoded_against_the_reference(self):
        """A dot where the parent agrees, the parent's base where it does not,
        so the designed change is the one column that stands out."""
        import re

        page = _page(parent=WT)
        encoded = re.search(r'var parent="([^"]*)"', page).group(1)
        assert len(encoded) == len(REF)
        differing = [i for i, c in enumerate(encoded) if c != "."]
        assert differing == [50]

    def test_the_row_is_labelled(self):
        assert "'Parent'" in _page(parent=WT)

    def test_a_parent_identical_to_the_reference_shows_no_difference(self):
        import re

        page = _page(parent=REF)
        encoded = re.search(r'var parent="([^"]*)"', page).group(1)
        assert set(encoded) == {"."}

    def test_the_wrong_length_is_refused(self):
        """A row that does not line up with the reference would put the change
        in the wrong column, which is worse than no row."""
        with pytest.raises(ValueError, match="parent is"):
            PileupGroup(name="v1", ref_seq=REF, rows=_rows(),
                        parent=REF[:-5])

    def test_groups_without_one_still_render(self):
        page = _page()
        assert page.lstrip().lower().startswith("<!doctype html")


class TestRegionTintExtent:
    """One labelled block names each span, and nothing else is tinted.

    The wash behind the features, the mismatch track and the ruler said the same
    thing a third time, in a colour those rows then had to be read through.
    """

    def test_no_tint_is_drawn_at_all(self):
        page = _page()
        assert "sv-region" not in page

    def test_the_band_still_names_the_spans(self):
        page = _page()
        assert 'class="sv-band"' in page
        assert "5\u2032 vector" in page or "5′ vector" in page

    def test_the_tracks_need_no_stacking_context_now(self):
        """The z-index stack existed only to sit above the tint."""
        page = _page()
        assert "z-index: 0" not in page

    def test_the_boundaries_still_run_the_full_height(self):
        """The dashed lines cost no contrast, so they stay: they are what marks
        the insert once the wash no longer does."""
        page = _page()
        assert "ctx.lineTo(bx, pH)" in page
