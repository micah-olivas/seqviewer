import pytest

from seqviewer import Feature, PileupGroup, PileupView, Reference


def _ref(length=100, insert=(20, 80), **kwargs):
    features = []
    if insert is not None:
        features.append(Feature("insert", insert[0], insert[1], **kwargs))
    return Reference(seq="A" * length, name="test", features=features)


def test_flank_lengths_from_annotated_insert():
    assert _ref().flank_lengths() == (20, 20)


def test_flank_lengths_are_asymmetric_when_the_insert_is():
    assert _ref(insert=(10, 90)).flank_lengths() == (10, 10)
    assert _ref(insert=(30, 90)).flank_lengths() == (30, 10)


def test_no_insert_means_no_flanks():
    assert _ref(insert=None).flank_lengths() is None


def test_insert_spanning_the_whole_reference_has_no_flanks():
    assert _ref(insert=(0, 100)).flank_lengths() is None


def test_origin_wrapping_insert_is_refused():
    """Its real extent is not [start, end), so deriving flanks from it would lie."""
    assert _ref(insert=(80, 20), wraps_origin=True).flank_lengths() is None


def test_find_returns_first_match_or_none():
    ref = _ref()
    assert ref.find("insert").start == 20
    assert ref.find("promoter") is None


def test_topology():
    assert not Reference("ACGT").is_circular
    assert Reference("ACGT", topology="circular").is_circular


def test_group_rejects_rows_that_do_not_match_the_reference_width():
    with pytest.raises(ValueError, match="3 cells wide"):
        PileupGroup("g", "ACGT", [[("A", True)] * 3])


def test_group_accepts_matching_rows():
    group = PileupGroup("g", "ACGT", [[("A", True)] * 4])
    assert len(group.rows) == 1


def test_view_from_reference_derives_flanks():
    view = PileupView.from_reference(
        "title", _ref(), [PileupGroup("g", "A" * 100, [])]
    )
    assert view.flanks == (20, 20)


def test_view_defaults_are_self_contained():
    view = PileupView(title="t")
    assert view.groups == []
    assert view.flanks is None
    assert view.theme.css_prefix == "cv"
