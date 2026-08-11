"""Reducing a pileup to depth, agreement, and called variants.

Pure arithmetic over grids, so nothing here needs a browser, an aligner, or
Biopython.
"""

import pytest

from seqviewer.construct import Feature, Reference
from seqviewer.pileup import PileupGroup, PileupView
from seqviewer.summary import (
    DEFAULT_FLAG_THRESHOLD, DEFAULT_MIN_COUNT, DEFAULT_MIN_FRACTION,
    GroupSummary, SummaryView, Variant, flagged_columns, summarize_group,
)

# A reference with a 9-base reading frame in the middle: ATG GTG AAA = M V K.
REF = "GGG" + "ATGGTGAAA" + "CCC"
FOCUS = (3, 12)


def _row(ref, spec):
    """One row from a spec string: ``.`` match, ``-`` gap, a letter mismatch."""
    assert len(spec) == len(ref)
    out = []
    for i, ch in enumerate(spec):
        if ch == ".":
            out.append((ref[i], True))
        elif ch == "-":
            out.append(("-", True))       # align writes deletions like this
        else:
            out.append((ch, ch == ref[i]))
    return out


def _group(ref, *specs, **kwargs):
    return PileupGroup(
        kwargs.pop("name", "g"), ref, [_row(ref, s) for s in specs], **kwargs
    )


def _view(group, flanks=(3, 3), translate=True, **kwargs):
    return PileupView(
        title="t", groups=[group], flanks=flanks, translate=translate, **kwargs
    )


def _summary(group, flanks=(3, 3), translate=True, **kwargs):
    return SummaryView.from_view(_view(group, flanks, translate), **kwargs)


# --- Per-position arrays --------------------------------------------------

def test_depth_counts_reads_reaching_a_position():
    group = _group("ACGT", "....", "..--", "----")
    s = summarize_group(group)
    # The all-gap read called nothing and contributes no coverage anywhere; the
    # second read's trailing gaps are absence of coverage, not deletions.
    assert s.depth == [2, 2, 1, 1]


def test_matches_counts_only_agreement_with_the_reference():
    s = summarize_group(_group("ACGT", "....", "T..."))
    assert s.matches == [1, 2, 2, 2]


def test_bases_called_at_a_position_are_depth_minus_deletions():
    s = summarize_group(_group("ACGTA", ".-...", "....."))
    assert s.depth == [2, 2, 2, 2, 2]
    assert s.deletions == [0, 1, 0, 0, 0]


def test_a_group_holds_no_reference_sequence_of_its_own():
    """The reference belongs to the view, which draws one coordinate system.

    A copy per group is how a page comes to serialise the same bases once for
    every subpopulation it shows.
    """
    s = summarize_group(_group("ACGT", "...."))
    assert not hasattr(s, "ref_seq")
    assert s.ref_len == 4


# --- Deletion versus no coverage -----------------------------------------

def test_leading_and_trailing_gaps_are_not_deletions():
    """A read that starts late has not deleted the bases before it."""
    s = summarize_group(_group("ACGTACGT", "--....--", "--....--"))
    assert s.deletions == [0] * 8
    assert s.depth == [0, 0, 2, 2, 2, 2, 0, 0]
    assert s.variants == []


def test_an_interior_gap_run_is_one_deletion_variant():
    group = _group("ACGTACGT", "..---...", "..---...", "........")
    s = summarize_group(group, min_depth=1, min_count=2, min_fraction=0.25)
    dels = [v for v in s.variants if v.kind == "del"]
    assert len(dels) == 1
    assert dels[0].pos == 2
    assert dels[0].ref == "GTA"
    assert dels[0].length == 3
    assert dels[0].count == 2
    assert dels[0].label == "Δ3 bp"


def test_a_deletion_reaching_the_end_of_coverage_is_not_called():
    """Only a run flanked by called bases on both sides is evidence."""
    s = summarize_group(_group("ACGTACGT", "....----", "....----"),
                        min_depth=1, min_count=1, min_fraction=0.1)
    assert [v for v in s.variants if v.kind == "del"] == []


# --- Calling thresholds ---------------------------------------------------

def test_a_single_supporting_read_is_never_a_variant():
    """One read disagreeing is the error rate, not an allele."""
    specs = ["T..."] + ["...."] * 3
    s = summarize_group(_group("ACGT", *specs), min_fraction=0.1, min_count=2)
    assert s.variants == []


def test_an_allele_below_the_fraction_floor_is_not_called():
    specs = ["T...", "T..."] + ["...."] * 18
    s = summarize_group(_group("ACGT", *specs), min_fraction=0.25, min_count=2)
    assert s.variants == []


def test_an_allele_clearing_both_floors_is_called():
    specs = ["T...", "T..."] + ["...."] * 2
    s = summarize_group(_group("ACGT", *specs), min_fraction=0.25, min_count=2)
    assert [(v.pos, v.ref, v.alt, v.count, v.depth) for v in s.variants] == [
        (0, "A", "T", 2, 4)
    ]
    assert s.variants[0].fraction == 0.5


def test_a_position_below_the_depth_floor_is_not_called_at_all():
    s = summarize_group(_group("ACGT", "T...", "T..."), min_depth=3, min_count=2)
    assert s.variants == []


def test_two_alternative_alleles_at_one_position_are_both_called():
    specs = ["T...", "T...", "G...", "G..."]
    s = summarize_group(_group("ACGT", *specs), min_fraction=0.25, min_count=2)
    assert sorted(v.alt for v in s.variants) == ["G", "T"]


def test_defaults_are_the_documented_ones():
    assert DEFAULT_MIN_FRACTION == 0.25
    assert DEFAULT_MIN_COUNT == 2


# --- Insertions -----------------------------------------------------------

def test_a_plain_grid_reports_no_insertions():
    """Not because there are none, but because a grid cannot carry one."""
    s = summarize_group(_group("ACGT", "....", "...."))
    assert [v for v in s.variants if v.kind == "ins"] == []


def test_sidecar_evidence_produces_insertion_variants():
    group = _group("ACGT", "....", "....", "....", "....")
    s = summarize_group(group, insertions={1: {"GG": 3}},
                        min_fraction=0.25, min_count=2)
    ins = [v for v in s.variants if v.kind == "ins"]
    assert len(ins) == 1
    assert ins[0].pos == 1
    assert ins[0].alt == "GG"
    assert ins[0].count == 3
    assert ins[0].label == "+2 bp"


def test_insertion_evidence_outside_the_reference_is_ignored():
    group = _group("ACGT", "....", "....", "....")
    s = summarize_group(group, insertions={99: {"G": 3}}, min_count=1)
    assert [v for v in s.variants if v.kind == "ins"] == []


# --- Consequences ---------------------------------------------------------

def _called(spec, **kwargs):
    """Summarize a two-read group over REF and return its variants."""
    view = _view(_group(REF, spec, spec), **kwargs)
    return SummaryView.from_view(view, min_depth=1, min_count=2,
                                 min_fraction=0.25).groups[0].variants


def test_a_substitution_that_keeps_the_residue_is_silent():
    # Third base of GTG -> GTA, still valine.
    variants = _called("..." + ".....A..." + "...")
    assert [(v.consequence, v.effect) for v in variants] == [
        ("silent", "silent (V2)")
    ]


def test_a_substitution_that_changes_the_residue_is_missense():
    # ATG -> GTG, methionine to valine.
    variants = _called("..." + "G........" + "...")
    assert [(v.consequence, v.effect) for v in variants] == [
        ("missense", "M1V")
    ]


def test_a_substitution_creating_a_stop_is_nonsense():
    # AAA -> TAA.
    variants = _called("..." + "......T.." + "...")
    assert [(v.consequence, v.effect) for v in variants] == [
        ("nonsense", "K3*")
    ]


def test_a_deletion_off_the_reading_frame_is_a_frameshift():
    variants = _called("..." + "...-....." + "...")
    assert variants[0].consequence == "frameshift"
    assert variants[0].effect == "frameshift at 2"


def test_a_deletion_of_whole_codons_is_an_in_frame_indel():
    variants = _called("..." + "...---..." + "...")
    assert variants[0].consequence == "inframe_indel"
    assert variants[0].effect == "Δ1 aa at 2"


def test_a_variant_outside_the_reading_frame_is_noncoding():
    variants = _called("A.." + "........." + "...")
    assert variants[0].consequence == "noncoding"
    assert variants[0].effect == ""


def test_no_focus_region_leaves_every_variant_unclassified():
    variants = _called("..." + "G........" + "...", flanks=None)
    assert [(v.consequence, v.effect) for v in variants] == [("", "")]


def test_asking_for_no_translation_leaves_variants_unclassified():
    variants = _called("..." + "G........" + "...", translate=False)
    assert variants[0].consequence == ""


# --- Verdicts -------------------------------------------------------------

def test_a_group_with_nothing_called_is_clean():
    s = _summary(_group(REF, "." * len(REF), "." * len(REF))).groups[0]
    assert s.verdict == "clean"


def test_a_verdict_is_the_worst_consequence_present():
    # A silent change and a frameshift in the same group.
    spec = "..." + "...GTA..." + "..."
    worse = "..." + "...-....." + "..."
    view = _view(_group(REF, spec, spec, worse, worse))
    s = SummaryView.from_view(view, min_depth=1, min_count=2,
                              min_fraction=0.25).groups[0]
    assert "frameshift" in [v.consequence for v in s.variants]
    assert s.verdict == "frameshift"


def test_a_called_variant_with_no_frame_reports_a_bare_verdict():
    spec = "..." + "G........" + "..."
    view = _view(_group(REF, spec, spec), flanks=None)
    s = SummaryView.from_view(view, min_depth=1, min_count=2,
                              min_fraction=0.25).groups[0]
    assert s.verdict == "variant"


# --- Derived figures ------------------------------------------------------

def test_identity_is_over_called_bases_not_over_the_reference():
    s = summarize_group(_group("ACGTA", "T....", "-...."))
    # Nine bases called, eight of them matching.
    assert s.identity == pytest.approx(8 / 9)


def test_identity_is_none_when_nothing_was_called():
    s = summarize_group(_group("ACGT", "----"))
    assert s.identity is None


def test_mean_depth_is_over_covered_positions_only():
    s = summarize_group(_group("ACGTACGT", "....----", "....----"))
    assert s.covered == 4
    assert s.mean_depth == 2.0
    assert s.max_depth == 2


def test_a_group_that_called_nothing_reports_zero_rather_than_failing():
    s = summarize_group(_group("ACGT"))
    assert s.depth == [0, 0, 0, 0]
    assert s.mean_depth == 0.0
    assert s.verdict == "clean"


# --- The view -------------------------------------------------------------

def test_from_view_carries_features_and_theme_across():
    reference = Reference(seq=REF, features=[Feature("insert", 3, 12, label="cds")])
    view = PileupView.from_reference(
        "t", reference, [_group(REF, "." * len(REF))],
    )
    s = SummaryView.from_view(view)
    assert [f.label for f in s.features] == ["cds"]
    assert s.focus == FOCUS
    assert s.ref_len == len(REF)


def test_a_view_with_no_groups_is_refused():
    with pytest.raises(ValueError, match="at least one group"):
        SummaryView.from_view(PileupView(title="t"))


def test_the_widest_group_sets_the_page_coordinate_system():
    """The pileup view allows groups over references of different lengths."""
    short = PileupGroup("short", "ACGT", [_row("ACGT", "....")])
    long_ = PileupGroup("long", "ACGTACGT", [_row("ACGTACGT", "........")])
    s = SummaryView.from_view(PileupView(title="t", groups=[short, long_]))
    assert s.ref_len == 8
    assert [g.ref_len for g in s.groups] == [4, 8]


def test_a_frame_too_short_to_hold_a_codon_is_no_frame_at_all():
    s = _summary(_group("ACGTAC", "......"), flanks=(2, 2))
    assert s.focus is None


def test_insertions_are_keyed_by_group_name():
    a = PileupGroup("a", "ACGT", [_row("ACGT", "....")] * 4)
    b = PileupGroup("b", "ACGT", [_row("ACGT", "....")] * 4)
    s = SummaryView.from_view(
        PileupView(title="t", groups=[a, b]),
        insertions={"b": {1: {"G": 4}}},
        min_count=2, min_fraction=0.25,
    )
    by_name = {g.name: g for g in s.groups}
    assert by_name["a"].variants == []
    assert [v.kind for v in by_name["b"].variants] == ["ins"]


def test_status_is_carried_through_but_the_verdict_is_computed():
    """A free-text status must never be what a page styles on."""
    group = _group(REF, "." * len(REF), "." * len(REF), status="Anything At All")
    s = _summary(group).groups[0]
    assert s.status == "Anything At All"
    assert s.verdict == "clean"


# --- The shared "worth a look" statistic ---------------------------------

def test_flagged_columns_marks_disagreement_above_the_threshold():
    ref = "ACGT"
    rows = [_row(ref, "T..."), _row(ref, "....")]
    assert flagged_columns(rows, ref) == [0]
    assert flagged_columns(rows, ref, threshold=0.6) == []


def test_flagged_columns_counts_a_deletion_as_disagreement():
    ref = "ACGTA"
    rows = [_row(ref, ".-..."), _row(ref, ".....")]
    assert 1 in flagged_columns(rows, ref)


def test_flagged_columns_measures_against_covering_reads_only():
    """A position two reads never reached is not a disagreement."""
    ref = "ACGTACGT"
    rows = [_row(ref, "T...----"), _row(ref, "....----")]
    assert flagged_columns(rows, ref) == [0]


def test_flagged_columns_is_looser_than_the_calling_threshold():
    """Marking a column and asserting an allele are different claims."""
    assert DEFAULT_FLAG_THRESHOLD < DEFAULT_MIN_FRACTION


def test_flagged_columns_on_a_clean_grid_is_empty():
    ref = "ACGT"
    assert flagged_columns([_row(ref, "....")] * 4, ref) == []


# --- Variant record ------------------------------------------------------

def test_fraction_of_an_uncovered_variant_does_not_divide_by_zero():
    assert Variant(pos=0, kind="snv", ref="A", alt="T", count=0, depth=0).fraction == 0.0


def test_variant_labels_read_as_changes():
    snv = Variant(pos=0, kind="snv", ref="A", alt="T", count=1, depth=2)
    deletion = Variant(pos=0, kind="del", ref="ACG", alt="", count=1, depth=2)
    insertion = Variant(pos=0, kind="ins", ref="", alt="AC", count=1, depth=2)
    assert (snv.label, deletion.label, insertion.label) == ("A→T", "Δ3 bp", "+2 bp")
    assert (snv.length, deletion.length, insertion.length) == (1, 3, 1)
