"""The annotated-file reader.

These build GenBank text rather than shipping fixture files, so the case each
test covers is readable next to the assertion.  Biopython is required, and the
whole module skips without it — importing :mod:`seqviewer` never needs it.
"""

import pytest

pytest.importorskip("Bio")

from seqviewer.genbank import (  # noqa: E402
    DEFAULT_SKIP_TYPES, feature_spans, load_reference,
)
from seqviewer.construct import Feature  # noqa: E402


def _feature(key, location, **qualifiers):
    """One GenBank feature, laid out on the columns the format requires.

    The key sits at column 6 and both the location and every qualifier at
    column 22.  Getting this wrong makes Biopython read the qualifier lines as
    locations, so it is worth building rather than hand-indenting.
    """
    lines = ["     " + key.ljust(16) + location]
    for name, value in qualifiers.items():
        name = name.rstrip("_")
        if value is True:
            lines.append(" " * 21 + f"/{name}")
        else:
            lines.append(" " * 21 + f'/{name}="{value}"')
    return "\n".join(lines)


def _record(*features, length=100, topology="circular"):
    """A minimal GenBank record carrying *features*."""
    seq = ("acgt" * ((length // 4) + 1))[:length]
    origin = []
    for i in range(0, length, 60):
        block = seq[i:i + 60]
        groups = " ".join(block[j:j + 10] for j in range(0, len(block), 10))
        origin.append(f"{i + 1:>9} {groups}")
    return (
        f"LOCUS       test  {length} bp    DNA     {topology} SYN 01-JAN-2020\n"
        "FEATURES             Location/Qualifiers\n"
        + "\n".join(features) + "\n"
        "ORIGIN\n" + "\n".join(origin) + "\n//\n"
    )


def _load(tmp_path, *features, length=100, topology="circular"):
    path = tmp_path / "t.gb"
    path.write_text(_record(*features, length=length, topology=topology))
    return load_reference(path)


def _by_label(reference):
    return {f.label: f for f in reference.features}


def test_the_builder_produces_a_record_biopython_can_read(tmp_path):
    """Guards the helper itself: a mis-indented qualifier reads as a location,
    which would make every other test in this file pass or fail for the wrong
    reason.
    """
    reference = _load(tmp_path, _feature("CDS", "10..40", label="x"))
    assert len(reference.features) == 1
    assert reference.features[0].type == "CDS"
    assert (reference.features[0].start, reference.features[0].end) == (9, 40)


# --- Origin wrapping ------------------------------------------------------

def test_a_spliced_feature_is_not_treated_as_wrapping(tmp_path):
    """`len(location.parts) > 1` is not the test: a spliced CDS has two parts
    and stays put.  The old rule called this a wrapping feature.
    """
    feature = _by_label(_load(
        tmp_path, _feature("CDS", "join(20..25,35..45)", label="spliced")
    ))["spliced"]
    assert feature.wraps_origin is False
    assert (feature.start, feature.end) == (19, 45)


def test_an_origin_crossing_feature_keeps_its_real_extent(tmp_path):
    """The hull of a wrapping location is the whole sequence, which says nothing.
    Storing start > end preserves the two pieces.
    """
    feature = _by_label(_load(
        tmp_path, _feature("misc_feature", "join(90..100,1..10)", label="wrapper")
    ))["wrapper"]
    assert feature.wraps_origin is True
    assert (feature.start, feature.end) == (89, 10)
    assert feature.start > feature.end


def test_a_minus_strand_wrapping_feature_sorts_its_pieces_by_position(tmp_path):
    """Minus-strand parts arrive in biological order, not reference order."""
    feature = _by_label(_load(
        tmp_path,
        _feature("misc_feature", "complement(join(98..100,1..3))", label="rev-wrap"),
    ))["rev-wrap"]
    assert feature.wraps_origin is True
    assert (feature.start, feature.end) == (97, 3)


def test_a_linear_record_never_wraps(tmp_path):
    feature = _by_label(_load(
        tmp_path,
        _feature("misc_feature", "join(90..100,1..10)", label="not-really"),
        topology="linear",
    ))["not-really"]
    assert feature.wraps_origin is False


def test_feature_spans_splits_a_wrapping_feature():
    wrapping = Feature("misc_feature", start=89, end=10, wraps_origin=True)
    assert feature_spans(wrapping, 100) == [(89, 100), (0, 10)]


def test_feature_spans_leaves_an_ordinary_feature_alone():
    assert feature_spans(Feature("CDS", 10, 40), 100) == [(10, 40)]


def test_feature_spans_clips_to_the_reference_and_drops_what_falls_outside():
    """Groups may have different reference lengths; a feature stated against a
    longer reference must not draw off the end of a shorter one.
    """
    assert feature_spans(Feature("CDS", 10, 400), 100) == [(10, 100)]
    assert feature_spans(Feature("CDS", 200, 400), 100) == []


def test_a_wrapping_feature_has_no_meaningful_len():
    """Guards the sizing trap: len() clamps to 0, so glyph width must come from
    feature_spans rather than from the dataclass.
    """
    assert len(Feature("misc_feature", 89, 10, wraps_origin=True)) == 0


# --- Labels ---------------------------------------------------------------

def test_label_prefers_the_explicit_label_qualifier(tmp_path):
    reference = _load(tmp_path, _feature(
        "CDS", "10..40", gene="gene-name", label="label-name"))
    assert "label-name" in _by_label(reference)


def test_label_falls_back_through_gene_and_product(tmp_path):
    reference = _load(tmp_path, _feature("CDS", "10..40", product="a product"))
    assert "a product" in _by_label(reference)


def test_a_snapgene_colour_note_is_not_used_as_a_label(tmp_path):
    """SnapGene round-trips styling through /note, so a note that is markup must
    not become a name even though notes are searched for one.
    """
    labels = _by_label(_load(tmp_path, _feature(
        "rep_origin", "10..40", note="color: #ffd281; direction: BOTH")))
    assert "rep_origin" in labels
    assert not any(label.startswith("color:") for label in labels)


def test_a_long_note_is_truncated_rather_than_used_whole(tmp_path):
    labels = _by_label(_load(tmp_path, _feature(
        "misc_feature", "10..40", note="z" * 120)))
    label = next(iter(labels))
    assert len(label) <= 40
    assert label.endswith("…")


def test_a_feature_with_no_usable_qualifier_is_labelled_by_type(tmp_path):
    assert "terminator" in _by_label(_load(tmp_path, _feature("terminator", "10..40")))


# --- Colour ---------------------------------------------------------------

def test_a_snapgene_colour_note_supplies_the_colour(tmp_path):
    reference = _load(tmp_path, _feature(
        "rep_origin", "10..40", label="ori",
        note="color: #ffd281; direction: BOTH"))
    assert _by_label(reference)["ori"].color == "#ffd281"


def test_an_ape_colour_qualifier_supplies_the_colour(tmp_path):
    reference = _load(tmp_path, _feature(
        "CDS", "10..40", label="ape", ApEinfo_fwdcolor="#00ffff"))
    assert _by_label(reference)["ape"].color == "#00ffff"


def test_an_ape_colour_name_is_resolved_to_a_hex(tmp_path):
    reference = _load(tmp_path, _feature(
        "CDS", "10..40", label="named", ApEinfo_fwdcolor="cyan"))
    assert _by_label(reference)["named"].color == "#00ffff"


def test_a_minus_strand_feature_prefers_the_reverse_colour(tmp_path):
    feature = _by_label(_load(tmp_path, _feature(
        "CDS", "complement(10..40)", label="rev",
        ApEinfo_fwdcolor="#00ffff", ApEinfo_revcolor="#00ff00")))["rev"]
    assert feature.strand == -1
    assert feature.color == "#00ff00"


def test_a_feature_with_no_colour_qualifier_reports_none(tmp_path):
    reference = _load(tmp_path, _feature("CDS", "10..40", label="plain"))
    assert _by_label(reference)["plain"].color is None


def test_a_nonsense_colour_value_is_ignored(tmp_path):
    reference = _load(tmp_path, _feature(
        "CDS", "10..40", label="bad", color="not a colour"))
    assert _by_label(reference)["bad"].color is None


# --- Type filtering ------------------------------------------------------

def test_source_and_primers_are_dropped_by_default(tmp_path):
    reference = _load(
        tmp_path,
        _feature("source", "1..100", organism="synthetic"),
        _feature("primer_bind", "10..30", label="fwd primer"),
        _feature("CDS", "40..70", label="orf"),
    )
    assert set(_by_label(reference)) == {"orf"}
    assert "primer_bind" in DEFAULT_SKIP_TYPES


def test_skip_types_can_be_overridden_to_keep_everything(tmp_path):
    path = tmp_path / "t.gb"
    path.write_text(_record(_feature("primer_bind", "10..30", label="fwd primer")))
    assert "fwd primer" in _by_label(load_reference(path, skip_types=()))


def test_a_gene_duplicating_a_cds_is_dropped(tmp_path):
    """Plasmid files annotate both; drawing both puts two glyphs on one ORF."""
    reference = _load(
        tmp_path,
        _feature("gene", "10..40", label="ampR gene"),
        _feature("CDS", "10..40", label="ampR"),
    )
    assert set(_by_label(reference)) == {"ampR"}


def test_a_gene_with_no_matching_cds_is_kept(tmp_path):
    """A record that annotates genes and nothing else still has reading frames."""
    reference = _load(tmp_path, _feature("gene", "10..40", label="lonely gene"))
    assert "lonely gene" in _by_label(reference)


# --- Record level --------------------------------------------------------

def test_topology_and_length_come_from_the_record(tmp_path):
    reference = _load(tmp_path, _feature("CDS", "10..40"), length=100)
    assert reference.is_circular
    assert len(reference) == 100


def test_fasta_needs_no_annotations(tmp_path):
    path = tmp_path / "r.fasta"
    path.write_text(">my-ref some description\nACGTACGT\nACGT\n")
    reference = load_reference(path)
    assert reference.seq == "ACGTACGTACGT"
    assert reference.name == "my-ref"
    assert reference.features == []
