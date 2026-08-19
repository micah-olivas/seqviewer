"""Tests for the command line driver's read collection and its output paths.

These cover the part that decides *which* reads reach the aligner — finding the
FASTQs at a path and pooling or splitting them — and where the pages it writes
go.  All of it is pure, so none of it needs pysam, minimap2, or samtools; the one
end-to-end test that does need them skips when they are absent.
"""

import random
import shutil
from pathlib import Path

import pytest

from seqviewer.cli import (
    build_parser, collect_samples, fastq_paths, main, read_fastq, summary_path,
)
from seqviewer.summary import DEFAULT_MIN_COUNT, DEFAULT_MIN_FRACTION

#: Fixed, so the end-to-end fixture is the same sequence every run.
RNG = random.Random(5)

READS = [("r1", "ACGTACGT"), ("r2", "ACGTTTTT"), ("r3", "AC")]


def _fastq(path, reads=READS):
    """Write a FASTQ of (name, seq) pairs, quality all 'I'."""
    path.write_text("".join(
        f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in reads
    ))
    return path


def test_a_directory_yields_its_fastqs_sorted(tmp_path):
    _fastq(tmp_path / "b.fastq")
    _fastq(tmp_path / "a.fq")
    assert [p.name for p in fastq_paths(tmp_path)] == ["a.fq", "b.fastq"]


def test_non_fastq_and_hidden_files_are_skipped(tmp_path):
    """A cloud-synced share carries .DS_Store and ``._`` stubs; neither is reads."""
    _fastq(tmp_path / "reads.fastq")
    (tmp_path / ".DS_Store").write_text("junk")
    (tmp_path / "._reads.fastq").write_text("junk")
    (tmp_path / "notes.txt").write_text("junk")
    assert [p.name for p in fastq_paths(tmp_path)] == ["reads.fastq"]


def test_a_single_file_is_still_accepted(tmp_path):
    path = _fastq(tmp_path / "reads.fastq")
    assert fastq_paths(path) == [path]


def test_a_missing_path_yields_nothing(tmp_path):
    assert fastq_paths(tmp_path / "absent") == []


def test_pooling_puts_every_file_in_one_group(tmp_path):
    _fastq(tmp_path / "a.fastq")
    _fastq(tmp_path / "b.fastq")
    samples = collect_samples(fastq_paths(tmp_path))
    assert len(samples) == 1
    label, reads, seen = samples[0]
    assert label is None                      # named for the reference, not a file
    assert len(reads) == seen == 2 * len(READS)


def test_per_file_gives_one_group_per_file_labelled_by_stem(tmp_path):
    _fastq(tmp_path / "a.fastq")
    _fastq(tmp_path / "b.fastq")
    samples = collect_samples(fastq_paths(tmp_path), pooled=False)
    assert [label for label, _, _ in samples] == ["a", "b"]
    assert all(len(reads) == len(READS) for _, reads, _ in samples)


def test_a_cap_bounds_the_pool_not_each_file(tmp_path):
    """Pooled, --max is what the one group draws, across all the files."""
    _fastq(tmp_path / "a.fastq")
    _fastq(tmp_path / "b.fastq")
    (_, reads, seen), = collect_samples(fastq_paths(tmp_path), limit=4)
    assert len(reads) == 4
    assert seen == 2 * len(READS)             # the depth behind the sample


def test_a_cap_bounds_each_group_when_split(tmp_path):
    _fastq(tmp_path / "a.fastq")
    _fastq(tmp_path / "b.fastq")
    samples = collect_samples(fastq_paths(tmp_path), limit=2, pooled=False)
    assert [len(reads) for _, reads, _ in samples] == [2, 2]


def test_a_file_with_no_surviving_reads_drops_out_of_the_split(tmp_path):
    _fastq(tmp_path / "a.fastq")
    _fastq(tmp_path / "b.fastq", [("short", "AC")])
    samples = collect_samples(fastq_paths(tmp_path), min_len=4, pooled=False)
    assert [label for label, _, _ in samples] == ["a"]


def test_pooling_returns_nothing_when_no_read_survives(tmp_path):
    _fastq(tmp_path / "a.fastq", [("short", "AC")])
    assert collect_samples(fastq_paths(tmp_path), min_len=4) == []


def test_a_sample_is_drawn_from_the_whole_file_not_its_front(tmp_path):
    """The point of sampling: a cap must not degenerate into a head.

    Reads are numbered in file order, so a head of 20 would hold none above 19.
    A uniform sample of a thousand reaches the end of the file, and its mean sits
    near the middle rather than near the front.
    """
    _fastq(tmp_path / "a.fastq",
           [(str(i), "ACGT") for i in range(1000)])
    (_, reads, seen), = collect_samples(fastq_paths(tmp_path), limit=20)
    positions = sorted(int(r.name) for r in reads)
    assert len(reads) == 20 and seen == 1000
    assert positions != list(range(20))       # not the front of the file
    assert max(positions) > 900               # and it reaches the back
    assert 300 < sum(positions) / len(positions) < 700


def test_pooled_sampling_spans_every_file(tmp_path):
    """A pool sampled per file would draw one barcode; this draws from both."""
    _fastq(tmp_path / "a.fastq", [(f"a{i}", "ACGT") for i in range(500)])
    _fastq(tmp_path / "b.fastq", [(f"b{i}", "ACGT") for i in range(500)])
    (_, reads, _), = collect_samples(fastq_paths(tmp_path), limit=50)
    stems = {r.name[0] for r in reads}
    assert stems == {"a", "b"}


def test_a_sample_smaller_than_the_cap_keeps_every_read(tmp_path):
    _fastq(tmp_path / "a.fastq")
    (_, reads, seen), = collect_samples(fastq_paths(tmp_path), limit=100)
    assert len(reads) == seen == len(READS)


def test_sampling_the_same_directory_twice_gives_the_same_reads(tmp_path):
    """The seed is fixed, so a page can be regenerated as it was drawn."""
    _fastq(tmp_path / "a.fastq", [(str(i), "ACGT") for i in range(200)])
    first, second = (collect_samples(fastq_paths(tmp_path), limit=10)[0][1]
                     for _ in range(2))
    assert [r.name for r in first] == [r.name for r in second]


def test_a_falsy_cap_draws_everything(tmp_path):
    _fastq(tmp_path / "a.fastq", [(str(i), "ACGT") for i in range(50)])
    (_, reads, seen), = collect_samples(fastq_paths(tmp_path), limit=0)
    assert len(reads) == seen == 50


def test_reads_shorter_than_the_minimum_are_dropped(tmp_path):
    path = _fastq(tmp_path / "reads.fastq")
    assert [r.name for r in read_fastq(path, min_len=4)] == ["r1", "r2"]


def test_reads_can_be_filtered_by_name(tmp_path):
    path = _fastq(tmp_path / "reads.fastq")
    assert [r.name for r in read_fastq(path, name_contains="r2")] == ["r2"]


def test_gzipped_fastqs_are_read_and_found(tmp_path):
    import gzip

    path = tmp_path / "reads.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r1\nACGT\n+\nIIII\n")
    assert fastq_paths(tmp_path) == [path]
    assert [r.name for r in read_fastq(path)] == ["r1"]


# --- The summarized page --------------------------------------------------

def test_the_summary_page_sits_beside_the_pileup():
    """A directory of runs sorts each pair together."""
    assert summary_path("1A12.html").name == "1A12.summary.html"


def test_the_summary_path_is_derived_from_the_stem_not_appended():
    """Not 1A12.html.summary.html, which sorts away from its own pileup."""
    assert summary_path(Path("runs/1A12.html")) == Path("runs/1A12.summary.html")


def test_a_path_with_dots_in_its_name_keeps_them():
    assert summary_path("plate.1.A12.html").name == "plate.1.A12.summary.html"


def test_the_cli_calling_thresholds_do_not_drift_from_the_reducer():
    """The flags default to the module's own floors rather than restating them."""
    defaults = {a.dest: a.default for a in build_parser()._actions}
    assert defaults["variant_freq"] == DEFAULT_MIN_FRACTION
    assert defaults["variant_reads"] == DEFAULT_MIN_COUNT
    assert defaults["summary"] is False


@pytest.mark.skipif(shutil.which("minimap2") is None
                    or shutil.which("samtools") is None,
                    reason="needs minimap2 and samtools on PATH")
def test_summary_writes_a_second_page_reporting_the_planted_variant(tmp_path):
    """End to end through the real aligner, not a synthetic grid.

    Also pins the two limits the reducer documents: a deletion is recovered from
    the alignment, and a planted insertion is not, because a grid is exactly as
    wide as the reference and has nowhere to put one.
    """
    ref = "".join(RNG.choice("ACGT") for _ in range(900))
    cds_start = 150
    snv = cds_start + 149            # inside the reading frame
    deletion = cds_start + 300

    (tmp_path / "ref.fasta").write_text(">pTEST\n" + ref + "\n")
    reads = []
    for i in range(20):
        reads.append((f"clean{i}", ref))
    for i in range(20):
        alt = "A" if ref[snv] != "A" else "C"
        reads.append((f"snv{i}", ref[:snv] + alt + ref[snv + 1:]))
    for i in range(20):
        reads.append((f"del{i}", ref[:deletion] + ref[deletion + 3:]))
    for i in range(20):
        reads.append((f"ins{i}", ref[:600] + "GGG" + ref[600:]))
    _fastq(tmp_path / "reads.fastq", reads)

    out = tmp_path / "page"
    code = main([str(tmp_path / "reads.fastq"), str(tmp_path / "ref.fasta"),
                 str(out), "--summary", "--variant-freq", "0.1"])
    assert code == 0

    pileup = out.with_suffix(".html")
    summary = summary_path(pileup)
    assert pileup.exists()
    assert summary.exists()

    html = summary.read_text()
    assert html.startswith("<!DOCTYPE html>")
    assert "SNV" in html
    assert "Deletion" in html
    # The grid cannot carry an insertion, so none is reported however many
    # reads carry one.  This is the documented limit, pinned.
    assert "Insertion" not in html


def test_no_summary_page_is_written_unless_asked_for(tmp_path):
    """The flag is opt-in; a plain run leaves no stray file."""
    assert not summary_path(tmp_path / "page.html").exists()
