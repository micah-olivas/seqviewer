"""Read-length distributions and their terminal histogram.

Pure arithmetic and text, so none of this needs a terminal or a real run.
"""

import gzip

from seqviewer.lengths import (
    DEFAULT_BINS, DEFAULT_BULK, Bin, Binning, bin_lengths, distribution,
    histogram, read_lengths, summarise, summary_lines,
)


def _fastq(path, lengths, prefix="r"):
    with open(path, "w") as handle:
        for i, n in enumerate(lengths):
            handle.write(f"@{prefix}{i}\n{'A' * n}\n+\n{'I' * n}\n")
    return path


# --- Reading --------------------------------------------------------------

def test_lengths_come_back_in_file_order(tmp_path):
    path = _fastq(tmp_path / "a.fastq", [10, 30, 20])
    assert list(read_lengths([path])) == [10, 30, 20]


def test_several_files_are_read_as_one_run(tmp_path):
    a = _fastq(tmp_path / "a.fastq", [10, 20])
    b = _fastq(tmp_path / "b.fastq", [30])
    assert sorted(read_lengths([a, b])) == [10, 20, 30]


def test_a_gzipped_file_is_read_directly(tmp_path):
    path = tmp_path / "a.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r0\nAAAA\n+\nIIII\n")
    assert list(read_lengths([path])) == [4]


def test_a_truncated_final_record_is_dropped(tmp_path):
    """The state a file interrupted mid-write is left in."""
    path = tmp_path / "a.fastq"
    path.write_text("@r0\nAAAA\n+\nIIII\n@r1\nAAAAAA\n+\n")
    assert list(read_lengths([path])) == [4]


def test_reading_holds_no_sequences(tmp_path):
    """A generator, so a deep run does not have to fit in memory."""
    import types

    path = _fastq(tmp_path / "a.fastq", [10])
    assert isinstance(read_lengths([path]), types.GeneratorType)


# --- Summary --------------------------------------------------------------

def test_an_empty_run_summarises_to_zero():
    s = summarise([])
    assert s.empty and s.reads == 0 and s.n50 == 0


def test_the_median_of_an_even_count_is_the_midpoint():
    assert summarise([10, 20, 30, 40]).median == 25


def test_the_median_of_an_odd_count_is_the_middle_value():
    assert summarise([10, 20, 100]).median == 20


def test_n50_is_the_length_holding_half_the_bases():
    """One 100 and ten 10s: 200 bases, so half is 100 and the single long read
    reaches it on its own.
    """
    s = summarise([100] + [10] * 10)
    assert s.bases == 200
    assert s.n50 == 100


def test_n50_and_mean_differ_on_a_mixed_run():
    """The distinction the two figures exist to show.  Twenty 1000s hold 20,000
    of 30,000 bases, so the N50 sits among them while the mean sits near the 50s.
    """
    s = summarise([1000] * 20 + [50] * 200)
    assert s.bases == 30_000
    assert s.mean < 150
    assert s.n50 == 1000


def test_the_summary_covers_every_read_not_only_the_axis():
    lengths = [100] * 500 + [9000]
    summary, binning = distribution(lengths, bulk=99)
    assert summary.longest == 9000
    assert binning.high < 9000
    assert binning.above == 1


# --- Binning --------------------------------------------------------------

def test_bins_are_equal_width():
    binning = bin_lengths(list(range(100, 200)), count=5, bulk=100)
    widths = {b.high - b.low for b in binning.bins}
    assert len(widths) == 1


def test_every_read_lands_in_exactly_one_bin():
    lengths = [10, 10, 25, 40, 55, 70]
    binning = bin_lengths(lengths, count=4, bulk=100)
    assert sum(b.count for b in binning.bins) == len(lengths)


def test_reads_of_one_length_give_one_bin():
    binning = bin_lengths([500] * 20, count=DEFAULT_BINS)
    assert len(binning.bins) == 1
    assert binning.bins[0].count == 20


def test_an_empty_run_gives_no_bins():
    assert bin_lengths([]).bins == []


def test_the_axis_excludes_the_long_tail_by_default():
    """The case the clipping exists for: a few concatemers among a tight run."""
    lengths = [800] * 1000 + [4000, 4200, 8000]
    binning = bin_lengths(lengths)
    assert binning.above == 3
    assert binning.high < 4000
    assert binning.longest == 8000


def test_clipping_keeps_the_product_resolved():
    """The comparison the default is chosen on, and the figures the README
    quotes.  A product spread over 780-820 bp with three concatemers: over the
    full range the product occupies one bin of twenty-four, and over the default
    axis it spreads across all twenty-one.
    """
    lengths = [780 + (i % 41) for i in range(1000)] + [1640, 2460, 4100]

    full = bin_lengths(lengths, count=24, bulk=100)
    assert len(full.bins) == 24
    assert sum(1 for b in full.bins if b.count) == 4
    assert max(b.count for b in full.bins) == 1000

    clipped = bin_lengths(lengths, count=24, bulk=99)
    assert len(clipped.bins) == 21
    assert all(b.count for b in clipped.bins)
    assert clipped.above == 3


def test_clipped_reads_are_counted_not_discarded():
    lengths = list(range(1, 1001))
    binning = bin_lengths(lengths, bulk=90)
    inside = sum(b.count for b in binning.bins)
    assert inside + binning.below + binning.above == len(lengths)
    assert binning.below and binning.above


def test_full_bulk_clips_nothing():
    lengths = [800] * 100 + [9000]
    binning = bin_lengths(lengths, bulk=100)
    assert not binning.clipped
    assert binning.high >= 9000


def test_the_default_bulk_is_the_documented_one():
    assert DEFAULT_BULK == 99.0
    assert DEFAULT_BINS == 24


def test_clipping_that_would_empty_the_axis_falls_back_to_the_full_range():
    """Two reads and a 99% window can leave nothing inside."""
    binning = bin_lengths([10, 10_000], bulk=1)
    assert sum(b.count for b in binning.bins) == 2
    assert not binning.clipped


# --- Histogram ------------------------------------------------------------

def test_an_empty_run_says_so():
    assert histogram(Binning()) == ["no reads"]


def test_a_bin_holding_reads_is_never_blank():
    """Rounding to nothing would read as an empty bin, a different fact."""
    lengths = [800] * 10_000 + [400]
    binning = bin_lengths(lengths, bulk=100)
    lines = histogram(binning, width=80)
    drawn = [ln for ln in lines[1:] if ln.split()[1:2]]
    assert all("█" in ln or "▏" in ln or ln.rstrip().endswith("0")
               for ln in drawn)


def test_the_tails_get_their_own_labelled_rows():
    lengths = [800] * 1000 + [80, 6000]
    text = "\n".join(histogram(bin_lengths(lengths), width=90))
    assert "shorter, down to 80" in text
    assert "longer, up to 6,000" in text


def test_the_tails_do_not_set_the_bar_scale():
    """The tallest binned count sets the scale.  A flat run clipped hard puts
    far more reads in each tail than in any bin, and the tail bars still stop at
    the axis width rather than running past it.
    """
    binning = bin_lengths(list(range(1, 1001)), bulk=50)
    tallest = max(b.count for b in binning.bins)
    assert binning.below > tallest and binning.above > tallest

    lines = histogram(binning, width=70)
    assert max(len(line) for line in lines) <= 70


def test_no_line_is_wider_than_asked_for():
    lengths = [800] * 500 + [80, 9000]
    for width in (60, 80, 120):
        lines = histogram(bin_lengths(lengths), width=width)
        assert max(len(line) for line in lines) <= width


def test_the_log_scale_lifts_the_small_bins():
    lengths = [800] * 5000 + [400] * 3
    binning = bin_lengths(lengths, bulk=100)
    linear = histogram(binning, width=80)
    logged = histogram(binning, width=80, log=True)
    assert sum(ln.count("█") for ln in logged) > \
           sum(ln.count("█") for ln in linear)


def test_the_axis_labels_stay_linear_under_a_log_scale():
    """Only bar length changes; a reader is never misled about where a bin is."""
    binning = bin_lengths([800] * 100 + [400] * 3, bulk=100)
    labels = lambda lines: [ln.split()[0] for ln in lines[1:]]
    assert labels(histogram(binning, 80)) == \
           labels(histogram(binning, 80, log=True))


# --- Summary lines --------------------------------------------------------

def test_the_figures_are_reported_over_every_read():
    lengths = [800] * 1000 + [9000]
    summary, binning = distribution(lengths)
    text = "\n".join(summary_lines(summary, binning))
    assert "9,000" in text            # the true maximum, past the axis
    assert "1,001 reads" in text


def test_a_clipped_axis_says_what_it_left_out():
    lengths = [800] * 1000 + [9000]
    summary, binning = distribution(lengths)
    text = "\n".join(summary_lines(summary, binning))
    assert "--bulk 100" in text
    assert "holding 1,000 of 1,001 reads" in text


def test_an_unclipped_axis_says_nothing_extra():
    summary, binning = distribution([800] * 10, bulk=100)
    assert len(summary_lines(summary, binning)) == 2


def test_an_empty_run_reports_no_reads():
    assert summary_lines(summarise([])) == ["no reads"]


# --- The command line -----------------------------------------------------

def test_the_command_is_wired_up():
    from seqviewer.cli import lengths_main

    assert callable(lengths_main)


def test_it_reports_a_directory_with_no_fastqs(tmp_path, capsys):
    from seqviewer.cli import lengths_main

    assert lengths_main([str(tmp_path)]) == 1


def test_it_draws_a_histogram_for_a_run(tmp_path, capsys):
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [800] * 200 + [90, 5000])
    assert lengths_main([str(tmp_path), "--width", "80"]) == 0
    out = capsys.readouterr().out
    assert "reads" in out and "█" in out
    assert "202 reads" in out
    assert "5,000" in out                      # the tail is reported


def test_per_file_draws_one_histogram_each(tmp_path, capsys):
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [800] * 20)
    _fastq(tmp_path / "b.fastq", [400] * 20)
    assert lengths_main([str(tmp_path), "--per-file", "--width", "80"]) == 0
    out = capsys.readouterr().out
    assert "a.fastq" in out and "b.fastq" in out
    assert out.count("20 reads") == 2
