"""Read-length distributions and their terminal histogram.

Pure arithmetic and text, so none of this needs a terminal or a real run.
"""

import gzip
import re

import pytest

from seqviewer import lengths
from seqviewer.lengths import (
    DEFAULT_BINS, DEFAULT_BULK, Bin, Binning, bin_lengths, distribution,
    histogram, read_lengths, summarise, summarise_counts, summary_lines,
)

_ANSI = re.compile(r"\x1b\[[0-9;]*[a-zA-Z]")   # colours, and the erase code


def _strip(text):
    """Return *text* without its ANSI codes."""
    return _ANSI.sub("", text)


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


# --- The two scanners -----------------------------------------------------

SCANNERS = [False] + ([True] if lengths.HAVE_NUMPY else [])


def _tally(path, fast):
    return lengths.count_lengths([path], fast=fast).items()


@pytest.mark.parametrize("body,label", [
    (b"", "an empty file"),
    (b"@r\nAAAA\n+\nIIII\n", "one record"),
    (b"@r\nAAAA\n+\nIIII", "one record, no trailing newline"),
    (b"@r\nAAAA\n+\nIIII\n@s\nAA\n+\n", "a record truncated after the plus"),
    (b"@r\nAAAA\n+\nIIII\n@s\nAA\n", "a record truncated after the sequence"),
    (b"@r\nAAAA\n+\nIIII\n@s\n", "a record truncated after the header"),
    (b"@r\n\n+\n\n", "a zero-length read"),
    (b"@r\r\nAAAA\r\n+\r\nIIII\r\n", "CRLF line endings"),
    (b"@r\r\nAAAA\r\n+\r\nIIII", "CRLF with no trailing newline"),
    (b"@r\nAAAA\n+\nIIII\n" * 500, "many identical records"),
])
def test_the_scanners_agree(tmp_path, body, label):
    """The array scanner and the pure-Python one must return the same tally,
    including where a file stops part way through a record.
    """
    path = tmp_path / "a.fastq"
    path.write_bytes(body)
    tallies = {fast: _tally(path, fast) for fast in SCANNERS}
    assert len(set(map(tuple, tallies.values()))) == 1, (label, tallies)


def test_the_scanners_agree_across_block_boundaries(tmp_path):
    """A record split across two reads is counted once.  The block size is
    lowered so the boundary is crossed by a small file.
    """
    path = tmp_path / "a.fastq"
    with open(path, "wb") as handle:
        for i in range(4000):
            n = 40 + (i % 700)
            handle.write(b"@r%d\n%s\n+\n%s\n"
                         % (i, b"A" * n, b"I" * n))

    original = lengths._BLOCK
    try:
        for block in (64, 1000, 4096, 1 << 16):
            lengths._BLOCK = block
            tallies = [_tally(path, fast) for fast in SCANNERS]
            assert len(set(map(tuple, tallies))) == 1, block
            assert sum(c for _, c in tallies[0]) == 4000
    finally:
        lengths._BLOCK = original


def test_the_scanners_agree_on_a_record_longer_than_a_block(tmp_path):
    """One read longer than the block size, so the carry path is taken."""
    path = tmp_path / "a.fastq"
    n = 5000
    path.write_bytes(b"@r\n%s\n+\n%s\n" % (b"A" * n, b"I" * n))
    original = lengths._BLOCK
    try:
        lengths._BLOCK = 512
        tallies = [_tally(path, fast) for fast in SCANNERS]
        assert len(set(map(tuple, tallies))) == 1
        assert tallies[0] == [(n, 1)]
    finally:
        lengths._BLOCK = original


def test_the_scanners_agree_on_gzip(tmp_path):
    path = tmp_path / "a.fastq.gz"
    with gzip.open(path, "wt") as handle:
        for i in range(300):
            n = 100 + i
            handle.write(f"@r{i}\n{'A' * n}\n+\n{'I' * n}\n")
    tallies = [_tally(path, fast) for fast in SCANNERS]
    assert len(set(map(tuple, tallies))) == 1
    assert sum(c for _, c in tallies[0]) == 300


def test_the_tally_matches_reading_every_length(tmp_path):
    """count_lengths and read_lengths describe the same file."""
    path = _fastq(tmp_path / "a.fastq", [800] * 50 + [90, 4000, 4000])
    streamed = sorted(lengths.read_lengths([path]))
    for fast in SCANNERS:
        counts = lengths.count_lengths([path], fast=fast)
        assert counts.total == len(streamed)
        assert counts.bases == sum(streamed)
        assert counts.shortest == streamed[0]
        assert counts.longest == streamed[-1]


def test_requiring_the_array_scanner_without_numpy_is_an_error(monkeypatch,
                                                              tmp_path):
    path = _fastq(tmp_path / "a.fastq", [10])
    monkeypatch.setattr(lengths, "HAVE_NUMPY", False)
    with pytest.raises(RuntimeError):
        lengths.count_lengths([path], fast=True)


# --- The tally ------------------------------------------------------------

def test_the_tally_is_bounded_by_the_length_range_not_the_read_count():
    """Why a multi-gigabyte file can be scanned in fixed memory."""
    few = lengths.LengthCounts.from_lengths([800] * 10)
    many = lengths.LengthCounts.from_lengths([800] * 1_000_000)
    assert len(few.items()) == len(many.items()) == 1
    assert many.total == 1_000_000


def test_the_tally_reports_the_same_figures_as_a_list_of_lengths():
    values = [780 + (i * 7) % 41 for i in range(1000)] + [4000]
    counts = lengths.LengthCounts.from_lengths(values)
    direct = summarise(values)
    assert summarise_counts(counts) == direct
    assert counts.quantile(0.5) == bin_lengths(values).low


def test_a_quantile_of_an_empty_tally_is_zero():
    counts = lengths.LengthCounts()
    assert counts.total == 0 and counts.quantile(50) == 0
    assert counts.median == 0 and counts.n50 == 0 and counts.mean == 0.0


# --- Colour ---------------------------------------------------------------

def test_a_palette_does_not_change_the_visible_width():
    """Codes are added after the columns are laid out."""
    binning = bin_lengths([800] * 500 + [80, 9000])
    plain = histogram(binning, width=80)
    painted = histogram(binning, width=80, palette=lengths.PALETTE)
    assert [_strip(line) for line in painted] == plain


def test_a_palette_marks_the_peak_apart_from_the_other_bars():
    binning = bin_lengths([800] * 500 + [400] * 5, bulk=100)
    painted = histogram(binning, width=80, palette=lengths.PALETTE)
    assert sum(lengths.PALETTE.peak in line for line in painted) == 1


def test_a_palette_dims_the_clipped_tails():
    binning = bin_lengths([800] * 1000 + [80, 9000])
    painted = histogram(binning, width=80, palette=lengths.PALETTE)
    assert sum(lengths.PALETTE.tail in line for line in painted) == 3  # + header


def test_no_codes_are_written_without_a_palette():
    lines = histogram(bin_lengths([800] * 20 + [90, 9000]), width=80)
    assert not any("\x1b" in line for line in lines)

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


# --- Progress and formatting ----------------------------------------------

class _Terminal:
    """A stream that reports itself a terminal and keeps what was written."""

    def __init__(self):
        self.chunks = []

    def isatty(self):
        return True

    def write(self, text):
        self.chunks.append(text)

    def flush(self):
        pass

    @property
    def text(self):
        return "".join(self.chunks)


def test_a_small_run_draws_no_progress():
    """Below the threshold the scan is over before a bar could be read."""
    from seqviewer.cli import PROGRESS_AFTER_BYTES, ScanProgress

    term = _Terminal()
    bar = ScanProgress(80, stream=term)
    bar.update(1 << 10, PROGRESS_AFTER_BYTES - 1, 4)
    assert term.text == ""


def test_a_large_run_draws_progress():
    from seqviewer.cli import PROGRESS_AFTER_BYTES, ScanProgress

    term = _Terminal()
    bar = ScanProgress(80, stream=term)
    bar.update(PROGRESS_AFTER_BYTES, PROGRESS_AFTER_BYTES * 4, 12_345)
    assert "25%" in term.text
    assert "12k reads" in term.text
    assert "█" in term.text and "░" in term.text


def test_progress_is_recorded_even_when_it_is_not_drawn():
    """The scan is reported once it finishes whether or not a bar was shown."""
    from seqviewer.cli import ScanProgress

    term = _Terminal()
    bar = ScanProgress(80, stream=term, enabled=False)
    bar.update(5 << 30, 5 << 30, 900)
    assert term.text == ""
    assert bar.bytes == 5 << 30 and bar.reads == 900


def test_progress_redraws_are_throttled():
    from seqviewer.cli import PROGRESS_AFTER_BYTES, ScanProgress

    total = PROGRESS_AFTER_BYTES * 10
    term = _Terminal()
    bar = ScanProgress(80, stream=term)
    for i in range(1, 40):
        bar.update(i * (total // 100), total, i * 100)
    assert 1 <= term.text.count("\r") <= 4      # not once per call


def test_progress_always_draws_the_last_position():
    from seqviewer.cli import PROGRESS_AFTER_BYTES, ScanProgress

    total = PROGRESS_AFTER_BYTES * 10
    term = _Terminal()
    bar = ScanProgress(80, stream=term)
    bar.update(total // 2, total, 10)
    bar.update(total, total, 20)                # not throttled away
    assert "100%" in term.text


def test_clearing_progress_erases_the_line():
    from seqviewer.cli import PROGRESS_AFTER_BYTES, ScanProgress

    term = _Terminal()
    bar = ScanProgress(80, stream=term)
    bar.update(PROGRESS_AFTER_BYTES, PROGRESS_AFTER_BYTES, 5)
    bar.clear()
    assert term.text.endswith("\r\x1b[2K")


def test_a_progress_line_fits_the_width():
    from seqviewer.cli import PROGRESS_AFTER_BYTES, ScanProgress

    for width in (32, 60, 80, 200):
        term = _Terminal()
        bar = ScanProgress(width, stream=term)
        bar.update(PROGRESS_AFTER_BYTES, PROGRESS_AFTER_BYTES * 3, 1_500_000)
        drawn = _strip(term.text).replace("\r", "")
        assert len(drawn) <= max(32, width)


def test_a_quick_scan_is_not_reported():
    from seqviewer.cli import ScanProgress

    bar = ScanProgress(80, stream=_Terminal())
    bar.update(1 << 30, 1 << 30, 10)            # took no measurable time
    assert bar.note() is None


def test_byte_and_read_counts_are_abbreviated():
    from seqviewer.cli import _duration, _si_bytes, _si_reads

    assert _si_bytes(512) == "512 B"
    assert _si_bytes(1 << 20) == "1.0 MiB"
    assert _si_bytes(3 << 30) == "3.0 GiB"
    assert _si_reads(42) == "42"
    assert _si_reads(4_200) == "4k"
    assert _si_reads(4_200_000) == "4.2M"
    assert _duration(0.44) == "0.4s"
    assert _duration(42) == "42s"
    assert _duration(125) == "2m05s"


# --- The new flags --------------------------------------------------------

def test_colour_is_left_out_when_stdout_is_not_a_terminal(tmp_path, capsys):
    """pytest captures stdout, so this is the path a pipe takes."""
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [800] * 40)
    lengths_main([str(tmp_path), "--width", "80"])
    assert "\x1b" not in capsys.readouterr().out


def test_no_color_is_accepted(tmp_path, capsys):
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [800] * 40)
    assert lengths_main([str(tmp_path), "--no-color"]) == 0
    assert "\x1b" not in capsys.readouterr().out


def test_the_slow_scanner_gives_the_same_output(tmp_path, capsys):
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [780 + (i * 7) % 41 for i in range(400)]
           + [90, 4000])
    lengths_main([str(tmp_path), "--width", "80"])
    fast = capsys.readouterr().out
    lengths_main([str(tmp_path), "--width", "80", "--slow"])
    assert capsys.readouterr().out == fast


def test_no_progress_is_accepted(tmp_path, capsys):
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [800] * 20)
    assert lengths_main([str(tmp_path), "--no-progress"]) == 0
