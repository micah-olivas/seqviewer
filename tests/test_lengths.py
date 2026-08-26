"""Read-length distributions and their terminal histogram.

Pure arithmetic and text, so none of this needs a terminal or a real run.
"""

import gzip
import os
import re
import time

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


# --- The live histogram ---------------------------------------------------

def _live(width=80, bins=6, rows=40, **kw):
    from seqviewer.cli import LiveView

    term = _Terminal()
    view = LiveView(width, bins, lengths.DEFAULT_BULK, False,
                    stream=term, rows=rows, **kw)
    return view, term


def test_the_live_view_needs_a_terminal():
    from seqviewer.cli import LiveView

    class Piped(_Terminal):
        def isatty(self):
            return False

    view = LiveView(80, 6, 99.0, False, stream=Piped(), rows=40)
    assert not view.enabled


def test_the_live_view_needs_a_window_tall_enough_for_the_frame():
    """Redrawing a frame taller than the window would scroll it apart."""
    tall, _ = _live(bins=40, rows=20)
    short, _ = _live(bins=6, rows=20)
    assert not tall.enabled
    assert short.enabled


def test_the_live_view_stays_quiet_on_a_small_run():
    from seqviewer.cli import PROGRESS_AFTER_BYTES

    view, term = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 100)
    view.update(1 << 10, PROGRESS_AFTER_BYTES - 1, 100, lambda: counts)
    assert term.text == ""


def test_the_live_view_draws_the_distribution_so_far():
    from seqviewer.cli import PROGRESS_AFTER_BYTES

    view, term = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 50 + [400] * 10)
    view.update(PROGRESS_AFTER_BYTES, PROGRESS_AFTER_BYTES * 4, 60,
                lambda: counts)
    drawn = _strip(term.text)
    assert "bp" in drawn and "reads" in drawn
    assert "60 reads" in drawn                  # the figures so far
    assert "25%" in drawn                       # and the progress line
    assert "█" in drawn


def test_a_redraw_moves_the_cursor_back_over_the_last_frame():
    from seqviewer.cli import PROGRESS_AFTER_BYTES, LIVE_INTERVAL

    view, term = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 50)
    view.update(PROGRESS_AFTER_BYTES, PROGRESS_AFTER_BYTES * 4, 50,
                lambda: counts)
    height = view.height
    assert height > 1
    term.chunks.clear()
    view.last -= LIVE_INTERVAL * 2              # let the next frame be due
    view.update(PROGRESS_AFTER_BYTES * 2, PROGRESS_AFTER_BYTES * 4, 90,
                lambda: counts)
    assert term.text.startswith(f"\x1b[{height}A")


def test_every_redrawn_line_is_erased_first():
    """A shorter line must not leave the end of the last one behind."""
    from seqviewer.cli import PROGRESS_AFTER_BYTES

    view, term = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 50)
    view.update(PROGRESS_AFTER_BYTES, PROGRESS_AFTER_BYTES * 2, 50,
                lambda: counts)
    body = term.text.split("A", 1)[-1]
    assert body.count("\x1b[2K") == view.height


def test_the_final_frame_leaves_the_cursor_below_the_histogram():
    """Otherwise the next group, or the shell prompt, prints into it."""
    from seqviewer.cli import PROGRESS_AFTER_BYTES

    view, term = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 50 + [90, 9000])
    view.update(PROGRESS_AFTER_BYTES, PROGRESS_AFTER_BYTES * 2, 52,
                lambda: counts)
    term.chunks.clear()
    view.finish(counts)
    assert not re.search(r"\x1b\[\d+A$", term.text)
    assert term.text.endswith("\n")


def test_the_final_frame_drops_the_progress_line():
    from seqviewer.cli import PROGRESS_AFTER_BYTES

    view, term = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 50)
    view.update(PROGRESS_AFTER_BYTES, PROGRESS_AFTER_BYTES * 2, 50,
                lambda: counts)
    term.chunks.clear()
    view.finish(counts)
    assert "%" not in _strip(term.text)


def test_a_view_that_never_drew_prints_the_histogram_once(capsys):
    """The path a pipe takes: one write, and no partial frames in it."""
    view, _ = _live(enabled=False)
    counts = lengths.LengthCounts.from_lengths([800] * 50 + [90, 9000])
    view.finish(counts)
    out = capsys.readouterr().out
    assert "\x1b[" not in out
    assert out.count("N50") == 1                # the frame, written once
    assert "9,000" in out                       # the tail, past the axis


def test_the_settled_frame_holds_the_histogram_and_the_figures():
    view, _ = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 50 + [90, 9000])
    lines = view.frame(counts, ["tail line"])
    text = "\n".join(lines)
    assert "bp" in text and "N50" in text and "52 reads" in text
    assert lines[-1] == "tail line"


def test_a_trailing_string_is_taken_whole_not_per_letter():
    view, _ = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 10)
    assert view.frame(counts, "one line")[-1] == "one line"


def test_a_live_frame_leaves_the_figures_out():
    """They churn under the histogram while the scan runs."""
    view, _ = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 50 + [90, 9000])
    live = "\n".join(view.frame(counts, ["bar"], stats=False))
    assert "N50" not in live and "bases" not in live
    assert "bp" in live                         # the histogram is still there
    assert "N50" in "\n".join(view.frame(counts))


def test_an_empty_frame_says_no_reads_once():
    view, _ = _live()
    assert view.frame(lengths.LengthCounts()) == ["", "no reads"]


def test_the_live_view_is_off_when_asked(tmp_path, capsys):
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [800] * 30)
    assert lengths_main([str(tmp_path), "--no-live"]) == 0
    assert "\x1b[" not in capsys.readouterr().out


def test_the_snapshot_grows_as_the_scan_proceeds(tmp_path):
    """What the live view draws from: each callback sees more than the last."""
    path = _fastq(tmp_path / "a.fastq", [500 + i % 200 for i in range(20000)])
    seen = []

    def progress(done, total, reads, snapshot=None):
        if snapshot is not None:
            seen.append(snapshot().total)

    original = lengths._BLOCK
    try:
        lengths._BLOCK = 4096                   # many blocks from a small file
        counts = lengths.count_lengths([path], progress=progress)
    finally:
        lengths._BLOCK = original

    assert len(seen) > 3
    assert seen == sorted(seen)                 # never goes backwards
    assert seen[-1] == counts.total == 20000


# --- The seqview dispatcher -----------------------------------------------

def test_seqview_lists_its_commands(capsys):
    from seqviewer.cli import seqview_main

    assert seqview_main(["--help"]) == 0
    out = capsys.readouterr().out
    assert "pileup" in out and "lengths" in out


def test_seqview_with_no_command_explains_itself(capsys):
    """Usage goes to stderr with a non-zero status, as a misuse."""
    from seqviewer.cli import seqview_main

    assert seqview_main([]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "usage: seqview" in captured.err


def test_seqview_rejects_an_unknown_command(capsys):
    from seqviewer.cli import seqview_main

    assert seqview_main(["pileip"]) == 2
    assert "unknown command 'pileip'" in capsys.readouterr().err


def test_seqview_lengths_runs_the_lengths_command(tmp_path, capsys):
    from seqviewer.cli import seqview_main

    _fastq(tmp_path / "a.fastq", [800] * 40 + [90, 5000])
    assert seqview_main(["lengths", str(tmp_path), "--width", "80"]) == 0
    out = capsys.readouterr().out
    assert "42 reads" in out and "█" in out


def test_seqview_lengths_and_the_hyphenated_command_agree(tmp_path, capsys):
    from seqviewer.cli import lengths_main, seqview_main

    _fastq(tmp_path / "a.fastq", [800] * 40 + [90, 5000])
    seqview_main(["lengths", str(tmp_path), "--width", "80"])
    through_seqview = capsys.readouterr().out
    lengths_main([str(tmp_path), "--width", "80"])
    assert capsys.readouterr().out == through_seqview


def test_a_subcommand_names_itself_in_its_usage(capsys):
    from seqviewer.cli import seqview_main

    with pytest.raises(SystemExit):
        seqview_main(["lengths", "--help"])
    assert "seqview lengths" in capsys.readouterr().out


def test_the_hyphenated_command_still_names_itself(capsys):
    """Scripts calling the old name should see the old name in errors."""
    from seqviewer.cli import lengths_main

    with pytest.raises(SystemExit):
        lengths_main(["--help"])
    assert "seqviewer-lengths" in capsys.readouterr().out


def test_a_blank_separator_survives_the_tail(): 
    """A blank line before the progress line, so it does not butt the last bar."""
    view, _ = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 20)
    lines = view.frame(counts, ["", "progress"], stats=False)
    assert lines[-2:] == ["", "progress"]


def test_a_missing_tail_line_is_left_out():
    view, _ = _live()
    counts = lengths.LengthCounts.from_lengths([800] * 20)
    assert view.frame(counts, [None])[-1] != None
    assert None not in view.frame(counts, [None])


# --- Stopping a scan ------------------------------------------------------

def test_a_scan_stops_when_asked(tmp_path):
    """The tally then covers the blocks counted, not the whole file."""
    path = _fastq(tmp_path / "a.fastq", [800] * 40000)
    calls = []

    def stop():
        calls.append(1)
        return len(calls) > 2                   # let two blocks through

    original = lengths._BLOCK
    try:
        lengths._BLOCK = 8192
        partial = lengths.count_lengths([path], stop=stop)
        whole = lengths.count_lengths([path])
    finally:
        lengths._BLOCK = original

    assert 0 < partial.total < whole.total
    assert whole.total == 40000


def test_stopping_leaves_the_remaining_files_unopened(tmp_path):
    a = _fastq(tmp_path / "a.fastq", [800] * 20000)
    b = _fastq(tmp_path / "b.fastq", [400] * 20000)
    original = lengths._BLOCK
    try:
        lengths._BLOCK = 8192
        counts = lengths.count_lengths([a, b], stop=lambda: True)
    finally:
        lengths._BLOCK = original
    # Only the first file's reads, and none of the second file's 400s.
    assert counts.total < 40000
    assert all(length != 400 for length, _ in counts.items())


def test_both_scanners_stop(tmp_path):
    path = _fastq(tmp_path / "a.fastq", [800] * 20000)
    original = lengths._BLOCK
    try:
        lengths._BLOCK = 8192
        totals = {fast: lengths.count_lengths([path], stop=lambda: True,
                                              fast=fast).total
                  for fast in SCANNERS}
    finally:
        lengths._BLOCK = original
    assert all(0 < t < 20000 for t in totals.values())


def test_not_asking_to_stop_reads_everything(tmp_path):
    path = _fastq(tmp_path / "a.fastq", [800] * 5000)
    assert lengths.count_lengths([path], stop=lambda: False).total == 5000
    assert lengths.count_lengths([path]).total == 5000


def test_a_stopped_scan_says_what_the_figures_cover():
    view, _ = _live()
    view.bytes, view.total, view.reads = 77 << 20, 463 << 20, 175_545
    view.stopped = True
    note = _strip(view.note())
    assert "stopped after" in note
    assert "175,545 reads counted" in note
    assert "463" in note                        # and of how much


def test_a_finished_scan_reports_its_rate_instead():
    view, _ = _live()
    view.bytes, view.total, view.took = 2 << 30, 2 << 30, 4.0
    note = _strip(view.note())
    assert note.startswith("scanned") and "/s" in note


# --- Watching for a keypress ---------------------------------------------

def test_key_watching_needs_a_terminal():
    from seqviewer.cli import KeyWatch

    class Piped:
        def isatty(self):
            return False

    watch = KeyWatch(stream=Piped())
    assert not watch.armed
    assert watch.pressed() is False
    watch.close()


def test_key_watching_can_be_turned_off():
    from seqviewer.cli import KeyWatch

    watch = KeyWatch(enabled=False)
    assert not watch.armed
    watch.close()


def test_a_stream_that_only_looks_like_a_terminal_does_not_arm():
    """isatty can lie; fileno then fails and the watch must stand down."""
    from seqviewer.cli import KeyWatch

    class Liar:
        def isatty(self):
            return True

        def fileno(self):
            raise OSError("not a real descriptor")

    watch = KeyWatch(stream=Liar())
    assert not watch.armed
    assert watch.pressed() is False
    watch.close()


def test_the_terminal_is_put_back_as_it_was(tmp_path):
    """Leaving cbreak mode set would break the shell the command returns to."""
    import pty

    from seqviewer.cli import KeyWatch, termios

    if termios is None:
        pytest.skip("no termios on this platform")
    leader, follower = pty.openpty()
    try:
        stream = open(follower, "rb", buffering=0)

        def line_mode():
            """Whether the terminal buffers lines and echoes them.

            Only these two flags are compared: the whole attribute list also
            carries status bits the kernel sets on its own, such as PENDIN.
            """
            flags = termios.tcgetattr(follower)[3]
            return bool(flags & termios.ICANON), bool(flags & termios.ECHO)

        assert line_mode() == (True, True)
        with KeyWatch(stream=stream) as watch:
            assert watch.armed
            assert line_mode() == (False, False)       # cbreak, so keys arrive
        assert line_mode() == (True, True)             # and put back after
        stream.close()
    finally:
        os.close(leader)


def test_a_keypress_is_seen(tmp_path):
    """What stops the scan: a byte waiting on the terminal."""
    import pty

    from seqviewer.cli import KeyWatch, termios

    if termios is None:
        pytest.skip("no termios on this platform")
    leader, follower = pty.openpty()
    try:
        stream = open(follower, "rb", buffering=0)
        with KeyWatch(stream=stream) as watch:
            assert watch.armed
            assert not watch.pressed()
            os.write(leader, b"q")
            time.sleep(0.05)
            assert watch.pressed()
            assert watch.pressed()              # and it stays pressed
        stream.close()
    finally:
        os.close(leader)


# --- Naming and opening the PNG ------------------------------------------

def test_a_run_stem_drops_the_fastq_suffix():
    from seqviewer.cli import _run_stem

    assert _run_stem("6YB866_1_sample_1.fastq.gz") == "6YB866_1_sample_1"
    assert _run_stem("/a/b/reads.fq") == "reads"
    assert _run_stem("/a/b/reads/") == "reads"
    assert _run_stem("plain_name") == "plain_name"


def test_a_png_is_named_after_the_run_and_put_in_downloads():
    from pathlib import Path

    from seqviewer.cli import png_path

    where = png_path("sample_1.fastq.gz")
    assert where.name == "sample_1-lengths.png"
    assert where.parent in (Path.home() / "Downloads", Path.cwd())


def test_a_given_png_path_is_used_as_it_stands(tmp_path):
    from seqviewer.cli import png_path

    given = tmp_path / "somewhere" / "mine.png"
    assert png_path("ignored.fastq", str(given)) == given


def test_opening_reports_when_there_is_nothing_to_open(monkeypatch):
    import seqviewer.cli as cli

    monkeypatch.setattr(cli.sys, "platform", "sunos5")
    assert cli.open_in_viewer("/tmp/whatever.png") is False


def test_opening_reports_when_the_opener_is_absent(monkeypatch):
    import seqviewer.cli as cli

    monkeypatch.setattr(cli.shutil, "which", lambda _name: None)
    assert cli.open_in_viewer("/tmp/whatever.png") is False


# --- Drawing the PNG -----------------------------------------------------

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def test_a_png_is_written(tmp_path):
    pytest.importorskip("matplotlib")
    from seqviewer import plot

    counts = lengths.LengthCounts.from_lengths(
        [780 + (i * 7) % 41 for i in range(500)] + [90, 9000])
    binning = lengths.bin_counts(counts)
    written = plot.write_png(tmp_path / "a.png", binning,
                             lengths.summarise_counts(counts), title="a run")
    assert written.read_bytes()[:8] == PNG_MAGIC
    assert written.stat().st_size > 5_000        # a real figure, not a stub


def test_the_png_directory_is_made_if_missing(tmp_path):
    pytest.importorskip("matplotlib")
    from seqviewer import plot

    counts = lengths.LengthCounts.from_lengths([800] * 100)
    target = tmp_path / "new" / "deeper" / "a.png"
    plot.write_png(target, lengths.bin_counts(counts),
                   lengths.summarise_counts(counts))
    assert target.exists()


def test_a_log_png_differs_from_a_linear_one(tmp_path):
    pytest.importorskip("matplotlib")
    from seqviewer import plot

    counts = lengths.LengthCounts.from_lengths([800] * 5000 + [400] * 3)
    binning = lengths.bin_counts(counts, bulk=100)
    summary = lengths.summarise_counts(counts)
    flat = plot.write_png(tmp_path / "flat.png", binning, summary)
    tall = plot.write_png(tmp_path / "log.png", binning, summary, log=True)
    assert flat.read_bytes() != tall.read_bytes()


def test_drawing_a_png_needs_no_import_at_module_level():
    """Importing seqviewer must not pull matplotlib in."""
    import ast
    import pathlib

    from seqviewer import plot

    tree = ast.parse(pathlib.Path(plot.__file__).read_text())
    top = set()
    for node in tree.body:                       # module level only
        if isinstance(node, ast.Import):
            top.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            top.add(node.module.split(".")[0])
    assert "matplotlib" not in top


# --- The --png flag ------------------------------------------------------

def test_png_writes_and_does_not_open_when_told_not_to(tmp_path, capsys):
    pytest.importorskip("matplotlib")
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [800] * 300 + [90, 9000])
    target = tmp_path / "out.png"
    assert lengths_main([str(tmp_path), "--width", "80", "--no-live",
                         "--png", str(target), "--no-open"]) == 0
    assert target.read_bytes()[:8] == PNG_MAGIC
    assert str(target) in capsys.readouterr().out


def test_png_is_not_drawn_without_the_flag(tmp_path, capsys):
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [800] * 50)
    lengths_main([str(tmp_path), "--width", "80"])
    assert "wrote" not in capsys.readouterr().out


def test_a_missing_matplotlib_is_reported_not_raised(tmp_path, capsys,
                                                    monkeypatch):
    from seqviewer import plot
    from seqviewer.cli import lengths_main

    def unavailable(*args, **kwargs):
        raise plot.PngUnavailable("drawing a PNG needs matplotlib")

    monkeypatch.setattr(plot, "write_png", unavailable)
    _fastq(tmp_path / "a.fastq", [800] * 50)
    code = lengths_main([str(tmp_path), "--no-live", "--png",
                         str(tmp_path / "x.png"), "--no-open"])
    assert code == 1
    assert "needs matplotlib" in capsys.readouterr().err


def test_a_path_that_cannot_be_written_is_reported(tmp_path, capsys):
    pytest.importorskip("matplotlib")
    from seqviewer.cli import lengths_main

    _fastq(tmp_path / "a.fastq", [800] * 50)
    blocker = tmp_path / "blocked"
    blocker.write_text("not a directory")
    code = lengths_main([str(tmp_path), "--no-live", "--png",
                         str(blocker / "x.png"), "--no-open"])
    assert code == 1
    assert "could not write" in capsys.readouterr().err


# --- Where the figure goes ------------------------------------------------

def test_a_directory_is_taken_as_the_folder_to_write_into(tmp_path):
    """--png ~/Downloads/ asks for a file in Downloads, not one called that.

    matplotlib appends the extension when a name has none, so passing the folder
    through would have written Downloads.png beside the folder.
    """
    from seqviewer.cli import png_path

    where = png_path("sample_1.fastq.gz", str(tmp_path))
    assert where.parent == tmp_path
    assert where.name == "sample_1-lengths.png"


def test_a_trailing_slash_means_a_folder_even_if_it_is_absent(tmp_path):
    from seqviewer.cli import png_path

    absent = tmp_path / "not_yet"
    where = png_path("a.fastq", f"{absent}/")
    assert where == absent / "a-lengths.png"


def test_a_name_without_a_suffix_gains_png(tmp_path):
    from seqviewer.cli import png_path

    assert png_path("a.fastq", str(tmp_path / "chart")).name == "chart.png"


def test_another_suffix_is_left_alone(tmp_path):
    """The suffix picks the format, so a PDF stays a PDF."""
    from seqviewer.cli import png_path

    for name in ("chart.pdf", "chart.svg", "chart.png"):
        assert png_path("a.fastq", str(tmp_path / name)).name == name


def test_writing_to_a_directory_is_refused(tmp_path):
    pytest.importorskip("matplotlib")
    from seqviewer import plot

    counts = lengths.LengthCounts.from_lengths([800] * 50)
    with pytest.raises(IsADirectoryError):
        plot.write_png(tmp_path, lengths.bin_counts(counts),
                       lengths.summarise_counts(counts))


def test_write_png_adds_the_suffix_itself(tmp_path):
    """A caller of the module, not just the command, gets a named file."""
    pytest.importorskip("matplotlib")
    from seqviewer import plot

    counts = lengths.LengthCounts.from_lengths([800] * 50)
    written = plot.write_png(tmp_path / "chart", lengths.bin_counts(counts),
                             lengths.summarise_counts(counts))
    assert written.name == "chart.png"
    assert written.read_bytes()[:8] == PNG_MAGIC


def test_the_command_writes_into_a_directory_given_one(tmp_path, capsys):
    pytest.importorskip("matplotlib")
    from seqviewer.cli import lengths_main

    reads = tmp_path / "reads"
    reads.mkdir()
    _fastq(reads / "sample_1.fastq", [800] * 300 + [90, 9000])
    out = tmp_path / "figures"
    out.mkdir()
    assert lengths_main([str(reads), "--no-live", "--no-open",
                         "--png", str(out)]) == 0
    written = list(out.glob("*.png"))
    assert len(written) == 1
    assert written[0].read_bytes()[:8] == PNG_MAGIC
    assert not (tmp_path / "figures.png").exists()   # not beside the folder


# --- How the figure is styled -------------------------------------------

def test_the_figure_is_drawn_at_three_hundred_dpi():
    from seqviewer import plot

    assert plot.DPI == 300


def test_the_panel_is_the_width_of_a_two_column_figure():
    from seqviewer import plot

    assert plot.FIGSIZE == (7.2, 3.6)               # 183 mm across


def test_arial_is_asked_for_first():
    from seqviewer.plot import _STYLE

    assert _STYLE["font.family"] == "sans-serif"
    assert _STYLE["font.sans-serif"][0] == "Arial"
    assert "DejaVu Sans" in _STYLE["font.sans-serif"]   # and a fallback


def test_the_panel_is_closed_on_all_four_sides():
    from seqviewer.plot import _STYLE

    assert _STYLE["axes.spines.top"] is True
    assert _STYLE["axes.spines.right"] is True
    assert "axes.grid" not in _STYLE                # gridlines stay off


def test_the_style_does_not_leak_into_the_callers_settings(tmp_path):
    """Applied through rc_context, so a caller keeps its own rcParams."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    from seqviewer import plot

    before = dict(plt.rcParams)
    counts = lengths.LengthCounts.from_lengths([800] * 50)
    plot.write_png(tmp_path / "a.png", lengths.bin_counts(counts),
                   lengths.summarise_counts(counts))
    assert plt.rcParams["font.size"] == before["font.size"]
    assert plt.rcParams["axes.spines.top"] == before["axes.spines.top"]


def test_the_figures_line_is_one_line_of_the_essentials():
    from seqviewer.plot import _figures

    counts = lengths.LengthCounts.from_lengths([800] * 1000 + [90, 9000])
    line = _figures(lengths.summarise_counts(counts))
    assert "\n" not in line
    assert "n = 1,002 reads" in line
    assert "median" in line and "N50" in line
    assert "bases" not in line                      # left off, as extraneous


def test_the_clipped_note_names_what_is_missing_and_the_true_range():
    from seqviewer.plot import _clipped_note

    counts = lengths.LengthCounts.from_lengths([800] * 1000 + [90, 9000])
    binning = lengths.bin_counts(counts)
    note = _clipped_note(binning, lengths.summarise_counts(counts))
    assert "not shown" in note
    assert "90–9,000 bp overall" in note            # the real extremes
    assert "%" in note


def test_there_is_no_clipped_note_when_nothing_is_clipped():
    from seqviewer.plot import _clipped_note

    counts = lengths.LengthCounts.from_lengths([800] * 100)
    assert _clipped_note(lengths.bin_counts(counts, bulk=100),
                         lengths.summarise_counts(counts)) is None
