"""Base-quality distributions and their terminal histogram.

Pure arithmetic and text, so none of this needs a terminal or a real run.
"""

import gzip
import io
import re
from collections import Counter

import pytest

from seqviewer import quality
from seqviewer.cli import quality_main, seqview_main
from seqviewer.quality import (
    QualityCounts, bin_counts, count_qualities, histogram, summarise,
    summary_lines,
)

needs_numpy = pytest.mark.skipif(not quality.HAVE_NUMPY,
                                 reason="the array scanner needs numpy")

_ANSI = re.compile(r"\x1b\[[0-9;]*[a-zA-Z]")


def _scores(text):
    """The Phred+33 scores a quality string stands for."""
    return [ord(c) - 33 for c in text]


def _fastq(path, quals, sep="\n"):
    """Write a FASTQ whose reads have the quality strings in *quals*."""
    with open(path, "w", newline="") as handle:
        for i, q in enumerate(quals):
            handle.write(f"@r{i}{sep}{'A' * len(q)}{sep}+{sep}{q}{sep}")
    return path


def _counts(path, **kwargs):
    return dict(count_qualities([path], **kwargs).items())


QUALS = ["IIII", "I5#!", "?"]
WANT = Counter(s for q in QUALS for s in _scores(q))


# --- Reading --------------------------------------------------------------

def test_every_base_is_counted_at_its_score(tmp_path):
    path = _fastq(tmp_path / "a.fastq", QUALS)
    assert _counts(path) == dict(WANT)


def test_reads_are_counted_apart_from_bases(tmp_path):
    path = _fastq(tmp_path / "a.fastq", QUALS)
    counts = count_qualities([path])
    assert counts.reads == 3
    assert counts.total == 9


def test_only_the_fourth_line_is_read(tmp_path):
    """A sequence made of quality characters must not be counted as scores."""
    path = tmp_path / "a.fastq"
    path.write_text("@IIII\nIIII\n+IIII\n5555\n")
    assert _counts(path) == {20: 4}


def test_several_files_are_read_as_one_run(tmp_path):
    a = _fastq(tmp_path / "a.fastq", ["II"])
    b = _fastq(tmp_path / "b.fastq", ["55"])
    counts = count_qualities([a, b])
    assert dict(counts.items()) == {40: 2, 20: 2}
    assert counts.reads == 2


def test_a_gzipped_file_is_read_directly(tmp_path):
    path = tmp_path / "a.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r0\nAAAA\n+\nII55\n")
    assert _counts(path) == {40: 2, 20: 2}


def test_a_truncated_final_record_is_dropped(tmp_path):
    """The state a file interrupted mid-write is left in."""
    path = tmp_path / "a.fastq"
    path.write_text("@r0\nAAAA\n+\nIIII\n@r1\nAAAAAA\n+\n")
    counts = count_qualities([path])
    assert counts.reads == 1
    assert dict(counts.items()) == {40: 4}


def test_a_last_line_without_a_newline_is_counted(tmp_path):
    path = tmp_path / "a.fastq"
    path.write_text("@r0\nAAAA\n+\nIIII")
    assert _counts(path) == {40: 4}


def test_an_empty_file_has_no_bases(tmp_path):
    path = tmp_path / "a.fastq"
    path.write_text("")
    assert count_qualities([path]).empty


@pytest.mark.parametrize("fast", [False, pytest.param(True, marks=needs_numpy)])
def test_windows_line_endings_leave_no_stray_scores(tmp_path, fast):
    path = _fastq(tmp_path / "a.fastq", QUALS, sep="\r\n")
    counts = count_qualities([path], fast=fast)
    assert dict(counts.items()) == dict(WANT)
    assert counts.skipped == 0


def test_the_offset_is_taken_off_each_byte(tmp_path):
    path = tmp_path / "a.fastq"
    path.write_text("@r\nAA\n+\nhh\n")           # 'h' is 104, so Q40 at 64
    assert _counts(path, offset=64) == {40: 2}


def test_bytes_below_the_offset_are_counted_apart(tmp_path):
    """Usually the wrong offset, so it is reported rather than dropped."""
    path = tmp_path / "a.fastq"
    path.write_text("@r\nAAAA\n+\n5555\n")        # '5' is 53, under 64
    counts = count_qualities([path], offset=64)
    assert counts.empty
    assert counts.skipped == 4


# --- The two scanners -----------------------------------------------------

def _stream(quals, sep="\n"):
    text = "".join(f"@r{i}{sep}{'A' * len(q)}{sep}+{sep}{q}{sep}"
                   for i, q in enumerate(quals))
    return io.BytesIO(text.encode())


def _random_quals(n=300, seed=3):
    import random
    rng = random.Random(seed)
    return ["".join(chr(33 + rng.randint(0, 41))
                    for _ in range(rng.randint(1, 120))) for _ in range(n)]


@needs_numpy
@pytest.mark.parametrize("sep", ["\n", "\r\n"])
def test_the_scanners_return_the_same_tally(sep):
    quals = _random_quals()
    slow, slow_reads = quality._tally_python(_stream(quals, sep))
    fast, fast_reads = quality._tally_numpy(_stream(quals, sep))
    a = QualityCounts.from_tally(slow, slow_reads)
    b = QualityCounts.from_tally(fast, fast_reads)
    assert a.items() == b.items()
    assert slow_reads == fast_reads == len(quals)


@pytest.mark.parametrize("fast", [False, pytest.param(True, marks=needs_numpy)])
def test_records_split_across_blocks_are_counted_once(monkeypatch, fast):
    """A block boundary can fall anywhere, including inside a record."""
    quals = _random_quals(n=200)
    monkeypatch.setattr(quality, "_BLOCK", 97)
    scan = quality._tally_numpy if fast else quality._tally_python
    tally, reads = scan(_stream(quals))
    counts = QualityCounts.from_tally(tally, reads)
    want = Counter(s for q in quals for s in _scores(q))
    assert dict(counts.items()) == dict(want)
    assert reads == len(quals)


@pytest.mark.parametrize("fast", [False, pytest.param(True, marks=needs_numpy)])
def test_progress_reports_every_block_and_ends_at_the_full_size(tmp_path, fast):
    path = _fastq(tmp_path / "a.fastq", _random_quals(n=50))
    seen = []
    count_qualities([path], progress=lambda *a: seen.append(a), fast=fast)
    done, total, reads = seen[-1]
    assert done == total == path.stat().st_size
    assert reads == 50


def test_the_array_scanner_is_refused_without_numpy(monkeypatch, tmp_path):
    monkeypatch.setattr(quality, "HAVE_NUMPY", False)
    path = _fastq(tmp_path / "a.fastq", ["I"])
    with pytest.raises(RuntimeError):
        count_qualities([path], fast=True)


# --- Figures --------------------------------------------------------------

def _tally(pairs, reads=1):
    return QualityCounts(dict(pairs), reads=reads)


def test_the_mean_is_weighted_by_count():
    assert _tally({10: 1, 40: 3}).mean == pytest.approx(32.5)


def test_the_median_is_the_score_half_the_bases_reach():
    assert _tally({10: 1, 20: 1, 30: 1}).median() == 20
    assert _tally({10: 2, 30: 2}).median() == 10     # half are at or below 10


def test_a_share_counts_the_threshold_itself():
    counts = _tally({19: 1, 20: 1, 30: 2})
    assert counts.share_at_least(20) == pytest.approx(0.75)
    assert counts.share_at_least(30) == pytest.approx(0.5)


def test_the_summary_reports_both_thresholds():
    summary = summarise(_tally({15: 1, 25: 1, 35: 2}, reads=2))
    assert summary.bases == 4
    assert summary.lowest == 15 and summary.highest == 35
    assert dict(summary.shares) == {20: pytest.approx(0.75),
                                    30: pytest.approx(0.5)}


def test_an_empty_tally_summarises_as_empty():
    assert summarise(QualityCounts()).empty
    assert summary_lines(summarise(QualityCounts())) == ["no bases"]


def test_the_summary_says_when_bases_fell_below_the_offset():
    lines = summary_lines(summarise(QualityCounts({30: 5}, reads=1,
                                                  skipped=2)))
    assert any("below the offset" in line for line in lines)


def test_a_clean_summary_says_nothing_about_the_offset():
    lines = summary_lines(summarise(_tally({30: 5})))
    assert not any("offset" in line for line in lines)


# --- Binning --------------------------------------------------------------

def test_every_base_lands_in_exactly_one_bin():
    counts = _tally({q: q + 1 for q in range(0, 42)})
    for bins in (1, 5, 14, 21, 60):
        binning = bin_counts(counts, bins)
        assert sum(b.count for b in binning.bins) == counts.total


def test_bins_are_whole_scores_wide_and_aligned_to_the_width():
    binning = bin_counts(_tally({41: 1, 5: 1}), 14)
    assert binning.width == 3
    assert all(b.high - b.low == 3 for b in binning.bins)
    assert all(b.low % 3 == 0 for b in binning.bins)


def test_the_axis_opens_on_the_bin_holding_the_lowest_score():
    binning = bin_counts(_tally({27: 5, 40: 5}), 14)
    first = binning.bins[0]
    assert first.low <= 27 < first.high
    assert first.count == 5                      # no empty rows before it


def test_a_run_reaching_score_zero_starts_there():
    assert bin_counts(_tally({0: 1, 40: 1}), 14).bins[0].low == 0


def test_the_axis_reaches_the_highest_score_seen():
    binning = bin_counts(_tally({38: 1}), 21)
    assert binning.bins[-1].low <= 38 < binning.bins[-1].high


def test_more_bins_than_scores_gives_one_score_a_bin():
    binning = bin_counts(_tally({0: 1, 4: 1}), 100)
    assert binning.width == 1
    assert len(binning.bins) == 5


def test_nothing_to_bin_is_empty():
    assert bin_counts(QualityCounts()).empty


# --- Histogram ------------------------------------------------------------

def test_no_line_is_wider_than_asked():
    binning = bin_counts(_tally({q: 1000 + q for q in range(42)}), 14)
    for width in (40, 60, 80):
        assert all(len(line) <= width for line in histogram(binning, width))


def test_the_tallest_bin_gets_the_longest_bar():
    binning = bin_counts(_tally({10: 1, 30: 100}), 42)
    rows = histogram(binning, 60)[1:]
    bars = [len(re.findall("[█▏▎▍▌▋▊▉]", r)) for r in rows]
    assert bars.index(max(bars)) == 20       # the axis opens at Q10


def test_a_bin_with_bases_is_never_blank():
    binning = bin_counts(_tally({10: 1, 30: 10_000_000}), 42)
    row = histogram(binning, 60)[1]            # the axis opens at Q10
    assert re.search("[▏▎▍▌▋▊▉█]", row)


def test_an_empty_bin_has_no_bar():
    binning = bin_counts(_tally({0: 5, 30: 5}), 42)
    assert not re.search("[▏▎▍▌▋▊▉█]", histogram(binning, 60)[1 + 15])


def test_log_scale_lifts_the_small_bins():
    binning = bin_counts(_tally({10: 10, 30: 100_000}), 42)

    def bar(lines, i):
        return len(re.findall("[█▏▎▍▌▋▊▉]", lines[1 + i]))

    assert bar(histogram(binning, 60, log=True), 0) > bar(
        histogram(binning, 60), 0)             # the axis opens at Q10


def test_labels_are_exact_ranges():
    lines = histogram(bin_counts(_tally({0: 1, 8: 1}), 3), 60)
    assert "Q0–2" in lines[1]
    assert "Q6–8" in lines[3]
    wide = histogram(bin_counts(_tally({0: 1, 2: 1}), 10), 60)
    assert "Q1 " in wide[2]                      # one score a bin is a bare Q


def test_a_palette_colours_the_bars_without_moving_them():
    binning = bin_counts(_tally({10: 1, 30: 9}), 14)
    plain = histogram(binning, 60)
    coloured = histogram(binning, 60, palette=quality.PALETTE)
    assert coloured != plain
    assert [_ANSI.sub("", line) for line in coloured] == plain


def test_an_empty_binning_says_so():
    assert histogram(bin_counts(QualityCounts()), 60) == ["no bases"]


# --- The command ----------------------------------------------------------

def test_the_command_prints_the_histogram_and_the_figures(tmp_path, capsys):
    path = _fastq(tmp_path / "a.fastq", ["IIII", "5555"])
    assert quality_main([str(path), "--no-progress"]) == 0
    out = capsys.readouterr().out
    assert "a.fastq" in out
    assert "bases" in out
    assert "2 reads · 8 bases" in out
    assert "50.0% of bases at Q30 or above" in out


def test_a_directory_is_read_as_one_run(tmp_path, capsys):
    _fastq(tmp_path / "a.fastq", ["IIII"])
    _fastq(tmp_path / "b.fastq", ["IIII"])
    assert quality_main([str(tmp_path), "--no-progress"]) == 0
    out = capsys.readouterr().out
    assert "2 files" in out
    assert "2 reads · 8 bases" in out


def test_per_file_draws_one_histogram_each(tmp_path, capsys):
    _fastq(tmp_path / "a.fastq", ["IIII"])
    _fastq(tmp_path / "b.fastq", ["5555"])
    assert quality_main([str(tmp_path), "--per-file", "--no-progress"]) == 0
    out = capsys.readouterr().out
    assert out.count("1 file") == 2
    assert out.count("1 reads · 4 bases") == 2


def test_output_to_a_pipe_has_no_colour(tmp_path, capsys):
    path = _fastq(tmp_path / "a.fastq", ["IIII"])
    quality_main([str(path), "--no-progress"])
    assert "\x1b" not in capsys.readouterr().out


def test_the_slow_flag_reports_the_same_figures(tmp_path, capsys):
    path = _fastq(tmp_path / "a.fastq", _random_quals(n=40))
    quality_main([str(path), "--no-progress"])
    fast = capsys.readouterr().out
    quality_main([str(path), "--no-progress", "--slow"])
    assert capsys.readouterr().out == fast


def test_a_missing_path_is_an_error(tmp_path, capsys):
    assert quality_main([str(tmp_path / "none.fastq")]) == 1
    assert "no FASTQ files" in capsys.readouterr().err


def test_seqview_dispatches_to_quality(tmp_path, capsys):
    path = _fastq(tmp_path / "a.fastq", ["IIII"])
    assert seqview_main(["quality", str(path), "--no-progress"]) == 0
    assert "4 bases" in capsys.readouterr().out


def test_seqview_lists_quality_among_its_commands(capsys):
    seqview_main([])
    assert "quality" in capsys.readouterr().err


# --- Scoring each read by its mean ----------------------------------------

def _by_read(path, **kwargs):
    counts = count_qualities([path], by="read", **kwargs)
    return dict(counts.items()), counts


def test_each_read_is_one_count_at_its_mean(tmp_path):
    path = _fastq(tmp_path / "a.fastq", ["IIII", "5555", "I5"])
    tally, counts = _by_read(path, mean="phred")
    assert tally == {40: 1, 20: 1, 30: 1}
    assert counts.reads == 3
    assert counts.total == 3                     # reads, not bases
    assert counts.bases == 10
    assert counts.unit == "reads"


def test_a_mean_on_a_half_rounds_up(tmp_path):
    path = _fastq(tmp_path / "a.fastq", ["I4", "I5", "55I"])
    tally, _ = _by_read(path, mean="phred")
    # (40 + 19) / 2 = 29.5 -> 30;  (40 + 20) / 2 = 30;  (20 + 20 + 40) / 3 = 26.67 -> 27
    assert tally == {30: 2, 27: 1}


def test_a_read_with_no_bases_has_no_mean(tmp_path):
    path = tmp_path / "a.fastq"
    path.write_text("@r0\n\n+\n\n@r1\nAA\n+\nII\n")
    tally, counts = _by_read(path)
    assert tally == {40: 1}
    assert counts.skipped == 0


def test_a_read_holding_a_byte_below_the_offset_is_set_aside(tmp_path):
    """Its mean would mean nothing, and it is reported rather than dropped."""
    path = tmp_path / "a.fastq"
    path.write_text("@r0\nAAA\n+\nII\x1f\n@r1\nAA\n+\nII\n")
    tally, counts = _by_read(path)
    assert tally == {40: 1}
    assert counts.skipped == 1
    assert counts.reads == 1


@pytest.mark.parametrize("sep", ["\n", "\r\n"])
@pytest.mark.parametrize("fast", [False, pytest.param(True, marks=needs_numpy)])
def test_both_scanners_find_the_same_means(tmp_path, sep, fast):
    quals = _random_quals(n=250)
    path = _fastq(tmp_path / "a.fastq", quals, sep=sep)
    tally, counts = _by_read(path, fast=fast, mean="phred")
    want = Counter()
    for q in quals:
        scores = _scores(q)
        want[(2 * sum(scores) + len(scores)) // (2 * len(scores))] += 1
    assert tally == dict(want)
    assert counts.reads == len(quals)
    assert counts.bases == sum(len(q) for q in quals)


@pytest.mark.parametrize("fast", [False, pytest.param(True, marks=needs_numpy)])
def test_reads_split_across_blocks_are_scored_once(monkeypatch, tmp_path, fast):
    quals = _random_quals(n=150)
    path = _fastq(tmp_path / "a.fastq", quals)
    monkeypatch.setattr(quality, "_BLOCK", 101)
    tally, counts = _by_read(path, fast=fast)
    assert counts.reads == len(quals)
    assert sum(tally.values()) == len(quals)


def test_a_gzipped_file_is_scored_by_read(tmp_path):
    path = tmp_path / "a.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r0\nAAAA\n+\nII55\n")
    tally, _ = _by_read(path, mean="phred")
    assert tally == {30: 1}


def test_several_files_are_pooled_by_read(tmp_path):
    a = _fastq(tmp_path / "a.fastq", ["II"])
    b = _fastq(tmp_path / "b.fastq", ["55", "55"])
    counts = count_qualities([a, b], by="read", mean="phred")
    assert dict(counts.items()) == {40: 1, 20: 2}
    assert counts.reads == 3 and counts.bases == 6


def test_the_offset_applies_to_means_too(tmp_path):
    path = tmp_path / "a.fastq"
    path.write_text("@r\nAA\n+\nhh\n")           # 'h' is 104: Q40 at offset 64
    tally, _ = _by_read(path, offset=64)
    assert tally == {40: 1}


def test_an_unknown_way_of_scoring_is_refused(tmp_path):
    path = _fastq(tmp_path / "a.fastq", ["I"])
    with pytest.raises(ValueError):
        count_qualities([path], by="run")


def test_counts_of_bases_and_of_reads_do_not_mix():
    with pytest.raises(ValueError):
        QualityCounts(by="base").merge(QualityCounts(by="read"))


def test_a_read_histogram_is_headed_reads():
    counts = QualityCounts.from_read_tally({30: 2, 20: 1}, 3, 12, 0)
    lines = histogram(bin_counts(counts, 14), 60)
    assert lines[0].split()[-1] == "reads"
    assert bin_counts(counts).unit == "reads"


def test_the_summary_names_the_unit_and_how_it_was_scored():
    counts = QualityCounts.from_read_tally({30: 2, 20: 1}, 3, 12, 0)
    lines = summary_lines(summarise(counts))
    assert "3 reads · 12 bases · scored by each read's mean error rate" in lines[0]
    assert "of reads at Q20 or above" in lines[2]


def test_skipped_reads_are_explained_as_reads():
    counts = QualityCounts.from_read_tally({30: 2}, 2, 8, 5)
    text = " ".join(summary_lines(summarise(counts)))
    assert "5 reads held a base scored below the offset" in text


def test_the_command_scores_by_read_on_request(tmp_path, capsys):
    path = _fastq(tmp_path / "a.fastq", ["IIII", "5555"])
    assert quality_main([str(path), "--no-progress", "--by", "read"]) == 0
    out = capsys.readouterr().out
    assert "scored by each read's mean error rate" in out
    assert "of reads at Q30 or above" in out
    assert out.splitlines()[2].split()[-1] == "reads"


def test_the_command_scores_by_base_unless_told_otherwise(tmp_path, capsys):
    path = _fastq(tmp_path / "a.fastq", ["IIII", "5555"])
    quality_main([str(path), "--no-progress"])
    assert "mean error rate" not in capsys.readouterr().out


def test_the_command_refuses_an_unknown_way_of_scoring(tmp_path, capsys):
    path = _fastq(tmp_path / "a.fastq", ["I"])
    with pytest.raises(SystemExit):
        quality_main([str(path), "--by", "run"])


# --- How a read's mean is taken -------------------------------------------

def _expected_error_score(text):
    import math
    errors = [10 ** (-q / 10) for q in _scores(text)]
    return int(math.floor(-10 * math.log10(sum(errors) / len(errors)) + 0.5))


def test_a_reads_mean_is_taken_over_error_rates_by_default(tmp_path):
    """Q40 and Q20 average to Q23, not Q30: the poor bases dominate."""
    path = _fastq(tmp_path / "a.fastq", ["I5"])
    tally, _ = _by_read(path)
    assert tally == {23: 1}
    assert _by_read(path, mean="phred")[0] == {30: 1}


def test_a_uniform_read_scores_the_same_either_way(tmp_path):
    path = _fastq(tmp_path / "a.fastq", ["IIII", "5555"])
    assert _by_read(path)[0] == _by_read(path, mean="phred")[0] == {40: 1, 20: 1}


def test_the_error_mean_is_never_above_the_plain_mean(tmp_path):
    quals = _random_quals(n=120)
    path = _fastq(tmp_path / "a.fastq", quals)
    error, _ = _by_read(path)
    plain, _ = _by_read(path, mean="phred")
    assert sum(q * n for q, n in error.items()) < sum(
        q * n for q, n in plain.items())


@pytest.mark.parametrize("sep", ["\n", "\r\n"])
@pytest.mark.parametrize("fast", [False, pytest.param(True, marks=needs_numpy)])
def test_both_scanners_find_the_same_error_means(tmp_path, sep, fast):
    quals = _random_quals(n=250)
    path = _fastq(tmp_path / "a.fastq", quals, sep=sep)
    tally, counts = _by_read(path, fast=fast)
    want = Counter(_expected_error_score(q) for q in quals)
    assert tally == dict(want)
    assert counts.reads == len(quals)


@pytest.mark.parametrize("fast", [False, pytest.param(True, marks=needs_numpy)])
def test_error_means_survive_block_boundaries(monkeypatch, tmp_path, fast):
    quals = _random_quals(n=150)
    path = _fastq(tmp_path / "a.fastq", quals)
    monkeypatch.setattr(quality, "_BLOCK", 101)
    tally, _ = _by_read(path, fast=fast)
    assert tally == dict(Counter(_expected_error_score(q) for q in quals))


def test_the_summary_says_how_the_mean_was_taken(tmp_path):
    path = _fastq(tmp_path / "a.fastq", ["I5"])
    error = summary_lines(summarise(count_qualities([path], by="read")))
    plain = summary_lines(summarise(count_qualities(
        [path], by="read", mean="phred")))
    assert "mean error rate" in error[0]
    assert "mean score" in plain[0]


def test_an_unknown_mean_is_refused(tmp_path):
    path = _fastq(tmp_path / "a.fastq", ["I"])
    with pytest.raises(ValueError):
        count_qualities([path], by="read", mean="median")


def test_reads_scored_differently_do_not_mix():
    with pytest.raises(ValueError):
        QualityCounts(by="read", read_mean="error").merge(
            QualityCounts(by="read", read_mean="phred"))


def test_the_command_takes_the_mean_as_asked(tmp_path, capsys):
    path = _fastq(tmp_path / "a.fastq", ["I5"])
    quality_main([str(path), "--no-progress", "--by", "read"])
    assert "Q23" in capsys.readouterr().out
    quality_main([str(path), "--no-progress", "--by", "read", "--mean", "phred"])
    assert "Q30" in capsys.readouterr().out
