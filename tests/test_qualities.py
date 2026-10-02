"""Per-read mean quality: the arithmetic, the two scanners, and the command."""
import gzip
import math
import random

import pytest

from seqviewer import lengths, qualities
from seqviewer.cli import build_qualities_parser, qualities_main, seqview_main


def _phred(*scores):
    """A quality line from Phred scores."""
    return bytes(q + 33 for q in scores)


def _fastq(path, quality_lines, crlf=False, final_newline=True):
    """Write one record per quality line, with a sequence the same length."""
    sep = "\r\n" if crlf else "\n"
    records = []
    for i, quality in enumerate(quality_lines):
        q = quality.decode()
        records.append(sep.join([f"@r{i}", "A" * len(q), "+", q]))
    text = sep.join(records) + (sep if final_newline else "")
    path.write_bytes(text.encode())
    return path


SCANNERS = [False] + ([True] if qualities.HAVE_NUMPY else [])


# --- the arithmetic -------------------------------------------------------

def test_a_uniform_read_scores_its_one_score():
    assert qualities.mean_quality(_phred(*[20] * 50)) == pytest.approx(20.0)


def test_the_mean_is_of_error_rates_not_of_scores():
    """One base at Q3 among ninety-nine at Q40, the case docs/index.html works.

    Its error rate is about one in two hundred, which is Q23; the mean of its
    scores is Q39.6.  The second figure is the one a naive average reports.
    """
    line = _phred(3, *[40] * 99)
    proper = qualities.mean_quality(line)
    naive = (3 + 40 * 99) / 100
    assert proper == pytest.approx(22.9, abs=0.05)
    assert naive == pytest.approx(39.63)
    expected_rate = (10 ** -0.3 + 99 * 1e-4) / 100
    assert proper == pytest.approx(-10 * math.log10(expected_rate))


def test_an_empty_quality_line_has_no_mean():
    assert qualities.mean_quality(b"") is None


def test_keys_are_floored_so_a_bar_holds_the_reads_inside_it():
    """Bar 19 is Q19.0 up to Q20.0; a read at Q19.96 is in it, not in bar 20."""
    assert qualities._key(19.96) == 199
    assert qualities._key(20.0) == 200
    assert qualities._key(19.0) == 190


def test_float_noise_on_an_exact_tenth_stays_on_that_tenth():
    """A uniform read lands on a tenth; both float paths must find the same."""
    assert qualities._key(20.0 - 1e-12) == 200
    assert qualities._key(20.0 + 1e-12) == 200


def test_scores_off_the_scale_are_clamped_not_raised():
    assert qualities._key(-5.0) == 0
    assert qualities._key(500.0) == qualities.MAX_Q * qualities.SCALE


# --- the scanners -----------------------------------------------------------

def _random_lines(n, seed=1):
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        length = rng.randint(1, 400)
        centre = rng.uniform(8, 38)
        out.append(_phred(*(max(0, min(60, int(rng.gauss(centre, 6))))
                            for _ in range(length))))
    return out


@pytest.mark.parametrize("fast", SCANNERS)
def test_a_scanner_tallies_each_read_once(tmp_path, fast):
    path = _fastq(tmp_path / "r.fq", _random_lines(300))
    tally = qualities.count_qualities([path], fast=fast)
    assert tally.counts.total == 300
    assert tally.unscored == 0


@pytest.mark.skipif(not qualities.HAVE_NUMPY, reason="needs numpy")
def test_the_scanners_agree(tmp_path):
    path = _fastq(tmp_path / "r.fq", _random_lines(500))
    slow = qualities.count_qualities([path], fast=False)
    fast = qualities.count_qualities([path], fast=True)
    assert slow.counts.items() == fast.counts.items()


@pytest.mark.skipif(not qualities.HAVE_NUMPY, reason="needs numpy")
@pytest.mark.parametrize("block", [512, 1031, 4096])
def test_the_scanners_agree_across_block_boundaries(tmp_path, block):
    """Records split between blocks are carried and scored once."""
    path = _fastq(tmp_path / "r.fq", _random_lines(200, seed=4))
    original = lengths._BLOCK
    try:
        lengths._BLOCK = block
        slow = qualities.count_qualities([path], fast=False)
        fast = qualities.count_qualities([path], fast=True)
    finally:
        lengths._BLOCK = original
    whole = qualities.count_qualities([path], fast=False)
    assert slow.counts.items() == fast.counts.items() == whole.counts.items()


@pytest.mark.parametrize("fast", SCANNERS)
def test_crlf_line_endings_do_not_score_the_cr(tmp_path, fast):
    """A CR is byte 13, below the offset, and would score as Q0 if counted."""
    lines = [_phred(*[30] * 20)] * 5
    unix = qualities.count_qualities([_fastq(tmp_path / "u.fq", lines)],
                                     fast=fast)
    dos = qualities.count_qualities(
        [_fastq(tmp_path / "d.fq", lines, crlf=True)], fast=fast)
    assert unix.counts.items() == dos.counts.items() == [(300, 5)]


@pytest.mark.parametrize("fast", SCANNERS)
def test_a_final_record_without_a_newline_is_scored(tmp_path, fast):
    lines = [_phred(*[25] * 10), _phred(*[35] * 10)]
    path = _fastq(tmp_path / "r.fq", lines, final_newline=False)
    tally = qualities.count_qualities([path], fast=fast)
    assert tally.counts.items() == [(250, 1), (350, 1)]


@pytest.mark.parametrize("fast", SCANNERS)
def test_an_empty_quality_line_is_counted_apart(tmp_path, fast):
    path = _fastq(tmp_path / "r.fq",
                  [_phred(*[20] * 10), b"", _phred(*[30] * 10), b""])
    tally = qualities.count_qualities([path], fast=fast)
    assert tally.counts.total == 2
    assert tally.unscored == 2


@pytest.mark.parametrize("fast", SCANNERS)
def test_a_gzipped_file_is_read_directly(tmp_path, fast):
    plain = _fastq(tmp_path / "r.fq", _random_lines(50, seed=9))
    packed = tmp_path / "r.fq.gz"
    with gzip.open(packed, "wb") as handle:
        handle.write(plain.read_bytes())
    a = qualities.count_qualities([plain], fast=fast)
    b = qualities.count_qualities([packed], fast=fast)
    assert a.counts.items() == b.counts.items()


def test_requiring_the_array_scanner_without_numpy_is_an_error(monkeypatch,
                                                              tmp_path):
    monkeypatch.setattr(qualities, "HAVE_NUMPY", False)
    path = _fastq(tmp_path / "r.fq", [_phred(20)])
    with pytest.raises(RuntimeError):
        qualities.count_qualities([path], fast=True)


# --- the figures ------------------------------------------------------------

def _tally(*scores_per_read):
    counts = lengths.LengthCounts()
    for q in scores_per_read:
        counts.add(qualities._key(q))
    return qualities.QualityTally(counts)


def test_a_threshold_counts_the_reads_at_or_above_it():
    tally = _tally(9.9, 10.0, 15.0, 19.96, 20.0, 31.0)
    summary = qualities.summarise_qualities(tally, (10, 20, 30))
    assert dict(summary.at_least) == {10: 5, 20: 2, 30: 1}


def test_the_figures_are_exact_to_a_tenth():
    summary = qualities.summarise_qualities(_tally(12.84, 24.25, 43.79))
    assert (summary.lowest, summary.median, summary.highest) == (12.8, 24.2,
                                                                 43.7)


def test_the_summary_names_reads_left_out():
    tally = _tally(20.0)
    tally.unscored = 3
    text = "\n".join(qualities.summary_lines(
        qualities.summarise_qualities(tally)))
    assert "3 records with an empty quality line are left out" in text


def test_a_run_with_no_scores_says_why():
    tally = qualities.QualityTally(lengths.LengthCounts(), unscored=4)
    lines = qualities.summary_lines(qualities.summarise_qualities(tally))
    assert lines == ["no reads with qualities; 4 records had an empty "
                     "quality line"]


# --- the drawing ------------------------------------------------------------

def test_bars_are_whole_q_and_the_header_names_the_unit():
    tally = _tally(*[12.3, 12.9, 13.1, 14.0, 14.5, 15.2])
    hist, _ = qualities.render(tally, bins=32, width=60)
    labels = [line.split()[0] for line in hist[1:]]
    assert hist[0].split()[0] == "Q"
    assert labels == ["12", "13", "14", "15"]


def test_fewer_bins_than_scores_merge_neighbouring_q():
    tally = _tally(*[q + 0.5 for q in range(10, 30)])
    hist, _ = qualities.render(tally, bins=10, width=60)
    assert hist[1].split()[0] == "10–11"


def test_nothing_is_clipped_since_scores_are_bounded():
    binning = qualities.bin_qualities(_tally(5.0, 20.0, 60.0))
    assert not binning.clipped


# --- the command -------------------------------------------------------------

def test_the_command_prints_the_histogram_and_figures(tmp_path, capsys):
    path = _fastq(tmp_path / "r.fq",
                  [_phred(*[20] * 30)] * 4 + [_phred(*[32] * 30)] * 2)
    assert qualities_main([str(path), "--no-live", "--width", "60"]) == 0
    out = capsys.readouterr().out
    assert "6 reads" in out
    assert "≥Q20 6 (100.0%)" in out
    assert "≥Q30 2 (33.3%)" in out


def test_thresholds_are_parsed_sorted_and_deduplicated():
    args = build_qualities_parser().parse_args(
        ["r.fq", "--thresholds", "30,15,15,20"])
    assert args.thresholds == (15, 20, 30)


@pytest.mark.parametrize("bad", ["ten", "10,x", "-1", "200"])
def test_bad_thresholds_are_refused(bad):
    with pytest.raises(SystemExit):
        build_qualities_parser().parse_args(["r.fq", "--thresholds", bad])


def test_a_missing_file_is_an_error(tmp_path, capsys):
    assert qualities_main([str(tmp_path / "absent.fq")]) == 1
    assert "no FASTQ files" in capsys.readouterr().err


def test_seqview_lists_and_dispatches_it(tmp_path, capsys):
    assert seqview_main(["--help"]) == 0
    assert "qualities" in capsys.readouterr().out
    path = _fastq(tmp_path / "r.fq", [_phred(*[25] * 10)])
    assert seqview_main(["qualities", str(path), "--no-live"]) == 0
