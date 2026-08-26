#!/usr/bin/env python
"""Align reads to a reference and write a seqviewer pileup page.

    seqviewer-pileup reads/ reference out.html [options]

Reads are a directory of FASTQs or a single file.  A sequencing run arrives as a
directory of per-barcode files, so a directory is the expected input; its files
are pooled into one pileup unless ``--per-file`` asks for a group each.

A page draws a few hundred reads by default, sampled from across the whole of
every file.  A deep run drawn whole is a page too large to open, and the reads
at the front of one file are not the run; ``--max 0`` draws all of them anyway.

Reading the reference is :mod:`seqviewer.genbank`'s job — FASTA, GenBank, ApE,
and SnapGene all arrive as a ``Reference`` with its topology and features.  What
is left here is the part the package does not cover: reading a FASTQ, and the
two adjustments a plasmid needs.  A circular reference is aligned against a
doubled copy so reads crossing the origin stay whole, and a named feature can
supply the flanks the page marks.

``--summary`` writes a second page beside the pileup, reduced from the same
view: the construct as an annotated map, one band per group with a lollipop per
called variant, and the variants as a table.  It answers "is anything wrong with
this clone" where the pileup answers "what does every read say".

seqviewer takes reads as ``Read(name, seq, qual)`` records, so FASTQ parsing is
the caller's job.  Everything after that is grid_from_reads plus render.
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import os
import random
import shutil
import sys
import textwrap
import time
from pathlib import Path

from . import lengths
from .align import Read, grid_from_reads
from .cluster import cluster_rows
from .genbank import load_reference
from .pileup import PileupGroup, PileupView
from .render import render
from .render_summary import render_summary
from .summary import DEFAULT_MIN_COUNT, DEFAULT_MIN_FRACTION, SummaryView

FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")

# A page is drawn to be read, and a deep run drawn whole is neither readable nor
# loadable: a row costs roughly 2 KB of HTML, so a 35,000-read pool renders to
# some 70 MB.  A few hundred reads is what a pileup can actually show, and is
# enough to see the subpopulations, so that is the default and the whole pile is
# what has to be asked for.
DEFAULT_MAX_READS = 500

# Fixed, so the same directory downsamples to the same page twice.
SAMPLE_SEED = 0


def summary_path(out):
    """Where the summarized page goes, given the pileup's *out* path.

    Beside the pileup rather than replacing it, and named from its stem, so a
    directory of runs sorts each pair together: ``1A12.html`` next to
    ``1A12.summary.html``.
    """
    out = Path(out)
    return out.with_name(out.stem + ".summary.html")


def fastq_paths(path):
    """The FASTQs at *path*: a directory's files, or the one file named.

    Naming a single file is the same code path with a one-element list.  Hidden
    files are skipped, which is what keeps the ``._`` stubs a cloud-synced share
    leaves behind from being read as reads.
    """
    path = Path(path)
    if not path.is_dir():
        return [path] if path.exists() else []
    return sorted(p for p in path.iterdir()
                  if p.is_file()
                  and not p.name.startswith(".")
                  and p.name.endswith(FASTQ_SUFFIXES))


def iter_fastq(path, name_contains=None, min_len=0):
    """Stream a FASTQ, plain or gzipped, as Read records."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            handle.readline()                      # the '+' line
            qual = handle.readline().strip()
            name = header[1:].strip().split()[0]
            if name_contains and name_contains not in name:
                continue
            if len(seq) < min_len:
                continue
            yield Read(name, seq, qual or None)


def read_fastq(path, name_contains=None, limit=None, min_len=0):
    """Read a FASTQ into Read records, taking the first *limit* if given."""
    reads = iter_fastq(path, name_contains, min_len)
    return list(reads if limit is None else itertools.islice(reads, limit))


def downsample(reads, k, seed=SAMPLE_SEED):
    """Take *k* of *reads* uniformly at random, returning ``(kept, seen)``.

    The reads are a stream of unknown length, so this is reservoir sampling: it
    holds only the k it keeps and makes one pass, which is what lets a run of any
    depth be capped without reading it into memory first.  Sampling rather than
    truncating matters for a directory — the first k reads of a pool are the
    first file's reads, and a page drawn from those is a page of one barcode.

    Every read is kept when there are no more than k of them, and the order the
    survivors come back in is not the order they were read; rows are sorted
    afterwards anyway.
    """
    kept = []
    rng = random.Random(seed)
    seen = 0
    for seen, read in enumerate(reads, start=1):
        if len(kept) < k:
            kept.append(read)
        else:
            j = rng.randrange(seen)
            if j < k:
                kept[j] = read
    return kept, seen


def _take(reads, limit):
    """Every read, or a sample of *limit* of them; ``(kept, seen)`` either way."""
    if not limit:
        kept = list(reads)
        return kept, len(kept)
    return downsample(reads, limit)


def collect_samples(paths, name_contains=None, limit=None, min_len=0, pooled=True):
    """Read *paths* into the ``(label, reads, seen)`` triples groups are built from.

    Pooled, every file's reads land in one list and the page draws a single
    pileup; per-file, each file becomes its own group so the samples can be read
    side by side.  A label of None means the group is named for the reference
    rather than for a file, which is what pooling leaves it with.  ``seen`` is
    how many reads the group was sampled from, which is the only place the depth
    behind a downsampled page survives.

    A cap samples the pool in the first case and each file in the second, so it
    bounds what one group draws either way.  A falsy cap draws everything.
    """
    if pooled:
        stream = itertools.chain.from_iterable(
            iter_fastq(path, name_contains, min_len) for path in paths)
        reads, seen = _take(stream, limit)
        return [(None, reads, seen)] if reads else []

    out = []
    for path in paths:
        reads, seen = _take(iter_fastq(path, name_contains, min_len), limit)
        if reads:
            out.append((path.stem, reads, seen))
        else:
            print(f"no reads matched in {path.name}", file=sys.stderr)
    return out


def write_fasta(reference, path, doubled=False):
    path = Path(path)
    seq = reference.seq * 2 if doubled else reference.seq
    path.write_text(f">{reference.name}\n{seq}\n")
    return str(path)


def write_log(out, args, lines):
    """Write a plain-text log of this run next to the HTML page it produced.

    Just the run's own stdout narration (reference, files, counts) plus the
    ordering function applied, so the page's provenance survives without
    re-running the command.
    """
    log_path = out.with_suffix(".log.txt")
    header = [f"seqviewer-pileup log: {out.stem}",
              f"ordering function: {args.order}", ""]
    log_path.write_text("\n".join(header + lines) + "\n")
    return log_path


def fold(row, length):
    """Collapse a row aligned against a doubled reference back to one copy.

    A read crossing the origin of a circular construct aligns contiguously in
    ``seq + seq``, so the two halves of the grid hold the two ends of that read.
    Merging them keeps the whole read on one line instead of discarding the
    segment minimap2 reported as supplementary.
    """
    merged = []
    for i in range(length):
        left, right = row[i], row[i + length]
        if left[0] == "-":
            merged.append(right)
        elif right[0] == "-":
            merged.append(left)
        else:
            # Both copies saw this base; a mismatch is the signal worth keeping.
            merged.append(left if not left[1] else right)
    return merged


def cluster(rows):
    """Order rows so reads sharing a mismatch pattern sit together."""
    return sorted(rows, key=lambda row: "".join(
        "." if is_match or base == "-" else base.upper() for base, is_match in row
    ))


def _extent(row):
    """Return (called bases, first covered position) for one row."""
    called = 0
    first = len(row)
    for i, (base, _) in enumerate(row):
        if base != "-":
            called += 1
            if i < first:
                first = i
    return called, first


def by_length(rows):
    """Longest read first, ties broken by where the read starts.

    Length is counted as bases actually called, not the aligned span, so a read
    with an internal deletion sorts by what it really contributed.  This replaces
    the mismatch clustering seqviewer applies by default, which groups
    subpopulations into blocks but interleaves long and short reads.
    """
    return sorted(rows, key=lambda row: (lambda c, f: (-c, f))(*_extent(row)))


def by_position(rows):
    """Leftmost read first, ties broken by decreasing length."""
    return sorted(rows, key=lambda row: (lambda c, f: (f, -c))(*_extent(row)))


#: Row orderings, by the name --order takes.  Each is called with the grid and
#: the reference; only the clustering one needs the second.
#:
#: "mismatch" and "cluster" both aim to put similar reads together and differ in
#: how: "mismatch" sorts the pattern as a string, which is cheap but is decided by
#: the leftmost difference, so one read with an early sequencing error sorts
#: between two halves of a subpopulation.  "cluster" builds the dendrogram.
ORDERINGS = {
    "length": lambda rows, ref: by_length(rows),
    "position": lambda rows, ref: by_position(rows),
    "mismatch": lambda rows, ref: cluster(rows),
    "cluster": lambda rows, ref: cluster_rows(rows, ref),
}


def focus_flanks(reference, label):
    """Derive ``(5' length, 3' length)`` around the feature matching *label*.

    The flanks are computed here rather than by re-typing the matched feature to
    ``"insert"`` and handing it to ``Reference.flank_lengths``: that added a
    second copy of one feature to the reference, which is invisible only for as
    long as nothing draws annotations.
    """
    if not label:
        return reference.flank_lengths()      # a feature already typed "insert"

    for feature in reference.features:
        if (label.lower() in (feature.label or "").lower()
                or label.lower() == feature.type.lower()):
            if feature.wraps_origin:
                print(f"feature {feature.label!r} wraps the origin; "
                      f"not using it for flanks", file=sys.stderr)
                return None
            five, three = feature.start, len(reference) - feature.end
            if five <= 0 and three <= 0:
                return None
            return (max(0, five), max(0, three))

    print(f"no feature matching {label!r}; flanks left unset", file=sys.stderr)
    return None


def build_parser(prog=None):
    """The command line, as a parser.

    Separated from :func:`main` so the flags and their defaults can be
    tested without running an alignment.  *prog* names the command in usage
    text, which differs between ``seqviewer-pileup`` and ``seqview pileup``.
    """
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("reads",
                        help="a directory of FASTQs, pooled into one pileup, or "
                             "a single FASTQ; .gz is read directly")
    parser.add_argument("reference",
                        help="FASTA, GenBank (.gb/.gbk/.ape), or SnapGene (.dna); "
                             "the annotated formats need biopython")
    parser.add_argument("--skip-types", metavar="TYPES",
                        help="comma-separated feature types to drop, replacing "
                             "the default (source,primer_bind).  Pass "
                             "'source' to keep primer sites, which is what "
                             "shows the Golden Gate junctions on an amplicon; "
                             "pass '' to keep every feature the file declares")
    parser.add_argument("out", help="HTML page to write; .html is appended "
                                     "if not already there")
    parser.add_argument("--per-file", action="store_true",
                        help="draw one group per FASTQ instead of pooling the "
                             "directory into a single pileup")
    parser.add_argument("--name", help="keep only reads whose name contains this")
    parser.add_argument("--max", type=int, default=DEFAULT_MAX_READS,
                        help="reads drawn per group, sampled uniformly from "
                             "across the whole of every file rather than taken "
                             f"from the front (default {DEFAULT_MAX_READS}); "
                             "0 draws every read, which for a deep run is a "
                             "page too large to open")
    parser.add_argument("--min-read-len", type=int, default=0,
                        help="drop reads shorter than this before aligning")
    parser.add_argument("--insert", help="feature label to derive flanks from")
    parser.add_argument("--min-overlap-pos", type=int, default=0,
                        help="drop reads not crossing this position; "
                             "0 = keep everything (default), "
                             "-1 = reference midpoint")
    parser.add_argument("--order", choices=sorted(ORDERINGS), default="length",
                        help="row order: length (longest first, the default), "
                             "position (leftmost first), mismatch (sorts the "
                             "mismatch pattern as a string), or cluster "
                             "(hierarchical, average linkage over what each "
                             "read disagrees about \u2014 groups a "
                             "subpopulation that mismatch splits when a read "
                             "carries an unrelated early error)")
    parser.add_argument("--no-circular", action="store_true",
                        help="align against one copy even if the reference is "
                             "circular, discarding origin-crossing segments")
    parser.add_argument("--mismatch-freq", action="store_true",
                        help="accepted and ignored: the mismatch track is now "
                             "always drawn, having replaced the row of flag "
                             "triangles. Kept so existing commands still run")
    parser.add_argument("--ref-name",
                        help="name for the reference, overriding what the file "
                             "carries (a GenBank/SnapGene file's declared name, "
                             "or a FASTA header, is often not human-friendly); "
                             "used as the pooled group name and in the default "
                             "title")
    parser.add_argument("--summary", action="store_true",
                        help="also write a summarized page beside the pileup: "
                             "the construct as an annotated map, one compact "
                             "band per group with a lollipop per called "
                             "variant, and the variants as a table. Named "
                             "<out>.summary.html")
    parser.add_argument("--variant-freq", type=float,
                        default=DEFAULT_MIN_FRACTION, metavar="F",
                        help="share of covering reads an allele needs before "
                             "the summary calls it (default "
                             f"{DEFAULT_MIN_FRACTION:g}). Lower it to see "
                             "events the default suppresses as sequencing "
                             "error; only meaningful with --summary")
    parser.add_argument("--variant-reads", type=int,
                        default=DEFAULT_MIN_COUNT, metavar="N",
                        help="reads that must support an allele whatever the "
                             f"share (default {DEFAULT_MIN_COUNT}). At shallow "
                             "depth a single read is the error rate, not a "
                             "variant; only meaningful with --summary")
    parser.add_argument("--title")
    return parser


def main(argv=None, prog=None):
    parser = build_parser(prog)
    args = parser.parse_args(argv)

    log_lines = []

    def log(msg):
        print(msg)
        log_lines.append(msg)

    skip_types = (None if args.skip_types is None
                  else tuple(t for t in args.skip_types.split(",") if t))
    reference = load_reference(args.reference, skip_types=skip_types)
    if args.ref_name:
        reference.name = args.ref_name
    log(f"reference {reference.name}: {len(reference)} bp, {reference.topology}, "
        f"{len(reference.features)} features")

    paths = fastq_paths(args.reads)
    if not paths:
        print(f"no FASTQ files at {args.reads}", file=sys.stderr)
        return 1
    log(f"{len(paths)} FASTQ file{'s' if len(paths) != 1 else ''}: "
        + ", ".join(p.name for p in paths))

    samples = collect_samples(paths, args.name, args.max, args.min_read_len,
                              pooled=not args.per_file)
    if not samples:
        print(f"no reads matched in {args.reads}", file=sys.stderr)
        return 1
    total_reads = sum(len(reads) for _, reads, _ in samples)
    total_seen = sum(seen for _, _, seen in samples)
    sampled = " sampled from {:,}".format(total_seen) if total_seen > total_reads else ""
    log(f"{total_reads:,} reads to align{sampled} in {len(samples)} "
        f"group{'s' if len(samples) != 1 else ''}, ordered by {args.order}")
    if sampled:
        print(f"pass --max 0 to draw all {total_seen:,}, at roughly "
              f"{total_seen * 2 // 1000:,} MB of HTML")

    # A circular reference is aligned against two copies of itself so that reads
    # crossing the origin stay in one piece, then folded back to one copy.
    circular = reference.is_circular and not args.no_circular
    out = Path(args.out)
    if out.suffix.lower() != ".html":
        out = out.with_name(out.name + ".html")
    fasta = write_fasta(reference, out.with_suffix(".ref.fasta"), doubled=circular)
    align_seq = reference.seq * 2 if circular else reference.seq
    if circular:
        log("circular: aligning against a doubled reference, then folding")

    groups = []
    for label, reads, seen in samples:
        name = label or reference.name
        rows = grid_from_reads(reads, fasta, align_seq,
                               min_overlap_pos=args.min_overlap_pos)
        if not rows:
            # One empty sample is not a failed run: the rest still draw.
            print(f"{name}: nothing aligned", file=sys.stderr)
            continue
        if circular:
            rows = [fold(row, len(reference)) for row in rows]
        # grid_from_reads has already clustered by mismatch pattern; reorder
        # unless that is what was asked for.
        rows = ORDERINGS[args.order](rows, reference.seq)
        covered = sum(1 for i in range(len(reference))
                      if any(row[i][0] != "-" for row in rows))
        of_seen = f", sampled from {seen:,}" if seen > len(reads) else ""
        log(f"{name}: {len(rows)} of {len(reads):,} reads drawn "
            f"({len(rows) / len(reads):.0%}){of_seen}; {covered} of "
            f"{len(reference)} positions covered "
            f"({covered / len(reference):.0%})")
        # No status: it is a consensus call, and this driver has no consensus to
        # report.  A constant string here would be styled as one and say nothing.
        groups.append(PileupGroup(name=name, ref_seq=reference.seq, rows=rows,
                                  n_reads=len(rows),
                                  fraction=len(rows) / total_reads))

    if not groups:
        print("nothing aligned", file=sys.stderr)
        return 1

    title = args.title or f"Pileup: {reference.name}"
    # Built directly rather than through PileupView.from_reference so that the
    # focus region comes from --insert by label, rather than from a feature the
    # file happens to have typed "insert".  Everything else from_reference would
    # carry across still has to be carried.
    view = PileupView(
        title=title, groups=groups, total_reads=total_reads,
        flanks=focus_flanks(reference, args.insert),
        features=reference.features,
        ref_len=len(reference),
    )
    log(f"flanks: {view.flanks} | {len(view.features)} features drawn")
    out.write_text(render(view))
    print(f"Wrote {out.resolve()}")

    if args.summary:
        # Reduced from the same view the pileup drew, so the two pages cannot
        # disagree about the reference, the focus region, or the features.
        summary = SummaryView.from_view(
            view,
            min_fraction=args.variant_freq,
            min_count=args.variant_reads,
        )
        log(f"summary: calling at >={args.variant_freq:.0%} of covering reads "
            f"and >={args.variant_reads} supporting")
        for group in summary.groups:
            called = len(group.variants)
            log(f"{group.name}: {called} variant{'s' if called != 1 else ''} "
                f"called, {group.verdict}; {group.mean_depth:.0f}x mean depth "
                f"over {group.covered} of {group.ref_len} positions")
        summary_out = summary_path(out)
        summary_out.write_text(render_summary(summary))
        print(f"Wrote {summary_out.resolve()}")

    log_path = write_log(out, args, log_lines)
    print(f"Wrote {log_path.resolve()}")
    return 0


#: Bytes a run must hold before a scan reports progress.  Below it the scan is
#: over before a bar would be read.
PROGRESS_AFTER_BYTES = 32 << 20

#: Seconds between redraws of the progress line.
PROGRESS_INTERVAL = 0.1

#: Seconds a scan must take before its rate is reported alongside the figures.
SCAN_NOTE_AFTER = 1.0

#: Seconds a scan must run before a rate is worth dividing out.
RATE_AFTER = 0.05


def _si_bytes(count):
    """Format a byte count in binary multiples."""
    for unit, size in (("GiB", 1 << 30), ("MiB", 1 << 20), ("KiB", 1 << 10)):
        if count >= size:
            return f"{count / size:.1f} {unit}"
    return f"{count} B"


def _si_reads(count):
    """Format a read count, abbreviated past a thousand."""
    if count >= 1_000_000:
        return f"{count / 1_000_000:.1f}M"
    if count >= 1_000:
        return f"{count / 1_000:.0f}k"
    return str(count)


def _duration(seconds):
    """Format a duration: to a tenth of a second under ten, then whole seconds,
    then minutes and seconds.  A fast scan otherwise reports as 0s or 1s.
    """
    if seconds < 10:
        return f"{seconds:.1f}s"
    seconds = int(round(seconds))
    if seconds < 60:
        return f"{seconds}s"
    return f"{seconds // 60}m{seconds % 60:02d}s"


#: Seconds between redraws of the live histogram.  Slower than the progress
#: line, because each frame reduces the whole tally to an axis.
LIVE_INTERVAL = 0.25

#: Rows a live histogram needs beyond its bins: a leading blank, the axis
#: header, two tail rows, a blank, three lines of figures and a progress line.
LIVE_EXTRA_ROWS = 9


def _progress_line(share, reads, rate, left, width, palette=None):
    """Return one line reporting how far a scan has got, at most *width* wide.

    The bar is dropped where the figures alone fill the width, since the figures
    are what the reader needs.
    """
    text = f"{share * 100:3.0f}%  {_si_reads(reads)} reads"
    if rate:
        text += f"  {_si_bytes(int(rate))}/s  {_duration(left)} left"
    room = width - len(text) - 4
    if room < 12:
        return text[:width]
    fill = int(round(room * share))
    done, todo = "█" * fill, "░" * (room - fill)
    if palette:
        done = f"{palette.bar}{done}{palette.reset}"
        todo = f"{palette.tail}{todo}{palette.reset}"
    return f"{done}{todo}  {text}"


class ScanProgress:
    """One rewritten line reporting the progress of a scan, on stderr.

    Drawn only where stderr is a terminal and the run holds at least
    :data:`PROGRESS_AFTER_BYTES`, so a small run and a redirected one print
    nothing.  Progress goes to stderr so that the histogram on stdout can be
    redirected to a file on its own.

    The bytes and reads seen are recorded whether or not the line is drawn, so
    the scan can be reported once it finishes.
    """

    def __init__(self, width, palette=None, enabled=True, stream=None):
        self.width = max(32, width)
        self.palette = palette
        self.stream = stream if stream is not None else sys.stderr
        self.enabled = enabled and self.stream.isatty()
        self.started = time.monotonic()
        self.took = 0.0
        self.last = 0.0
        self.drawn = False
        self.bytes = 0
        self.reads = 0

    def figures(self, done, total):
        """Return the share done, and the rate and time left once measurable.

        Under :data:`RATE_AFTER` the elapsed time is too short to divide by, and
        the rate it gives is wrong by orders of magnitude, so it is left out.
        """
        share = min(1.0, done / total) if total else 1.0
        if self.took < RATE_AFTER:
            return share, 0.0, 0.0
        rate = done / self.took
        return share, rate, (total - done) / rate if rate else 0.0

    def update(self, done, total, reads, snapshot=None):
        """Record progress, and redraw the line if it is due."""
        self.bytes, self.reads = done, reads
        self.took = time.monotonic() - self.started
        if not self.enabled or total < PROGRESS_AFTER_BYTES:
            return
        now = time.monotonic()
        if done < total and now - self.last < PROGRESS_INTERVAL:
            return
        self.last = now
        share, rate, left = self.figures(done, total)
        line = _progress_line(share, reads, rate, left, self.width,
                              self.palette)
        self.stream.write("\r\x1b[2K" + line)
        self.stream.flush()
        self.drawn = True

    def clear(self):
        """Erase the line, leaving the cursor where the next output starts."""
        if self.drawn:
            self.stream.write("\r\x1b[2K")
            self.stream.flush()
            self.drawn = False

    def note(self):
        """Return how much was scanned and how fast, or None if it was quick."""
        if self.took < SCAN_NOTE_AFTER or not self.bytes:
            return None
        rate = _si_bytes(int(self.bytes / self.took))
        return (f"scanned {_si_bytes(self.bytes)} in "
                f"{_duration(self.took)} at {rate}/s")


class LiveView:
    """The histogram, redrawn in place while the run is scanned.

    Used only where stdout is a terminal and the frame fits the window, since
    redrawing moves the cursor back over rows that must still be on screen.  A
    caller that is not drawing gets one write at the end instead.

    The axis is recomputed for each frame, so it moves until the last block is
    counted.  Frames are drawn at :data:`LIVE_INTERVAL` rather than once per
    block, since each one reduces the whole tally to an axis.
    """

    def __init__(self, width, bins, bulk, log, palette=None, rows=24,
                 enabled=True, stream=None):
        self.width = width
        self.bins = bins
        self.bulk = bulk
        self.log = log
        self.palette = palette
        self.stream = stream if stream is not None else sys.stdout
        self.enabled = (enabled and self.stream.isatty()
                        and rows > bins + LIVE_EXTRA_ROWS)
        self.started = time.monotonic()
        self.last = 0.0
        self.took = 0.0
        self.bytes = 0
        self.reads = 0
        self.height = 0
        self.drew = False

    def update(self, done, total, reads, snapshot=None):
        """Record progress, and redraw the histogram if a frame is due."""
        self.bytes, self.reads = done, reads
        self.took = time.monotonic() - self.started
        if not (self.enabled and snapshot) or total < PROGRESS_AFTER_BYTES:
            return
        now = time.monotonic()
        if now - self.last < LIVE_INTERVAL:
            return
        self.last = now
        share = min(1.0, done / total) if total else 1.0
        rate = done / self.took if self.took >= RATE_AFTER else 0.0
        left = (total - done) / rate if rate else 0.0
        self.draw(self.frame(snapshot(),
                             _progress_line(share, reads, rate, left,
                                            self.width, self.palette)))

    def frame(self, counts, last_line=None):
        """Return the lines of one frame: the histogram, then the figures."""
        if counts.empty:
            lines = ["", "no reads"]
        else:
            binning = lengths.bin_counts(counts, self.bins, self.bulk)
            summary = lengths.summarise_counts(counts)
            lines = [""]
            lines += lengths.histogram(binning, self.width, self.log,
                                       self.palette)
            lines.append("")
            for text in lengths.summary_lines(summary, binning):
                lines += textwrap.fill(text, self.width).split("\n")
        if last_line:
            lines.append(last_line)
        return lines

    def draw(self, lines, final=False):
        """Write *lines* over the frame already on screen.

        Each line is erased before it is written, so a shorter line does not
        leave the end of the last one behind, and a frame with fewer rows than
        the one before erases the rows it no longer fills.

        The cursor is left at the top of where the next frame goes, which is the
        row after the last line written.  On the *final* frame it is left below
        every row instead, so that whatever prints next -- another group, or the
        shell prompt -- does not land inside the histogram.
        """
        out = []
        if self.height:
            out.append(f"\x1b[{self.height}A")
        out += [f"\x1b[2K{line}\n" for line in lines]
        spare = max(0, self.height - len(lines))
        if spare:
            out += ["\x1b[2K\n"] * spare
            if not final:
                out.append(f"\x1b[{spare}A")
        self.stream.write("".join(out))
        self.stream.flush()
        self.height = len(lines)
        self.drew = True

    def note(self):
        """Return how much was scanned and how fast, or None if it was quick."""
        if self.took < SCAN_NOTE_AFTER or not self.bytes:
            return None
        rate = _si_bytes(int(self.bytes / self.took))
        text = (f"scanned {_si_bytes(self.bytes)} in {_duration(self.took)} "
                f"at {rate}/s")
        if self.palette:
            return f"{self.palette.tail}{text}{self.palette.reset}"
        return text

    def finish(self, counts):
        """Write the settled histogram, over the last frame or on its own."""
        lines = self.frame(counts, self.note())
        if self.drew:
            self.draw(lines, final=True)
        else:
            for line in lines:
                print(line)


def lengths_main(argv=None, prog="seqviewer-lengths"):
    """Print a read-length histogram for a directory of FASTQs, or one file."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Plot the read-length distribution of a sequencing run as "
                    "a histogram in the terminal. Takes a directory of FASTQs "
                    "or a single file; .gz is read directly.")
    parser.add_argument("reads",
                        help="a directory of FASTQs, or a single FASTQ")
    parser.add_argument("--bins", type=int, default=lengths.DEFAULT_BINS,
                        metavar="N",
                        help="bins across the axis (default "
                             f"{lengths.DEFAULT_BINS})")
    parser.add_argument("--bulk", type=float, default=lengths.DEFAULT_BULK,
                        metavar="PCT",
                        help="percent of reads the axis covers, centred "
                             f"(default {lengths.DEFAULT_BULK:g}). Concatemers "
                             "and other ultra-long artifacts otherwise reach "
                             "the top of the axis on their own. 100 spans the "
                             "full range")
    parser.add_argument("--log", action="store_true",
                        help="scale bar length by log(1 + count), which keeps "
                             "the smaller bins of a peaked distribution "
                             "distinguishable")
    parser.add_argument("--width", type=int, metavar="COLS",
                        help="output width (default: the terminal's, or 80)")
    parser.add_argument("--per-file", action="store_true",
                        help="one histogram per FASTQ instead of one for the "
                             "whole directory")
    parser.add_argument("--no-color", action="store_true",
                        help="write plain text. Colour is already left out "
                             "when stdout is not a terminal, or when NO_COLOR "
                             "is set")
    parser.add_argument("--no-progress", action="store_true",
                        help="do not report progress while scanning")
    parser.add_argument("--no-live", action="store_true",
                        help="draw the histogram once the scan finishes rather "
                             "than filling it in as reads are counted")
    parser.add_argument("--slow", action="store_true",
                        help="scan without numpy, which is slower on a large "
                             "file and reports the same figures")
    args = parser.parse_args(argv)

    paths = fastq_paths(args.reads)
    if not paths:
        print(f"no FASTQ files at {args.reads}", file=sys.stderr)
        return 1

    window = shutil.get_terminal_size((80, 24))
    width = args.width or window.columns
    palette = None
    if not (args.no_color or os.environ.get("NO_COLOR")
            or not sys.stdout.isatty()):
        palette = lengths.PALETTE

    def dim(text):
        return f"{palette.tail}{text}{palette.reset}" if palette else text

    groups = ([(p.name, [p]) for p in paths] if args.per_file
              else [(str(args.reads), paths)])

    for index, (label, group) in enumerate(groups):
        if index:
            print()
        files = f"{len(group)} file{'s' if len(group) != 1 else ''}"
        print(f"{label}  {dim('·')}  {dim(files)}")

        view = LiveView(width, args.bins, args.bulk, args.log, palette=palette,
                        rows=window.lines, enabled=not args.no_live)
        # A live view already reports its own progress, so the line on stderr
        # is only for a scan that is not being drawn.
        bar = ScanProgress(width, palette=palette,
                           enabled=not args.no_progress and not view.enabled)

        def report(done, total, reads, snapshot=None):
            view.update(done, total, reads, snapshot)
            bar.update(done, total, reads)

        counts = lengths.count_lengths(group, progress=report,
                                       fast=False if args.slow else None)
        bar.clear()
        view.took = view.took or bar.took
        view.bytes = view.bytes or bar.bytes
        view.finish(counts)
    return 0


def seqview_main(argv=None):
    """Dispatch ``seqview <command>`` to that command's own parser.

    Each command is installed under its own name as well -- ``seqviewer-pileup``
    and ``seqviewer-lengths`` -- and takes the same arguments either way.
    """
    commands = {
        "pileup": (main, "align reads to a reference and write a pileup page"),
        "lengths": (lengths_main,
                    "plot the read-length distribution in the terminal"),
    }
    argv = list(sys.argv[1:] if argv is None else argv)

    if argv and argv[0] in commands:
        name = argv[0]
        return commands[name][0](argv[1:], prog=f"seqview {name}")

    label = max(len(name) for name in commands)
    usage = ["usage: seqview <command> [options]", "", "commands:"]
    usage += [f"  {name:<{label}}  {blurb}"
              for name, (_, blurb) in commands.items()]
    usage += ["", "seqview <command> --help lists that command's options."]

    asked = argv[:1] in (["-h"], ["--help"])
    stream = sys.stdout if asked else sys.stderr
    if argv and not asked:
        print(f"seqview: unknown command {argv[0]!r}", file=stream)
    print("\n".join(usage), file=stream)
    return 0 if asked else 2


if __name__ == "__main__":
    raise SystemExit(seqview_main())
