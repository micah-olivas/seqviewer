#!/usr/bin/env python
"""Align a FASTQ to a reference and write a seqviewer pileup page.

    pileup.py reads.fastq reference out.html [options]

Reading the reference is :mod:`seqviewer.genbank`'s job — FASTA, GenBank, ApE,
and SnapGene all arrive as a ``Reference`` with its topology and features.  What
is left here is the part the package does not cover: reading a FASTQ, and the
two adjustments a plasmid needs.  A circular reference is aligned against a
doubled copy so reads crossing the origin stay whole, and a named feature can
supply the flanks the page marks.

seqviewer takes reads as ``Read(name, seq, qual)`` records, so FASTQ parsing is
the caller's job.  Everything after that is grid_from_reads plus render.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

from seqviewer import PileupGroup, PileupView, render
from seqviewer.align import Read, grid_from_reads
from seqviewer.genbank import load_reference


def read_fastq(path, name_contains=None, limit=None, min_len=0):
    """Read a FASTQ, plain or gzipped, into Read records."""
    opener = gzip.open if str(path).endswith(".gz") else open
    out = []
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
            out.append(Read(name, seq, qual or None))
            if limit and len(out) >= limit:
                break
    return out


def write_fasta(reference, path, doubled=False):
    path = Path(path)
    seq = reference.seq * 2 if doubled else reference.seq
    path.write_text(f">{reference.name}\n{seq}\n")
    return str(path)


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


ORDERINGS = {"length": by_length, "position": by_position, "mismatch": cluster}


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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq")
    parser.add_argument("reference",
                        help="FASTA, GenBank (.gb/.gbk/.ape), or SnapGene (.dna); "
                             "the annotated formats need biopython")
    parser.add_argument("--skip-types", metavar="TYPES",
                        help="comma-separated feature types to drop, replacing "
                             "the default (source,primer_bind).  Pass "
                             "'source' to keep primer sites, which is what "
                             "shows the Golden Gate junctions on an amplicon; "
                             "pass '' to keep every feature the file declares")
    parser.add_argument("out", help="HTML page to write")
    parser.add_argument("--name", help="keep only reads whose name contains this")
    parser.add_argument("--max", type=int, help="cap on reads drawn")
    parser.add_argument("--min-read-len", type=int, default=0,
                        help="drop reads shorter than this before aligning")
    parser.add_argument("--insert", help="feature label to derive flanks from")
    parser.add_argument("--min-overlap-pos", type=int, default=0,
                        help="drop reads not crossing this position; "
                             "0 = keep everything (default), "
                             "-1 = reference midpoint")
    parser.add_argument("--order", choices=sorted(ORDERINGS), default="length",
                        help="row order: length (longest first, the default), "
                             "position (leftmost first), or mismatch "
                             "(seqviewer's clustering, which blocks "
                             "subpopulations together)")
    parser.add_argument("--no-circular", action="store_true",
                        help="align against one copy even if the reference is "
                             "circular, discarding origin-crossing segments")
    parser.add_argument("--title")
    args = parser.parse_args(argv)

    skip_types = (None if args.skip_types is None
                  else tuple(t for t in args.skip_types.split(",") if t))
    reference = load_reference(args.reference, skip_types=skip_types)
    print(f"reference {reference.name}: {len(reference)} bp, {reference.topology}, "
          f"{len(reference.features)} features")

    reads = read_fastq(args.fastq, args.name, args.max, args.min_read_len)
    if not reads:
        print(f"no reads matched in {args.fastq}", file=sys.stderr)
        return 1
    print(f"{len(reads)} reads to align")

    # A circular reference is aligned against two copies of itself so that reads
    # crossing the origin stay in one piece, then folded back to one copy.
    circular = reference.is_circular and not args.no_circular
    out = Path(args.out)
    fasta = write_fasta(reference, out.with_suffix(".ref.fasta"), doubled=circular)
    align_seq = reference.seq * 2 if circular else reference.seq
    if circular:
        print("circular: aligning against a doubled reference, then folding")

    rows = grid_from_reads(reads, fasta, align_seq,
                           min_overlap_pos=args.min_overlap_pos)
    if not rows:
        print("nothing aligned", file=sys.stderr)
        return 1
    if circular:
        rows = [fold(row, len(reference)) for row in rows]
    # grid_from_reads has already clustered by mismatch pattern; reorder unless
    # that is what was asked for.
    rows = ORDERINGS[args.order](rows)
    print(f"{len(rows)} of {len(reads)} reads drawn "
          f"({len(rows) / len(reads):.0%}), ordered by {args.order}")
    covered = sum(1 for i in range(len(reference))
                  if any(row[i][0] != "-" for row in rows))
    print(f"{covered} of {len(reference)} reference positions covered "
          f"({covered / len(reference):.0%})")

    # No status: it is a consensus call, and this driver has no consensus to
    # report.  A constant string here would be styled as one and say nothing.
    group = PileupGroup(name=reference.name, ref_seq=reference.seq, rows=rows,
                        n_reads=len(rows), fraction=len(rows) / len(reads))
    title = args.title or f"Pileup: {reference.name}"
    view = PileupView(
        title=title, groups=[group], total_reads=len(reads),
        flanks=focus_flanks(reference, args.insert),
    )
    print(f"flanks: {view.flanks}")
    out.write_text(render(view))
    print(f"Wrote {out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
