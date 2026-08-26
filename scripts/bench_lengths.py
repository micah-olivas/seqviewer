"""Time the read-length scanners over a synthetic FASTQ.

Writes a file of the requested size, then scans it three ways and reports the
wall time and peak resident memory of each:

* ``list`` holds every length, then reduces, which is what a caller keeping the
  lengths themselves pays;
* ``python`` tallies by length, splitting a block and taking every fourth line;
* ``array`` tallies by length, locating the newlines in a block with numpy and
  taking the gaps between them.

The three agree on every figure, so the difference between them is in time and
memory only.  Each runs in its own process, because peak resident memory is a
high-water mark that one method would otherwise inherit from the last.

The figures quoted in the README come from this script.  Run it to recheck them,
and say which machine and whether the file was already in the page cache: a
first read from disk is bound by the disk rather than by the scanner.

    python scripts/bench_lengths.py --gigabytes 2
"""

from __future__ import annotations

import argparse
import json
import random
import resource
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from seqviewer import lengths  # noqa: E402  (after the path is set)

METHODS = ("list", "python", "array")


def write_run(path: Path, gigabytes: float, product: int = 820) -> int:
    """Write a FASTQ of about *gigabytes*: one product, with a few concatemers.

    One block of bases is reused and sliced at varying offsets, so the time goes
    into writing rather than into generating sequence.
    """
    rng = random.Random(7)
    block = "".join(rng.choice("ACGT") for _ in range(8192))
    quals = "I" * 8192
    target = int(gigabytes * (1 << 30))
    written = reads = 0
    with open(path, "w", buffering=1 << 22) as handle:
        while written < target:
            n = max(60, int(rng.gauss(product, 60)))
            if rng.random() < 0.002:
                n = product * rng.choice([2, 3, 5])
            start = (reads * 131) % (8192 - n)
            record = f"@read_{reads}\n{block[start:start + n]}\n+\n{quals[:n]}\n"
            handle.write(record)
            written += len(record)
            reads += 1
    return reads


def peak_mb() -> float:
    """Peak resident size of this process, in megabytes."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":
        peak *= 1024                    # Linux reports kilobytes, macOS bytes
    return peak / (1 << 20)


def run_one(method: str, path: Path) -> dict:
    """Scan *path* by one method and return its timing, memory and figures."""
    started = time.monotonic()
    if method == "list":
        summary, _ = lengths.distribution(list(lengths.read_lengths([path])))
    else:
        counts = lengths.count_lengths([path], fast=(method == "array"))
        summary = lengths.summarise_counts(counts)
    took = time.monotonic() - started
    return {
        "method": method,
        "seconds": took,
        "peak_mb": peak_mb(),
        "reads": summary.reads,
        "bases": summary.bases,
        "median": summary.median,
        "n50": summary.n50,
        "longest": summary.longest,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gigabytes", type=float, default=2.0,
                        help="size of the FASTQ to write (default 2)")
    parser.add_argument("--path", type=Path,
                        default=Path("bench_reads.fastq"),
                        help="where to write it (default ./bench_reads.fastq)")
    parser.add_argument("--keep", action="store_true",
                        help="leave the file in place afterwards")
    parser.add_argument("--method", choices=METHODS,
                        help="run one method in this process and print its "
                             "result as JSON. Used by the parent run")
    args = parser.parse_args(argv)

    if args.method:
        print(json.dumps(run_one(args.method, args.path)))
        return 0

    fresh = not args.path.exists()
    if fresh:
        print(f"writing {args.gigabytes:g} GiB to {args.path} ...")
        print(f"{write_run(args.path, args.gigabytes):,} reads")

    size = args.path.stat().st_size / (1 << 30)
    numpy = "installed" if lengths.HAVE_NUMPY else "absent"
    print(f"\n{args.path}, {size:.2f} GiB.  numpy {numpy}.  "
          f"python {sys.version.split()[0]} on {sys.platform}.\n")

    results = []
    for method in METHODS:
        if method == "array" and not lengths.HAVE_NUMPY:
            print(f"  {'array':<26}skipped, numpy is not installed")
            continue
        out = subprocess.run(
            [sys.executable, __file__, "--method", method,
             "--path", str(args.path)],
            capture_output=True, text=True, check=True)
        r = json.loads(out.stdout)
        results.append(r)
        print(f"  {method:<12}{r['seconds']:6.2f} s   "
              f"{size / r['seconds']:5.2f} GiB/s   "
              f"peak {r['peak_mb']:5.0f} MB")

    figures = {(r["reads"], r["bases"], r["median"], r["n50"], r["longest"])
               for r in results}
    print(f"\nthe methods agree on every figure: {len(figures) == 1}")
    first = results[0]
    print(f"{first['reads']:,} reads · {first['bases']:,} bases · "
          f"median {first['median']:,} · N50 {first['n50']:,} · "
          f"max {first['longest']:,}")

    if fresh and not args.keep:
        args.path.unlink()
        print(f"\nremoved {args.path}; pass --keep to reuse it")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
