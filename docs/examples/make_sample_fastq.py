"""Write a synthetic FASTQ for the lengths and quality examples on the docs.

    python docs/examples/make_sample_fastq.py docs/examples/sample_reads.fastq

The real run the page was first written against is not in the repository, so
the example is built from this instead.  The file is seeded, so it is the same
every time, and it is not committed: rerun this, then

    seqview lengths docs/examples/sample_reads.fastq --bins 14 \\
        > docs/examples/lengths.txt
    seqview qualities docs/examples/sample_reads.fastq --bins 14 \\
        > docs/examples/quality.txt
    seqview qualities docs/examples/sample_reads.fastq --bins 14 --by read \\
        > docs/examples/quality_by_read.txt

The base qualities start high and fall along the read, with a read-to-read
spread, a few poor reads and the odd poor base, as on an Illumina run.  The
lengths mimic a small amplicon library: a large population of short fragments,
a smaller one at the full construct length, scattered reads between, and a few
concatemers longer than the construct.
"""
from __future__ import annotations

import random
import sys

CONSTRUCT = 3549
SEED = 7


def lengths(rng: random.Random) -> list[int]:
    out = []
    out += [int(rng.gauss(160, 28)) for _ in range(4300)]           # fragments
    out += [rng.randint(373, 1352) for _ in range(640)]             # in between
    out += [rng.randint(1353, 3300) for _ in range(60)]
    out += [CONSTRUCT - rng.randint(0, 40) for _ in range(260)]     # full length
    out += [rng.randint(CONSTRUCT + 5, 4700) for _ in range(25)]   # concatemers
    out += [rng.randint(105, 127) for _ in range(24)]               # truncated
    out = [max(100, n) for n in out]
    rng.shuffle(out)
    return out


def qualities(n: int, rng: random.Random) -> str:
    """Phred+33 scores for a read of *n* bases, falling toward the 3' end."""
    if rng.random() < 0.07:                         # a poor read
        start, drop = rng.gauss(26, 4), rng.uniform(8, 16)
    else:
        start, drop = rng.gauss(37, 2), rng.uniform(2, 9)
    out = []
    for i in range(n):
        q = rng.gauss(start - drop * i / n, 2.2)
        if rng.random() < 0.004:
            q = rng.uniform(2, 12)                  # a poor base
        out.append(chr(33 + max(2, min(41, int(round(q))))))
    return "".join(out)


def main(path: str) -> None:
    rng = random.Random(SEED)
    with open(path, "w") as fh:
        for i, n in enumerate(lengths(rng), 1):
            seq = "".join(rng.choices("ACGT", k=n))
            fh.write(f"@read{i}\n{seq}\n+\n{qualities(n, rng)}\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "sample_reads.fastq")
