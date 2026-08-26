"""The read-length distribution as a PNG.

Draws the binning the terminal histogram draws, from the same
:class:`seqviewer.lengths.Binning`, so the figure and the terminal agree on
where the axis falls and which reads are outside it.

matplotlib is imported inside :func:`write_png` rather than at the top of the
module, so that importing :mod:`seqviewer` does not require it.  It is the only
part of the package that needs anything beyond the standard library to produce
output, and it is reached only when a PNG is asked for.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .lengths import Binning, Summary

__all__ = ["DPI", "FIGSIZE", "PngUnavailable", "write_png"]

#: Dots per inch.  At the figure size below this is about 1800 px across, which
#: is enough to paste into a slide without the text going soft.
DPI = 200

#: Figure size in inches, wide rather than tall: a length axis needs the width.
FIGSIZE = (9.0, 5.0)

_BAR = "#2f8f9d"
_TAIL = "#b8c4c9"
_GRID = "#dfe5e7"
_INK = "#22333b"


class PngUnavailable(RuntimeError):
    """Raised when matplotlib is not installed."""


def _bar_caption(binning: Binning, summary: Summary) -> str:
    """Return the lines describing what the axis leaves out and what was read."""
    lines = [
        f"{summary.reads:,} reads · {summary.bases:,} bases · "
        f"median {summary.median:,} · mean {summary.mean:,.0f} · "
        f"N50 {summary.n50:,} · min {summary.shortest:,} · "
        f"max {summary.longest:,} bp"
    ]
    if binning.clipped:
        inside = summary.reads - binning.below - binning.above
        share = 100 * inside / summary.reads if summary.reads else 0.0
        outside = []
        if binning.below:
            outside.append(f"{binning.below:,} shorter than {binning.low:,}")
        if binning.above:
            outside.append(f"{binning.above:,} longer than {binning.high:,}")
        lines.append(
            f"axis holds {inside:,} of {summary.reads:,} reads ({share:.1f}%); "
            + " and ".join(outside) + " are not drawn"
        )
    return "\n".join(lines)


def write_png(
    path,
    binning: Binning,
    summary: Summary,
    title: Optional[str] = None,
    log: bool = False,
    dpi: int = DPI,
) -> Path:
    """Draw *binning* to *path* as a PNG, and return the path written.

    Bars are placed at the middle of each bin and are as wide as the bin, so the
    length axis reads linearly whether or not *log* scales the counts.  Reads
    outside the axis are named in the caption rather than drawn, since they have
    no width on it.

    Raises :class:`PngUnavailable` where matplotlib is not installed.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")           # a file, so no display is needed
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
    except ImportError as exc:          # pragma: no cover - matplotlib is opt-in
        raise PngUnavailable(
            "drawing a PNG needs matplotlib: pip install 'seqviewer[plot]'"
        ) from exc

    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(figsize=FIGSIZE)
    if binning.bins:
        width = binning.bins[0].high - binning.bins[0].low
        middles = [b.low + (b.high - b.low) / 2 for b in binning.bins]
        axes.bar(middles, [b.count for b in binning.bins],
                 width=width, color=_BAR, edgecolor="white", linewidth=0.4)
        axes.set_xlim(binning.bins[0].low - width * 0.5,
                      binning.bins[-1].high + width * 0.5)
    if log:
        axes.set_yscale("log")

    axes.set_xlabel("read length (bp)", color=_INK)
    axes.set_ylabel("reads" + (" (log)" if log else ""), color=_INK)
    if title:
        axes.set_title(title, color=_INK, loc="left")

    axes.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
    axes.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
    axes.grid(axis="y", color=_GRID, linewidth=0.7)
    axes.set_axisbelow(True)
    for side in ("top", "right"):
        axes.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axes.spines[side].set_color(_GRID)
    axes.tick_params(colors=_INK, labelsize=9)

    # The caption is placed from the bottom edge and the axes are lifted clear
    # of it, so its second line does not open a gap when there is only one.
    caption = _bar_caption(binning, summary)
    figure.text(0.012, 0.025, caption, fontsize=8, color=_INK, va="bottom")
    figure.subplots_adjust(bottom=0.19 if binning.clipped else 0.155,
                           left=0.11, right=0.98, top=0.91)
    figure.savefig(path, dpi=dpi, facecolor="white")
    plt.close(figure)
    return path
