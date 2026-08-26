"""The read-length distribution as a figure.

Draws the binning the terminal histogram draws, from the same
:class:`seqviewer.lengths.Binning`, so the figure and the terminal agree on
where the axis falls and which reads are outside it.

The styling follows the conventions of a journal figure: a sans-serif face at
small sizes, thin rules, ticks outward, no gridlines, a full box, and the panel
sized to the width a two-column figure is printed at.  The file's suffix chooses
the format, so a PDF or an SVG keeps its text as text and can be edited in a
drawing program.

matplotlib is imported inside :func:`write_png` rather than at the top of the
module, so that importing :mod:`seqviewer` does not require it.  It is the only
part of the package that needs anything beyond the standard library to produce
output, and it is reached only when a figure is asked for.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .lengths import Binning, Summary

__all__ = ["DPI", "EDGE_UNTIL", "FIGSIZE", "PngUnavailable",
           "write_png"]

#: Dots per inch.  300 is the floor most journals accept for a raster figure.
DPI = 300

#: Figure size in inches.  183 mm is the width of a two-column figure, and the
#: height leaves the panel about twice as wide as it is tall.
FIGSIZE = (7.2, 3.6)

#: Bins up to which the bars are separated by a hairline.
EDGE_UNTIL = 80

#: The bars, tab10 blue.  One colour rather than a scale: a histogram of one
#: variable carries no second quantity for colour to stand for.
_BAR = "#1f77b4"

_INK = "#1a1a1a"
_NOTE = "#595959"
_FAINT = "#8a8a8a"

#: Rules the figure is drawn under.  Applied through ``rc_context`` so a caller
#: that also uses matplotlib keeps its own settings.
_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica",
                        "Liberation Sans", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "axes.edgecolor": _INK,
    "axes.labelcolor": _INK,
    # All four spines, so the panel is closed on every side.
    "axes.spines.top": True,
    "axes.spines.right": True,
    "text.color": _INK,
    "xtick.color": _INK,
    "ytick.color": _INK,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
    # Text stays text in a PDF or an EPS, rather than being drawn as outlines,
    # so a figure can still be edited after it is written.
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


class PngUnavailable(RuntimeError):
    """Raised when matplotlib is not installed."""


def _figures(summary: Summary) -> str:
    """Return the one line of figures that sits under the panel."""
    return (f"n = {summary.reads:,} reads · median {summary.median:,} bp · "
            f"mean {summary.mean:,.0f} bp · N50 {summary.n50:,} bp")


def _clipped_note(binning: Binning, summary: Summary) -> Optional[str]:
    """Return the line naming the reads the axis leaves out, or None.

    Without it the panel reads as the whole distribution, when the tails have
    been left off it deliberately.
    """
    if not binning.clipped:
        return None
    inside = summary.reads - binning.below - binning.above
    share = 100 * inside / summary.reads if summary.reads else 0.0
    parts = []
    if binning.below:
        parts.append(f"{binning.below:,} shorter than {binning.low:,} bp")
    if binning.above:
        parts.append(f"{binning.above:,} longer than {binning.high:,} bp")
    return (f"Axis covers {share:.1f}% of reads; " + " and ".join(parts)
            + f" not shown ({summary.shortest:,}–{summary.longest:,} bp overall).")


def write_png(
    path,
    binning: Binning,
    summary: Summary,
    title: Optional[str] = None,
    log: bool = False,
    dpi: int = DPI,
) -> Path:
    """Draw *binning* to *path* and return the path written.

    Bars are placed at the middle of each bin and are as wide as the bin, so the
    length axis reads linearly whether or not *log* scales the counts.  Reads
    outside the axis are named under the panel rather than drawn, since they have
    no width on it.

    A path with no suffix gains ``.png``.  Matplotlib would otherwise append the
    extension itself, which turns a directory such as ``~/Downloads/`` into a
    file beside it.

    Raises :class:`PngUnavailable` where matplotlib is not installed, and
    ``IsADirectoryError`` where *path* names a directory.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")           # a file, so no display is needed
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
    except ImportError as exc:          # pragma: no cover - matplotlib is opt-in
        raise PngUnavailable(
            "drawing a figure needs matplotlib: pip install 'seqviewer[plot]'"
        ) from exc

    path = Path(path).expanduser()
    if path.is_dir():
        raise IsADirectoryError(
            f"{path} is a directory; give the name of the file to write")
    if not path.suffix:
        path = path.with_name(path.name + ".png")
    path.parent.mkdir(parents=True, exist_ok=True)

    note = _clipped_note(binning, summary)
    with plt.rc_context(_STYLE):
        figure, axes = plt.subplots(figsize=FIGSIZE)

        if binning.bins:
            width = binning.bins[0].high - binning.bins[0].low
            middles = [b.low + (b.high - b.low) / 2 for b in binning.bins]
            # A hairline between bars separates them while they are wide enough
            # to read one at a time.  Past that it is a stripe over each bar, so
            # the bars are left to meet and the distribution reads as one shape.
            edge = 0.3 if len(binning.bins) <= EDGE_UNTIL else 0.0
            axes.bar(middles, [b.count for b in binning.bins], width=width,
                     color=_BAR, edgecolor="white", linewidth=edge)
            axes.set_xlim(binning.bins[0].low - width * 0.5,
                          binning.bins[-1].high + width * 0.5)
        if log:
            axes.set_yscale("log")

        axes.set_xlabel("Read length (bp)")
        axes.set_ylabel("Reads")
        if title:
            axes.set_title(title, loc="left", pad=6)

        thousands = FuncFormatter(lambda v, _: f"{int(v):,}")
        axes.xaxis.set_major_formatter(thousands)
        if not log:
            axes.yaxis.set_major_formatter(thousands)

        # The figures go under the panel, where they do not sit over a bar.
        figure.text(0.012, 0.075 if note else 0.045, _figures(summary),
                    fontsize=6.5, color=_NOTE, va="bottom")
        if note:
            figure.text(0.012, 0.018, note, fontsize=6, color=_FAINT,
                        va="bottom")
        figure.subplots_adjust(left=0.085, right=0.985, top=0.90,
                               bottom=0.245 if note else 0.20)
        figure.savefig(path, dpi=dpi, facecolor="white")
        plt.close(figure)
    return path
