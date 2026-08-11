"""One figure for a whole sample series: the observed patterns, stacked.

Every trace is rescaled to the same height and lifted above the one below,
so the series reads top to bottom as one picture instead of as six separate
figures. The refined fit, the background and the residual stay in the
per-sample figures the notebook writes; this one answers a different
question, which reflection appears, moves or goes away across the series.

Run it from the repository root, with the exports already in data/::

    python make_series.py

It writes output/pdf/ and output/png/ like the notebook does.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

import xrd_plotter as xp

# --- Settings ----------------------------------------------------------
# The series, bottom trace first. Listed by hand because the order that
# means something is the physical one, composition or temperature or time,
# and no sort of the file names recovers it. A name here that is not in
# DATA_FOLDER stops the run rather than leaving a gap in the series.
SERIES: list[str] = [
    "sample_1.csv",
    "sample_2.csv",
    "sample_3.csv",
]

# Vertical gap between traces, in the normalised units stack() produces,
# where every pattern spans 1.0 from baseline to tallest reflection. Above
# 1.0 leaves clear air between traces, below 1.0 lets them overlap, which
# saves height on a long series at the cost of reading the weak reflections.
OFFSET = 1.15

# One row of reflection ticks per phase, below the bottom trace, in the
# colours the per-sample figures use. Set False for traces alone.
SHOW_TICKS = True

# sqrt(I) on the drawn copy, as in the notebook, so the weak reflections
# survive beside a strong one. It also flattens the difference in contrast
# between samples, which matters more here than in a single figure.
USE_SQRT = True

# The 2theta window shared by every trace. None takes the widest measured
# range in the series. A window fixed here is what makes the traces
# comparable, so prefer setting both.
PLOT_X_MIN: float | None = None
PLOT_X_MAX: float | None = None

DATA_FOLDER = Path("data")
METADATA_FILE = Path("Samples_metadata.csv")
OUTPUT_FOLDER = Path("output")
OUTPUT_BASENAME = "series_XRD_stacked"

# Thinner than a single figure's fit line: six traces at 1.5 pt turn the
# figure into a solid block.
LINEWIDTH_TRACE = 0.9


def stack(patterns: list[np.ndarray],
          offset: float = OFFSET) -> list[np.ndarray]:
    """Rescale every pattern to a 0 to 1 span, then lift each by its index.

    Rescaling is per pattern, so a weakly scattering sample is drawn as tall
    as a strong one and no trace climbs into the one above it. Heights are
    therefore not comparable between traces, and the caption of the figure
    has to say so.

    >>> [t.tolist() for t in stack([np.array([0.0, 2.0]),
    ...                             np.array([1.0, 3.0])], 1.0)]
    [[0.0, 1.0], [1.0, 2.0]]
    >>> stack([np.array([5.0, 5.0])], 1.0)[0].tolist()
    [0.0, 0.0]
    """
    stacked = []
    for index, values in enumerate(patterns):
        low = float(np.nanmin(values))
        span = float(np.nanmax(values)) - low
        # A pattern with no span has nothing to normalise by. It is drawn on
        # its own baseline rather than divided by zero into NaN, which would
        # take the trace out of the figure without saying anything.
        scaled = (values - low) / span if span > 0 else np.zeros_like(values)
        stacked.append(scaled + index * offset)
    return stacked


def load_series(filenames: list[str], data_folder: Path, metadata_file: Path
                ) -> tuple[list[tuple[np.ndarray, np.ndarray, str]],
                           dict[str, np.ndarray]]:
    """Read the series in order: (2theta, observed, label) per file.

    The reflection positions come back separately, from the first file that
    has any, since one row of ticks per phase is drawn for the series rather
    than one per sample. A file that cannot be read stops the run: a series
    figure silently missing a member is worse than no figure.

    >>> load_series([], Path("data"), Path("nowhere.csv"))
    ([], {})
    """
    meta = xp.load_metadata(metadata_file)
    traces: list[tuple[np.ndarray, np.ndarray, str]] = []
    phases: dict[str, np.ndarray] = {}

    for filename in filenames:
        path = Path(data_folder) / filename
        # weighted=False asks the parser for obs and calc rather than for a
        # 'diff/sigma' column. No residual is drawn here, and requiring the
        # weighted one would refuse an export that has everything this
        # figure needs.
        data, phase_cols, error = xp.read_gsas2_csv(path, weighted=False)
        if error is not None:
            raise SystemExit(f"{filename}: {error}")
        theta, obs, _calc, _bkg, _resid, file_phases = xp.prepare_data(
            data, phase_cols, use_sqrt=USE_SQRT)
        name, _pct, _colors, _window = xp.sample_info(meta, filename,
                                                      path.stem)
        traces.append((theta, obs, name))
        if not phases:
            phases = file_phases
    return traces, phases


def plot_series(traces: list[tuple[np.ndarray, np.ndarray, str]],
                phases: dict[str, np.ndarray],
                offset: float = OFFSET) -> Figure:
    """Draw the stacked series on one set of axes; returns the Figure.

    Typical use, after load_series::

        traces, phases = load_series(SERIES, DATA_FOLDER, METADATA_FILE)
        fig = plot_series(traces, phases)
    """
    fig, ax = plt.subplots(figsize=(xp.FIGURE_WIDTH, xp.FIGURE_HEIGHT),
                           dpi=110)
    stacked = stack([obs for _theta, obs, _name in traces], offset)

    # Every trace in the same black. A colour per sample would be a second
    # scale to decode, and the vertical position already says which sample
    # is which. Line, not the scatter of the per-sample figure: at this
    # density the markers of six patterns merge into a block.
    for (theta, _obs, _name), values in zip(traces, stacked):
        ax.plot(theta, values, color=xp.COLOR_OBS, lw=LINEWIDTH_TRACE,
                zorder=2)

    # Labels sit outside the right border, at the baseline of their own
    # trace: the axes fraction places them horizontally, the data
    # coordinate vertically.
    outside_right = blended_transform_factory(ax.transAxes, ax.transData)
    for index, (_theta, _obs, name) in enumerate(traces):
        ax.text(1.01, index * offset, name, transform=outside_right,
                va="bottom", ha="left", fontsize=xp.FONT_SIZE_LEGEND)

    labels = {phase: xp.phase_label(phase) for phase in phases}
    ordered = sorted(phases, key=lambda phase: labels[phase])
    tick_colors = xp.phase_colors([labels[phase] for phase in ordered])
    rows = 0
    if SHOW_TICKS:
        for phase in ordered:
            # A reflection sits at a real angle, so a zero in a phase column
            # is padding and would draw a tick against the left border.
            locations = phases[phase][phases[phase] > 0.1]
            if len(locations):
                row_y = -0.10 * offset - rows * 0.07 * offset
                ax.vlines(locations, row_y, row_y + 0.05 * offset,
                          colors=tick_colors[labels[phase]],
                          lw=xp.LINEWIDTH_TICKS, zorder=3)
                rows += 1
        if rows:
            ax.legend(handles=[Line2D([0], [0], color=tick_colors[label],
                                      lw=3, label=label)
                               for label in dict.fromkeys(labels.values())],
                      loc="upper left", fontsize=xp.FONT_SIZE_LEGEND,
                      framealpha=1, edgecolor="black")

    ax.set_xlabel(r"$2\theta\:/\:^\circ$", fontsize=xp.FONT_SIZE_LABEL)
    ax.set_ylabel(r"$\sqrt{\mathrm{Intensity}}$ / a.u." if USE_SQRT
                  else r"Intensity / a.u.", fontsize=xp.FONT_SIZE_LABEL,
                  labelpad=10)
    widest = xp.plot_window(np.concatenate([t for t, _o, _n in traces]),
                            PLOT_X_MIN, PLOT_X_MAX)
    ax.set_xlim(*widest)
    ax.set_ylim(bottom=(-0.10 * offset - (rows - 0.5) * 0.07 * offset
                        if rows else -0.05 * offset),
                top=(len(traces) - 1) * offset + 1.0 + 0.08 * offset)
    # No numbers on the intensity axis, for the reason the module gives: the
    # unit is arbitrary, and here the normalisation makes a comparison of
    # heights between traces meaningless on top of that.
    ax.tick_params(direction="in", top=False, right=False, left=False,
                   labelleft=False, width=1.5, length=6,
                   labelsize=xp.FONT_SIZE_TICK)
    for spine in ax.spines.values():
        spine.set_linewidth(xp.LINEWIDTH_BORDER)
    return fig


def main() -> None:
    """Read SERIES, draw it and write the two files under OUTPUT_FOLDER."""
    traces, phases = load_series(SERIES, DATA_FOLDER, METADATA_FILE)
    if not traces:
        raise SystemExit("SERIES is empty, nothing to draw")
    fig = plot_series(traces, phases)
    base = xp.save_figure(fig, OUTPUT_FOLDER, OUTPUT_BASENAME)
    print(f"wrote {OUTPUT_FOLDER}/pdf/{base}.pdf "
          f"and {OUTPUT_FOLDER}/png/{base}.png")


if __name__ == "__main__":
    main()
