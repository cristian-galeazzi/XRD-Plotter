"""One figure for a whole sample series: the observed patterns, stacked.

Every trace is rescaled to the same height and lifted above the one below,
so the series reads top to bottom as one picture instead of as six separate
figures. The refined fit, the background and the residual stay in the
per-sample figures the notebook writes; this one answers a different
question, which reflection appears, moves or goes away across the series.

The tick rows and the guides are read from the first file in SERIES that
carries a reflection column, one row of ticks per phase for the whole
series rather than one per sample. In a series drawn precisely because
reflections move, remember that the ticks are anchored to that one file
and not to every trace on top of them. A first file whose reflection
column holds nothing but padding wins that choice all the same and leaves
the series with no ticks at all, so put a well refined pattern first.

Run it from the repository root, with the exports already in data/::

    python make_series.py

It writes output/pdf/ and output/png/ like the notebook does.
"""
from __future__ import annotations

import textwrap
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

# Legend text per file, winning over the 'formula' column of the metadata.
# Empty leaves the behaviour as it is: the metadata first, then the file
# stem. The key is the file name as it appears in SERIES.
SERIES_LABELS: dict[str, str] = {}

# Vertical gap between traces, in the normalised units stack() produces,
# where every pattern spans 1.0 from baseline to tallest reflection. The
# 0.35 above a full trace is the room its label needs: below about 1.25 the
# text lands on the trace above it. Raise it to spread a long series out.
# Past roughly a dozen traces, OFFSET buys label room only by squeezing the
# patterns themselves, since xp.FIGURE_HEIGHT is fixed and does not grow
# with the series: a long series wants that raised too.
OFFSET = 1.35

# Height of a trace label above its own baseline, where 1.0 is the top of
# the trace, since stack() gives every pattern a span of exactly 1.0. Below
# 1.0 the label drops towards the pattern, above it the label rises into the
# gap. Keep it under OFFSET, or the label lands on the trace above.
LABEL_HEIGHT = 0.90

# Where a trace label starts, in degrees 2theta, measured from its left
# edge since the text runs rightwards from there. None instead pins every
# label just inside the left border, whatever window is drawn, which is the
# only default that works without knowing the range:
#   LABEL_X = 14.0
LABEL_X: float | None = None

# Weight of the trace labels: 'normal', 'medium', 'semibold' or 'bold'.
# A heavier label separates itself from the pattern it sits over without a
# box behind it. STIXGeneral, the face these figures are set in, carries a
# real bold, so nothing is synthesised from the upright.
LABEL_WEIGHT = "normal"

# One row of reflection ticks per phase, below the bottom trace, in the
# colours the per-sample figures use. Set False for traces alone.
SHOW_TICKS = True

# Height of one reflection tick, in units of OFFSET. Taller ticks are easier
# to follow across a wide pattern. The rows sit one TICK_HEIGHT + 0.02
# apart, so raising this moves them apart instead of overlapping them.
TICK_HEIGHT = 0.10

# 2theta of the reflections to follow up through the stack, as dotted lines
# behind the traces. Empty draws none, which is the default. Each value
# snaps to the nearest reflection position, so a value read off the PDF by
# eye still lands exactly on its tick:
#   GUIDE_LINES = [30.9, 37.0]
GUIDE_LINES: list[float] = []

# How far a GUIDE_LINES value may sit from a reflection and still snap to
# it, in degrees. A value with nothing this close is reported and not drawn.
GUIDE_SNAP = 0.3

# sqrt(I) on the drawn copy, as in the notebook, so the weak reflections
# survive beside a strong one. It also flattens the difference in contrast
# between samples, which matters more here than in a single figure.
USE_SQRT = True

# The 2theta window shared by every trace. None takes the widest measured
# range in the series, as long as xp.PLOT_X_MIN/xp.PLOT_X_MAX are unset too,
# since xp.plot_window falls back to those before the data. A window fixed
# here is what makes the traces comparable, so prefer setting both.
PLOT_X_MIN: float | None = None
PLOT_X_MAX: float | None = None

DATA_FOLDER = Path("data")
METADATA_FILE = Path("Samples_metadata.csv")
OUTPUT_FOLDER = Path("output")
OUTPUT_BASENAME = "series_XRD_stacked"

# Thinner than a single figure's fit line: six traces at 1.5 pt turn the
# figure into a solid block.
LINEWIDTH_TRACE = 0.9


def stack(patterns: list[np.ndarray], offset: float = OFFSET,
          scopes: list[np.ndarray] | None = None) -> list[np.ndarray]:
    """Rescale every pattern to a 0 to 1 span, then lift each by its index.

    Rescaling is per pattern, so a weakly scattering sample is drawn as tall
    as a strong one and no trace climbs into the one above it. Heights are
    therefore not comparable between traces, and the caption of the figure
    has to say so.

    'scopes' names, per pattern, the points that set the 0 to 1 range: a
    boolean mask matching the 2theta window actually plotted, so a
    reflection outside that window does not set the span of a trace whose
    only visible part is a small one. Omitted, the whole pattern sets it, as
    before. A point outside the scope can land above 1.0 or below 0.0; that
    is correct, since set_xlim crops it before it is drawn.

    >>> [t.tolist() for t in stack([np.array([0.0, 2.0]),
    ...                             np.array([1.0, 3.0])], 1.0)]
    [[0.0, 1.0], [1.0, 2.0]]
    >>> stack([np.array([5.0, 5.0])], 1.0)[0].tolist()
    [0.0, 0.0]
    """
    stacked = []
    for index, values in enumerate(patterns):
        scope = values[scopes[index]] if scopes is not None else values
        low = float(np.nanmin(scope))
        span = float(np.nanmax(scope)) - low
        # A pattern with no span has nothing to normalise by. It is drawn on
        # its own baseline rather than divided by zero into NaN, which would
        # take the trace out of the figure without saying anything.
        scaled = (values - low) / span if span > 0 else np.zeros_like(values)
        stacked.append(scaled + index * offset)
    return stacked


def snap_to_reflection(value: float, positions: np.ndarray,
                       tolerance: float = GUIDE_SNAP) -> float | None:
    """The reflection position nearest to value, or None past the tolerance.

    A guide line is worth drawing only where a reflection actually is, so a
    value with nothing near it returns None for the caller to report rather
    than being drawn at the typed position.

    >>> snap_to_reflection(31.0, np.array([20.0, 30.95, 37.04]), 0.3)
    30.95
    >>> snap_to_reflection(31.0, np.array([20.0, 37.04]), 0.3) is None
    True
    """
    if not len(positions):
        return None
    nearest = float(positions[np.argmin(np.abs(positions - value))])
    return nearest if abs(nearest - value) <= tolerance else None


def snap_to_phase(value: float, rows: dict[str, np.ndarray],
                  tolerance: float = GUIDE_SNAP
                  ) -> tuple[float | None, str | None]:
    """The nearest reflection across the phases, and the phase it belongs to.

    Snapping per phase rather than against one merged list is what lets a
    guide be drawn in the colour of its own tick row. Two phases equally
    close to the same value are a real ambiguity, and the first in the
    order the rows were given wins it, so the figure is at least the same
    on every run.

    >>> rows = {"Phase 1": np.array([30.95]), "Phase 2": np.array([37.04])}
    >>> snap_to_phase(37.0, rows, 0.3)
    (37.04, 'Phase 2')
    >>> snap_to_phase(27.0, rows, 0.3)
    (None, None)
    """
    best: tuple[float | None, str | None] = (None, None)
    smallest = tolerance
    for label, positions in rows.items():
        snapped = snap_to_reflection(value, positions, tolerance)
        if snapped is not None and abs(snapped - value) <= smallest:
            if best[0] is None or abs(snapped - value) < smallest:
                best, smallest = (snapped, label), abs(snapped - value)
    return best


def reflections_in_window(positions: np.ndarray,
                          window: tuple[float, float] | None = None
                          ) -> np.ndarray:
    """The drawable reflection positions of one phase, inside the window.

    A reflection sits at a real angle, so a zero left in a phase column is
    padding rather than a position and is dropped. Outside the plotted
    window a position is neither drawn nor worth listing, so the window
    filters it out too. This is the one place that decides what counts as a
    reflection, so the ticks, the guide snapping and the printed list can
    never disagree.

    >>> reflections_in_window(np.array([0.0, 30.95, 90.0]), (13.0, 85.0))
    array([30.95])
    >>> reflections_in_window(np.array([0.0, 30.95, 90.0])).tolist()
    [30.95, 90.0]
    """
    drawable = positions[positions > 0.1]
    if window is None:
        return drawable
    low, high = window
    return drawable[(drawable >= low) & (drawable <= high)]


def format_reflections(rows: dict[str, np.ndarray]) -> str:
    """One line per phase, its reflections ready to paste into GUIDE_LINES.

    >>> print(format_reflections({"Phase 1": np.array([30.95, 37.04])}))
      Phase 1: 30.95, 37.04
    >>> format_reflections({})
    ''
    """
    lines = []
    for label, positions in rows.items():
        values = ", ".join(f"{p:.2f}" for p in positions)
        # Wrapped rather than run off the terminal: a phase can carry
        # dozens of reflections and the point of printing them is that they
        # can be read and picked from.
        lines.extend(textwrap.wrap(f"{label}: {values}", width=76,
                                   initial_indent="  ",
                                   subsequent_indent="    "))
    return "\n".join(lines)


def load_series(filenames: list[str], data_folder: Path, metadata_file: Path,
                labels: dict[str, str] | None = None
                ) -> tuple[list[tuple[np.ndarray, np.ndarray, str]],
                           dict[str, np.ndarray], dict[str, str]]:
    """Read the series in order: (2theta, observed, label) per file.

    The label is what SERIES_LABELS says for that file name, or the
    'formula' cell of its metadata row, or the file stem, first hit
    winning. 'labels' overrides the SERIES_LABELS constant for this call.

    The reflection positions and the phase colour overrides come back
    separately, from the first file that carries a reflection column, since
    one row of ticks per phase is drawn for the series rather than one per
    sample, and its colours are what the tick row is drawn in. That file
    wins the choice on having the column, not on the column holding a
    position worth drawing. A file that cannot be read stops the run: a
    series figure silently missing a member is worse than no figure.

    >>> load_series([], Path("data"), Path("nowhere.csv"))
    ([], {}, {})
    """
    overrides = SERIES_LABELS if labels is None else labels
    meta = xp.load_metadata(metadata_file)
    traces: list[tuple[np.ndarray, np.ndarray, str]] = []
    phases: dict[str, np.ndarray] = {}
    colors: dict[str, str] = {}

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
        name, _pct, file_colors, _window = xp.sample_info(meta, filename,
                                                          path.stem)
        # get with a default rather than 'or': a label deliberately set to
        # an empty string stays empty instead of falling back.
        traces.append((theta, obs, overrides.get(filename, name)))
        if not phases:
            phases = file_phases
            colors = file_colors
    return traces, phases, colors


def plot_series(traces: list[tuple[np.ndarray, np.ndarray, str]],
                phases: dict[str, np.ndarray],
                offset: float = OFFSET,
                colors: dict[str, str] | None = None) -> Figure:
    """Draw the stacked series on one set of axes; returns the Figure.

    Typical use, after load_series::

        traces, phases, colors = load_series(SERIES, DATA_FOLDER,
                                             METADATA_FILE)
        fig = plot_series(traces, phases, colors=colors)
    """
    fig, ax = plt.subplots(figsize=(xp.FIGURE_WIDTH, xp.FIGURE_HEIGHT),
                           dpi=110)
    # The window is set before stack() runs, so the 0 to 1 span of a trace
    # is set by the part of the pattern the window actually shows, as
    # xp.create_plot already does. A reflection outside the window would
    # otherwise set the span of a trace whose visible part is much smaller,
    # sinking its label and making two samples with different off-window
    # peaks scale by different amounts, exactly the comparison the
    # normalisation exists to prevent.
    widest = xp.plot_window(np.concatenate([t for t, _o, _n in traces]),
                            PLOT_X_MIN, PLOT_X_MAX)
    x_low, x_high = widest
    scopes = []
    for theta, _obs, _name in traces:
        visible = (theta >= x_low) & (theta <= x_high)
        # A trace with nothing inside the window keeps its own full range
        # rather than dividing by an empty slice.
        scopes.append(visible if visible.any() else np.ones(len(theta), bool))
    stacked = stack([obs for _theta, obs, _name in traces], offset, scopes)

    # Every trace in the same black. A colour per sample would be a second
    # scale to decode, and the vertical position already says which sample
    # is which. Line, not the scatter of the per-sample figure: at this
    # density the markers of six patterns merge into a block.
    #
    # Sliced to the scope, not the whole trace: a point outside the window
    # can land far above 1.0 once normalised, and set_xlim discards it
    # anyway. Plotting it regardless would draw the segment crossing the
    # window edge up to that off-scale height, a spike cutting across every
    # trace above it.
    for (theta, _obs, _name), values, scope in zip(traces, stacked, scopes):
        ax.plot(theta[scope], values[scope], color=xp.COLOR_OBS,
                lw=LINEWIDTH_TRACE, zorder=2)

    # Labels sit inside the frame, above their own trace. The y is always a
    # data coordinate, so it follows the trace: stack() normalises every
    # pattern to a span of 1.0, so index * offset + 1.0 is the top of that
    # trace and LABEL_HEIGHT moves the text down towards the pattern or up
    # into the gap. The x is a real 2theta when LABEL_X names one, and an
    # axes fraction otherwise, so the shipped default needs no window.
    if LABEL_X is None:
        placement = blended_transform_factory(ax.transAxes, ax.transData)
        label_x = 0.02
    else:
        placement = ax.transData
        label_x = LABEL_X
    for index, (_theta, _obs, name) in enumerate(traces):
        ax.text(label_x, index * offset + LABEL_HEIGHT, name,
                transform=placement, va="bottom", ha="left",
                fontsize=xp.FONT_SIZE_LEGEND, fontweight=LABEL_WEIGHT)

    labels = {phase: xp.phase_label(phase) for phase in phases}
    ordered = sorted(phases, key=lambda phase: labels[phase])
    tick_colors = xp.phase_colors([labels[phase] for phase in ordered],
                                  colors)
    rows = 0
    pitch = (TICK_HEIGHT + 0.02) * offset  # 0.02 of clearance between rows
    # The tick block hangs below the bottom trace by 0.05 of a trace height,
    # whatever TICK_HEIGHT is. Derived from it rather than written out again,
    # so a taller tick moves the block down instead of growing into the
    # pattern above it.
    block_top = -(TICK_HEIGHT + 0.05) * offset
    listed: dict[str, np.ndarray] = {}
    if SHOW_TICKS:
        for phase in ordered:
            locations = reflections_in_window(phases[phase], widest)
            if len(locations):
                row_y = block_top - rows * pitch
                ax.vlines(locations, row_y, row_y + TICK_HEIGHT * offset,
                          colors=tick_colors[labels[phase]],
                          lw=xp.LINEWIDTH_TICKS, zorder=3)
                # Two columns can carry the same phase name, and the legend
                # gives them one entry between them. Their positions join
                # instead of the second replacing the first, or the list and
                # the guides would know about half the phase.
                listed[labels[phase]] = (
                    np.concatenate([listed[labels[phase]], locations])
                    if labels[phase] in listed else locations)
                rows += 1

        # Printed so a GUIDE_LINES value is picked from the reflections that
        # are actually there rather than read off the drawn figure by eye.
        if listed:
            print(f"reflections between {widest[0]:g} and {widest[1]:g} deg, "
                  "for GUIDE_LINES:")
            print(format_reflections(listed))

        # Guides span the full height of the axes, behind the black traces,
        # so a tick can be followed up to the peak it belongs to. axvline
        # works in axes fractions vertically, so it does not depend on the
        # limits being set yet. Nested under SHOW_TICKS: a guide exists only
        # to connect a tick to a peak, so with the ticks off there is
        # nothing left for it to point at.
        for value in GUIDE_LINES:
            snapped, phase_label = snap_to_phase(value, listed, GUIDE_SNAP)
            if snapped is None:
                print(f"  ! no reflection within {GUIDE_SNAP:g} deg of "
                      f"{value:g}, no guide drawn")
                continue
            # Coloured like the tick it came from, so the guide says which
            # phase owns the peak it points at as well as where it is.
            ax.axvline(snapped, color=tick_colors[phase_label], ls=":",
                       lw=1.2, zorder=4)
        if rows:
            # Built from 'ordered', the order the tick rows are drawn in,
            # not from the column order of the export: a reader pairs a
            # legend entry to a row by position, and with three or more
            # phases the two orders need not agree.
            seen: set[str] = set()
            handles = []
            for phase in ordered:
                label = labels[phase]
                if label in seen:
                    continue
                seen.add(label)
                handles.append(Line2D([0], [0], color=tick_colors[label],
                                      lw=3, label=label))
            # Opposite corner from the trace labels, which start at the
            # left. A short handle: the colour is the whole message, and a
            # long rule beside a short phase name reads as a second scale.
            ax.legend(handles=handles, loc="upper right",
                      fontsize=xp.FONT_SIZE_LEGEND, framealpha=1,
                      handlelength=1, handletextpad=1)

    ax.set_xlabel(r"$2\theta\:/\:^\circ$", fontsize=xp.FONT_SIZE_LABEL)
    ax.set_ylabel(r"$\sqrt{\mathrm{Intensity}}$ / a.u." if USE_SQRT
                  else r"Intensity / a.u.", fontsize=xp.FONT_SIZE_LABEL,
                  labelpad=10)
    ax.set_xlim(*widest)
    ax.set_ylim(bottom=(block_top - (rows - 0.5) * pitch
                        if rows else -0.05 * offset),
                # 0.08 was clearance for a curve. A line of text needs the
                # same room the labels get between traces, or the topmost
                # one is cut by the frame.
                top=(len(traces) - 1) * offset + 1.0 + 0.30 * offset)
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
    traces, phases, colors = load_series(SERIES, DATA_FOLDER, METADATA_FILE)
    if not traces:
        raise SystemExit("SERIES is empty, nothing to draw")
    fig = plot_series(traces, phases, colors=colors)
    base = xp.save_figure(fig, OUTPUT_FOLDER, OUTPUT_BASENAME)
    print(f"wrote {OUTPUT_FOLDER}/pdf/{base}.pdf "
          f"and {OUTPUT_FOLDER}/png/{base}.png")


if __name__ == "__main__":
    main()
