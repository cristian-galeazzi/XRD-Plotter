"""One figure for a whole sample series: the observed patterns, stacked.

Every trace is rescaled to the same height and lifted above the one below,
so the series reads top to bottom as one picture instead of as several
separate figures. The refined fit, the background and the residual stay in
the per-sample figures the notebook writes; this one answers a different
question, which reflection appears, moves or goes away across the series.

The tick rows and the guides are read from the first file of the series
that carries a reflection column, one row of ticks per phase for the whole
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
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

import xrd_plotter as xp

# --- Settings ----------------------------------------------------------
# Which samples the series holds, in which order, and what each trace is
# called all come from Samples_metadata.csv, the private file beside the
# notebook that .gitignore excludes: a row with a number in its
# 'series_order' column joins the series at that position, and its
# 'series_label' cell names its trace. See docs/metadata.md.
#
# Everything below is the look of the figure. Section 5 of the notebook
# offers the same settings as controls, so a figure can be tuned there
# without editing this file; what you set here is what a plain
# 'python make_series.py' run draws, and what the notebook opens on.
# Each one is also an argument of plot_series, so nothing has to be
# written to this module to change one drawing.
#
# Two units appear throughout, and they are worth keeping apart:
#   degrees 2theta   a position on the horizontal axis
#   trace heights    stack() rescales every pattern to span exactly 1.0
#                    from its own baseline to its tallest reflection, so
#                    0.5 is half a pattern tall, wherever it sits
# Widths are in points, matplotlib's own unit for a line.

# --- The traces ---------------------------------------------------------
# Vertical distance from one trace's baseline to the next, in trace
# heights. At 1.0 the top of a pattern touches the baseline above it, so
# anything larger leaves a gap and anything smaller overlaps them on
# purpose. It also has to clear LABEL_HEIGHT, the room the label below
# needs, or that text lands on the trace above. Raise it to spread a short
# series out. Past roughly a dozen traces it buys label room only by
# squeezing the patterns themselves, since xp.FIGURE_HEIGHT does not grow
# with the series: a long series wants that raised too.
OFFSET = 1.35

# Width of a trace, in points. Thinner than a single figure's fit line on
# purpose: a whole series at 1.5 pt turns the figure into a solid block.
LINEWIDTH_TRACE = 0.9

# sqrt(I) on the drawn copy, so the weak reflections survive beside a
# strong one. It also flattens the difference in contrast between samples,
# which matters more here than in a single figure. False draws the
# intensity as it is. The axis label follows, and so does the output file
# name, so the two versions never overwrite each other.
USE_SQRT = True

# The 2theta window every trace is drawn in, in degrees. None at either end
# takes the widest measured range in the series, as long as
# xp.PLOT_X_MIN/xp.PLOT_X_MAX are unset too, since xp.plot_window falls
# back to those before the data. A window fixed here is part of what makes
# the traces comparable, so prefer setting both. It also sets what each
# trace is rescaled by, since a pattern is normalised over the part of it
# the window shows.
PLOT_X_MIN: float | None = None
PLOT_X_MAX: float | None = None

# --- The label on each trace --------------------------------------------
# Height of a trace's label above that trace's own baseline, in trace
# heights, so 1.0 sits exactly at the top of the pattern. Below 1.0 the
# label drops towards the pattern, above it the label rises into the gap.
# Keep it under OFFSET, or the label lands on the trace above.
LABEL_HEIGHT = 0.90

# Where a label starts, in degrees 2theta, measured from its left edge
# since the text runs rightwards from there. None instead pins every label
# just inside the left border, whatever window is drawn, which is the only
# default that works without knowing the range:
#   LABEL_X = 14.0
LABEL_X: float | None = None

# Weight of the trace labels: 'normal', 'medium', 'semibold' or 'bold'.
# A heavier label separates itself from the pattern it sits over without a
# box behind it. STIXGeneral, the face these figures are set in, carries a
# real bold, so nothing is synthesised from the upright.
LABEL_WEIGHT = "normal"

# --- The reflection ticks -----------------------------------------------
# One row of ticks per phase, below the bottom trace, in the colours the
# per-sample figures use. False draws the traces alone, and also silences
# the printed reflection list, the tool for picking guide positions, since
# with no ticks drawn there is nothing to list.
SHOW_TICKS = True

# Height of one tick, in trace heights, so 0.10 is a tenth of a pattern.
# Taller ticks are easier to follow across a wide pattern. The rows sit one
# TICK_HEIGHT + 0.02 apart, so raising this moves them apart instead of
# overlapping them.
TICK_HEIGHT = 0.10

# --- The guides -----------------------------------------------------
# 2theta of the reflections to follow up through the stack, as vertical
# lines drawn in front of the traces, so a guide is not hidden by the dense
# stack of patterns it is meant to be followed through. Empty draws none,
# which is the default. Each value snaps to the nearest reflection
# position, so a value read off the PDF by eye still lands exactly on its
# tick, and each guide is drawn in the colour of the phase that owns the
# reflection it landed on:
#   GUIDE_LINES = [30.9, 37.0]
#
# The one setting in this file that carries data rather than a preference: a
# reflection position is a measurement, and a few of them identify the
# phase. Empty it before committing, or set the guides in section 5 of the
# notebook, where they stay in the widget. See docs/privacy.md.
GUIDE_LINES: list[float] = []

# How far a GUIDE_LINES value may sit from a reflection and still snap to
# it, in degrees. A value with nothing this close is reported and not
# drawn, rather than being drawn at the typed position where no reflection
# is. Widen it for a value read off a printed figure, narrow it where two
# reflections sit close together and the wrong one keeps winning.
GUIDE_SNAP = 0.3

# Line style of the guides, in matplotlib's notation: ':' dotted, '--'
# dashed, '-.' dash-dot, '-' solid. Dotted is the default because a guide
# crosses a figure already full of vertical lines, and a dotted line reads
# as an annotation rather than as one more reflection. Dashed carries
# further across a tall stack; solid is worth it only for one guide alone.
GUIDE_STYLE = ":"

# Width of a guide, in points. Raise it when a guide disappears into a
# dense stack. Past the width of the traces it crosses, a guide starts
# hiding the peaks it is there to point at.
GUIDE_WIDTH = 1.2

# --- Where the files are ------------------------------------------------
DATA_FOLDER = Path("data")
METADATA_FILE = Path("Samples_metadata.csv")
OUTPUT_FOLDER = Path("output")
OUTPUT_BASENAME = "series_XRD_stacked"


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
        if snapped is not None and (best[0] is None
                                    or abs(snapped - value) < smallest):
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


def series_from_metadata(meta: pd.DataFrame) -> list[tuple[str, str | None]]:
    """The series as (file name, label) pairs, in the order given.

    A row joins the series by carrying a number in 'series_order', and the
    number is the position rather than the row's place in the file, so the
    series is reordered in a spreadsheet without touching any code. A row
    with no number is left out in silence, since a blank cell is a choice;
    a row with something that is not a number is left out and says so. The
    label is the 'series_label' cell, or None to let load_series fall back
    to the formula and then the file name.

    >>> meta = pd.DataFrame({"series_order": ["2", "1"],
    ...                      "series_label": ["B", ""]},
    ...                     index=pd.Index(["b.csv", "a.csv"],
    ...                                    name="filename"))
    >>> series_from_metadata(meta)
    [('a.csv', None), ('b.csv', 'B')]
    >>> series_from_metadata(pd.DataFrame())
    []
    """
    if "series_order" not in getattr(meta, "columns", []):
        return []
    ordered = []
    for filename in meta.index:
        row = meta.loc[filename]
        cell = row.get("series_order")
        # pd.isna as well as the string test: load_metadata reads every
        # cell as a string, and an empty one arrives as NaN rather than as
        # '', which to_number would then report as unreadable.
        if cell is None or pd.isna(cell) or not str(cell).strip():
            continue
        position = xp.to_number(cell)
        if position is None:
            # Reported rather than dropped in silence: a typo in one cell
            # would otherwise remove a sample from the figure invisibly.
            print(f"  ! {filename}: 'series_order' is not a number "
                  f"({cell!r}), not drawn in the series")
            continue
        label = row.get("series_label")
        if not isinstance(label, str) or not label.strip():
            label = None
        ordered.append((position, str(filename), label))
    # The file name breaks a tie, so two rows sharing a position give the
    # same figure on every run instead of following the file's row order.
    ordered.sort(key=lambda row: (row[0], row[1]))
    return [(name, label) for _position, name, label in ordered]


def load_series(entries: list[tuple[str, str | None]], data_folder: Path,
                meta: pd.DataFrame, use_sqrt: bool | None = None
                ) -> tuple[list[tuple[np.ndarray, np.ndarray, str]],
                           dict[str, np.ndarray], dict[str, str]]:
    """Read the series in order: (2theta, observed, label) per file.

    'entries' is what series_from_metadata returns, the file names in the
    order they are drawn with the label each carries. A label of None falls
    back to the 'formula' cell of that sample's metadata row, and then to
    the file name without its extension. The metadata arrives already
    loaded, so a run reads it once.

    The square root is applied here rather than in plot_series, because it
    transforms the data and not the drawing. A caller offering it as a
    control has to pass the same value to both, or the axis would name a
    transform nobody applied. None takes USE_SQRT.

    The reflection positions and the phase colour overrides come back
    separately, from the first file that carries a reflection column, since
    one row of ticks per phase is drawn for the series rather than one per
    sample, and its colours are what the tick row is drawn in. That file
    wins the choice on having the column, not on the column holding a
    position worth drawing. A file that cannot be read stops the run: a
    series figure silently missing a member is worse than no figure.

    >>> load_series([], Path("data"), pd.DataFrame())
    ([], {}, {})
    """
    use_sqrt = USE_SQRT if use_sqrt is None else use_sqrt
    traces: list[tuple[np.ndarray, np.ndarray, str]] = []
    phases: dict[str, np.ndarray] = {}
    colors: dict[str, str] = {}

    for filename, override in entries:
        path = Path(data_folder) / filename
        # weighted=False asks the parser for obs and calc rather than for a
        # 'diff/sigma' column. No residual is drawn here, and requiring the
        # weighted one would refuse an export that has everything this
        # figure needs.
        data, phase_cols, error = xp.read_gsas2_csv(path, weighted=False)
        if error is not None:
            raise SystemExit(f"{filename}: {error}")
        theta, obs, _calc, _bkg, _resid, file_phases = xp.prepare_data(
            data, phase_cols, use_sqrt=use_sqrt)
        name, _pct, file_colors, _window = xp.sample_info(meta, filename,
                                                          path.stem)
        traces.append((theta, obs,
                       override if override is not None else name))
        if not phases:
            phases = file_phases
            colors = file_colors
    return traces, phases, colors


def plot_series(traces: list[tuple[np.ndarray, np.ndarray, str]],
                phases: dict[str, np.ndarray], *,
                colors: dict[str, str] | None = None,
                offset: float | None = None,
                label_height: float | None = None,
                label_x: float | None = None,
                label_weight: str | None = None,
                tick_height: float | None = None,
                show_ticks: bool | None = None,
                guide_lines: list[float] | None = None,
                guide_snap: float | None = None,
                guide_style: str | None = None,
                guide_width: float | None = None,
                linewidth: float | None = None,
                use_sqrt: bool | None = None,
                window: tuple[float | None, float | None] | None = None
                ) -> Figure:
    """Draw the stacked series on one set of axes; returns the Figure.

    Every appearance setting can be given per call, and each one left as
    None takes the module setting of the same meaning, so the notebook can
    drive one redraw without writing to module state. 'label_x' is the one
    exception: None means take LABEL_X, and a caller that wants the labels
    pinned to the left border passes float('nan').

    'use_sqrt' only names the axis. The transform itself belongs to
    load_series, which reads the file, so a caller offering it as a control
    passes the same value to both.

    Typical use, after load_series::

        meta = xp.load_metadata(METADATA_FILE)
        traces, phases, colors = load_series(series_from_metadata(meta),
                                             DATA_FOLDER, meta)
        fig = plot_series(traces, phases, colors=colors)
    """
    # One place where an argument falls back to its setting, so the rest of
    # the body reads one name per value and never both.
    offset = OFFSET if offset is None else offset
    label_height = LABEL_HEIGHT if label_height is None else label_height
    label_weight = LABEL_WEIGHT if label_weight is None else label_weight
    tick_height = TICK_HEIGHT if tick_height is None else tick_height
    show_ticks = SHOW_TICKS if show_ticks is None else show_ticks
    guide_lines = GUIDE_LINES if guide_lines is None else guide_lines
    guide_snap = GUIDE_SNAP if guide_snap is None else guide_snap
    guide_style = GUIDE_STYLE if guide_style is None else guide_style
    guide_width = GUIDE_WIDTH if guide_width is None else guide_width
    linewidth = LINEWIDTH_TRACE if linewidth is None else linewidth
    use_sqrt = USE_SQRT if use_sqrt is None else use_sqrt
    window = (PLOT_X_MIN, PLOT_X_MAX) if window is None else window
    # NaN is the caller's way of asking for the left border, since None
    # already means 'take the setting'. It is never a 2theta.
    if label_x is None:
        label_x = LABEL_X
    left_border = label_x is None or label_x != label_x

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
                            *window)
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
    # density the markers of a full series merge into a block.
    #
    # Sliced to the scope, not the whole trace: a point outside the window
    # can land far above 1.0 once normalised, and set_xlim discards it
    # anyway. Plotting it regardless would draw the segment crossing the
    # window edge up to that off-scale height, a spike cutting across every
    # trace above it.
    for (theta, _obs, _name), values, scope in zip(traces, stacked, scopes):
        ax.plot(theta[scope], values[scope], color=xp.COLOR_OBS,
                lw=linewidth, zorder=2)

    # Labels sit inside the frame, above their own trace. The y is always a
    # data coordinate, so it follows the trace: stack() normalises every
    # pattern to a span of 1.0, so index * offset + 1.0 is the top of that
    # trace and label_height moves the text down towards the pattern or up
    # into the gap. The x is a real 2theta when one is given, and an axes
    # fraction otherwise, so the shipped default needs no window.
    if left_border:
        placement = blended_transform_factory(ax.transAxes, ax.transData)
        label_x = 0.02
    else:
        placement = ax.transData
    for index, (_theta, _obs, name) in enumerate(traces):
        ax.text(label_x, index * offset + label_height, name,
                transform=placement, va="bottom", ha="left",
                fontsize=xp.FONT_SIZE_LEGEND, fontweight=label_weight)

    labels = {phase: xp.phase_label(phase) for phase in phases}
    ordered = sorted(phases, key=lambda phase: labels[phase])
    tick_colors = xp.phase_colors([labels[phase] for phase in ordered],
                                  colors)
    rows = 0
    pitch = (tick_height + 0.02) * offset  # 0.02 of clearance between rows
    # The tick block hangs below the bottom trace by 0.05 of a trace height,
    # whatever tick_height is. Derived from it rather than written out again,
    # so a taller tick moves the block down instead of growing into the
    # pattern above it.
    block_top = -(tick_height + 0.05) * offset
    listed: dict[str, np.ndarray] = {}
    if show_ticks:
        for phase in ordered:
            locations = reflections_in_window(phases[phase], widest)
            if len(locations):
                row_y = block_top - rows * pitch
                ax.vlines(locations, row_y, row_y + tick_height * offset,
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

        # Printed so a guide value is picked from the reflections that are
        # actually there rather than read off the drawn figure by eye.
        if listed:
            print(f"reflections between {widest[0]:g} and {widest[1]:g} deg, "
                  "to pick the guides from:")
            print(format_reflections(listed))

        # Guides span the full height of the axes, in front of the black
        # traces, drawn at zorder=4 above their zorder=2, so a guide behind
        # a dense stack of patterns is not hidden by them and a tick can
        # still be followed up to the peak it belongs to. axvline works in
        # axes fractions vertically, so it does not depend on the limits
        # being set yet. Nested under show_ticks: a guide exists only to
        # connect a tick to a peak, so with the ticks off there is nothing
        # left for it to point at.
        for value in guide_lines:
            snapped, owner = snap_to_phase(value, listed, guide_snap)
            if snapped is None:
                print(f"  ! no reflection within {guide_snap:g} deg of "
                      f"{value:g}, no guide drawn")
                continue
            # Coloured like the tick it came from, so the guide says which
            # phase owns the peak it points at as well as where it is.
            ax.axvline(snapped, color=tick_colors[owner], ls=guide_style,
                       lw=guide_width, zorder=4)
        if rows:
            # Built from 'listed', not 'ordered': a phase whose reflections
            # all fall outside the window gets no tick row, and must get no
            # legend entry either, or the legend would name a phase that is
            # nowhere in the figure. 'listed' is already keyed by label, one
            # entry per phase, in the order the tick rows were drawn, which
            # is what pairs a legend entry to a row by position.
            handles = [Line2D([0], [0], color=tick_colors[label], lw=3,
                              label=label) for label in listed]
            # Opposite corner from the trace labels, which start at the
            # left. A short handle: the colour is the whole message, and a
            # long rule beside a short phase name reads as a second scale.
            ax.legend(handles=handles, loc="upper right",
                      fontsize=xp.FONT_SIZE_LEGEND, framealpha=1,
                      handlelength=1, handletextpad=1)

    ax.set_xlabel(r"$2\theta\:/\:^\circ$", fontsize=xp.FONT_SIZE_LABEL)
    ax.set_ylabel(r"$\sqrt{\mathrm{Intensity}}$ / a.u." if use_sqrt
                  else r"Intensity / a.u.", fontsize=xp.FONT_SIZE_LABEL,
                  labelpad=10)
    ax.set_xlim(*widest)
    ax.set_ylim(bottom=(block_top - (rows - 0.5) * pitch
                        if rows else -0.05 * offset),
                # 0.08 was clearance for a curve. A line of text needs the
                # same room the labels get between traces, or the topmost
                # one is cut by the frame. max(1.0, label_height): the top
                # of the trace is at 1.0, but label_height may sit above
                # it, and the frame has to clear whichever is higher.
                top=(len(traces) - 1) * offset + max(1.0, label_height)
                + 0.30 * offset)
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
    """Draw the series the metadata names and write the two output files."""
    meta = xp.load_metadata(METADATA_FILE)
    entries = series_from_metadata(meta)
    if not entries:
        raise SystemExit(
            f"no sample carries a 'series_order' in {METADATA_FILE}, so "
            "there is no series to draw. See docs/metadata.md")
    traces, phases, colors = load_series(entries, DATA_FOLDER, meta)
    fig = plot_series(traces, phases, colors=colors)
    # weighted=True: this figure draws no residual, so the unweighted
    # suffix output_basename adds for weighted=False would be meaningless
    # here. USE_SQRT still tells it sqrt from linear, so the two settings
    # write to different names instead of overwriting each other.
    name = xp.output_basename(OUTPUT_BASENAME, USE_SQRT, weighted=True)
    base = xp.save_figure(fig, OUTPUT_FOLDER, name)
    print(f"wrote {OUTPUT_FOLDER}/pdf/{base}.pdf "
          f"and {OUTPUT_FOLDER}/png/{base}.png")


if __name__ == "__main__":
    main()
