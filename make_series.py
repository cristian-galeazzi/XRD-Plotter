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
from matplotlib.axes import Axes
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
#
# One warning runs through the whole block. A setting in degrees 2theta is
# not a preference: it is a position you read off your own patterns, and a
# few of them give a lattice spacing through Bragg's law and identify the
# phase. That is PLOT_X_MIN, PLOT_X_MAX, LABEL_X and GUIDE_LINES, each
# marked again where it stands. Leave them as shipped and set them in
# section 5 of the notebook, where they stay in the widget, or per sample in
# the metadata file. See docs/privacy.md.

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
#
# In degrees, so a window fixed here is measured data: a pair framing the
# reflection the series is about says where that reflection is. Set it in
# section 5, or per sample through the 'x_min'/'x_max' metadata columns.
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
#
# In degrees, so a value here is chosen against your own pattern: it names
# an empty stretch of your figure, which says where the reflections are not.
# Weaker than a window, but it is still a position. Set it in section 5.
LABEL_X: float | None = None

# Weight of the trace labels: 'normal', 'medium', 'semibold' or 'bold'.
# A heavier label separates itself from the pattern it sits over without a
# box behind it. STIXGeneral, the face these figures are set in, carries a
# real bold, so nothing is synthesised from the upright.
LABEL_WEIGHT = "normal"

# --- The reflection ticks -----------------------------------------------
# One row of ticks per phase, below the bottom trace, in the colours the
# per-sample figures use. False hides those rows and the legend that names
# them, and nothing else: the reflection list is still printed and the
# guides are still drawn, since a guide carries a reflection up through the
# stack whether or not a tick marks it at the bottom.
SHOW_TICKS = True

# Height of one tick, in trace heights, so 0.10 is a tenth of a pattern and
# does not change when OFFSET does. Taller ticks are easier to follow across
# a wide pattern. The rows sit one TICK_HEIGHT + 0.02 apart, so raising this
# moves them apart instead of overlapping them.
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
# The sharpest of the four settings in degrees: a guide sits exactly on a
# reflection, so a filled list is a measured peak table. Empty it before
# committing, or set the guides in section 5 of the notebook, where they
# stay in the widget. See docs/privacy.md.
GUIDE_LINES: list[float] = []

# The name of the reflection each guide marks, one per GUIDE_LINES value and
# in the same order, written above the stack in the colour of the phase that
# owns it. An empty string leaves that guide unnamed, and a list shorter than
# GUIDE_LINES leaves the rest of them unnamed:
#   GUIDE_LABELS = ["(111)", "(110)"]
#
# The text is drawn exactly as written, so an overbar is typed as mathtext,
# r"$(\bar{1}11)$", in the same STIX face as the axis labels. A Miller index
# is not a position and gives no lattice spacing on its own, but it names a
# reflection out of your own pattern and it pairs with GUIDE_LINES, so it
# ships empty for the same reason. Set it in section 5 of the notebook, where
# the position and the name are typed together as '30.95=(111)'. See
# docs/privacy.md.
GUIDE_LABELS: list[str] = []

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

# Weight of the guide names: 'normal', 'medium', 'semibold' or 'bold'. Its
# own setting rather than the trace labels' LABEL_WEIGHT, because the two
# sit differently: a trace label lies over a pattern, a guide name stands
# alone in the space above the stack, where bold keeps a short name legible
# once the figure is reduced to a column width. STIXGeneral carries a real
# bold, so nothing is synthesised from the upright.
GUIDE_LABEL_WEIGHT = "bold"

# The room a name keeps over the pattern beneath it, in trace heights, so
# 0.20 is a fifth of a pattern of clearance and does not change when the
# trace spacing does. At 0.0 a name rests on what it stands over.
#
# Measured over the whole width the name covers on the topmost trace, not at
# the one angle its guide sits at. A name is some two degrees wide on a
# seventy degree axis, and measured at a point it would clear its own
# reflection and still land on a taller neighbour, which is what a reader
# sees as a name stuck to a peak. So a name over a small reflection standing
# on its own sits low, and one over a crowded stretch rides above all of it.
# This is the clearance each name gets, not one height for all of them.
GUIDE_LABEL_HEIGHT = 0.20

# Angle of the names, in degrees anticlockwise from upright. 0 writes them
# across, 90 on end, reading upwards from just above the peak.
#
# It is the one setting that decides whether crowded names can be read at
# all. Upright, a Miller index covers some two degrees of a seventy degree
# axis and reaches its neighbours, so a group of reflections close together
# ends up as a staircase of names climbing away from the peaks they belong
# to. On end the same name covers about the width of one line of text, and
# the same group stands side by side, each over its own peak. Nothing else
# changes: the names are measured as they are drawn, so the clearance and
# the lifting follow the turn on their own.
GUIDE_LABEL_ROTATION = 0

# --- Where the files are ------------------------------------------------
DATA_FOLDER = Path("data")
METADATA_FILE = Path("Samples_metadata.csv")
OUTPUT_FOLDER = Path("output")
# Start of the output file name. xp.output_basename appends the rest, so
# this ships as 'series_stacked_XRD_analysis_sqrt.pdf' and its .png.
OUTPUT_BASENAME = "series_stacked"


def stack(patterns: list[np.ndarray], offset: float | None = None,
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
    # Resolved here, not in the signature: a default evaluated at import
    # would freeze the setting as it stood then, so a caller that raised
    # OFFSET afterwards would still be stacked at the shipped one.
    offset = OFFSET if offset is None else offset
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
                       tolerance: float | None = None) -> float | None:
    """The reflection position nearest to value, or None past the tolerance.

    A guide line is worth drawing only where a reflection actually is, so a
    value with nothing near it returns None for the caller to report rather
    than being drawn at the typed position.

    >>> snap_to_reflection(31.0, np.array([20.0, 30.95, 37.04]), 0.3)
    30.95
    >>> snap_to_reflection(31.0, np.array([20.0, 37.04]), 0.3) is None
    True
    """
    tolerance = GUIDE_SNAP if tolerance is None else tolerance
    if not len(positions):
        return None
    nearest = float(positions[np.argmin(np.abs(positions - value))])
    return nearest if abs(nearest - value) <= tolerance else None


def snap_to_phase(value: float, rows: dict[str, np.ndarray],
                  tolerance: float | None = None
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
    tolerance = GUIDE_SNAP if tolerance is None else tolerance
    best: tuple[float | None, str | None] = (None, None)
    smallest = tolerance
    for label, positions in rows.items():
        snapped = snap_to_reflection(value, positions, tolerance)
        if snapped is not None and (best[0] is None
                                    or abs(snapped - value) < smallest):
            best, smallest = (snapped, label), abs(snapped - value)
    return best


def parse_guides(text: str) -> tuple[list[float], list[str], list[str]]:
    """A typed guide list as positions, the name of each, and what would not read.

    One entry per guide, separated by a comma or a semicolon. An entry is a
    2theta on its own, or a 2theta, an '=' and the name to write over that
    guide. Typing the two together is what keeps them paired: two boxes, one
    of positions and one of names, would fall out of step the first time a
    guide was inserted in the middle of the list. The name is taken exactly
    as typed, so an overbar can be written as mathtext.

    The comma separates, it is not a decimal point: '28,4' is two guides.

    >>> parse_guides("30.95=(111), 37.04")
    ([30.95, 37.04], ['(111)', ''], [])
    >>> parse_guides("here=(110), 30.95")
    ([30.95], [''], ['here=(110)'])
    """
    positions: list[float] = []
    labels: list[str] = []
    unreadable: list[str] = []
    for part in text.replace(";", ",").split(","):
        if not part.strip():
            continue
        head, _sign, name = part.partition("=")
        value = xp.to_number(head)
        if value is None:
            # Named rather than dropped in silence: a typo would otherwise
            # take a guide out of the figure while its text sat in the box,
            # read back afterwards as if it had been drawn.
            unreadable.append(part.strip())
            continue
        positions.append(value)
        labels.append(name.strip())
    return positions, labels, unreadable


def peak_height(theta: np.ndarray, values: np.ndarray, position: float,
                window: float | None = None) -> float:
    """How high the pattern rises at a reflection, over a window around it.

    The name of a reflection belongs just above the peak it marks, so it has
    to know how tall that peak is on the trace it is written over. Read over
    a window rather than at the one angle: the ticks and the guides come
    from the first sample of the series, and in a series drawn because a
    reflection moves, the peak on the topmost trace has moved off it.
    'window' defaults to GUIDE_SNAP, the same tolerance the guide itself
    snapped with.

    A reflection with no pattern inside the window answers 0.0, the
    baseline of that trace, which is where a phase that has gone leaves its
    name.

    >>> theta, values = np.array([19.8, 20.0, 20.2]), np.array([0.1, 0.9, 0.2])
    >>> peak_height(theta, values, 20.0, 0.3)
    0.9
    >>> peak_height(theta, values, 40.0, 0.3)
    0.0
    """
    window = GUIDE_SNAP if window is None else window
    near = np.abs(theta - position) <= window
    return float(np.max(values[near])) if near.any() else 0.0


def stagger(boxes: list[tuple[float, float, float, float]],
            obstacles: list[tuple[float, float, float, float]],
            step: float, gap: float = 4.0) -> list[float]:
    """How far each box has to be lifted to clear the ones already there.

    Every box is (left, bottom, right, top) in one unit, pixels in the
    figure, and comes back with the distance to raise it, 0.0 for one that
    already stands clear. Boxes settle from the lowest up, so a name on a
    small reflection keeps the place its own peak gives it and a name that
    would land on top of another is the one that moves. 'obstacles' are
    boxes that never move, the trace labels among them.

    Sideways the boxes are kept 'gap' apart; vertically touching is enough,
    since 'step' is a line of text with air already in it.

    >>> stagger([(0.0, 0.0, 10.0, 5.0), (20.0, 0.0, 30.0, 5.0)], [], 6.0)
    [0.0, 0.0]
    >>> stagger([(0.0, 0.0, 10.0, 5.0), (5.0, 0.0, 15.0, 5.0)], [], 6.0)
    [0.0, 6.0]
    """
    lifts = [0.0] * len(boxes)
    placed = list(obstacles)
    for index in sorted(range(len(boxes)),
                        key=lambda i: (boxes[i][1], boxes[i][0])):
        left, low, right, high = boxes[index]
        lift = 0.0
        # Bounded by the number of boxes there are to clear: a lift that has
        # passed every one of them is clear by construction, so a wrong
        # overlap test cannot turn this into an endless climb.
        for _ in range(len(placed) + 1):
            if not any(right + gap > other_left and other_right + gap > left
                       and high + lift > other_low and other_high > low + lift
                       for other_left, other_low, other_right, other_high
                       in placed):
                break
            lift += step
        lifts[index] = lift
        placed.append((left, low + lift, right, high + lift))
    return lifts


def draw_guide_labels(ax: Axes, labels: list[tuple[float, str, str]],
                      pattern: tuple[np.ndarray, np.ndarray],
                      clearance: float | None = None,
                      snap: float | None = None,
                      weight: str | None = None,
                      rotation: float | None = None) -> None:
    """Write one name per guide, each clear of the pattern under it.

    'labels' is one (2theta, text, colour) per guide that was drawn, and
    'pattern' the (2theta, height) of the topmost trace, the one the names
    are written over. The limits of the axes have to be set before this
    runs, since a string has no width until it is rendered and what the
    renderer answers depends on the scale it is drawn at.

    'clearance' is the room each name keeps over that pattern, in trace
    heights, and it is measured over the whole width the name covers, not at
    the one angle its guide sits at. A name is about two degrees wide on a
    seventy degree axis, so measuring at a point would let it land on a
    taller neighbour while formally clearing its own peak. 'snap' widens
    that reach to at least the tolerance the guide itself snapped with, so a
    reflection that has moved on the topmost trace is still covered.

    A name over a small reflection therefore sits low and one over a tall
    reflection sits high, which is what keeps them apart without a rule: two
    names side by side only collide when the pattern under them is of a
    height. Where they do, the lower one keeps its place and the other is
    lifted a line at a time until it is clear, of the names already placed
    and of the trace labels, which were on the axes first.

    Whatever room is missing above them is taken by growing the figure, not
    by widening the axis inside a figure of the same size. The second would
    leave every pattern with less of the figure than it had, and the small
    reflections a series is drawn to show are the first thing that costs.
    The phase legend hangs from the top of the frame, so a name standing
    under it asks for that much more; a name that clears it sideways, or
    that sits low enough to pass beneath it, asks for nothing.

    Typical use, from inside plot_series::

        draw_guide_labels(ax, [(30.95, "(111)", "#EE8031")], (theta, values))
    """
    clearance = GUIDE_LABEL_HEIGHT if clearance is None else clearance
    rotation = GUIDE_LABEL_ROTATION if rotation is None else rotation
    snap = GUIDE_SNAP if snap is None else snap
    weight = GUIDE_LABEL_WEIGHT if weight is None else weight
    figure = ax.figure
    theta, values = pattern
    # Read before the names are added, so they are the obstacles and not
    # each other: the trace labels are the only text already on the axes.
    figure.canvas.draw()  # a string has no width until it is rendered
    obstacles = [tuple(text.get_window_extent().extents) for text in ax.texts]
    # Placed at the height of their own reflection to begin with. The width
    # they cover is not known until they are drawn, and the height that
    # width asks for is settled in the pass below.
    texts = [ax.text(position, peak_height(theta, values, position, snap)
                     + clearance, name, color=color, ha="center", va="bottom",
                     # No rotation_mode='anchor': that turns the name about
                     # its anchor and leaves the turned box hanging to one
                     # side of the guide. The default aligns after the turn,
                     # so a name on end still stands centred on its own line.
                     fontsize=xp.FONT_SIZE_LEGEND, fontweight=weight,
                     rotation=rotation, zorder=5)
             for position, name, color in labels]
    figure.canvas.draw()
    frame = ax.get_window_extent()
    inverse = ax.transData.inverted()
    bottom, top = ax.get_ylim()
    scale = frame.height / (top - bottom)  # pixels per trace height
    boxes, line = [], 0.0
    for text in texts:
        box = text.get_window_extent()
        # A name centred on a guide near the border would hang outside the
        # frame, and bbox_inches='tight' would then widen the saved figure
        # around it. Nudged in by the overhang, it sits a little off its own
        # guide, which is still the guide nearest to it.
        shift = max(0.0, frame.x0 - box.x0) - max(0.0, box.x1 - frame.x1)
        left = inverse.transform((box.x0 + shift, 0.0))[0]
        right = inverse.transform((box.x1 + shift, 0.0))[0]
        # Clear of everything the name covers, its own reflection included:
        # the wider of the two reaches wins.
        height = peak_height(theta, values, (left + right) / 2,
                             max(snap, (right - left) / 2)) + clearance
        # Moved by the overhang, not to the middle of its own box: the two
        # are the same only while a name is centred on its guide, and a
        # name that is not would walk further off it on every pass.
        moved = inverse.transform((shift, 0.0))[0] - inverse.transform(
            (0.0, 0.0))[0]
        text.set_position((text.get_position()[0] + moved, height))
        # va='bottom' anchors the box at that height, so the box the packing
        # works on follows from it without a second render.
        low = ax.transData.transform((0.0, height))[1]
        boxes.append((box.x0 + shift, low, box.x1 + shift, low + box.height))
        line = max(line, box.height)

    lifts = stagger(boxes, obstacles, line * 1.25)
    for text, lift in zip(texts, lifts):
        # In the scale the figure has now, which growing it below is built to
        # preserve, so a name keeps the distance to the pattern either way.
        text.set_y(text.get_position()[1] + lift / scale)

    # What the names need above them, in pixels: to stay inside the frame,
    # and to stay under the legend where one of them stands in its way. The
    # legend is pinned to the top of the frame, so the two are kept apart by
    # making room rather than by moving either of them.
    legend = ax.get_legend()
    legend_box = None if legend is None else legend.get_window_extent()
    wanted = 0.0
    for (left_px, _low, right_px, high), lift in zip(boxes, lifts):
        wanted = max(wanted, high + lift + 0.25 * line - frame.y1)
        if (legend_box is not None and right_px >= legend_box.x0
                and left_px <= legend_box.x1):
            wanted = max(wanted, high + lift + 0.25 * line - legend_box.y0)
    # Doubling the figure is as far as this goes. Past it the names are worth
    # less than the patterns under them, so the run says what it could not do
    # rather than growing without end or cutting the top names in silence.
    deficit = min(max(0.0, wanted), frame.height)
    if deficit < wanted:
        print(f"  ! {len(labels)} names need more room than the figure can "
              "grow to hold; the top ones are cut. Name fewer guides, or "
              "widen the window so fewer of them are lifted")
    if deficit:
        # The room is taken from a taller figure, not from the stack, and
        # the axes keeps its share of a figure whose margins are relative,
        # so the figure has to grow by more than the axes needs.
        width, height = figure.get_size_inches()
        figure.set_size_inches(
            width, height + deficit / (frame.height / figure.bbox.height)
            / figure.dpi)
        ax.set_ylim(bottom, top + deficit / scale)


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
                guide_labels: list[str] | None = None,
                guide_label_weight: str | None = None,
                guide_label_height: float | None = None,
                guide_label_rotation: float | None = None,
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

    'guide_labels' pairs with 'guide_lines' by position: the name of the
    reflection that guide marks, or '' for none. Each name is written
    'guide_label_height' above its own peak on the topmost trace, in the
    colour of the phase its guide landed on, so a name over a small
    reflection sits low and one over a tall reflection sits high.

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
    guide_labels = GUIDE_LABELS if guide_labels is None else guide_labels
    guide_label_weight = (GUIDE_LABEL_WEIGHT if guide_label_weight is None
                          else guide_label_weight)
    guide_label_height = (GUIDE_LABEL_HEIGHT if guide_label_height is None
                          else guide_label_height)
    guide_label_rotation = (GUIDE_LABEL_ROTATION
                            if guide_label_rotation is None
                            else guide_label_rotation)
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

    # Guarded here rather than left to np.concatenate below, which answers an
    # empty series with 'need at least one array to concatenate'. This is a
    # public entry point, so it says what the caller actually got wrong.
    if not traces:
        raise ValueError("no traces to draw: the series is empty")

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
    pitch = tick_height + 0.02  # 0.02 of a trace height between rows
    # The tick block hangs below the bottom trace by 0.05 of a trace height,
    # whatever tick_height is. Derived from it rather than written out again,
    # so a taller tick moves the block down instead of growing into the
    # pattern above it.
    block_top = -(tick_height + 0.05)
    listed: dict[str, np.ndarray] = {}
    # Collected whether or not the rows are drawn. The printed list is how a
    # guide value is picked and the guides are what carry a reflection up
    # through the stack, and neither needs a tick under the bottom trace to
    # make sense: show_ticks hides the rows and the legend naming them, and
    # nothing else.
    for phase in ordered:
        locations = reflections_in_window(phases[phase], widest)
        if not len(locations):
            continue
        # Two columns can carry the same phase name, and the legend gives
        # them one entry between them. Their positions join instead of the
        # second replacing the first, or the list and the guides would know
        # about half the phase.
        listed[labels[phase]] = (
            np.concatenate([listed[labels[phase]], locations])
            if labels[phase] in listed else locations)
        if show_ticks:
            row_y = block_top - rows * pitch
            ax.vlines(locations, row_y, row_y + tick_height,
                      colors=tick_colors[labels[phase]],
                      lw=xp.LINEWIDTH_TICKS, zorder=3)
            rows += 1

    # Printed so a guide value is picked from the reflections that are
    # actually there rather than read off the drawn figure by eye.
    if listed:
        print(f"reflections between {widest[0]:g} and {widest[1]:g} deg, "
              "to pick the guides from:")
        print(format_reflections(listed))

    # Guides span the full height of the axes, in front of the black traces,
    # drawn at zorder=4 above their zorder=2, so a guide behind a dense stack
    # of patterns is not hidden by them and a reflection can still be
    # followed up to the peak it belongs to. axvline works in axes fractions
    # vertically, so it does not depend on the limits being set yet.
    # Reported rather than ignored: a name typed without the position it
    # belongs to would otherwise vanish, and the reader would look for it on
    # a guide that was never asked for.
    if len(guide_labels) > len(guide_lines):
        print(f"  ! {len(guide_labels)} names for {len(guide_lines)} guides, "
              "the ones past the end of the guide list are not drawn")
    named: list[tuple[float, str, str]] = []
    for index, value in enumerate(guide_lines):
        snapped, owner = snap_to_phase(value, listed, guide_snap)
        if snapped is None:
            print(f"  ! no reflection within {guide_snap:g} deg of "
                  f"{value:g}, no guide drawn")
            continue
        # Coloured like the tick it came from, so the guide says which phase
        # owns the peak it points at as well as where it is.
        ax.axvline(snapped, color=tick_colors[owner], ls=guide_style,
                   lw=guide_width, zorder=4)
        # By position in the list, so a guide dropped for having no
        # reflection near it does not shift the next one's name onto itself.
        name = guide_labels[index] if index < len(guide_labels) else ""
        if name.strip():
            named.append((snapped, name.strip(), tick_colors[owner]))
    if rows:
        # Built from 'listed', not 'ordered': a phase whose reflections all
        # fall outside the window gets no tick row, and must get no legend
        # entry either, or the legend would name a phase that is nowhere in
        # the figure. 'listed' is already keyed by label, one entry per
        # phase, in the order the tick rows were drawn, which is what pairs
        # a legend entry to a row by position.
        handles = [Line2D([0], [0], color=tick_colors[label], lw=3,
                          label=label) for label in listed]
        # Opposite corner from the trace labels, which start at the left. A
        # short handle: the colour is the whole message, and a long rule
        # beside a short phase name reads as a second scale.
        ax.legend(handles=handles, loc="upper right",
                  fontsize=xp.FONT_SIZE_LEGEND, framealpha=1,
                  handlelength=1, handletextpad=1)

    ax.set_xlabel(r"$2\theta\:/\:^\circ$", fontsize=xp.FONT_SIZE_LABEL)
    ax.set_ylabel(r"$\sqrt{\mathrm{Intensity}}$ / a.u." if use_sqrt
                  else r"Intensity / a.u.", fontsize=xp.FONT_SIZE_LABEL,
                  labelpad=10)
    ax.set_xlim(*widest)
    # 0.08 was clearance for a curve. A line of text needs the same room the
    # labels get between traces, or the topmost one is cut by the frame.
    # max(1.0, label_height): the top of the trace is at 1.0, but
    # label_height may sit above it, and the frame has to clear whichever is
    # higher.
    stack_top = (len(traces) - 1) * offset + max(1.0, label_height)
    ax.set_ylim(bottom=(block_top - (rows - 0.5) * pitch
                        if rows else -0.05),
                top=stack_top + 0.30 * offset)
    if named:
        # After the limits and after the legend, both of which it measures.
        # It makes whatever room the names need itself; a figure with no name
        # keeps the limit set above, untouched. The topmost trace is the one
        # the names are written over, cropped to the window, so a name is
        # placed against what is actually drawn under it.
        top_theta, _top_obs, _top_name = traces[-1]
        draw_guide_labels(ax, named,
                          (top_theta[scopes[-1]], stacked[-1][scopes[-1]]),
                          guide_label_height, guide_snap,
                          guide_label_weight, guide_label_rotation)
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
