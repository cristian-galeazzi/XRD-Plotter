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
