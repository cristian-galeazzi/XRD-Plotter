"""Rietveld plots from the CSV a GSAS-II Rietveld plot saves.

Parsing, data preparation, plotting and the batch driver. The notebook and
the tests both import this module, so a figure drawn by hand matches the one
the batch writes.

Appearance is set by the constants below. Override them on the module, not
on a copy, so the functions see the change:

    import xrd_plotter as xp
    xp.PLOT_X_MIN, xp.PLOT_X_MAX = 13, 85

Importing this module also sets the serif typography on the matplotlib
defaults, for the whole process rather than for one figure, so any other
figure drawn in the same kernel is set in that face too.
"""
import io
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

# --- Plot appearance --------------------------------------------------
# 2theta window (deg). None takes the measured range of each file, so a
# pattern is never cropped without being asked. Set both to fix every
# figure to the same window, or give x_min and x_max per sample in the
# metadata file, which wins over these.
PLOT_X_MIN, PLOT_X_MAX = None, None
FIGURE_WIDTH, FIGURE_HEIGHT = 12, 10   # inches
DPI_EXPORT = 600                       # PNG export resolution
PREVIEW_WIDTH_PX = 640                 # on-screen width of the section-3 inline
                                       # preview; the saved PNG and PDF are
                                       # untouched, so lowering this only shrinks
                                       # what scrolls past during a run

# Lower panel: True draws diff/sigma, the residual divided by the standard
# deviation of the point, so a well-fitted pattern stays inside a band of a
# few units. False draws the raw diff, where the tall reflections dominate.
# Files carry the suffix '_unweighted' when this is False.
WEIGHTED_RESIDUALS = True
# Smallest half height of the weighted residual panel, in units of sigma. The
# panel is never narrower than this, so every well-fitted pattern is drawn on
# the same residual scale and two figures can be read side by side. A fit that
# leaves a residual above it widens the panel to fit, since clipping the misfit
# is the one thing this panel must not do. Raise it to flatten the trace,
# lower it to magnify. It does not apply to the raw diff, which has no unit
# that means the same thing twice.
RESIDUAL_SPAN = 5.0

COLOR_OBS = "#000000"
COLOR_CALC = "#D62728"
COLOR_BKG = "#7f7f7f"
COLOR_RESIDUALS = "#696969"

LINEWIDTH_BORDER = 2.0
LINEWIDTH_CALC = 1.5
LINEWIDTH_RESIDUALS = 0.8
LINEWIDTH_TICKS = 2.0
MARKER_SIZE_OBS = 15

FONT_SIZE_LABEL = 20
FONT_SIZE_LEGEND = 16
FONT_SIZE_TICK = 16

# Publication typography: serif text, STIX maths. STIXGeneral is the text
# companion of that maths font and ships inside matplotlib, so a machine
# without Times New Roman still draws words that match the symbols beside
# them. DejaVu Serif, the matplotlib default, is far wider than Times and
# would not, which is why it is last rather than second.
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
})

# --- Phase columns ----------------------------------------------------
# Headers GSAS-II writes itself, plus the fit statistics people keep beside
# the pattern. Everything else that holds only a few values is taken to be
# the reflection positions of a phase, so phases are found by their own
# names instead of a keyword list. Entries are compared folded, so 'rw'
# covers 'Rw', 'Rw%' and 'Rw / %' without an entry for each spelling. Add
# your own header here to keep a column of your own out of the legend.
NON_PHASE_COLUMNS = frozenset({
    "used", "obs", "calc", "bkg", "diff", "diff/sigma", "weight", "weights",
    "sig", "sigma", "tick-pos", "tick pos", "axis-limits", "axis limits",
    "excluded", "gof", "rw", "rwp", "rexp", "chi2",
})
# 'rp' is deliberately absent: GSAS-II writes Rwp, and 'RP' is how a
# phase is labelled. Every entry is a name a phase could
# carry, 'sigma' for the intermetallic among them, and a phase whose header
# matches one is not drawn. The 'phases:' line each run prints is the check.
# A reflection list is short: at most this share of the pattern's rows. A
# column above the ceiling is reported and left undrawn, never dropped in
# silence.
PHASE_MAX_FILL = 0.5
# Legend name for a phase whose column header is not the name to print.
# Keys are matched case-insensitively as a fragment of the header, and the
# longest matching key wins:
#   PHASE_LABELS = {"phase 2": "Phase 2, high temperature"}
PHASE_LABELS = {}
# Tick colour tied to a phase name rather than to a row, so a phase keeps
# its colour whether or not the others are present in that sample. Keys are
# matched as a fragment of the legend name, the name PHASE_LABELS produced:
#   PHASE_COLORS = {"phase 1": "#1f77b4"}
# The same is available per sample, and privately, through the
# '<phase>_color' columns of the metadata file, which win over this.
PHASE_COLORS = {}
# Colours for every phase with no colour of its own, in legend order.
PHASE_COLOR_CYCLE = ("#EE8031", "#1f77b4", "#FF1493", "#2CA02C", "#9467BD",
                     "#8C564B")


def fmt_pct(p: float) -> str:
    """Phase fraction label: one decimal, trailing '.0' dropped (62.0 -> 62).

    >>> fmt_pct(62.0)
    '62'
    >>> fmt_pct(62.3)
    '62.3'
    """
    s = f"{p:.1f}"
    return s[:-2] if s.endswith(".0") else s


def residual_column(weighted: bool | None = None) -> str:
    """Residual the lower panel is drawn from, per WEIGHTED_RESIDUALS.

    'diff/sigma' is a column read from the file. 'diff' names the raw
    residual, which is computed as obs - calc rather than read, because the
    column of that name in a GSAS-II export is the plotted curve and carries
    the plot's offset.

    >>> residual_column()
    'diff/sigma'
    >>> residual_column(False)
    'diff'
    """
    if weighted is None:
        weighted = WEIGHTED_RESIDUALS
    return "diff/sigma" if weighted else "diff"


def residual_limits(resid: np.ndarray,
                    weighted: bool | None = None) -> tuple[float, float]:
    """Symmetric limits for the lower panel, so zero sits at its middle.

    The half height is the largest residual drawn, with a margin above it. In
    the weighted panel it is never smaller than RESIDUAL_SPAN, because
    diff/sigma is measured in standard deviations of the point and so means
    the same thing in every sample: a floor puts every well-fitted pattern on
    one scale, and only a fit that genuinely misses widens the panel. The raw
    diff is in the intensities of one measurement, comparable across nothing,
    so it gets no floor and only the centring.

    >>> residual_limits(np.array([-1.0, 2.0]), weighted=True)
    (-5.0, 5.0)
    >>> residual_limits(np.array([-2.0, 1.0]), weighted=False)
    (-2.1, 2.1)
    """
    finite = resid[~np.isnan(resid)]
    half = 1.05 * float(np.max(np.abs(finite))) if len(finite) else 0.0
    if residual_column(weighted) == "diff/sigma":
        half = max(half, RESIDUAL_SPAN)
    # A residual that is flat, empty or all NaN leaves no scale to take, and
    # matplotlib refuses limits that are equal. RESIDUAL_SPAN is an arbitrary
    # but harmless choice here: the panel is empty either way.
    return (-half, half) if half > 0.0 else (-RESIDUAL_SPAN, RESIDUAL_SPAN)


def fold(header: str) -> str:
    """Header reduced to the characters that carry its meaning.

    Case, the Greek theta and the punctuation exporters vary are dropped, so
    '2 theta', '2-theta' and '2θ' become one key, and so do 'Rw', 'Rw%' and
    'Rw / %'. Used for the angle header and for NON_PHASE_COLUMNS, which is
    why one entry there covers a column's spellings.

    >>> fold("2-theta")
    '2theta'
    >>> fold("Rw%")
    'rw'
    """
    text = header.lower().replace("θ", "theta")
    for ch in " _-/%()":
        text = text.replace(ch, "")
    return text


def is_theta_header(header: str) -> bool:
    """True when a column header names the diffraction angle.

    Matched on the folded header, so '2theta', '2 theta', 'two-theta' and
    '2θ' all count. A header whose first field is 'x' counts as well, which
    covers both 'x' alone and 'x, 2theta (deg)'.

    >>> is_theta_header("2-theta (deg)")
    True
    >>> is_theta_header("Yobs")
    False
    """
    folded = fold(header)
    return ("2theta" in folded or "twotheta" in folded
            or header.lower().split(",")[0].strip() == "x")


def is_monotonic(values: np.ndarray, tolerance: float = 0.01) -> bool:
    """True when the finite values run one way from end to end.

    A 2theta axis climbs, or on some instruments descends, from the first row
    to the last. An intensity column read under the 2theta header wanders.

    The two are told apart by distance travelled, not by counting turns: a
    column that scores one backward step out of five hundred passes either
    way, but a single backward jump the width of the scan does not, and that
    is what two exports pasted into one file looks like. Repeated points
    travel nowhere and so count for neither direction, which leaves a column
    of one repeated value with no distance at all: not an axis, and refused.

    >>> is_monotonic(np.array([1.0, 2.0, 3.0]))
    True
    >>> is_monotonic(np.array([1.0, 5.0, 2.0, 4.0, 3.0]))
    False
    """
    finite = values[~np.isnan(values)]
    if len(finite) < 3:
        return True  # too short to judge, the other checks still apply
    steps = np.diff(finite)
    rise = float(steps[steps > 0].sum())
    fall = float(-steps[steps < 0].sum())
    if rise + fall == 0.0:
        return False  # a constant column, such as a mask of ones
    return min(rise, fall) <= tolerance * (rise + fall)


def plot_window(theta: np.ndarray, x_min: float | None = None,
                x_max: float | None = None) -> tuple[float, float]:
    """2theta limits: the per-sample pair, then the constants, then the data.

    >>> plot_window(np.array([10.0, 20.0, 30.0]))
    (10.0, 30.0)
    >>> plot_window(np.array([10.0, 20.0, 30.0]), x_min=12.0, x_max=28.0)
    (12.0, 28.0)
    """
    low = x_min if x_min is not None else PLOT_X_MIN
    high = x_max if x_max is not None else PLOT_X_MAX
    return (float(np.nanmin(theta)) if low is None else float(low),
            float(np.nanmax(theta)) if high is None else float(high))


def phase_label(column: str) -> str:
    """Legend name for a phase column: 'Phase 1 hkl' -> 'Phase 1'.

    >>> phase_label("Phase 1 hkl")
    'Phase 1'
    >>> phase_label("phase2")
    'Phase2'
    """
    label = column.strip()
    if label.lower().endswith("hkl"):
        label = label[:-3].strip()
    override = longest_match(PHASE_LABELS, label)
    return override if override else label[:1].upper() + label[1:]


def longest_match(mapping: dict[str, str] | None, label: str) -> str | None:
    """Value whose key is the longest fragment of label, or None.

    The longest key wins so that a specific column ('phase 1_color') beats a
    general one ('phase_color') instead of resolving on dictionary order.

    >>> longest_match({"phase": "gray", "phase 1": "red"}, "Phase 1")
    'red'
    >>> longest_match({"Phase 1": "red"}, "phase 1")
    'red'
    >>> longest_match({"phase": "gray"}, "Other")
    """
    # Both sides folded: a key typed as it appears in the legend ('Phase 1')
    # has to match as readily as the lowercase keys the metadata file yields.
    keys = [k for k in (mapping or {}) if k and k.lower() in label.lower()]
    return mapping[max(keys, key=len)] if keys else None


def phase_colors(ordered_labels: list[str],
                 overrides: dict[str, str] | None = None) -> dict[str, str]:
    """Tick colour per phase: the metadata, then PHASE_COLORS, then the cycle.

    Both mappings are keyed by a fragment of the legend name. Repeated
    labels share one colour and consume one slot of the cycle, matching the
    single legend entry they get.

    >>> phase_colors(["Quartz", "Quartz", "Calcite"])
    {'Quartz': '#EE8031', 'Calcite': '#1f77b4'}
    """
    colors, cycled = {}, 0
    for label in dict.fromkeys(ordered_labels):
        color = longest_match(overrides, label) or longest_match(PHASE_COLORS,
                                                                 label)
        if color is None:
            color = PHASE_COLOR_CYCLE[cycled % len(PHASE_COLOR_CYCLE)]
            cycled += 1
        colors[label] = color
    return colors


def phase_fraction(pct: dict[str, float], label: str) -> float:
    """Percentage for one phase, matched on the metadata key as a fragment.

    >>> phase_fraction({"phase 1": 62.0}, "Phase 1")
    62.0
    >>> phase_fraction({"a": 1.0, "ab": 2.0}, "ab")
      ! metadata columns a_pct, ab_pct all match 'ab', no percentage printed.
    0.0
    """
    hits = sorted(k for k in pct if k in label.lower())
    if len(hits) > 1:
        print(f"  ! metadata columns {', '.join(h + '_pct' for h in hits)} "
              f"all match '{label}', no percentage printed.")
        return 0.0
    return pct[hits[0]] if hits else 0.0


# --- Parsing -----------------------------------------------------------
def read_gsas2_csv(csv_path: str | Path, weighted: bool | None = None
                   ) -> tuple[dict[str, np.ndarray] | None,
                              list[str] | None, str | None]:
    """Read one CSV saved from a GSAS-II Rietveld plot.

    Returns (data, phase_cols, error): on success error is None, on failure
    data is None and error holds the reason. Non-fatal issues are printed
    as warnings. 'weighted' overrides WEIGHTED_RESIDUALS for this call.

    >>> import tempfile
    >>> from pathlib import Path
    >>> rows = ["2theta;Obs;Calc;Bkg;diff/sigma", "10.0;5;4;1;0.5",
    ...         "10.5;6;5;1;0.4", "11.0;7;6;1;-0.3"]
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = Path(d) / "sample.csv"
    ...     with open(p, "w") as fh:
    ...         for row in rows:
    ...             print(row, file=fh)
    ...     data, phase_cols, error = read_gsas2_csv(p)
    >>> error is None
    True
    >>> sorted(data)
    ['bkg', 'calc', 'obs', 'resid', 'x']
    """
    csv_path = Path(csv_path)
    warnings = []
    try:
        # Columns are read as strings: pandas' fast C float parser is not
        # correctly rounded (up to 1 ULP off), Python's float() is.
        # GSAS-II uses ';' or ',' as separator depending on locale.
        try:
            df = pd.read_csv(csv_path, sep=";", encoding="utf-8", dtype=str)
            if len(df.columns) <= 1:
                df = pd.read_csv(csv_path, sep=",", encoding="utf-8", dtype=str)
        except (pd.errors.ParserError, UnicodeDecodeError):
            # Ragged or mis-encoded file. The python engine hands each bad
            # row to a callback rather than giving up on the whole file,
            # and dtype=str keeps the surviving rows correctly rounded.
            bad = []
            for sep in (";", ","):
                bad.clear()
                df = pd.read_csv(csv_path, sep=sep, encoding="latin-1",
                                 dtype=str, engine="python",
                                 on_bad_lines=lambda row: bad.append(row))
                if len(df.columns) > 1:
                    break
            if bad:
                warnings.append(f"{len(bad)} malformed row(s) skipped")

        df.columns = df.columns.str.strip()
        if df.empty:
            return None, None, "empty CSV file"

        col_lower_map = {c.lower(): c for c in df.columns}

        def clean(s):
            """Correctly-rounded str -> float64, ',' decimal mark accepted."""
            def to_float(v):
                try:
                    return float(v.replace(",", "."))
                except (AttributeError, ValueError, TypeError):
                    return np.nan  # missing or non-numeric entry
            return np.array([to_float(v) for v in s], dtype=np.float64)

        data = {}

        # 2theta column, across the spellings exporters use. Among the
        # candidates take the first one filled like a data column rather than
        # like a reflection list, so a phase column whose header happens to
        # match never becomes the axis. Half is a property of an axis, not the
        # PHASE_MAX_FILL knob, which the reader is free to move.
        candidates = [col for col in df.columns if is_theta_header(col)]
        if not candidates:
            return None, None, "2theta column not found"
        theta = next((c for c in candidates
                      if np.count_nonzero(~np.isnan(clean(df[c])))
                      > 0.5 * len(df)), candidates[0])
        if len(candidates) > 1:
            warnings.append(f"{len(candidates)} columns could be the 2theta "
                            f"axis ({', '.join(candidates)}), '{theta}' used")
        data["x"] = clean(df[theta])
        if np.all(np.isnan(data["x"])):
            return None, None, "no valid data in 2theta column"
        n_nan = int(np.isnan(data["x"]).sum())
        if n_nan:
            warnings.append(f"column '{theta}': {n_nan} non-numeric values "
                            "-> NaN, those points are not drawn")
        if not is_monotonic(data["x"]):
            return None, None, (
                f"the '{theta}' column does not run from one end of the scan "
                "to the other, so it holds something other than angles and "
                "the header names do not line up with the columns they sit "
                "above. See docs/input-format.md")
        # Degrees, so the axis lives inside 0 to 180. A shift can land a
        # smooth monotone column such as the background under the 2theta
        # header, which the test above cannot see but this one can.
        low, high = float(np.nanmin(data["x"])), float(np.nanmax(data["x"]))
        if low < 0.0 or high > 180.0:
            return None, None, (
                f"the '{theta}' column runs from {low:g} to {high:g}, outside "
                "the 0 to 180 degrees a 2theta axis occupies, so it holds "
                "something other than angles. See docs/input-format.md")

        # Main pattern columns, in the fixed GSAS-II layout.
        missing = []
        for k in ("obs", "calc", "bkg"):
            col = col_lower_map.get(k)
            if col:
                data[k] = clean(df[col])
                n_nan = int(np.isnan(data[k]).sum())
                if n_nan:
                    warnings.append(f"column '{k}': {n_nan} non-numeric values -> NaN")
            else:
                data[k] = np.zeros_like(data["x"])
                missing.append(k)
        if missing:
            warnings.append(f"missing columns (using zeros): {', '.join(missing)}")
        if np.all(np.isnan(data["obs"])):
            return None, None, "no valid data in the obs column"

        # Residuals are not filled with zeros when absent: a flat lower
        # panel reads as a perfect fit, which is the one lie the figure
        # must not tell.
        wanted = residual_column(weighted)
        if wanted == "diff":
            # GSAS-II copies its 'diff' column off the plotted line, which the
            # plot holds below the pattern by delOffset: 2% of the largest
            # observed intensity by default, anything at all once the curve is
            # dragged, and the square roots of both arrays in a square-root
            # plot. Subtracting the two columns already read is exact and does
            # not depend on how the export was drawn. 'diff/sigma' needs no
            # such care: GSAS-II recomputes that one from the data.
            absent = [k for k in ("obs", "calc") if k in missing]
            if absent:
                # Both are named when both are gone, so the reader fixes the
                # file once instead of learning about the second on a rerun.
                names = " and ".join(f"'{k}'" for k in absent)
                was = "column was" if len(absent) == 1 else "columns were"
                return None, None, ("the raw residual is computed as obs - "
                                    f"calc, and the {names} {was} not found")
            data["resid"] = data["obs"] - data["calc"]
            # Named for the message below: sending the reader to a 'diff'
            # column would send them to one this branch never opened.
            source = "obs - calc"
        else:
            resid_col = col_lower_map.get(wanted)
            if resid_col is None:
                return None, None, f"residual column '{wanted}' not found"
            data["resid"] = clean(df[resid_col])
            source = f"the '{wanted}' column"
        if np.all(np.isnan(data["resid"])):
            return None, None, f"no valid data in {source}"
        # Only the read path counts NaN here. On the computed path the NaN
        # came from obs or calc, and both have already reported their own.
        n_nan = int(np.isnan(data["resid"]).sum()) if wanted != "diff" else 0
        if n_nan:
            warnings.append(f"column '{wanted}': {n_nan} non-numeric values -> NaN")

        # Per-phase reflection-position columns: whatever is left over and
        # sparse enough to be a reflection list rather than a data column.
        phase_cols = []
        # Folded here, not at import, so a header added to NON_PHASE_COLUMNS
        # on the module takes effect on the next call.
        blocked = {fold(name) for name in NON_PHASE_COLUMNS}
        for cl, col in col_lower_map.items():
            # pandas renames a repeated header 'tick-pos' to 'tick-pos.1',
            # which the blocklist would otherwise miss.
            head, _, tail = cl.rpartition(".")
            base = head if head and tail.isdigit() else cl
            if col == theta or fold(base) in blocked:
                continue
            vals = clean(df[col])
            n_values = int(np.count_nonzero(~np.isnan(vals)))
            if not n_values:
                continue
            if n_values > PHASE_MAX_FILL * len(vals):
                warnings.append(f"column '{col}': {n_values} values in "
                                f"{len(vals)} rows, too many for a reflection "
                                "list, not drawn")
                continue
            phase_cols.append(col)
            data[col] = vals

        return data, phase_cols, None

    except Exception as e:  # isolate unreadable files, keep the batch running
        return None, None, f"read error: {e}"
    finally:
        # On every way out, so a file that is refused still shows what was
        # noticed before the refusal. Two candidate 2theta columns, or a
        # count of skipped rows, is usually the explanation of the refusal.
        if warnings:
            print("\n".join(f"  ! {w}" for w in warnings))


def load_metadata(metadata_path: str | Path) -> pd.DataFrame:
    """Read the (private) sample metadata CSV, indexed by filename.

    Expected columns: 'filename' (or 'file'), optional 'formula', one
    '<phase>_pct' and one '<phase>_color' column per phase, and an
    optional 'x_min'/'x_max' pair. Returns an empty DataFrame when the
    file is absent. The contents are never displayed by this notebook.

    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = Path(d) / "meta.csv"
    ...     with open(p, "w") as fh:
    ...         print("filename;formula", file=fh)
    ...         print("sample.csv;NaAlSi3O8", file=fh)
    ...     df = load_metadata(p)
    Metadata loaded for 1 sample(s).
    >>> list(df.index)
    ['sample.csv']
    """
    metadata_path = Path(metadata_path)
    if not metadata_path.is_file():
        return pd.DataFrame()
    try:
        df = pd.read_csv(metadata_path, sep=";", encoding="utf-8", dtype=str)
        if len(df.columns) <= 1:
            df = pd.read_csv(metadata_path, sep=",", encoding="utf-8", dtype=str)
    except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError):
        try:
            df = pd.read_csv(metadata_path, sep=",", encoding="latin-1", dtype=str)
        except Exception as e:
            # A file with no header row at all reaches here. It is read
            # outside the per-file guard of the batch, so raising would end
            # the whole run over a metadata file someone has just created.
            print(f"  ! metadata could not be read ({e}), ignored.")
            return pd.DataFrame()

    df.columns = df.columns.str.strip().str.lower()
    fname_col = next((c for c in df.columns if c in ("filename", "file")), None)
    if fname_col is None:
        print("  ! metadata needs a 'filename' or 'file' column, ignored.")
        return pd.DataFrame()
    df["filename"] = df[fname_col].astype(str).str.strip()
    df = df.drop_duplicates("filename", keep="first").set_index("filename")
    print(f"Metadata loaded for {len(df)} sample(s).")
    return df


def to_number(value: object) -> float | None:
    """One metadata cell as a float, or None when it is blank or not a number.

    >>> to_number("12,5")
    12.5
    >>> to_number("abc")
    """
    try:
        number = float(str(value).replace(",", "."))  # decimal comma
    except (ValueError, TypeError):
        return None
    return None if np.isnan(number) else number


def to_color(value: object) -> str | None:
    """One metadata cell as a colour, or None when blank or not a colour.

    >>> to_color("red")
    'red'
    >>> to_color("notacolor")
      ! 'notacolor' is not a colour, the cycle is used instead.
    """
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    if not is_color_like(text):
        print(f"  ! '{text}' is not a colour, the cycle is used instead.")
        return None
    return text


def metadata_keys(columns: Iterable[str], suffix: str) -> dict[str, str]:
    """Column names ending in suffix, keyed by the phase fragment before it.

    A column named exactly like the suffix carries no phase name and is
    skipped, since its empty key would match every phase.

    >>> metadata_keys(["phase_1_pct", "sample_id"], "_pct")
    {'phase 1': 'phase_1_pct'}
    >>> metadata_keys(["_pct"], "_pct")
    {}
    """
    return {c[:-len(suffix)].replace("_", " ").strip(): c for c in columns
            if c.endswith(suffix) and len(c) > len(suffix)}


def sample_info(meta_df: pd.DataFrame, filename: str, default_name: str
                ) -> tuple[str, dict[str, float], dict[str, str],
                           tuple[float | None, float | None]]:
    """Name, fractions, colours and 2theta window for one file (empty-safe).

    Every '<phase>_pct' and '<phase>_color' column becomes one entry, keyed
    by the part before the suffix with underscores as spaces, so
    'phase_1_color' reaches the phase named 'Phase 1' in the figure.
    'x_min' and 'x_max' are returned as given, for plot_window to apply.

    >>> import pandas as pd
    >>> df = pd.DataFrame({"formula": ["NaAlSi3O8"],
    ...                    "phase1_pct": ["62.0"]},
    ...                   index=pd.Index(["a.csv"], name="filename"))
    >>> name, pct, colors, window = sample_info(df, "a.csv", "a")
    >>> name, pct
    ('NaAlSi3O8', {'phase1': 62.0})
    """
    name, pct, colors, window = default_name, {}, {}, (None, None)
    if filename not in meta_df.index:
        return name, pct, colors, window

    row = meta_df.loc[filename]
    formula = row.get("formula")
    if isinstance(formula, str) and formula.strip():
        name = formula.strip()  # an empty cell is NaN, which would print as 'nan'

    for key, col in metadata_keys(meta_df.columns, "_pct").items():
        value = to_number(row[col])
        if value is not None:
            pct[key] = value
    for key, col in metadata_keys(meta_df.columns, "_color").items():
        value = to_color(row[col])
        if value is not None:
            colors[key] = value

    window = (to_number(row.get("x_min")), to_number(row.get("x_max")))
    return name, pct, colors, window


# --- Data preparation ---------------------------------------------------
def prepare_data(data: dict[str, np.ndarray], phase_cols: list[str],
                 use_sqrt: bool = True
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray,
                            np.ndarray, dict[str, np.ndarray]]:
    """Mask invalid points; optionally apply the display-only sqrt transform.

    >>> data = {"x": np.array([1.0, np.nan]), "obs": np.array([4.0, 5.0]),
    ...         "calc": np.array([9.0, 1.0]), "bkg": np.array([0.0, 0.0]),
    ...         "resid": np.array([0.0, 0.0])}
    >>> x, obs, calc, bkg, resid, phases = prepare_data(data, [])
    >>> len(x)
    1
    >>> obs
    array([2.])
    """
    mask = ~(np.isnan(data["x"]) | np.isnan(data["obs"]))
    if not mask.any():
        raise ValueError("no point has both a 2theta and an obs value")
    x = data["x"][mask]

    if use_sqrt:
        # Clipped, not made positive: a negative intensity has no square root,
        # and a background-subtracted export full of them would be mirrored
        # into peaks that were never measured. Clipping puts them on the axis.
        obs, calc, bkg = (np.sqrt(np.clip(data[k][mask], 0.0, None))
                          for k in ("obs", "calc", "bkg"))
    else:
        obs, calc, bkg = (data[k][mask] for k in ("obs", "calc", "bkg"))

    resid = data["resid"][mask]

    # Phase columns hold independent reflection positions (different length
    # from the pattern): drop their NaN padding individually.
    phases = {k: data[k][~np.isnan(data[k])] for k in phase_cols}
    return x, obs, calc, bkg, resid, phases


# --- Plotting -------------------------------------------------------------
def create_plot(theta: np.ndarray, obs: np.ndarray, calc: np.ndarray,
                bkg: np.ndarray, resid: np.ndarray,
                phases: dict[str, np.ndarray], name: str,
                pct: dict[str, float], use_sqrt: bool = True,
                xlim: tuple[float | None, float | None] = (None, None),
                ylim: tuple[float | None, float | None] = (None, None),
                colors: dict[str, str] | None = None,
                weighted: bool | None = None) -> Figure:
    """Two-panel Rietveld plot; returns the matplotlib Figure.

    Typical use, after read_gsas2_csv and prepare_data::

        data, phase_cols, error = read_gsas2_csv("sample.csv")
        if error is None:
            fig = create_plot(*prepare_data(data, phase_cols),
                              "Sample", pct={}, colors={})
    """
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(FIGURE_WIDTH, FIGURE_HEIGHT), dpi=110,
        gridspec_kw={"height_ratios": [4, 1]}, sharex=True)

    ax1.plot(theta, bkg, color=COLOR_BKG, ls="--", lw=1.5,
             label="_nolegend_", zorder=1)
    ax1.plot(theta, calc, color=COLOR_CALC, lw=LINEWIDTH_CALC,
             label="Calculated Fit", zorder=3)
    ax1.scatter(theta, obs, color=COLOR_OBS, s=MARKER_SIZE_OBS,
                label="Observed", zorder=2, lw=0)

    # Reflection tick rows below the pattern, one row per phase, in the
    # order of the legend. The intensity axis is scaled on the part of the
    # pattern the window actually shows, so cropping in 2theta does not
    # leave the figure scaled by a peak nobody sees.
    x_low, x_high = plot_window(theta, *xlim)
    visible = (theta >= x_low) & (theta <= x_high)
    shown = obs[visible] if visible.any() else obs
    ymin, ymax = shown.min(), shown.max()
    yrange = ymax - ymin
    base_y, step_y, tick_h = (ymin - 0.05 * yrange, 0.06 * yrange,
                              0.04 * yrange)
    # The upper margin clears the calculated curve as well: where the fit
    # overshoots, its peak rises above the observed scatter that set ymax,
    # and scaling on obs alone would clip it against the top border.
    calc_shown = calc[visible] if visible.any() else calc
    peak = max(ymax, float(np.nanmax(calc_shown))) if len(calc_shown) else ymax

    labels = {ph: phase_label(ph) for ph in phases}
    ordered = sorted(phases, key=lambda ph: labels[ph])
    tick_colors = phase_colors([labels[ph] for ph in ordered], colors)
    rows = 0
    for ph in ordered:
        # A reflection sits at a real angle, so a zero left in a phase column
        # is padding rather than a position and would draw a tick against the
        # left border.
        locs = phases[ph][phases[ph] > 0.1]
        if len(locs):
            y = base_y - rows * step_y
            ax1.vlines(locs, y, y + tick_h, colors=tick_colors[labels[ph]],
                       lw=LINEWIDTH_TICKS, zorder=4)
            rows += 1

    # 'quantity / unit' per IUPAC: the slash divides by the unit, so the axis
    # would carry pure numbers. Parentheses are needed only for a compound
    # unit, where 'I / counts s^-1' would read as (I/counts)*s^-1.
    # \mathrm keeps the name upright inside the radical. Italic is for symbols,
    # and an italic word is read as a product of variables.
    ylabel = (r"$\sqrt{\mathrm{Intensity}}$ / a.u." if use_sqrt
              else r"Intensity / a.u.")
    ax1.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL, labelpad=10)
    # No tick marks and no numbers on either intensity axis: labelleft does
    # both ends of that in one call, and unlike set_yticklabels([]) it does
    # not depend on a later tick_params resetting the tick length.
    ax1.tick_params(direction="in", top=False, right=False, left=False,
                    labelleft=False, width=1.5, length=6,
                    labelsize=FONT_SIZE_TICK)
    for spine in ax1.spines.values():
        spine.set_linewidth(LINEWIDTH_BORDER)
    ax1.spines["bottom"].set_visible(False)
    ax1.tick_params(bottom=False)
    # Default limits leave room for the tick rows below and a margin above;
    # both ends can be given explicitly, in the units actually drawn (sqrt
    # intensity when USE_SQRT).
    y_low, y_high = ylim
    ax1.set_ylim(bottom=(base_y - (max(rows, 1) - 0.5) * step_y
                         if y_low is None else float(y_low)),
                 top=(peak + 0.08 * yrange if y_high is None
                      else float(y_high)))

    # Legend: pattern entries, then one entry per detected phase.
    handles = [
        Line2D([0], [0], color=COLOR_OBS, marker="o", ls="", label=name),
        Line2D([0], [0], color=COLOR_CALC, lw=2, label="Calculated Fit"),
        Line2D([0], [0], color=COLOR_BKG, lw=2, ls="--", label="Background"),
    ]
    seen = set()
    for ph in ordered:
        lbl = labels[ph]
        if lbl in seen:
            continue  # two columns of the same phase share one entry
        seen.add(lbl)
        frac = phase_fraction(pct, lbl)
        text = f"{lbl} ({fmt_pct(frac)}%)" if frac > 0 else lbl
        handles.append(Line2D([0], [0], color=tick_colors[lbl], lw=3,
                              label=text))
    ax1.legend(handles=handles, loc="upper right", fontsize=FONT_SIZE_LEGEND,
               framealpha=1, edgecolor="black")

    orphans = sorted(k for k in pct
                     if not any(k in lbl.lower() for lbl in seen))
    if orphans:
        print("  ! no phase in this file matches metadata column(s): "
              + ", ".join(o + "_pct" for o in orphans))

    ax2.plot(theta, resid, color=COLOR_RESIDUALS, lw=LINEWIDTH_RESIDUALS)
    ax2.set_xlabel(r"$2\theta\:/\:^\circ$", fontsize=FONT_SIZE_LABEL)
    # Both panels carry a delta in I: the residual in intensity units, and
    # the same difference standardised by the uncertainty of the point.
    ax2.set_ylabel(r"$\Delta I\:/\:\sigma$"
                   if residual_column(weighted) == "diff/sigma"
                   else r"$\Delta I$ / a.u.", fontsize=FONT_SIZE_LABEL)
    ax2.set_xlim(x_low, x_high)
    # The residual axis loses its numbers too, and no line is drawn at zero.
    # Zero is held at the middle of the panel instead, by limits symmetric
    # about it: under autoscale it landed wherever the misfit happened to be
    # asymmetric, so the trace sat at a different height in every sample and
    # no two figures could be read side by side.
    ax2.set_ylim(*residual_limits(resid[visible] if visible.any() else resid,
                                  weighted))
    ax2.tick_params(direction="in", right=False, left=False, labelleft=False,
                    bottom=True, width=1.5, length=6,
                    labelsize=FONT_SIZE_TICK)
    for spine in ax2.spines.values():
        spine.set_linewidth(LINEWIDTH_BORDER)
    ax2.spines["top"].set_visible(False)
    ax2.tick_params(top=False)
    plt.subplots_adjust(hspace=0)
    return fig


def report_colors(phase_cols: list[str],
                  colors: dict[str, str] | None) -> None:
    """Print each detected phase with the colour it was drawn in.

    >>> report_colors(["Quartz hkl"], None)
      phases: Quartz #EE8031
    """
    labels = list(dict.fromkeys(sorted(phase_label(c) for c in phase_cols)))
    if not labels:
        print("  phases: none detected")
        return
    drawn = phase_colors(labels, colors)
    print("  phases: " + ", ".join(f"{lbl} {drawn[lbl]}" for lbl in labels))


def replot_file(csv_path: str | Path, metadata_file: str | Path,
                x_min: float | None = None, x_max: float | None = None,
                y_min: float | None = None, y_max: float | None = None,
                use_sqrt: bool = True,
                weighted: bool | None = None) -> tuple[Figure, str]:
    """Draw one file with the given window and toggles, without saving it.

    weighted defaults to the WEIGHTED_RESIDUALS setting, as everywhere else,
    so a caller who set it on the module gets the panel it asks for.
    Returns (figure, metadata_line): the line is the row to paste into the
    metadata file so that the batch run reproduces this 2theta window.
    Raises ValueError when the file cannot be read or drawn.

    Typical use, in section 4 of the notebook::

        fig, line = replot_file("sample.csv", "Samples_metadata.csv",
                                x_min=20, x_max=60)
        show_inline(fig)
        print(line)  # paste into the metadata file to keep this window
    """
    csv_path = Path(csv_path)
    data, phase_cols, error = read_gsas2_csv(csv_path, weighted=weighted)
    if error:
        raise ValueError(f"{csv_path.name}: {error}")
    name, pct, colors, _ = sample_info(load_metadata(metadata_file),
                                       csv_path.name, csv_path.stem)
    report_colors(phase_cols, colors)
    fig = create_plot(*prepare_data(data, phase_cols, use_sqrt=use_sqrt),
                      name, pct, use_sqrt=use_sqrt, xlim=(x_min, x_max),
                      ylim=(y_min, y_max), colors=colors, weighted=weighted)

    cells = [csv_path.name, name,
             "" if x_min is None else f"{float(x_min):g}",
             "" if x_max is None else f"{float(x_max):g}"]
    return fig, "filename;formula;x_min;x_max\n" + ";".join(cells)


# --- Batch driver ---------------------------------------------------------
def output_basename(stem: str, use_sqrt: bool = True,
                    weighted: bool | None = None) -> str:
    """Output file name for a sample stem, identical to the batch run.

    weighted defaults to the WEIGHTED_RESIDUALS setting; pass it explicitly
    (as section 4 does) so a figure saved with the other residual mode does
    not overwrite the batch file under the same name.

    >>> output_basename("sample1")
    'sample1_XRD_analysis_sqrt'
    >>> output_basename("sample1", weighted=False)
    'sample1_XRD_analysis_sqrt_unweighted'
    """
    weighted = WEIGHTED_RESIDUALS if weighted is None else weighted
    return (f"{stem}_XRD_analysis{'_sqrt' if use_sqrt else '_linear'}"
            f"{'' if weighted else '_unweighted'}")


def save_figure(fig: Figure, output_folder: str | Path, base: str) -> str:
    """Write base.pdf and base.png into output_folder/pdf and .../png.

    Sorting by extension keeps a large output folder tidy. Returns the base
    name, which the caller pairs with the 'pdf/' and 'png/' subfolders.

    >>> import tempfile, os
    >>> from matplotlib.figure import Figure
    >>> fig = Figure()
    >>> with tempfile.TemporaryDirectory() as d:
    ...     base = save_figure(fig, d, "sample")
    ...     sorted(os.listdir(d))
    ['pdf', 'png']
    """
    out_dir = Path(output_folder)
    (out_dir / "pdf").mkdir(parents=True, exist_ok=True)
    (out_dir / "png").mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "pdf" / f"{base}.pdf", bbox_inches="tight",
                facecolor="white")
    fig.savefig(out_dir / "png" / f"{base}.png", dpi=DPI_EXPORT,
                bbox_inches="tight", facecolor="white")
    return base


def sample_window(metadata_file: str | Path, filename: str
                  ) -> tuple[float | None, float | None]:
    """The (x_min, x_max) a sample's metadata row sets, each None if unset.

    Section 4 prefills its 2theta boxes with this, so a file opens on the
    same window the batch would draw instead of the full measured range.

    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as d:
    ...     p = Path(d) / "meta.csv"
    ...     with open(p, "w") as fh:
    ...         print("filename;x_min;x_max", file=fh)
    ...         print("a.csv;20;40", file=fh)
    ...     sample_window(p, "a.csv")
    Metadata loaded for 1 sample(s).
    (20.0, 40.0)
    """
    _, _, _, window = sample_info(load_metadata(metadata_file), filename,
                                  filename)
    return window


def show_inline(fig: Figure) -> None:
    """Render fig under the running cell at PREVIEW_WIDTH_PX; no-op off-kernel.

    A saved PNG shown at an explicit width renders the same in classic
    Jupyter, JupyterLab and VS Code, where relying on the active backend to
    honour ``plt.show()`` does not. Outside a kernel (a plain script, the
    pytest run) there is no display, so this returns without drawing.

    Typical use, inside a notebook cell, after create_plot::

        fig = create_plot(theta, obs, calc, bkg, resid, phases, name, pct)
        show_inline(fig)
    """
    try:
        from IPython import get_ipython
        from IPython.display import Image, display
    except ImportError:
        return
    if get_ipython() is None:
        return
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight",
                facecolor="white")
    display(Image(data=buf.getvalue(), width=PREVIEW_WIDTH_PX))


def process_folder(data_folder: str | Path,
                   metadata_file: str | Path,
                   output_folder: str | Path,
                   use_sqrt: bool = True, show: bool = True
                   ) -> list[tuple[str, str, str]]:
    """Plot every GSAS-II CSV export in data_folder; save PDF + PNG.

    Typical use, section 3 of the notebook::

        results = process_folder("data", "Samples_metadata.csv", "output")
        failed = [name for name, status, _ in results if status == "error"]
    """
    data_dir = Path(data_folder)
    data_dir.mkdir(parents=True, exist_ok=True)  # first run: created empty, ready for your files

    meta_name = Path(metadata_file).name
    files = sorted(f for f in data_dir.glob("*.csv") if f.name != meta_name)
    if not files:
        print(f"No CSV files found in '{data_dir}': add your GSAS-II exports and re-run.")
        return []

    meta_df = load_metadata(metadata_file)
    out_dir = Path(output_folder)
    out_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for f in files:
        print(f"\n{f.name}")
        # One bad file must not end the batch, whether it fails in the
        # parser, in the drawing or on the way to disk.
        try:
            data, phase_cols, error = read_gsas2_csv(f)
            if error:
                raise ValueError(error)

            name, pct, colors, window = sample_info(meta_df, f.name, f.stem)
            report_colors(phase_cols, colors)
            theta, obs, calc, bkg, resid, phases = prepare_data(
                data, phase_cols, use_sqrt=use_sqrt)
            # Show the 2theta window actually drawn and where it came from, so
            # a limit that did not land (name mismatch, blank cell) is visible
            # at a glance rather than read off the figure.
            low, high = plot_window(theta, *window)
            source = ("metadata" if window != (None, None)
                      else "settings" if (PLOT_X_MIN, PLOT_X_MAX) != (None, None)
                      else "auto")
            print(f"  2theta window: {low:g} to {high:g} ({source})")
            fig = create_plot(theta, obs, calc, bkg, resid, phases, name, pct,
                              use_sqrt=use_sqrt, xlim=window, colors=colors)

            # The residual suffix is added only when it is not the default,
            # so figures made before the setting existed keep their names.
            base = output_basename(f.stem, use_sqrt, WEIGHTED_RESIDUALS)
            save_figure(fig, out_dir, base)
            results.append((f.name, "ok", base))

            # The figure goes out under the header that names its file and
            # above the line that names its output, so the notebook reads as
            # one block per sample while the run proceeds.
            if show:
                show_inline(fig)
            print(f"  saved pdf/{base}.pdf and png/{base}.png")
            plt.close(fig)  # release memory between files
        except Exception as e:
            print(f"  FAILED: {e}")
            results.append((f.name, "error", str(e)))
            plt.close("all")

    # The summary repeats every failure. Between six figures a single line in
    # the middle of the run is easy to scroll past.
    failed = [name for name, status, _ in results if status == "error"]
    print(f"\n{'=' * 60}")
    print(f"{len(files) - len(failed)} of {len(files)} file(s) plotted, "
          f"output in '{out_dir}'.")
    if failed:
        print(f"FAILED ({len(failed)}): " + ", ".join(failed))
    return results
