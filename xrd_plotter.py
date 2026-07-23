"""Rietveld plots from the CSV a GSAS-II Rietveld plot saves.

Parsing, data preparation, plotting and the batch driver. The notebook and
the tests both import this module, so a figure drawn by hand matches the one
the batch writes.

Appearance is set by the constants below. Override them on the module, not
on a copy, so the functions see the change:

    import xrd_plotter as xp
    xp.PLOT_X_MIN, xp.PLOT_X_MAX = 13, 85
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D

# --- Plot appearance --------------------------------------------------
# 2theta window (deg). None takes the measured range of each file, so a
# pattern is never cropped without being asked. Set both to fix every
# figure to the same window, or give x_min and x_max per sample in the
# metadata file, which wins over these.
PLOT_X_MIN, PLOT_X_MAX = None, None
FIGURE_WIDTH, FIGURE_HEIGHT = 12, 10   # inches
DPI_EXPORT = 600                       # PNG export resolution

# Lower panel: True draws diff/sigma, the residual divided by the standard
# deviation of the point, so a well-fitted pattern stays inside a band of a
# few units. False draws the raw diff in counts, where the tall reflections
# dominate. Files carry the suffix '_counts' when this is False.
WEIGHTED_RESIDUALS = True

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

# --- Phase columns ----------------------------------------------------
# Headers GSAS-II writes itself. Everything else that holds only a few
# values is taken to be the reflection positions of a phase, so phases are
# found by their own names instead of a keyword list.
NON_PHASE_COLUMNS = frozenset({
    "used", "obs", "calc", "bkg", "diff", "diff/sigma", "weight", "weights",
    "sig", "sigma", "tick-pos", "tick pos", "axis-limits", "axis limits",
    "excluded", "gof",
})
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


def fmt_pct(p):
    """Phase fraction label: one decimal, trailing '.0' dropped (62.0 -> 62)."""
    s = f"{p:.1f}"
    return s[:-2] if s.endswith(".0") else s


def residual_column(weighted=None):
    """Header the lower panel is drawn from, per WEIGHTED_RESIDUALS."""
    if weighted is None:
        weighted = WEIGHTED_RESIDUALS
    return "diff/sigma" if weighted else "diff"


def plot_window(theta, x_min=None, x_max=None):
    """2theta limits: the per-sample pair, then the constants, then the data."""
    low = x_min if x_min is not None else PLOT_X_MIN
    high = x_max if x_max is not None else PLOT_X_MAX
    return (float(np.nanmin(theta)) if low is None else float(low),
            float(np.nanmax(theta)) if high is None else float(high))


def phase_label(column):
    """Legend name for a phase column: 'Phase 1 hkl' -> 'Phase 1'."""
    label = column.strip()
    if label.lower().endswith("hkl"):
        label = label[:-3].strip()
    override = longest_match(PHASE_LABELS, label)
    return override if override else label[:1].upper() + label[1:]


def longest_match(mapping, label):
    """Value whose key is the longest fragment of label, or None.

    The longest key wins so that a specific column ('phase 1_color') beats a
    general one ('phase_color') instead of resolving on dictionary order.
    """
    keys = [k for k in (mapping or {}) if k and k in label.lower()]
    return mapping[max(keys, key=len)] if keys else None


def phase_colors(ordered_labels, overrides=None):
    """Tick colour per phase: the metadata, then PHASE_COLORS, then the cycle.

    Both mappings are keyed by a fragment of the legend name. Repeated
    labels share one colour and consume one slot of the cycle, matching the
    single legend entry they get.
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


def phase_fraction(pct, label):
    """Percentage for one phase, matched on the metadata key as a fragment."""
    hits = sorted(k for k in pct if k in label.lower())
    if len(hits) > 1:
        print(f"  ! metadata columns {', '.join(h + '_pct' for h in hits)} "
              f"all match '{label}', no percentage printed.")
        return 0.0
    return pct[hits[0]] if hits else 0.0


# --- Parsing -----------------------------------------------------------
def read_gsas2_csv(csv_path, weighted=None):
    """Read one CSV saved from a GSAS-II Rietveld plot.

    Returns (data, phase_cols, error): on success error is None, on failure
    data is None and error holds the reason. Non-fatal issues are printed
    as warnings. 'weighted' overrides WEIGHTED_RESIDUALS for this call.
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

        # 2theta column: header contains '2theta' or starts like 'x,'.
        theta = next((col for cl, col in col_lower_map.items()
                      if "2theta" in cl or "x," in cl), None)
        if theta is None:
            return None, None, "2theta column not found"
        data["x"] = clean(df[theta])
        if np.all(np.isnan(data["x"])):
            return None, None, "no valid data in 2theta column"

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
        resid_col = col_lower_map.get(wanted)
        if resid_col is None:
            return None, None, f"residual column '{wanted}' not found"
        data["resid"] = clean(df[resid_col])
        n_nan = int(np.isnan(data["resid"]).sum())
        if n_nan:
            warnings.append(f"column '{wanted}': {n_nan} non-numeric values -> NaN")

        # Per-phase reflection-position columns: whatever is left over and
        # sparse enough to be a reflection list rather than a data column.
        phase_cols = []
        for cl, col in col_lower_map.items():
            # pandas renames a repeated header 'tick-pos' to 'tick-pos.1',
            # which the blocklist would otherwise miss.
            head, _, tail = cl.rpartition(".")
            base = head if head and tail.isdigit() else cl
            if col == theta or base in NON_PHASE_COLUMNS:
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

        if warnings:
            print("\n".join(f"  ! {w}" for w in warnings))
        return data, phase_cols, None

    except Exception as e:  # isolate unreadable files, keep the batch running
        return None, None, f"read error: {e}"


def load_metadata(metadata_path):
    """Read the (private) sample metadata CSV, indexed by filename.

    Expected columns: 'filename' (or 'file'), optional 'formula', one
    '<phase>_pct' and one '<phase>_color' column per phase, and an
    optional 'x_min'/'x_max' pair. Returns an empty DataFrame when the
    file is absent. The contents are never displayed by this notebook.
    """
    metadata_path = Path(metadata_path)
    if not metadata_path.is_file():
        return pd.DataFrame()
    try:
        df = pd.read_csv(metadata_path, sep=";", encoding="utf-8", dtype=str)
        if len(df.columns) <= 1:
            df = pd.read_csv(metadata_path, sep=",", encoding="utf-8", dtype=str)
    except (pd.errors.ParserError, UnicodeDecodeError):
        df = pd.read_csv(metadata_path, sep=",", encoding="latin-1", dtype=str)

    df.columns = df.columns.str.strip().str.lower()
    fname_col = next((c for c in df.columns if c in ("filename", "file")), None)
    if fname_col is None:
        print("  ! metadata needs a 'filename' or 'file' column, ignored.")
        return pd.DataFrame()
    df["filename"] = df[fname_col].astype(str).str.strip()
    df = df.drop_duplicates("filename", keep="first").set_index("filename")
    print(f"Metadata loaded for {len(df)} sample(s).")
    return df


def to_number(value):
    """One metadata cell as a float, or None when it is blank or not a number."""
    try:
        number = float(str(value).replace(",", "."))  # decimal comma
    except (ValueError, TypeError):
        return None
    return None if np.isnan(number) else number


def to_color(value):
    """One metadata cell as a colour, or None when blank or not a colour."""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    if not is_color_like(text):
        print(f"  ! '{text}' is not a colour, the cycle is used instead.")
        return None
    return text


def metadata_keys(columns, suffix):
    """Column names ending in suffix, keyed by the phase fragment before it.

    A column named exactly like the suffix carries no phase name and is
    skipped, since its empty key would match every phase.
    """
    return {c[:-len(suffix)].replace("_", " ").strip(): c for c in columns
            if c.endswith(suffix) and len(c) > len(suffix)}


def sample_info(meta_df, filename, default_name):
    """Name, fractions, colours and 2theta window for one file (empty-safe).

    Every '<phase>_pct' and '<phase>_color' column becomes one entry, keyed
    by the part before the suffix with underscores as spaces, so
    'phase_1_color' reaches the phase named 'Phase 1' in the figure.
    'x_min' and 'x_max' are returned as given, for plot_window to apply.
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
def prepare_data(data, phase_cols, use_sqrt=True):
    """Mask invalid points; optionally apply the display-only sqrt transform."""
    mask = ~(np.isnan(data["x"]) | np.isnan(data["obs"]))
    if not mask.any():
        raise ValueError("no point has both a 2theta and an obs value")
    x = data["x"][mask]

    if use_sqrt:
        obs = np.sqrt(np.abs(data["obs"][mask]))
        calc = np.sqrt(np.abs(data["calc"][mask]))
        bkg = np.sqrt(np.abs(data["bkg"][mask]))
    else:
        obs, calc, bkg = (data[k][mask] for k in ("obs", "calc", "bkg"))

    resid = data["resid"][mask]

    # Phase columns hold independent reflection positions (different length
    # from the pattern): drop their NaN padding individually.
    phases = {k: data[k][~np.isnan(data[k])] for k in phase_cols}
    return x, obs, calc, bkg, resid, phases


# --- Plotting -------------------------------------------------------------
def create_plot(theta, obs, calc, bkg, resid, phases, name, pct,
                use_sqrt=True, xlim=(None, None), ylim=(None, None),
                colors=None, weighted=None):
    """Two-panel Rietveld plot; returns the matplotlib Figure."""
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

    labels = {ph: phase_label(ph) for ph in phases}
    ordered = sorted(phases, key=lambda ph: labels[ph])
    tick_colors = phase_colors([labels[ph] for ph in ordered], colors)
    rows = 0
    for ph in ordered:
        locs = phases[ph][phases[ph] > 0.1]
        if len(locs):
            y = base_y - rows * step_y
            ax1.vlines(locs, y, y + tick_h, colors=tick_colors[labels[ph]],
                       lw=LINEWIDTH_TICKS, zorder=4)
            rows += 1

    ylabel = r"$\sqrt{Counts}$ / (a.u.)" if use_sqrt else r"Counts / (a.u.)"
    ax1.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL, labelpad=10)
    ax1.set_yticklabels([])
    ax1.tick_params(axis="y", length=0)
    ax1.tick_params(direction="in", top=False, right=False, left=True,
                    width=1.5, length=6, labelsize=FONT_SIZE_TICK)
    for spine in ax1.spines.values():
        spine.set_linewidth(LINEWIDTH_BORDER)
    ax1.spines["bottom"].set_visible(False)
    ax1.tick_params(bottom=False)
    # Default limits leave room for the tick rows below and a margin above;
    # both ends can be given explicitly, in the units actually drawn (sqrt
    # counts when USE_SQRT).
    y_low, y_high = ylim
    ax1.set_ylim(bottom=(base_y - (max(rows, 1) - 0.5) * step_y
                         if y_low is None else float(y_low)),
                 top=(ymax + 0.05 * yrange if y_high is None
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
    ax2.axhline(0, color="black", lw=1.5)
    ax2.set_xlabel(r"2$\theta$ / ($^\circ$)", fontsize=FONT_SIZE_LABEL)
    ax2.set_ylabel(r"diff/$\sigma$" if residual_column(weighted) == "diff/sigma"
                   else r"diff / (a.u.)", fontsize=FONT_SIZE_LABEL)
    ax2.set_xlim(x_low, x_high)
    ax2.tick_params(direction="in", right=False, left=True, bottom=True,
                    width=1.5, length=6, labelsize=FONT_SIZE_TICK)
    for spine in ax2.spines.values():
        spine.set_linewidth(LINEWIDTH_BORDER)
    ax2.spines["top"].set_visible(False)
    ax2.tick_params(top=False)
    plt.subplots_adjust(hspace=0)
    return fig


def report_colors(phase_cols, colors):
    """Print each detected phase with the colour it was drawn in."""
    labels = list(dict.fromkeys(sorted(phase_label(c) for c in phase_cols)))
    if not labels:
        print("  phases: none detected")
        return
    drawn = phase_colors(labels, colors)
    print("  phases: " + ", ".join(f"{lbl} {drawn[lbl]}" for lbl in labels))


def replot_file(csv_path, metadata_file, x_min=None, x_max=None,
                y_min=None, y_max=None, use_sqrt=True, weighted=True):
    """Draw one file with the given window and toggles, without saving it.

    Returns (figure, metadata_line): the line is the row to paste into the
    metadata file so that the batch run reproduces this 2theta window.
    Raises ValueError when the file cannot be read or drawn.
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
def process_folder(data_folder, metadata_file, output_folder,
                   use_sqrt=True, show=True):
    """Plot every GSAS-II CSV export in data_folder; save PDF + PNG."""
    data_dir = Path(data_folder)
    data_dir.mkdir(exist_ok=True)  # first run: created empty, ready for your files

    meta_name = Path(metadata_file).name
    files = sorted(f for f in data_dir.glob("*.csv") if f.name != meta_name)
    if not files:
        print(f"No CSV files found in '{data_dir}': add your GSAS-II exports and re-run.")
        return []

    meta_df = load_metadata(metadata_file)
    out_dir = Path(output_folder)
    out_dir.mkdir(exist_ok=True)

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
            fig = create_plot(theta, obs, calc, bkg, resid, phases, name, pct,
                              use_sqrt=use_sqrt, xlim=window, colors=colors)

            # The residual suffix is added only when it is not the default,
            # so figures made before the setting existed keep their names.
            base = (f"{f.stem}_XRD_analysis{'_sqrt' if use_sqrt else '_linear'}"
                    f"{'' if WEIGHTED_RESIDUALS else '_counts'}")
            fig.savefig(out_dir / f"{base}.pdf", bbox_inches="tight",
                        facecolor="white")
            fig.savefig(out_dir / f"{base}.png", dpi=DPI_EXPORT,
                        bbox_inches="tight", facecolor="white")
            results.append((f.name, "ok", base))

            # The figure goes out under the header that names its file and
            # above the line that names its output, so the notebook reads as
            # one block per sample while the run proceeds. A file backend has
            # no window, so it is skipped rather than warned about.
            if show and plt.get_backend().lower() != "agg":
                plt.show()
            print(f"  saved {out_dir / base}.pdf and .png")
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
