"""Validation suite for xrd_plotter.

Every input is built here from an analytic pattern with a fixed seed, so the
suite runs on any machine, including a fresh clone with an empty data folder,
and no measurement is ever read. Run it with `pytest -q`, or through section 2
of the notebook, which calls pytest for you.
"""
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")  # a headless run must not need a display

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import xrd_plotter as xp

N = 500
CONFIG = ("PLOT_X_MIN", "PLOT_X_MAX", "WEIGHTED_RESIDUALS", "PHASE_LABELS",
          "PHASE_COLORS", "PHASE_MAX_FILL", "NON_PHASE_COLUMNS", "DPI_EXPORT")


@pytest.fixture(autouse=True)
def config_restored():
    """A test that changes a module constant must not leak it to the next."""
    saved = {name: getattr(xp, name) for name in CONFIG}
    # The export resolution is not what any test checks, and 600 dpi rasters
    # cost more than everything else here put together.
    xp.DPI_EXPORT = 100
    yield
    for name, value in saved.items():
        setattr(xp, name, value)
    plt.close("all")


@pytest.fixture(scope="module")
def synth(tmp_path_factory):
    """The analytic pattern and every CSV the suite reads."""
    tmp = tmp_path_factory.mktemp("exports")
    rng = np.random.default_rng(0)

    x = np.linspace(10.0, 90.0, N)
    calc = 1000.0 * np.exp(-((x - 30.0) ** 2) / 2.0) + 50.0
    obs = calc + rng.normal(0.0, 5.0, N)
    bkg = np.full(N, 50.0)
    ds = (obs - calc) / 5.0

    df = pd.DataFrame({"2theta": x, "Obs": obs, "Calc": calc, "Bkg": bkg,
                       "diff/sigma": ds})
    df["Phase 1 hkl"] = pd.Series([30.0, 35.0, 50.0, 60.0])
    df["Phase 2 hkl"] = pd.Series([32.0, 38.0, 52.0])

    # Variant A: semicolon separator + decimal commas (locale export).
    a = tmp / "sample_A.csv"
    df.to_csv(a, sep=";", index=False, decimal=",")
    # Variant B: comma separator, decimal points, 'x, deg' header.
    b = tmp / "sample_B.csv"
    df.rename(columns={"2theta": "x, deg"}).to_csv(b, index=False)

    (tmp / "corrupt.csv").write_bytes(b"\x00\x01\x02 not a csv \xff")
    (tmp / "no_theta.csv").write_text("A;B\n1;2\n")
    df[["2theta", "Obs"]].to_csv(tmp / "only_obs.csv", sep=";", index=False)
    df[["2theta", "Obs", "diff/sigma"]].to_csv(tmp / "no_calc.csv", sep=";",
                                               index=False)
    df.assign(Obs="n/a").to_csv(tmp / "obs_all_text.csv", sep=";", index=False)
    df.assign(**{"diff/sigma": "n/a"}).to_csv(tmp / "resid_all_text.csv",
                                              sep=";", index=False)

    # One row two fields too long: the C parser cannot tokenise the file.
    lines = a.read_text().splitlines()
    lines[5] += ";99;99"
    (tmp / "ragged.csv").write_text("\n".join(lines) + "\n")

    # A complete export, in the column order and with the header GSAS-II
    # writes, and with phase names carrying no keyword.
    full = pd.DataFrame({"used": np.ones(N), "x, 2theta (deg)": x, "obs": obs,
                         "calc": calc, "bkg": bkg, "diff": obs - calc})
    full["Alpha hkl"] = pd.Series([30.0, 35.0, 50.0, 60.0])
    full["Beta hkl"] = pd.Series([32.0, 38.0, 52.0])
    full["tick-pos"] = pd.Series([-0.5])
    full["diff/sigma"] = ds
    full["Axis-limits"] = pd.Series([10.0, 90.0])
    full.to_csv(tmp / "full_export.csv", sep=";", index=False)

    (tmp / "repeated_header.csv").write_text(
        "2theta;Obs;diff/sigma;tick-pos;tick-pos\n"
        "10.0;1;0;-0.5;-0.5\n20.0;1;0;;\n30.0;1;0;;\n")

    dense = df.copy()
    dense["Wide hkl"] = pd.Series(np.linspace(10.0, 90.0, N - 100))
    dense.to_csv(tmp / "dense_column.csv", sep=";", index=False)

    (tmp / "Samples_metadata.csv").write_text(
        "filename;formula;phase_1_pct;phase_2_pct\n"
        "sample_A.csv;Sample A (synthetic);60;40\n")
    (tmp / "metadata_colour.csv").write_text(
        "filename;phase_2_color;phase_1_color\n"
        "sample_A.csv;#123456;not a colour\n")
    (tmp / "metadata_variants.csv").write_text(
        "filename;formula;phase_1_pct;phase_9_pct\nsample_A.csv;;60,5;10\n")
    (tmp / "metadata_window.csv").write_text(
        "filename;x_min;x_max\nsample_A.csv;20;60,5\n")

    return SimpleNamespace(tmp=tmp, x=x, calc=calc, obs=obs, bkg=bkg, ds=ds,
                           df=df, a=a, b=b, meta=tmp / "Samples_metadata.csv",
                           full=tmp / "full_export.csv")


@pytest.fixture
def drawn(synth):
    """sample_A parsed and prepared, with its metadata, ready to plot."""
    data, phase_cols, error = xp.read_gsas2_csv(synth.a)
    assert error is None
    name, pct, colors, _ = xp.sample_info(xp.load_metadata(synth.meta),
                                          "sample_A.csv", "sample_A")
    theta, o, c, b, r, phases = xp.prepare_data(data, phase_cols, use_sqrt=True)
    return SimpleNamespace(args=(theta, o, c, b, r, phases), theta=theta,
                           phases=phases, name=name, pct=pct, colors=colors)


def legend_of(fig):
    """(texts, colours) of the phase entries, the ones after the three fixed."""
    legend = fig.axes[0].get_legend()
    return ([t.get_text() for t in legend.get_texts()][3:],
            [h.get_color() for h in legend.legend_handles[3:]])


# 1 and 2. Bit-exact parsing of both locale variants, and phase detection.
@pytest.mark.parametrize("variant", ["a", "b"])
def test_parse_is_bit_exact(synth, variant):
    data, phase_cols, error = xp.read_gsas2_csv(getattr(synth, variant))
    assert error is None, error
    assert np.array_equal(data["x"], synth.x), "2theta round trip not exact"
    assert np.array_equal(data["obs"], synth.obs), "Obs round trip not exact"
    assert np.array_equal(data["resid"], synth.ds), "diff/sigma not exact"
    assert sorted(phase_cols) == ["Phase 1 hkl", "Phase 2 hkl"]


# Degraded inputs must be isolated with a reason, not raise.
@pytest.mark.parametrize("name, reason", [
    ("corrupt.csv", "empty CSV file"),
    ("no_theta.csv", "2theta column not found"),
    ("only_obs.csv", "residual column 'diff/sigma' not found"),
    ("obs_all_text.csv", "no valid data in the obs column"),
    ("resid_all_text.csv", "no valid data in the 'diff/sigma' column"),
])
def test_degraded_files_are_isolated(synth, name, reason):
    data, _, error = xp.read_gsas2_csv(synth.tmp / name)
    assert data is None and error == reason, error


# 2b. The native GSAS-II export, whose header is one name too long.
def test_shifted_export_is_refused(tmp_path):
    # GSAS-II writes 11 header names and 10 fields per row, so every column
    # is read under its left neighbour's name: 'used' takes the angles and
    # the 2theta header takes the intensities. The figure would still draw.
    x = np.linspace(10.0, 90.0, N)
    obs = 1000.0 * np.exp(-((x - 30.0) ** 2) / 2.0) + 50.0
    head = ("used;x, 2theta (deg);obs;calc;bkg;diff;Alpha;Beta;tick-pos;"
            "diff/sigma;Axis-limits")
    rows = [f"{x[i]:.6f};{obs[i]:.6f};{obs[i]:.6f};50;{obs[i] - 50:.6f};"
            f"{'12.5' if not i else ''};{'30.0' if not i else ''};"
            f"{'-0.5' if not i else ''};0.1;{'10' if not i else ''}"
            for i in range(N)]
    shifted = tmp_path / "native_export.csv"
    shifted.write_text(head + "\n" + "\n".join(rows) + "\n")

    data, _, error = xp.read_gsas2_csv(shifted)
    assert data is None, "a shifted export must not be drawn"
    assert "do not line up" in error, error


def test_an_axis_outside_0_to_180_degrees_is_refused(tmp_path):
    # A shift can land a smooth monotone column, the background for one,
    # under the 2theta header. Monotonicity cannot see that; degrees can.
    x = np.linspace(400.0, 40.0, 300)  # a falling background, not angles
    pd.DataFrame({"x, 2theta (deg)": x, "Obs": np.linspace(1.0, 2.0, 300),
                  "diff/sigma": np.zeros(300)}).to_csv(tmp_path / "wide.csv",
                                                       sep=";", index=False)
    data, _, error = xp.read_gsas2_csv(tmp_path / "wide.csv")
    assert data is None, "a column of intensities was drawn as an axis"
    assert "0 to 180 degrees" in error, error


def test_a_descending_2theta_axis_is_still_accepted(tmp_path):
    # The guard rejects a column that turns, not one that runs backwards.
    x = np.linspace(90.0, 10.0, N)
    df = pd.DataFrame({"2theta": x, "Obs": np.linspace(1.0, 2.0, N),
                       "diff/sigma": 0.0})
    csv = tmp_path / "descending.csv"
    df.to_csv(csv, sep=";", index=False)
    data, _, error = xp.read_gsas2_csv(csv)
    assert error is None, error
    assert np.array_equal(data["x"], x)


def test_header_matching(tmp_path):
    """The angle takes any spelling; every other column is matched whole."""
    for header in ("x, 2theta (deg)", "2theta", "2Theta (deg)", "X, 2THETA",
                   "x, deg", "angle, 2theta", "2theta_deg", "2θ",
                   "2 theta (deg)", "two-theta", "x"):
        assert xp.is_theta_header(header), header
    for header in ("theta", "angle", "obs", "intensity", "Max, phase",
                   "tick-pos"):
        assert not xp.is_theta_header(header), header

    # End to end: 'obs' is matched whole, so a renamed one is drawn flat, and
    # a phase whose header merely holds a comma never becomes the axis.
    x = np.linspace(10.0, 90.0, 300)
    obs = np.linspace(100.0, 200.0, 300)
    for header, read in ((" OBS ", True), ("observed", False)):
        csv = tmp_path / f"{header.strip()}.csv"
        pd.DataFrame({"Max, phase": pd.Series([30.0]).reindex(range(300)),
                      "x, 2theta (deg)": x, header: obs,
                      "diff/sigma": np.zeros(300)}).to_csv(csv, sep=";",
                                                           index=False)
        data, phase_cols, error = xp.read_gsas2_csv(csv)
        assert error is None, error
        assert np.array_equal(data["x"], x), "a phase column was read as 2theta"
        assert np.any(data["obs"] != 0) is np.True_ or not read, header
        assert phase_cols == ["Max, phase"], phase_cols


def test_is_monotonic():
    """Distance travelled, not the number of turns.

    Counting turns accepted two files that draw a plausible figure from the
    wrong column: one backward jump scores 1, which fits inside any tolerance
    on a long file, and a constant column scores 0.
    """
    rising = np.arange(500.0)
    assert xp.is_monotonic(rising)
    assert xp.is_monotonic(rising[::-1]), "a descending axis is still an axis"
    assert xp.is_monotonic(np.repeat(np.linspace(10.0, 90.0, 250), 2)), (
        "repeated points travel nowhere and must not count as a turn")
    rising[100], rising[101] = rising[101], rising[100]  # one swapped pair
    assert xp.is_monotonic(rising), "one glitch must not reject a real axis"

    assert not xp.is_monotonic(np.tile([1.0, 2.0], 250)), "intensities accepted"
    scan = np.linspace(10.0, 90.0, 500)
    assert not xp.is_monotonic(np.concatenate([scan, scan])), (
        "two scans pasted into one file accepted")
    assert not xp.is_monotonic(np.ones(500)), (
        "a constant column, such as a mask of ones, accepted as an axis")


def test_missing_curves_are_filled_with_zeros(synth):
    data, _, error = xp.read_gsas2_csv(synth.tmp / "no_calc.csv")
    assert error is None
    assert np.all(data["calc"] == 0.0) and np.all(data["bkg"] == 0.0)


# 3. The fallback parser keeps the rows it can read.
def test_ragged_file_loses_one_row_only(synth):
    data, phase_cols, error = xp.read_gsas2_csv(synth.tmp / "ragged.csv")
    assert error is None, error
    assert np.array_equal(data["obs"], np.delete(synth.obs, 4)), "not exact"
    assert sorted(phase_cols) == ["Phase 1 hkl", "Phase 2 hkl"]


# 4. Phase detection against a complete export.
def test_only_phases_are_detected_in_a_full_export(synth):
    data, phase_cols, error = xp.read_gsas2_csv(synth.full)
    assert error is None, error
    assert sorted(phase_cols) == ["Alpha hkl", "Beta hkl"], phase_cols


def test_repeated_header_is_not_a_phase(synth):
    _, phase_cols, error = xp.read_gsas2_csv(synth.tmp / "repeated_header.csv")
    assert error is None and phase_cols == [], phase_cols


def test_fit_statistics_are_not_read_as_phases(tmp_path):
    """A statistic beside the pattern is sparse, so it would get a tick row.

    Blocklist entries are folded, so 'rw' covers 'Rw', 'Rw%' and 'Rw / %'.
    'RP' names a phase and is deliberately not blocked, and
    a header added to the blocklist on the module takes effect at once, which
    is why the folded set cannot be built at import.
    """
    x = np.linspace(10.0, 90.0, 300)

    def phases_of(header):
        csv = tmp_path / "stat.csv"
        pd.DataFrame({"x, 2theta (deg)": x, "obs": np.linspace(1.0, 2.0, 300),
                      "diff/sigma": np.zeros(300),
                      header: pd.Series([30.0]).reindex(range(300))}).to_csv(
            csv, sep=";", index=False)
        _, phase_cols, error = xp.read_gsas2_csv(csv)
        assert error is None, error
        return phase_cols

    for header in ("Rw / %", "Rw%", "rw", "Rwp / %", "chi2", "GOF", "Tick-Pos"):
        assert phases_of(header) == [], header
    for header in ("Rp", "Rw phase", "Rutile"):
        assert phases_of(header) == [header], header

    assert phases_of("Scan note") == ["Scan note"]
    xp.NON_PHASE_COLUMNS = xp.NON_PHASE_COLUMNS | {"scan note"}
    assert phases_of("Scan note") == [], "the blocklist addition was ignored"


def test_column_too_full_is_reported_not_drawn(synth, capsys):
    _, phase_cols, error = xp.read_gsas2_csv(synth.tmp / "dense_column.csv")
    assert error is None and "Wide hkl" not in phase_cols, phase_cols
    assert "too many for a reflection list" in capsys.readouterr().out


# 5. Metadata reaches the legend.
def test_metadata_reaches_the_legend(synth, drawn):
    assert drawn.name == "Sample A (synthetic)"
    assert drawn.pct == {"phase 1": 60.0, "phase 2": 40.0}, drawn.pct
    assert drawn.colors == {}
    fig = xp.create_plot(*drawn.args, drawn.name, drawn.pct, use_sqrt=True)
    texts, _ = legend_of(fig)
    assert texts == ["Phase 1 (60%)", "Phase 2 (40%)"], texts


# 6. Order, colours and rows.
def test_phase_order_and_colour_cycle(drawn):
    phases = dict(drawn.phases, **{"Phase 3 hkl": np.array([40.0, 45.0])})
    fig = xp.create_plot(*drawn.args[:5], phases, drawn.name, drawn.pct)
    texts, colours = legend_of(fig)
    assert texts == ["Phase 1 (60%)", "Phase 2 (40%)", "Phase 3"], texts
    assert colours == list(xp.PHASE_COLOR_CYCLE[:3])
    rows = [c for c in fig.axes[0].collections if hasattr(c, "get_segments")]
    assert len(rows) == 3, "one tick row per phase"


def test_repeated_label_takes_one_entry_and_one_colour(drawn):
    twin = dict(drawn.phases, **{"Phase 1": drawn.phases["Phase 1 hkl"]})
    fig = xp.create_plot(*drawn.args[:5], twin, drawn.name, drawn.pct)
    texts, colours = legend_of(fig)
    assert texts == ["Phase 1 (60%)", "Phase 2 (40%)"], texts
    assert colours == list(xp.PHASE_COLOR_CYCLE[:2]), "a label ate a colour"


def test_phase_label_capitalises_and_strips_hkl():
    assert xp.phase_label("Phase 1 hkl") == "Phase 1"
    assert xp.phase_label("alpha") == "Alpha"


# 7. Colours pinned from the private metadata file.
def test_metadata_colour_follows_its_phase(synth, drawn):
    _, _, colors, _ = xp.sample_info(
        xp.load_metadata(synth.tmp / "metadata_colour.csv"), "sample_A.csv",
        "sample_A")
    assert colors == {"phase 2": "#123456"}, "the invalid cell must be dropped"

    solo = {"Phase 2 hkl": drawn.phases["Phase 2 hkl"]}
    fig = xp.create_plot(*drawn.args[:5], solo, drawn.name, drawn.pct,
                         colors=colors)
    assert legend_of(fig)[1] == ["#123456"], "a lone phase lost its colour"

    fig = xp.create_plot(*drawn.args, drawn.name, drawn.pct, colors=colors)
    assert legend_of(fig)[1] == [xp.PHASE_COLOR_CYCLE[0], "#123456"]


def test_longest_key_wins_and_nameless_column_is_skipped():
    assert xp.longest_match({"phase": "#000000", "phase 1": "#ffffff"},
                            "Phase 1") == "#ffffff"
    assert xp.metadata_keys(["_pct", "phase_1_pct"], "_pct") == {
        "phase 1": "phase_1_pct"}


# 8. The two dictionaries in the module.
def test_module_label_and_colour_dictionaries(drawn):
    xp.PHASE_LABELS = {"phase 1": "Alpha"}
    xp.PHASE_COLORS = {"alpha": "#654321"}
    assert xp.phase_label("Phase 1 hkl") == "Alpha"
    fig = xp.create_plot(*drawn.args, drawn.name, drawn.pct)
    texts, colours = legend_of(fig)
    assert texts == ["Alpha", "Phase 2 (40%)"], texts
    assert colours == ["#654321", xp.PHASE_COLOR_CYCLE[0]], "the pin missed"


# 9. Metadata quirks.
def test_metadata_variants(synth, drawn):
    name, pct, _, _ = xp.sample_info(
        xp.load_metadata(synth.tmp / "metadata_variants.csv"), "sample_A.csv",
        "sample_A")
    assert name == "sample_A", "an empty formula cell must not print as nan"
    assert pct == {"phase 1": 60.5, "phase 9": 10.0}, pct
    fig = xp.create_plot(*drawn.args, name, pct)
    texts, _ = legend_of(fig)
    assert "Phase 1 (60.5%)" in texts, texts
    assert fig.axes[1].get_ylabel() == r"diff/$\sigma$"


def test_colliding_percentage_columns_print_nothing():
    assert xp.phase_fraction({"phase": 40.0, "phase 1": 45.0}, "Phase 1") == 0.0


# 10. The 2theta window and the intensity axis.
def test_window_precedence(synth, drawn):
    theta = drawn.theta
    assert xp.plot_window(theta) == (10.0, 90.0)
    xp.PLOT_X_MIN, xp.PLOT_X_MAX = 13, 85
    assert xp.plot_window(theta) == (13.0, 85.0), "the constants must win"
    assert xp.plot_window(theta, 20.0, 60.0) == (20.0, 60.0), "metadata first"


def test_window_from_metadata_rescales_the_intensity_axis(synth, drawn):
    _, _, _, window = xp.sample_info(
        xp.load_metadata(synth.tmp / "metadata_window.csv"), "sample_A.csv",
        "sample_A")
    assert window == (20.0, 60.5), window
    cropped = xp.create_plot(*drawn.args, drawn.name, drawn.pct, xlim=window)
    assert cropped.axes[1].get_xlim() == (20.0, 60.5)
    whole = xp.create_plot(*drawn.args, drawn.name, drawn.pct)
    assert cropped.axes[0].get_ylim()[1] < whole.axes[0].get_ylim()[1]


def test_top_margin_clears_a_calc_peak_above_the_obs(tmp_path):
    # A fit that overshoots draws a calc peak taller than the obs scatter;
    # the top margin is set on the obs max, so it must still clear calc.
    x = np.linspace(5, 90, 800)
    obs = 50 + 100 * np.exp(-((x - 30) ** 2) / 0.4)
    calc = 50 + 150 * np.exp(-((x - 30) ** 2) / 0.4)
    df = pd.DataFrame({"2theta": x, "Obs": obs, "Calc": calc, "bkg": 50.0,
                       "diff/sigma": 0.0,
                       "Ph": np.where((x > 29) & (x < 31), 1.0, np.nan)})
    csv = tmp_path / "overshoot.csv"
    df.to_csv(csv, sep=";", index=False)
    data, phase_cols, error = xp.read_gsas2_csv(csv)
    assert error is None, error
    theta, o, c, b, r, ph = xp.prepare_data(data, phase_cols, use_sqrt=True)
    fig = xp.create_plot(theta, o, c, b, r, ph, "n", {}, use_sqrt=True)
    assert c.max() < fig.axes[0].get_ylim()[1], "the calc peak is clipped"
    plt.close(fig)


# 10b. Only the 2theta axis is numbered.
@pytest.mark.parametrize("use_sqrt, ylabel", [
    (True, r"$\sqrt{Intensity}$ / a.u."),
    (False, r"Intensity / a.u."),
])
def test_axis_labels_carry_no_parentheses(drawn, use_sqrt, ylabel):
    fig = xp.create_plot(*drawn.args, drawn.name, drawn.pct, use_sqrt=use_sqrt)
    assert fig.axes[0].get_ylabel() == ylabel
    assert fig.axes[1].get_xlabel() == r"2$\theta$ / $^\circ$"


def test_neither_intensity_axis_is_ticked_or_numbered(drawn):
    fig = xp.create_plot(*drawn.args, drawn.name, drawn.pct)
    upper, lower = fig.axes
    fig.canvas.draw()  # ticks are laid out lazily
    for ax, panel in ((upper, "upper"), (lower, "residual")):
        assert ax.get_yticklabels() == [], f"the {panel} panel kept its numbers"
        assert not any(t.get_visible() for t in ax.yaxis.get_ticklines()), (
            f"the {panel} panel kept its y tick marks")
    # The 2theta axis keeps both, and only the lower panel draws them.
    assert [t.get_text() for t in lower.get_xticklabels()], "2theta lost its numbers"
    assert any(t.get_visible() for t in lower.xaxis.get_ticklines())
    assert upper.get_xticklabels() == [], "the upper panel repeated the 2theta numbers"


# 11. The unweighted panel.
def test_unweighted_residuals(synth, drawn):
    data, phase_cols, error = xp.read_gsas2_csv(synth.full, weighted=False)
    assert error is None, error
    assert np.array_equal(data["resid"], synth.obs - synth.calc), "not exact"
    fig = xp.create_plot(*xp.prepare_data(data, phase_cols), drawn.name,
                         drawn.pct, weighted=False)
    assert fig.axes[1].get_ylabel() == r"diff / a.u."

    _, _, error = xp.read_gsas2_csv(synth.a, weighted=False)
    assert error == "residual column 'diff' not found", error
    assert xp.WEIGHTED_RESIDUALS is True, "the constant must not be mutated"


@pytest.mark.parametrize("weighted", [True, False])
def test_residual_panel_carries_the_trace_alone(synth, drawn, weighted):
    """Nothing is drawn on the panel besides the residual, and it is not clipped.

    The panel carries no numbers, so a line at zero used to be the only
    reference on it. The line goes and the limits hold zero at the middle
    instead, which is a reference the eye reads without anything drawn on it.
    """
    data, phase_cols, error = xp.read_gsas2_csv(synth.full, weighted=weighted)
    assert error is None, error
    fig = xp.create_plot(*xp.prepare_data(data, phase_cols), drawn.name,
                         drawn.pct, weighted=weighted)
    lower = fig.axes[1]
    assert len(lower.lines) == 1, "the residual panel drew more than the trace"
    low, high = lower.get_ylim()
    assert low <= np.nanmin(data["resid"]) and high >= np.nanmax(data["resid"]), (
        f"the trace is clipped: panel {low} to {high}")


def test_the_weighted_panel_keeps_one_scale_across_samples(synth, drawn):
    """Zero at the middle, and the weighted panel no narrower than the span.

    Autoscale put the trace at a different height in every sample, since the
    misfit is not symmetric, and no two figures could be read side by side.
    diff/sigma is in standard deviations of the point, so a floor in those
    units is the same scale everywhere. The raw diff is in intensities and has no
    such unit, so it is centred and nothing more.
    """
    data, phase_cols, error = xp.read_gsas2_csv(synth.full)
    assert error is None, error
    theta, o, c, b, resid, phases = xp.prepare_data(data, phase_cols)

    fig = xp.create_plot(theta, o, c, b, resid, phases, drawn.name, drawn.pct)
    assert fig.axes[1].get_ylim() == (-xp.RESIDUAL_SPAN, xp.RESIDUAL_SPAN)

    # A residual past the floor widens the panel; clipping the misfit is the
    # one thing this panel must not do.
    spiked = resid.copy()
    spiked[0] = -3.0 * xp.RESIDUAL_SPAN
    fig = xp.create_plot(theta, o, c, b, spiked, phases, drawn.name, drawn.pct)
    low, high = fig.axes[1].get_ylim()
    assert low == -high and low <= spiked.min(), (low, high)

    # No floor on the raw diff: a small one keeps a small panel.
    fig = xp.create_plot(theta, o, c, b, np.full_like(resid, 0.01), phases,
                         drawn.name, drawn.pct, weighted=False)
    low, high = fig.axes[1].get_ylim()
    assert low == -high and high < xp.RESIDUAL_SPAN, (low, high)


# 12. The function behind the interactive panel.
def test_replot_file(synth):
    fig, line = xp.replot_file(synth.a, synth.meta, x_min=20, x_max=60.5,
                               y_max=40)
    assert fig.axes[1].get_xlim() == (20.0, 60.5)
    assert fig.axes[0].get_ylim()[1] == 40.0
    assert line.splitlines()[1] == "sample_A.csv;Sample A (synthetic);20;60.5"

    fig, _ = xp.replot_file(synth.full, synth.meta, weighted=False)
    assert fig.axes[1].get_ylabel() == r"diff / a.u."
    assert xp.WEIGHTED_RESIDUALS is True, "the constant must not be mutated"

    with pytest.raises(ValueError, match="residual column"):
        xp.replot_file(synth.tmp / "only_obs.csv", synth.meta)


# 12b. The section 4 Save button reuses the batch naming and writer.
def test_output_basename_matches_the_batch():
    assert xp.output_basename("s", use_sqrt=True) == "s_XRD_analysis_sqrt"
    assert xp.output_basename("s", use_sqrt=False) == "s_XRD_analysis_linear"
    # weighted=False adds _unweighted, so a save with the other residual mode
    # does not overwrite the weighted file under the same name.
    assert xp.output_basename("s", weighted=False) == "s_XRD_analysis_sqrt_unweighted"
    # weighted left None follows the constant, which the suite keeps at True.
    assert xp.output_basename("s") == "s_XRD_analysis_sqrt"


def test_sample_window_reads_the_metadata_row(synth, tmp_path):
    meta = tmp_path / "meta.csv"
    meta.write_text("filename;formula;x_min;x_max\n"
                    "sample_A.csv;F;20;70\nsample_B.csv;G;;\n")
    assert xp.sample_window(meta, "sample_A.csv") == (20.0, 70.0)
    # blank cells and an absent file both give the widen-to-full sentinel.
    assert xp.sample_window(meta, "sample_B.csv") == (None, None)
    assert xp.sample_window(meta, "not_listed.csv") == (None, None)


def test_save_figure_writes_both_formats(tmp_path):
    fig, _ = plt.subplots()
    base = xp.save_figure(fig, tmp_path / "out", "sample_XRD_analysis_sqrt")
    plt.close(fig)
    assert base == "sample_XRD_analysis_sqrt"
    # Files are sorted into pdf/ and png/ subfolders, not the output root.
    assert (tmp_path / "out" / "pdf" / "sample_XRD_analysis_sqrt.pdf").is_file()
    assert (tmp_path / "out" / "png" / "sample_XRD_analysis_sqrt.png").is_file()
    assert not (tmp_path / "out" / "sample_XRD_analysis_sqrt.pdf").exists()


# 13. One unusable file must not end the batch, and the report must say so.
def test_the_batch_isolates_a_bad_file_and_reports_every_block(synth, tmp_path,
                                                              capsys):
    """One run over a good, a broken and a shifted file.

    Covers what the reader depends on: the two good files are written, the
    other two are isolated with a reason, each block prints its name, phases,
    2theta window and output in that order, and every failure is repeated in
    the summary at the end.
    """
    folder, out = tmp_path / "in", tmp_path / "out"
    folder.mkdir()
    synth.df.to_csv(folder / "a_good.csv", sep=";", index=False)
    synth.df.assign(Obs="n/a").to_csv(folder / "b_broken.csv", sep=";",
                                      index=False)
    # A shifted export: the angles land under 'used', the intensities under the
    # 2theta header. It must be isolated, not drawn.
    rows = [f"{10 + 0.2 * i:.4f};{100 + 50 * (i % 2)};1;1;1" for i in range(300)]
    (folder / "c_shifted.csv").write_text(
        "used;x, 2theta (deg);obs;calc;diff/sigma\n" + "\n".join(rows) + "\n")
    synth.df.to_csv(folder / "d_good.csv", sep=";", index=False)

    outcome = xp.process_folder(folder, synth.meta, out, show=False)
    assert [s for _, s, _ in outcome] == ["ok", "error", "error", "ok"], outcome
    assert len(list((out / "png").glob("*.png"))) == 2
    assert len(list((out / "pdf").glob("*.pdf"))) == 2

    text = capsys.readouterr().out
    assert "\na_good.csv\n" in text, "the file name must open its own block"
    assert "\n  phases: Phase 1 #" in text, text
    saved = next(l for l in text.splitlines() if l.startswith("  saved "))
    assert saved.endswith("pdf/a_good_XRD_analysis_sqrt.pdf and "
                          "png/a_good_XRD_analysis_sqrt.png"), saved
    assert (text.index("\na_good.csv\n") < text.index("\n  phases: ")
            < text.index("\n  2theta window: ")
            < text.index("\n  saved ")), "the block lost its order"
    assert "(auto)" in text, "the window source must be named"
    assert "  FAILED: " in text, "the failing file must say so in its block"
    tail = text.rsplit("=" * 60, 1)[1]
    assert "2 of 4 file(s) plotted" in tail, tail
    assert "b_broken.csv" in tail and "c_shifted.csv" in tail, (
        "every failure must survive to the summary: " + tail)


def test_the_batch_window_line_names_the_metadata_source(synth, tmp_path,
                                                         capsys):
    folder, out = tmp_path / "in", tmp_path / "out"
    folder.mkdir()
    synth.df.to_csv(folder / "sample_A.csv", sep=";", index=False)
    meta = folder / "meta.csv"
    meta.write_text("filename;formula;x_min;x_max\nsample_A.csv;F;20;70\n")

    xp.process_folder(folder, meta, out, show=False)
    assert "  2theta window: 20 to 70 (metadata)" in capsys.readouterr().out
