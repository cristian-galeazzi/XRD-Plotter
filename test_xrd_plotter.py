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
          "PHASE_COLORS", "PHASE_MAX_FILL")


@pytest.fixture(autouse=True)
def config_restored():
    """A test that changes a module constant must not leak it to the next."""
    saved = {name: getattr(xp, name) for name in CONFIG}
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
])
def test_degraded_files_are_isolated(synth, name, reason):
    data, _, error = xp.read_gsas2_csv(synth.tmp / name)
    assert data is None and error == reason, error


# 2b. The native GSAS-II export, whose header is one name too long.
def test_shifted_export_is_refused(tmp_path):
    # GSAS-II writes 11 header names and 10 fields per row, so every column
    # is read under its left neighbour's name: 'used' takes the angles and
    # the 2theta header takes the counts. The figure would still draw.
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


@pytest.mark.parametrize("header", [
    "x, 2theta (deg)", "2theta", "2Theta (deg)", "X, 2THETA", "x, deg",
    "angle, 2theta", "2theta_deg", "2θ", "2 theta (deg)", "two-theta", "x",
])
def test_theta_header_spellings_accepted(header):
    assert xp.is_theta_header(header), header


@pytest.mark.parametrize("header", [
    "theta", "angle", "obs", "counts", "intensity", "Max, phase", "tick-pos",
])
def test_theta_header_spellings_refused(header):
    assert not xp.is_theta_header(header), header


def test_a_phase_name_containing_x_comma_is_not_taken_as_the_axis(tmp_path):
    # 'x,' sits inside ordinary words, so 'Max, phase' matches the same rule
    # as 'x, 2theta (deg)'. The short column must not become the axis.
    x = np.linspace(10.0, 90.0, 300)
    df = pd.DataFrame({"Max, phase": pd.Series([30.0, 35.0]).reindex(range(300)),
                       "x, 2theta (deg)": x,
                       "obs": np.linspace(1.0, 2.0, 300),
                       "diff/sigma": np.zeros(300)})
    csv = tmp_path / "collision.csv"
    df.to_csv(csv, sep=";", index=False)
    data, phase_cols, error = xp.read_gsas2_csv(csv)
    assert error is None, error
    assert np.array_equal(data["x"], x), "a phase column was read as 2theta"
    assert phase_cols == ["Max, phase"], phase_cols


def test_column_order_does_not_change_the_result(tmp_path):
    x = np.linspace(10.0, 90.0, 300)
    obs = np.linspace(1.0, 2.0, 300)
    base = {"x, 2theta (deg)": x, "obs": obs, "calc": obs, "bkg": obs,
            "diff": np.zeros(300), "diff/sigma": np.zeros(300),
            "GOF": pd.Series([1.8]).reindex(range(300)),
            "Rutile": pd.Series([30.0, 35.0]).reindex(range(300)),
            "Anatase": pd.Series([32.0]).reindex(range(300))}
    seen = set()
    for i, order in enumerate([list(base), list(base)[::-1],
                               sorted(base), sorted(base, reverse=True)]):
        csv = tmp_path / f"order_{i}.csv"
        pd.DataFrame({k: base[k] for k in order}).to_csv(csv, sep=";",
                                                        index=False)
        data, phase_cols, error = xp.read_gsas2_csv(csv)
        assert error is None, error
        seen.add((data["x"].tobytes(), data["obs"].tobytes(),
                  tuple(sorted(phase_cols))))
    assert len(seen) == 1, "the column order changed what was parsed"


def test_is_monotonic_tolerates_a_few_out_of_order_points():
    rising = np.arange(500.0)
    assert xp.is_monotonic(rising)
    rising[100], rising[101] = rising[101], rising[100]  # one swapped pair
    assert xp.is_monotonic(rising), "one glitch must not reject a real axis"
    assert not xp.is_monotonic(np.tile([1.0, 2.0], 250)), "counts accepted"


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
    # weighted=False adds _counts, so a save with the other residual mode
    # does not overwrite the weighted file under the same name.
    assert xp.output_basename("s", weighted=False) == "s_XRD_analysis_sqrt_counts"
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


# 13. One unusable file must not end the batch.
def test_batch_survives_a_broken_file(synth, tmp_path):
    folder, out = tmp_path / "in", tmp_path / "out"
    folder.mkdir()
    synth.df.to_csv(folder / "a_good.csv", sep=";", index=False)
    synth.df.assign(Obs="n/a").to_csv(folder / "b_broken.csv", sep=";",
                                      index=False)
    synth.df.to_csv(folder / "c_good.csv", sep=";", index=False)

    outcome = xp.process_folder(folder, synth.meta, out, show=False)
    assert [status for _, status, _ in outcome] == ["ok", "error", "ok"]
    assert len(list((out / "png").glob("*.png"))) == 2
    assert len(list((out / "pdf").glob("*.pdf"))) == 2


# 14. The printed report: one block per file, failures visible at the end.
def test_batch_prints_one_block_per_file(synth, tmp_path, capsys):
    folder, out = tmp_path / "in", tmp_path / "out"
    folder.mkdir()
    synth.df.to_csv(folder / "a_good.csv", sep=";", index=False)

    xp.process_folder(folder, synth.meta, out, show=False)
    text = capsys.readouterr().out
    assert "\na_good.csv\n" in text, "the file name must open its own block"
    assert "\n  phases: Phase 1 #" in text, text
    saved_line = next(line for line in text.splitlines()
                      if line.startswith("  saved "))
    assert saved_line.endswith("pdf/a_good_XRD_analysis_sqrt.pdf and "
                               "png/a_good_XRD_analysis_sqrt.png"), (
        "the saved line must name the file it wrote: " + saved_line)
    assert (text.index("\na_good.csv\n") < text.index("\n  phases: ")
            < text.index("\n  saved ")), "the block lost its order"
    # The window line is shown, marked auto when no metadata row set it.
    assert "\n  2theta window: " in text and "(auto)" in text, text
    assert (text.index("\n  phases: ") < text.index("\n  2theta window: ")
            < text.index("\n  saved ")), "the window line is out of place"


def test_batch_window_line_names_the_metadata_source(synth, tmp_path):
    folder, out = tmp_path / "in", tmp_path / "out"
    folder.mkdir()
    synth.df.to_csv(folder / "sample_A.csv", sep=";", index=False)
    meta = folder / "meta.csv"
    meta.write_text("filename;formula;x_min;x_max\nsample_A.csv;F;20;70\n")

    import io as _io
    import contextlib
    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf):
        xp.process_folder(folder, meta, out, show=False)
    assert "  2theta window: 20 to 70 (metadata)" in buf.getvalue()


def test_batch_repeats_failures_in_the_summary(synth, tmp_path, capsys):
    folder, out = tmp_path / "in", tmp_path / "out"
    folder.mkdir()
    synth.df.to_csv(folder / "a_good.csv", sep=";", index=False)
    synth.df.assign(Obs="n/a").to_csv(folder / "b_broken.csv", sep=";",
                                      index=False)
    synth.df.to_csv(folder / "c_good.csv", sep=";", index=False)

    xp.process_folder(folder, synth.meta, out, show=False)
    text = capsys.readouterr().out
    assert "  FAILED: " in text, "the failing file must say so in its block"
    assert "=" * 60 in text, "the summary lost its separator"
    tail = text.rsplit("=" * 60, 1)[1]
    assert "2 of 3 file(s) plotted" in tail, tail
    assert "FAILED (1): b_broken.csv" in tail, "a failure must survive the run"
