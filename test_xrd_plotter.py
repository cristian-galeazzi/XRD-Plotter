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


# 11. The unweighted panel.
def test_unweighted_residuals(synth, drawn):
    data, phase_cols, error = xp.read_gsas2_csv(synth.full, weighted=False)
    assert error is None, error
    assert np.array_equal(data["resid"], synth.obs - synth.calc), "not exact"
    fig = xp.create_plot(*xp.prepare_data(data, phase_cols), drawn.name,
                         drawn.pct, weighted=False)
    assert fig.axes[1].get_ylabel() == r"diff / (a.u.)"

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
    assert fig.axes[1].get_ylabel() == r"diff / (a.u.)"
    assert xp.WEIGHTED_RESIDUALS is True, "the constant must not be mutated"

    with pytest.raises(ValueError, match="residual column"):
        xp.replot_file(synth.tmp / "only_obs.csv", synth.meta)


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
    assert len(list(out.glob("*.png"))) == 2
