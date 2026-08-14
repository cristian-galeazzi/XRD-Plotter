"""Tests for the stacked series figure script."""
import matplotlib

matplotlib.use("Agg")  # no display in CI, same as test_xrd_plotter.py

import ast
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.legend import Legend
from matplotlib.pyplot import close as plt_close

import make_series as ms
import xrd_plotter as xp


@pytest.fixture(autouse=True)
def shipped_defaults(monkeypatch):
    """Every test runs against the settings the repository ships.

    The copy of the script an owner runs carries their own tuning, so a
    test that read a setting off the module would pass in CI and fail on
    their machine.
    """
    for name, value in (("GUIDE_LINES", []), ("GUIDE_LABELS", []),
                        ("GUIDE_LABEL_WEIGHT", "bold"),
                        ("GUIDE_LABEL_HEIGHT", 0.20), ("LABEL_X", None),
                        ("LABEL_WEIGHT", "normal"), ("USE_SQRT", True),
                        ("PLOT_X_MIN", None), ("PLOT_X_MAX", None),
                        ("LABEL_HEIGHT", 0.90), ("SHOW_TICKS", True),
                        ("OFFSET", 1.35), ("TICK_HEIGHT", 0.10),
                        ("GUIDE_SNAP", 0.3), ("LINEWIDTH_TRACE", 0.9),
                        ("GUIDE_STYLE", ":"), ("GUIDE_WIDTH", 1.2)):
        monkeypatch.setattr(ms, name, value)


def test_each_pattern_is_normalised_to_its_own_span():
    low_contrast = np.array([10.0, 12.0, 10.0])
    high_contrast = np.array([0.0, 100.0, 0.0])
    first, second = ms.stack([low_contrast, high_contrast], offset=0.0)
    assert np.array_equal(first, np.array([0.0, 1.0, 0.0]))
    assert np.array_equal(second, np.array([0.0, 1.0, 0.0]))


def test_each_trace_is_raised_by_its_index_times_the_offset():
    flat_top = np.array([0.0, 1.0])
    traces = ms.stack([flat_top, flat_top, flat_top], offset=1.5)
    assert [float(t[0]) for t in traces] == [0.0, 1.5, 3.0]


def test_a_flat_pattern_sits_on_its_baseline_instead_of_becoming_nan():
    traces = ms.stack([np.array([5.0, 5.0, 5.0])], offset=1.0)
    assert np.array_equal(traces[0], np.zeros(3))


def test_a_scope_sets_the_range_instead_of_the_whole_pattern():
    # Only the first two points set the range; the third, a tall value the
    # window would crop, must not flatten the visible part of the trace.
    pattern = np.array([0.0, 2.0, 100.0])
    scope = np.array([True, True, False])
    stacked = ms.stack([pattern], offset=0.0, scopes=[scope])[0]
    assert stacked[:2].tolist() == [0.0, 1.0]
    assert stacked[2] > 1.0  # outside the scope, left to be cropped by xlim


def no_labels(*names):
    """The (file name, label) pairs load_series takes, with no label set.

    Most tests care about the figure rather than about the names on it, and
    a label of None is what a metadata row without a 'series_label' cell
    produces.
    """
    return [(name, None) for name in names]


HEADER = "2theta;Obs;Calc;Bkg;diff/sigma;Phase 1"


def write_export(path, peak_at, height, second_phase=None):
    """A synthetic GSAS-II publication export: one peak on a flat background.

    The reflection column holds one position and is padded, so the parser
    reads it as a phase rather than as a data column. 'second_phase' adds a
    second reflection column when given, for tests that need a second row.
    """
    angles = np.arange(10.0, 50.0, 0.5)
    obs = 100.0 + height * np.exp(-((angles - peak_at) ** 2) / 0.5)
    header = HEADER + (";Phase 2" if second_phase is not None else "")
    lines = [header]
    for i, (angle, value) in enumerate(zip(angles, obs)):
        reflection = f"{peak_at:.2f}" if i == 0 else ""
        row = f"{angle:.2f};{value:.4f};{value:.4f};100.0;0.1;{reflection}"
        if second_phase is not None:
            row += f";{second_phase:.2f}" if i == 0 else ";"
        lines.append(row)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_export_two_peaks(path, near_peak, near_height, far_peak,
                           far_height):
    """A synthetic export with two peaks, one meant to sit outside a window.

    The reflection column snaps to near_peak, so a guide or tick test can
    still use it; far_peak exists only to distort a normalisation that looks
    at the whole pattern instead of the plotted window.
    """
    angles = np.arange(10.0, 50.0, 0.5)
    obs = (100.0
          + near_height * np.exp(-((angles - near_peak) ** 2) / 0.5)
          + far_height * np.exp(-((angles - far_peak) ** 2) / 0.5))
    lines = [HEADER]
    for i, (angle, value) in enumerate(zip(angles, obs)):
        reflection = f"{near_peak:.2f}" if i == 0 else ""
        lines.append(f"{angle:.2f};{value:.4f};{value:.4f};100.0;0.1;"
                     f"{reflection}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@pytest.fixture
def series(tmp_path):
    """Three synthetic exports whose peak walks across the pattern."""
    folder = tmp_path / "data"
    folder.mkdir()
    for index, (peak, height) in enumerate([(20.0, 900.0), (25.0, 400.0),
                                            (30.0, 8000.0)], start=1):
        write_export(folder / f"s{index}.csv", peak, height)
    return folder


def test_every_named_file_becomes_one_trace(series, tmp_path):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv", "s3.csv"), series, pd.DataFrame())
    assert len(traces) == 3
    assert phases, "the reflection column should be read as a phase"


def test_a_missing_file_stops_the_run_and_names_itself(series, tmp_path):
    with pytest.raises(SystemExit, match="s9.csv"):
        ms.load_series(no_labels("s9.csv"), series, pd.DataFrame())


def test_the_trace_label_falls_back_to_the_file_stem(series, tmp_path):
    traces, _, _colors = ms.load_series(
            no_labels("s2.csv"), series, pd.DataFrame())
    assert traces[0][2] == "s2"


def test_the_figure_carries_one_line_and_one_label_per_sample(series,
                                                              tmp_path):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv", "s3.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    ax = fig.axes[0]
    assert len(ax.lines) == 3
    assert sorted(t.get_text() for t in ax.texts) == ["s1", "s2", "s3"]
    plt_close(fig)


def test_the_strongest_sample_does_not_overrun_the_trace_above(series,
                                                              tmp_path):
    """s3 scatters many times harder than s1; normalisation must hide that."""
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv", "s3.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, offset=1.5)
    tops = [float(line.get_ydata().max()) for line in fig.axes[0].lines]
    assert tops[0] < 1.5 and tops[1] < 3.0
    plt_close(fig)


def test_the_intensity_axis_carries_no_numbers(series, tmp_path):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    ax = fig.axes[0]
    assert not any(label.get_visible() for label in ax.get_yticklabels())
    plt_close(fig)


def write_metadata(path, filename, formula):
    """A one-row private metadata file, the format load_metadata expects."""
    path.write_text(f"filename;formula\n{filename};{formula}\n",
                    encoding="utf-8")


def test_a_written_label_wins_over_the_metadata_formula(series, tmp_path):
    meta = tmp_path / "meta.csv"
    write_metadata(meta, "s1.csv", "Formula From Metadata")
    traces, _, _colors = ms.load_series(
            [("s1.csv", "Label I Wrote")], series, xp.load_metadata(meta))
    assert traces[0][2] == "Label I Wrote"


def test_without_a_written_label_the_metadata_formula_is_used(series,
                                                              tmp_path):
    meta = tmp_path / "meta.csv"
    write_metadata(meta, "s1.csv", "Formula From Metadata")
    traces, _, _colors = ms.load_series(
            [("s1.csv", None)], series, xp.load_metadata(meta))
    assert traces[0][2] == "Formula From Metadata"


def test_each_label_sits_at_the_top_of_its_own_trace(series, tmp_path):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv", "s3.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, offset=1.35)
    placed = {text.get_text(): text.get_position()[1]
              for text in fig.axes[0].texts}
    # stack() normalises every pattern to a span of 1.0, so the top of the
    # trace at index i sits at i * offset + 1.0 and LABEL_HEIGHT is measured
    # from the same baseline. One label per trace, one offset apart, none of
    # them on the trace above.
    assert placed["s1"] == ms.LABEL_HEIGHT
    assert round(placed["s2"] - placed["s1"], 10) == 1.35
    assert round(placed["s3"] - placed["s2"], 10) == 1.35
    assert ms.LABEL_HEIGHT < 1.35, "a label at OFFSET hits the trace above"
    plt_close(fig)


def test_by_default_the_labels_start_just_inside_the_left_border(series,
                                                                 tmp_path):
    # LABEL_X is pinned to None by the shipped_defaults fixture: the copy
    # of the script an owner runs carries their own tuning, and this is a
    # claim about the default the repository ships.
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    label = fig.axes[0].texts[0]
    assert label.get_ha() == "left"
    # An axes fraction while LABEL_X is unset, so it holds for any window.
    assert 0.0 <= label.get_position()[0] < 0.1
    plt_close(fig)


def test_a_label_x_in_degrees_puts_the_label_at_that_angle(series, tmp_path,
                                                           monkeypatch):
    monkeypatch.setattr(ms, "LABEL_X", 18.0)
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    ax = fig.axes[0]
    label = ax.texts[0]
    assert label.get_position()[0] == 18.0
    # In data coordinates now, so 18 means 18 degrees rather than 18 widths
    # of the axes: the transform has to have changed with the setting.
    assert label.get_transform() is ax.transData
    plt_close(fig)


def test_the_topmost_label_is_not_cut_off_by_the_frame(series, tmp_path):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv", "s3.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, offset=1.35)
    ax = fig.axes[0]
    highest_label = max(text.get_position()[1] for text in ax.texts)
    assert ax.get_ylim()[1] > highest_label + 0.2, (
        "no room above the top label for the text itself")
    plt_close(fig)


def test_the_topmost_label_is_not_cut_off_above_a_full_trace(series,
                                                              tmp_path,
                                                              monkeypatch):
    """LABEL_HEIGHT may sit above the top of the trace, past 1.0."""
    monkeypatch.setattr(ms, "LABEL_HEIGHT", 1.20)
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv", "s3.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, offset=1.35)
    ax = fig.axes[0]
    highest_label = max(text.get_position()[1] for text in ax.texts)
    assert ax.get_ylim()[1] > highest_label + 0.2, (
        "no room above the top label for the text itself")
    plt_close(fig)


def test_the_phase_legend_sits_opposite_the_trace_labels(series, tmp_path):
    """The labels start at the left, so the legend keeps to the right."""
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    legend = fig.axes[0].get_legend()
    assert legend._get_loc() == Legend.codes["upper right"]
    # An opaque frame, so a trace passing behind it does not show through
    # the phase names.
    assert legend.get_frame_on() is True
    assert legend.get_frame().get_alpha() == 1
    plt_close(fig)


def test_snap_takes_the_nearest_reflection_inside_the_tolerance():
    positions = np.array([20.0, 30.95, 37.04])
    assert ms.snap_to_reflection(31.0, positions, tolerance=0.3) == 30.95


def test_snap_returns_none_when_nothing_is_near_enough():
    positions = np.array([20.0, 37.04])
    assert ms.snap_to_reflection(31.0, positions, tolerance=0.3) is None


def test_snap_includes_a_value_exactly_at_the_tolerance_boundary():
    # abs(nearest - value) <= tolerance: a value exactly tolerance away must
    # still snap, pinning the inclusive bound against a future '<'. 0.5 is
    # exact in binary floating point, so the difference lands on the
    # boundary exactly rather than a hair past it.
    positions = np.array([20.0])
    assert ms.snap_to_reflection(20.5, positions, tolerance=0.5) == 20.0


def test_snap_on_an_empty_reflection_list_returns_none():
    assert ms.snap_to_reflection(31.0, np.array([]), tolerance=0.3) is None


def test_a_guide_is_drawn_at_the_reflection_not_at_the_typed_value(
        series, tmp_path, monkeypatch):
    # The first file's only reflection sits at 20.0; 20.2 is inside the
    # default 0.3 tolerance, so the guide must land on 20.0.
    monkeypatch.setattr(ms, "GUIDE_LINES", [20.2])
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    guides = [line for line in fig.axes[0].lines
              if line.get_linestyle() in (":", "dotted")]
    assert len(guides) == 1
    assert float(guides[0].get_xdata()[0]) == 20.0
    plt_close(fig)


def test_a_guide_with_no_reflection_near_it_is_reported_and_not_drawn(
        series, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(ms, "GUIDE_LINES", [24.0])
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    guides = [line for line in fig.axes[0].lines
              if line.get_linestyle() in (":", "dotted")]
    assert guides == []
    assert "24" in capsys.readouterr().out
    plt_close(fig)


def test_no_guides_are_drawn_by_default(series, tmp_path):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    assert [line for line in fig.axes[0].lines
            if line.get_linestyle() in (":", "dotted")] == []
    plt_close(fig)


def test_the_ticks_are_taller_than_the_traces_are_apart_is_not_the_case(
        series, tmp_path):
    """A tick row must fit under the bottom trace without reaching it."""
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    offset = 1.35
    fig = ms.plot_series(traces, phases, offset=offset)
    ax = fig.axes[0]
    tick_collection = ax.collections[0]
    tops = [segment[:, 1].max() for segment in
            [np.array(s) for s in tick_collection.get_segments()]]
    # The tick row's own top sits at block_top + TICK_HEIGHT, where
    # block_top is -(TICK_HEIGHT + 0.05): the two cancel to -0.05, the
    # clearance plot_series intends between the bottom trace's baseline at
    # y=0 and the nearest tick. In trace heights, so it does not move with
    # the offset.
    expected_top = -0.05
    assert max(tops) == pytest.approx(expected_top), (
        "the ticks reach into the bottom trace")
    assert ax.get_ylim()[0] < min(
        np.array(s)[:, 1].min()
        for s in tick_collection.get_segments()), "the ticks are clipped"
    plt_close(fig)


def test_the_window_not_the_off_window_peak_sets_the_normalisation(
        tmp_path, monkeypatch):
    """A tall peak outside the window must not flatten the trace inside it."""
    folder = tmp_path / "data"
    folder.mkdir()
    write_export(folder / "plain.csv", peak_at=20.0, height=50.0)
    write_export_two_peaks(folder / "loaded.csv", near_peak=20.0,
                           near_height=50.0, far_peak=45.0, far_height=8000.0)
    monkeypatch.setattr(ms, "PLOT_X_MIN", 15.0)
    monkeypatch.setattr(ms, "PLOT_X_MAX", 25.0)
    traces, phases, _colors = ms.load_series(
            no_labels("plain.csv", "loaded.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    tops = []
    for line in fig.axes[0].lines:
        theta, values = line.get_xdata(), line.get_ydata()
        visible = (theta >= 15.0) & (theta <= 25.0)
        tops.append(float(values[visible].max() - values[visible].min()))
    # Without the fix 'loaded' is stretched flat by its off-window peak at
    # 45 deg and its visible span collapses far below 1.0.
    assert tops[1] > 0.9, "the visible peak should reach close to full span"
    assert tops[0] == pytest.approx(tops[1], abs=0.05), (
        "an off-window peak must not change how tall the two samples "
        "compare inside the window")
    plt_close(fig)


def test_the_drawn_line_does_not_spike_at_the_window_edge(tmp_path,
                                                           monkeypatch):
    """A grid point just outside the window must not draw a spike.

    The window edge sits between two grid points, 25.0 inside and 25.5
    outside. Plotting the whole trace would draw the segment between them
    across the frame, interpolated to a height set by the tall off-window
    point; plotting only the scope slice must not.
    """
    folder = tmp_path / "data"
    folder.mkdir()
    write_export_two_peaks(folder / "spiked.csv", near_peak=20.0,
                           near_height=5.0, far_peak=27.0, far_height=9000.0)
    monkeypatch.setattr(ms, "PLOT_X_MIN", 15.0)
    monkeypatch.setattr(ms, "PLOT_X_MAX", 25.25)
    traces, phases, _colors = ms.load_series(
            no_labels("spiked.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    values = fig.axes[0].lines[0].get_ydata()
    assert float(np.max(values)) < 1.5, (
        "the drawn line must not spike far above the trace's own slot")
    plt_close(fig)


def test_a_trace_entirely_outside_the_window_still_draws(tmp_path,
                                                          monkeypatch):
    """A member whose whole measured range misses the window must not raise."""
    folder = tmp_path / "data"
    folder.mkdir()
    write_export(folder / "inside.csv", peak_at=20.0, height=900.0)
    angles = np.arange(60.0, 70.0, 0.5)
    obs = 100.0 + 900.0 * np.exp(-((angles - 65.0) ** 2) / 0.5)
    lines = [HEADER]
    for i, (angle, value) in enumerate(zip(angles, obs)):
        reflection = "65.00" if i == 0 else ""
        lines.append(f"{angle:.2f};{value:.4f};{value:.4f};100.0;0.1;"
                     f"{reflection}")
    (folder / "outside.csv").write_text("\n".join(lines) + "\n",
                                        encoding="utf-8")
    monkeypatch.setattr(ms, "PLOT_X_MIN", 10.0)
    monkeypatch.setattr(ms, "PLOT_X_MAX", 50.0)
    traces, phases, _colors = ms.load_series(
            no_labels("inside.csv", "outside.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    assert isinstance(fig, Figure)
    plt_close(fig)


def test_a_metadata_colour_reaches_the_tick_row(series, tmp_path):
    meta = tmp_path / "meta.csv"
    meta.write_text("filename;phase 1_color\ns1.csv;#123456\n",
                    encoding="utf-8")
    traces, phases, colors = ms.load_series(
            [("s1.csv", None)], series, xp.load_metadata(meta))
    fig = ms.plot_series(traces, phases, colors=colors)
    tick_collection = fig.axes[0].collections[0]
    assert tuple(tick_collection.get_color()[0]) == to_rgba("#123456")
    plt_close(fig)


def test_the_legend_order_matches_the_tick_row_order_not_the_column_order(
        tmp_path):
    """The export lists Phase 2 before Phase 1; the legend must not follow."""
    folder = tmp_path / "data"
    folder.mkdir()
    path = folder / "reversed.csv"
    angles = np.arange(10.0, 50.0, 0.5)
    obs = 100.0 + 900.0 * np.exp(-((angles - 20.0) ** 2) / 0.5)
    lines = ["2theta;Obs;Calc;Bkg;diff/sigma;Phase 2;Phase 1"]
    for i, (angle, value) in enumerate(zip(angles, obs)):
        phase2 = "30.00" if i == 0 else ""
        phase1 = "20.00" if i == 0 else ""
        lines.append(f"{angle:.2f};{value:.4f};{value:.4f};100.0;0.1;"
                     f"{phase2};{phase1}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    traces, phases, colors = ms.load_series(
            no_labels("reversed.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases, colors=colors)
    legend_labels = [t.get_text()
                     for t in fig.axes[0].get_legend().get_texts()]
    assert legend_labels == ["Phase 1", "Phase 2"]
    plt_close(fig)


def test_a_phase_with_no_reflections_in_window_gets_no_legend_entry(
        tmp_path, monkeypatch):
    """A phase entirely outside the window must draw no tick row either."""
    folder = tmp_path / "data"
    folder.mkdir()
    write_export(folder / "two.csv", peak_at=20.0, height=900.0,
                second_phase=60.0)
    monkeypatch.setattr(ms, "PLOT_X_MIN", 10.0)
    monkeypatch.setattr(ms, "PLOT_X_MAX", 50.0)
    traces, phases, colors = ms.load_series(
            no_labels("two.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases, colors=colors)
    ax = fig.axes[0]
    legend_labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert legend_labels == ["Phase 1"]
    assert len(ax.collections) == 1, "only one phase should draw a tick row"
    plt_close(fig)


def test_the_second_tick_row_sits_one_pitch_below_the_first(tmp_path):
    """A two-phase export must exercise the second and later tick rows."""
    folder = tmp_path / "data"
    folder.mkdir()
    write_export(folder / "two.csv", 20.0, 900.0, second_phase=30.0)
    traces, phases, _colors = ms.load_series(
            no_labels("two.csv"), folder, pd.DataFrame())
    offset = 1.35
    fig = ms.plot_series(traces, phases, offset=offset)
    ax = fig.axes[0]
    assert len(ax.collections) == 2, "both phase columns should draw a row"
    bottoms = [np.array(c.get_segments()[0])[:, 1].min()
              for c in ax.collections]
    tops = [np.array(c.get_segments()[0])[:, 1].max()
           for c in ax.collections]
    pitch = ms.TICK_HEIGHT + 0.02
    assert bottoms[0] - bottoms[1] == pytest.approx(pitch)
    assert tops[1] < bottoms[0], "the second row overlaps the first"
    assert ax.get_ylim()[0] < bottoms[1], "the second row is clipped"
    plt_close(fig)


def test_padding_and_out_of_window_positions_are_not_reflections():
    positions = np.array([0.0, 12.0, 30.95, 90.0])
    inside = ms.reflections_in_window(positions, (13.0, 85.0))
    assert inside.tolist() == [30.95]


def test_without_a_window_only_the_padding_is_dropped():
    positions = np.array([0.0, 12.0, 30.95, 90.0])
    assert ms.reflections_in_window(positions).tolist() == [12.0, 30.95, 90.0]


def test_the_printed_list_wraps_instead_of_running_off_the_terminal():
    many = np.arange(20.0, 60.0, 0.5)
    lines = ms.format_reflections({"Phase 1": many}).splitlines()
    assert len(lines) > 1
    assert all(len(line) <= 76 for line in lines)
    assert lines[0].startswith("  Phase 1: 20.00")


def test_the_run_prints_the_reflections_it_drew(series, tmp_path, capsys):
    traces, phases, colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, colors=colors)
    printed = capsys.readouterr().out
    # The fixture's only reflection sits at 20.0, and it is what the tick
    # row was drawn from, so it is what the list must offer.
    assert "to pick the guides from:" in printed
    assert "20.00" in printed
    plt_close(fig)


def test_snap_to_phase_names_the_phase_the_reflection_belongs_to():
    rows = {"Phase 1": np.array([30.95]), "Phase 2": np.array([37.04])}
    assert ms.snap_to_phase(37.0, rows, 0.3) == (37.04, "Phase 2")
    assert ms.snap_to_phase(30.9, rows, 0.3) == (30.95, "Phase 1")


def test_snap_to_phase_returns_nothing_when_no_phase_is_near_enough():
    rows = {"Phase 1": np.array([30.95]), "Phase 2": np.array([37.04])}
    assert ms.snap_to_phase(27.0, rows, 0.3) == (None, None)


def test_the_nearer_phase_wins_and_a_tie_goes_to_the_first_row():
    rows = {"Phase 1": np.array([20.1]), "Phase 2": np.array([20.4])}
    assert ms.snap_to_phase(20.0, rows, 0.5)[1] == "Phase 1"
    assert ms.snap_to_phase(20.5, rows, 0.5)[1] == "Phase 2"
    tie = {"Phase 1": np.array([19.8]), "Phase 2": np.array([20.2])}
    assert ms.snap_to_phase(20.0, tie, 0.5)[1] == "Phase 1"


def test_a_guide_takes_the_colour_of_its_own_phase(tmp_path, monkeypatch):
    folder = tmp_path / "data"
    folder.mkdir()
    write_export(folder / "two.csv", 20.0, 900.0, second_phase=30.0)
    monkeypatch.setattr(ms, "GUIDE_LINES", [30.1])
    traces, phases, colors = ms.load_series(
            no_labels("two.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases, colors=colors)
    ax = fig.axes[0]
    guides = [line for line in ax.lines
              if line.get_linestyle() in (":", "dotted")]
    assert len(guides) == 1
    # The tick rows are drawn as collections, in the same order as the
    # legend; the second phase owns the reflection at 30.0.
    second_row = ax.collections[1]
    assert to_rgba(guides[0].get_color()) == to_rgba(
        second_row.get_colors()[0])
    plt_close(fig)


def test_two_columns_of_one_phase_are_listed_together(tmp_path, capsys):
    folder = tmp_path / "data"
    folder.mkdir()
    write_export(folder / "two.csv", 20.0, 900.0, second_phase=30.0)
    traces, phases, colors = ms.load_series(
            no_labels("two.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases, colors=colors)
    printed = capsys.readouterr().out
    assert "20.00" in printed and "30.00" in printed
    plt_close(fig)


def test_the_label_weight_reaches_the_drawn_text(series, tmp_path,
                                                 monkeypatch):
    monkeypatch.setattr(ms, "LABEL_WEIGHT", "bold")
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    assert fig.axes[0].texts[0].get_fontweight() == "bold"
    plt_close(fig)


def test_the_sqrt_and_linear_variants_do_not_share_an_output_name(
        monkeypatch, series):
    """main() must pick the output name up from xp.output_basename.

    A fixed OUTPUT_BASENAME passed straight to xp.save_figure would write
    both USE_SQRT settings to the same PDF and PNG, the sqrt run silently
    overwriting the linear one or the other way round. Driven through main()
    rather than by calling xp.output_basename twice here, which would assert
    something about xrd_plotter and nothing about this script.
    """
    meta = pd.DataFrame({"series_order": ["1"]},
                        index=pd.Index(["s1.csv"], name="filename"))
    monkeypatch.setattr(xp, "load_metadata", lambda _path: meta)
    monkeypatch.setattr(ms, "DATA_FOLDER", series)
    written = []

    def record(fig, _folder, base):
        written.append(base)
        plt_close(fig)  # main() hands its figure to save_figure and forgets
        return base

    monkeypatch.setattr(xp, "save_figure", record)
    for use_sqrt in (True, False):
        monkeypatch.setattr(ms, "USE_SQRT", use_sqrt)
        ms.main()
    assert len(written) == 2
    assert written[0] != written[1], (
        "both USE_SQRT settings write to the same file")


def write_series_metadata(path, rows):
    """A metadata file with the series columns; rows are dicts of cells."""
    columns = ["filename", "formula", "series_order", "series_label"]
    lines = [";".join(columns)]
    for row in rows:
        lines.append(";".join(str(row.get(c, "")) for c in columns))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_only_the_rows_with_an_order_join_the_series(tmp_path):
    meta_file = tmp_path / "meta.csv"
    write_series_metadata(meta_file, [
        {"filename": "s1.csv", "series_order": 2},
        {"filename": "s2.csv"},
        {"filename": "s3.csv", "series_order": 1},
    ])
    entries = ms.series_from_metadata(xp.load_metadata(meta_file))
    assert [name for name, _label in entries] == ["s3.csv", "s1.csv"]


def test_the_order_is_the_number_not_the_row_position(tmp_path):
    meta_file = tmp_path / "meta.csv"
    write_series_metadata(meta_file, [
        {"filename": "s1.csv", "series_order": 30},
        {"filename": "s2.csv", "series_order": 4},
        {"filename": "s3.csv", "series_order": 100},
    ])
    entries = ms.series_from_metadata(xp.load_metadata(meta_file))
    assert [name for name, _label in entries] == ["s2.csv", "s1.csv", "s3.csv"]


def test_a_series_label_is_carried_and_a_blank_one_is_not(tmp_path):
    meta_file = tmp_path / "meta.csv"
    write_series_metadata(meta_file, [
        {"filename": "s1.csv", "series_order": 1, "series_label": "First"},
        {"filename": "s2.csv", "series_order": 2},
    ])
    entries = ms.series_from_metadata(xp.load_metadata(meta_file))
    assert entries == [("s1.csv", "First"), ("s2.csv", None)]


def test_an_unreadable_order_is_reported_and_the_row_left_out(tmp_path,
                                                              capsys):
    meta_file = tmp_path / "meta.csv"
    write_series_metadata(meta_file, [
        {"filename": "s1.csv", "series_order": "first"},
        {"filename": "s2.csv", "series_order": 1},
    ])
    entries = ms.series_from_metadata(xp.load_metadata(meta_file))
    assert [name for name, _label in entries] == ["s2.csv"]
    assert "s1.csv" in capsys.readouterr().out


def test_a_blank_order_is_left_out_in_silence(tmp_path, capsys):
    """A blank cell is a choice, not a typo, and must not be reported."""
    meta_file = tmp_path / "meta.csv"
    write_series_metadata(meta_file, [
        {"filename": "s1.csv", "series_order": 1},
        {"filename": "s2.csv"},
    ])
    ms.series_from_metadata(xp.load_metadata(meta_file))
    assert "s2.csv" not in capsys.readouterr().out


def test_a_metadata_without_the_column_yields_no_series(tmp_path):
    meta_file = tmp_path / "meta.csv"
    meta_file.write_text("filename;formula\ns1.csv;Something\n",
                         encoding="utf-8")
    assert ms.series_from_metadata(xp.load_metadata(meta_file)) == []


def test_an_absent_metadata_yields_no_series():
    assert ms.series_from_metadata(pd.DataFrame()) == []


def test_the_series_label_wins_over_the_formula(series, tmp_path):
    meta_file = tmp_path / "meta.csv"
    write_series_metadata(meta_file, [
        {"filename": "s1.csv", "formula": "From Formula",
         "series_order": 1, "series_label": "From Label"},
    ])
    meta = xp.load_metadata(meta_file)
    traces, _phases, _colors = ms.load_series(
            ms.series_from_metadata(meta), series, meta)
    assert traces[0][2] == "From Label"


def test_without_a_series_label_the_formula_names_the_trace(series,
                                                            tmp_path):
    meta_file = tmp_path / "meta.csv"
    write_series_metadata(meta_file, [
        {"filename": "s1.csv", "formula": "From Formula", "series_order": 1},
    ])
    meta = xp.load_metadata(meta_file)
    traces, _phases, _colors = ms.load_series(
            ms.series_from_metadata(meta), series, meta)
    assert traces[0][2] == "From Formula"


def test_the_module_holds_no_list_of_file_names():
    """The point of this task: no setting an owner fills with real names."""
    assert not hasattr(ms, "SERIES")
    assert not hasattr(ms, "SERIES_LABELS")


def test_the_arguments_win_over_the_module_settings(series):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, offset=2.0, label_height=1.5)
    placed = sorted(t.get_position()[1] for t in fig.axes[0].texts)
    assert placed == [1.5, 3.5]
    plt_close(fig)


def test_passing_nothing_draws_what_the_settings_say(series):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    placed = sorted(t.get_position()[1] for t in fig.axes[0].texts)
    assert placed == [0.90, 2.25]
    plt_close(fig)


def test_a_nan_label_x_pins_the_labels_to_the_left_border(series,
                                                          monkeypatch):
    monkeypatch.setattr(ms, "LABEL_X", 18.0)
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, label_x=float("nan"))
    label = fig.axes[0].texts[0]
    assert 0.0 <= label.get_position()[0] < 0.1
    assert label.get_transform() is not fig.axes[0].transData
    plt_close(fig)


def test_a_label_x_argument_wins_over_the_setting(series, monkeypatch):
    monkeypatch.setattr(ms, "LABEL_X", None)
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, label_x=18.0)
    assert fig.axes[0].texts[0].get_position()[0] == 18.0
    plt_close(fig)


def test_the_window_argument_wins_over_the_module_window(series,
                                                          monkeypatch):
    monkeypatch.setattr(ms, "PLOT_X_MIN", None)
    monkeypatch.setattr(ms, "PLOT_X_MAX", None)
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, window=(15.0, 40.0))
    assert fig.axes[0].get_xlim() == (15.0, 40.0)
    plt_close(fig)


def test_the_guides_can_be_given_per_call(series, monkeypatch):
    monkeypatch.setattr(ms, "GUIDE_LINES", [])
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.1])
    guides = [line for line in fig.axes[0].lines
              if line.get_linestyle() in (":", "dotted")]
    assert len(guides) == 1
    plt_close(fig)


def test_the_ticks_can_be_turned_off_per_call(series):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, show_ticks=False)
    assert len(fig.axes[0].collections) == 0
    assert fig.axes[0].get_legend() is None
    plt_close(fig)


def test_the_axis_label_follows_the_sqrt_argument(series, monkeypatch):
    monkeypatch.setattr(ms, "USE_SQRT", True)
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, use_sqrt=False)
    assert fig.axes[0].get_ylabel() == r"Intensity / a.u."
    plt_close(fig)


def test_the_sqrt_argument_reaches_the_data_and_not_just_the_label(series):
    """A caller toggling sqrt must get the transform, not only the label.

    The square root is applied where the file is read, so a plot_series
    that relabelled the axis on its own would name a transform nobody
    applied.
    """
    linear, _phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame(), use_sqrt=False)
    rooted, _phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame(), use_sqrt=True)
    assert np.allclose(np.sqrt(linear[0][1]), rooted[0][1])


def test_the_guide_style_reaches_the_drawn_line(series):
    """A guide must be able to be dashed, not only dotted."""
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.1],
                         guide_style="--")
    dashed = [line for line in fig.axes[0].lines
              if line.get_linestyle() in ("--", "dashed")]
    assert len(dashed) == 1
    plt_close(fig)


def test_the_guide_width_reaches_the_drawn_line(series):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.1], guide_width=3.0)
    guides = [line for line in fig.axes[0].lines
              if line.get_linestyle() in (":", "dotted")]
    assert guides[0].get_linewidth() == 3.0
    plt_close(fig)


def test_the_guide_style_and_width_fall_back_to_the_settings(series,
                                                              monkeypatch):
    monkeypatch.setattr(ms, "GUIDE_STYLE", "-.")
    monkeypatch.setattr(ms, "GUIDE_WIDTH", 2.5)
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.1])
    drawn = [line for line in fig.axes[0].lines
             if line.get_linestyle() in ("-.", "dashdot")]
    assert len(drawn) == 1
    assert drawn[0].get_linewidth() == 2.5
    plt_close(fig)


def test_the_trace_width_reaches_every_drawn_trace(series):
    """The traces, not the guides: a series too thick reads as a block."""
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, linewidth=2.0)
    solid = [line for line in fig.axes[0].lines
             if line.get_linestyle() in ("-", "solid")]
    assert len(solid) == 2
    assert all(line.get_linewidth() == 2.0 for line in solid)
    plt_close(fig)


def test_the_guides_survive_the_ticks_being_turned_off(series, capsys):
    """SHOW_TICKS hides the rows and the legend, and nothing else.

    A guide carries a reflection up through the stack, which is worth doing
    whether or not a tick marks it at the bottom. Dropping the guides along
    with the ticks would take away, in silence, a control the reader set on
    purpose.
    """
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, show_ticks=False,
                         guide_lines=[20.1], guide_labels=["(110)"],
                         guide_style=":")
    ax = fig.axes[0]
    assert len(ax.collections) == 0, "no tick row should be drawn"
    assert ax.get_legend() is None, "no legend without the rows it names"
    # The name goes with the guide, and the room made for it is measured
    # against a legend that is not there.
    assert [text.get_text() for text in ax.texts
            if text.get_ha() == "center"] == ["(110)"]
    guides = [line for line in ax.lines
              if line.get_linestyle() in (":", "dotted")]
    assert len(guides) == 1, "the guide went away with the ticks"
    # The list is what a guide value is picked from, so it has to survive too.
    assert "to pick the guides from" in capsys.readouterr().out
    plt_close(fig)


def guide_names(ax):
    """The names written over the guides, told apart from the trace labels.

    A trace label starts at its own x and runs rightwards; a guide name is
    centred on the guide it belongs to, and that is enough to find them.
    """
    return [text for text in ax.texts if text.get_ha() == "center"]


def test_a_guide_name_is_written_over_its_own_guide(series):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.2],
                         guide_labels=["(110)"])
    names = guide_names(fig.axes[0])
    assert [text.get_text() for text in names] == ["(110)"]
    # On the reflection, not on the typed value: the guide snapped and the
    # name has to travel with it.
    assert names[0].get_position()[0] == pytest.approx(20.0)
    plt_close(fig)


def test_a_guide_name_takes_the_colour_of_its_phase(tmp_path):
    folder = tmp_path / "data"
    folder.mkdir()
    write_export(folder / "two.csv", 20.0, 900.0, second_phase=30.0)
    traces, phases, colors = ms.load_series(
            no_labels("two.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases, colors=colors, guide_lines=[30.1],
                         guide_labels=["(110)"])
    ax = fig.axes[0]
    # The tick rows are drawn as collections, in legend order; the second
    # phase owns the reflection at 30.0, so its row is the second one.
    assert to_rgba(guide_names(ax)[0].get_color()) == to_rgba(
        ax.collections[1].get_colors()[0])
    plt_close(fig)


def test_a_guide_with_no_name_gets_no_text(series):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.2],
                         guide_labels=[""])
    assert guide_names(fig.axes[0]) == []
    plt_close(fig)


def test_a_guide_that_does_not_snap_takes_its_name_with_it(series, capsys):
    """The name that survives is the one of the guide that survived.

    Names pair with positions by their place in the list, so a guide dropped
    for having no reflection near it must not shift the name of the next one
    onto itself.
    """
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[24.0, 20.2],
                         guide_labels=["(999)", "(110)"])
    assert [text.get_text()
            for text in guide_names(fig.axes[0])] == ["(110)"]
    assert "no reflection within" in capsys.readouterr().out
    plt_close(fig)


def test_two_names_too_close_together_are_stacked_not_overlapped(tmp_path):
    folder = tmp_path / "data"
    folder.mkdir()
    # Two reflections a fifth of a degree apart: their names cannot sit side
    # by side in a forty degree window.
    write_export(folder / "two.csv", 20.0, 900.0, second_phase=20.2)
    traces, phases, colors = ms.load_series(
            no_labels("two.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases, colors=colors,
                         guide_lines=[20.0, 20.2],
                         guide_labels=["(111)", "(110)"])
    fig.canvas.draw()
    first, second = (text.get_window_extent()
                     for text in guide_names(fig.axes[0]))
    assert not first.overlaps(second)
    plt_close(fig)


def test_a_name_sits_above_its_own_peak(series):
    """The top trace's peak here reaches 1.0, since stack() normalises it."""
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.2],
                         guide_labels=["(110)"], guide_label_height=0.10)
    assert guide_names(fig.axes[0])[0].get_position()[1] == pytest.approx(1.10)
    plt_close(fig)


def label_width_in_degrees(fig, text):
    """How much of the 2theta axis one name covers, as it is drawn."""
    fig.canvas.draw()
    box = text.get_window_extent()
    inverse = fig.axes[0].transData.inverted()
    return (inverse.transform((box.x1, 0.0))[0]
            - inverse.transform((box.x0, 0.0))[0])


def test_a_turned_name_covers_far_less_of_the_axis(series):
    """Turning a name is what stops it colliding, so this is the mechanism.

    Upright, a Miller index is some two degrees wide on a forty degree axis
    and reaches its neighbours; on end it is about as wide as one line of
    text is tall.
    """
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    flat = ms.plot_series(traces, phases, guide_lines=[20.2],
                          guide_labels=["(889)"])
    turned = ms.plot_series(traces, phases, guide_lines=[20.2],
                            guide_labels=["(889)"], guide_label_rotation=90)
    wide = label_width_in_degrees(flat, guide_names(flat.axes[0])[0])
    narrow = label_width_in_degrees(turned, guide_names(turned.axes[0])[0])
    assert narrow < wide / 2
    plt_close(flat)
    plt_close(turned)


def test_two_names_that_collide_upright_stand_side_by_side_turned(tmp_path):
    """The crowded case the turn exists for.

    Two reflections a degree and a half apart, on an axis where a name is
    two and a half degrees wide upright and one degree on end: upright the
    two have to be stacked, turned they stand side by side.
    """
    folder = tmp_path / "data"
    folder.mkdir()
    write_export(folder / "two.csv", 20.0, 900.0, second_phase=21.5)
    traces, phases, colors = ms.load_series(
            no_labels("two.csv"), folder, pd.DataFrame())
    upright = ms.plot_series(traces, phases, colors=colors,
                             guide_lines=[20.0, 21.5],
                             guide_labels=["(889)", "(123)"])
    turned = ms.plot_series(traces, phases, colors=colors,
                            guide_lines=[20.0, 21.5],
                            guide_labels=["(889)", "(123)"],
                            guide_label_rotation=90)

    def spans(fig):
        fig.canvas.draw()
        inverse = fig.axes[0].transData.inverted()
        return sorted(
            (inverse.transform((text.get_window_extent().x0, 0.0))[0],
             inverse.transform((text.get_window_extent().x1, 0.0))[0])
            for text in guide_names(fig.axes[0]))

    (_first_left, first_right), (second_left, _second_right) = spans(upright)
    assert first_right > second_left, "upright, the two names reach each other"
    (_first_left, first_right), (second_left, _second_right) = spans(turned)
    assert first_right < second_left, "turned, they stand clear side by side"
    plt_close(upright)
    plt_close(turned)


def test_the_turn_falls_back_to_its_setting(series, monkeypatch):
    monkeypatch.setattr(ms, "GUIDE_LABEL_ROTATION", 90)
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.2],
                         guide_labels=["(110)"])
    assert guide_names(fig.axes[0])[0].get_rotation() == pytest.approx(90.0)
    plt_close(fig)


def test_a_taller_peak_under_the_name_lifts_it(tmp_path):
    """The clearance is over everything the name covers, not over one point.

    A name is about two degrees wide on a forty degree axis. Measured at its
    own reflection alone, it would clear that peak and land on the taller
    one beside it, which is what a reader sees as a name stuck to a peak.
    """
    folder = tmp_path / "data"
    folder.mkdir()
    angles = np.arange(10.0, 50.0, 0.5)
    # A small peak at the reflection, a tall one under the same name.
    obs = (100.0 + 200.0 * np.exp(-((angles - 20.0) ** 2) / 0.5)
           + 4000.0 * np.exp(-((angles - 20.8) ** 2) / 0.5))
    lines = [HEADER]
    for index, (angle, value) in enumerate(zip(angles, obs)):
        reflection = "20.00" if index == 0 else ""
        lines.append(f"{angle:.2f};{value:.4f};{value:.4f};100.0;0.1;"
                     f"{reflection}")
    (folder / "pair.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    traces, phases, _colors = ms.load_series(
            no_labels("pair.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.0],
                         guide_labels=["(110)"], guide_label_height=0.10)
    ax = fig.axes[0]
    # The tall peak reaches 1.0, since stack() normalises the trace to it.
    assert guide_names(ax)[0].get_position()[1] == pytest.approx(1.10)
    plt_close(fig)


def test_a_name_on_a_reflection_the_top_trace_lost_sits_at_its_baseline(
        series):
    """A phase that went away must not have its name floating over nothing.

    The second sample's peak is elsewhere, so at the first sample's
    reflection the top trace is flat, and the name belongs down there where
    the peak is not rather than up at the height of the tallest one.
    """
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, offset=1.35, guide_lines=[20.2],
                         guide_labels=["(110)"], guide_label_height=0.10)
    baseline_of_the_top_trace = 1.35
    assert guide_names(fig.axes[0])[0].get_position()[1] == pytest.approx(
        baseline_of_the_top_trace + 0.10)
    plt_close(fig)


@pytest.fixture
def edge_series(tmp_path):
    """Two series of three, one reflection at the left, one at the right.

    Three traces because the taller the stack, the less of the figure the
    clearance above it is worth, and the right-hand one because a name there
    lands under the phase legend.
    """
    folder = tmp_path / "data"
    folder.mkdir()
    for index in range(1, 4):
        write_export(folder / f"left{index}.csv", 12.0, 900.0)
        write_export(folder / f"right{index}.csv", 48.0, 900.0)
    return folder


def test_the_frame_is_raised_to_hold_a_name_under_the_legend(edge_series):
    """A full stack has no spare room above it, so the frame has to grow."""
    traces, phases, _colors = ms.load_series(
            no_labels("right1.csv", "right2.csv", "right3.csv"), edge_series,
            pd.DataFrame())
    plain = ms.plot_series(traces, phases, guide_lines=[48.0])
    named = ms.plot_series(traces, phases, guide_lines=[48.0],
                           guide_labels=["(110)"])
    assert named.axes[0].get_ylim()[1] > plain.axes[0].get_ylim()[1]
    plt_close(plain)
    plt_close(named)


def test_the_frame_is_not_raised_when_no_guide_is_named(series):
    """A figure without names is the figure this repository already drew."""
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv", "s2.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, offset=1.35, label_height=0.90,
                         guide_lines=[20.2])
    assert fig.axes[0].get_ylim()[1] == pytest.approx(1.35 + 1.0 + 0.30 * 1.35)
    plt_close(fig)


def test_a_name_under_the_legend_never_reaches_it(tmp_path):
    """The legend hangs from the top of the frame, and the names rise to it.

    A reflection at the right of the window puts its name directly under the
    legend, which is the collision the extra room exists to prevent, so it
    is asserted on the drawn boxes rather than on the setting meant to
    avoid it.
    """
    folder = tmp_path / "data"
    folder.mkdir()
    write_export(folder / "right.csv", 48.0, 900.0)
    traces, phases, _colors = ms.load_series(
            no_labels("right.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[48.0],
                         guide_labels=["(110)"])
    ax = fig.axes[0]
    fig.canvas.draw()
    legend = ax.get_legend().get_window_extent()
    names = guide_names(ax)
    assert names, "the name was not drawn at all"
    assert not any(text.get_window_extent().overlaps(legend)
                   for text in names)
    plt_close(fig)


def test_the_room_comes_from_a_taller_figure_and_not_from_the_stack(
        edge_series):
    """Making room must not cost the patterns any of the figure.

    A series is drawn to show the small reflections, and they are the first
    thing lost when the stack is squeezed to fit something above it. The
    check is the scale itself, pixels per trace height, which has to come
    out the same whether or not the guides were named.
    """
    def scale_of(fig):
        fig.canvas.draw()
        ax = fig.axes[0]
        low, high = ax.get_ylim()
        return ax.get_window_extent().height / (high - low)

    traces, phases, _colors = ms.load_series(
            no_labels("right1.csv", "right2.csv", "right3.csv"), edge_series,
            pd.DataFrame())
    plain = ms.plot_series(traces, phases, guide_lines=[48.0])
    named = ms.plot_series(traces, phases, guide_lines=[48.0],
                           guide_labels=["(110)"])
    assert scale_of(named) == pytest.approx(scale_of(plain), rel=0.01)
    assert (named.get_size_inches()[1] > plain.get_size_inches()[1]), (
        "the figure did not grow, so the room came out of the stack")
    assert named.get_size_inches()[0] == plain.get_size_inches()[0], (
        "the width must not move: it is what the figure is scaled to on a "
        "page, and the patterns would shrink with it")
    plt_close(plain)
    plt_close(named)


def test_a_name_clear_of_the_legend_does_not_pay_for_it(edge_series):
    """Room for the legend is made only where a name stands under it.

    Reserving its height for a name at the other end of the axis opened a
    band of empty figure above every stack, which is what this refuses.
    """
    left, left_phases, _colors = ms.load_series(
            no_labels("left1.csv", "left2.csv", "left3.csv"), edge_series,
            pd.DataFrame())
    right, right_phases, _colors = ms.load_series(
            no_labels("right1.csv", "right2.csv", "right3.csv"), edge_series,
            pd.DataFrame())
    away = ms.plot_series(left, left_phases, guide_lines=[12.0],
                          guide_labels=["(110)"])
    under = ms.plot_series(right, right_phases, guide_lines=[48.0],
                           guide_labels=["(110)"])
    assert away.axes[0].get_ylim()[1] < under.axes[0].get_ylim()[1]
    plt_close(away)
    plt_close(under)


def test_a_name_on_a_guide_at_the_border_stays_inside_the_frame(series):
    """Centred on a guide at the edge, a name would hang outside the axes.

    bbox_inches='tight' would then widen the saved figure around it, so the
    name is nudged in far enough to fit and sits a little off its own guide.
    """
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, window=(20.0, 40.0),
                         guide_lines=[20.2], guide_labels=["(110)"])
    ax = fig.axes[0]
    fig.canvas.draw()
    frame = ax.get_window_extent()
    box = guide_names(ax)[0].get_window_extent()
    assert box.x0 >= frame.x0 and box.x1 <= frame.x1
    plt_close(fig)


def test_more_names_than_guides_is_reported(series, capsys):
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.2],
                         guide_labels=["(110)", "(200)"])
    assert "2 names for 1 guides" in capsys.readouterr().out
    assert len(guide_names(fig.axes[0])) == 1
    plt_close(fig)


def test_names_past_what_the_figure_can_grow_to_hold_are_reported(tmp_path,
                                                                  capsys):
    """More rows than a doubled figure can carry must not fail in silence.

    Forty-five reflections a hundredth of a degree apart cannot share a row
    at any window, so their names ask for more rows than the figure is
    allowed to grow into. The run has to say so: the alternative is a saved
    PDF with its top rows cut off and nothing to explain why.
    """
    folder = tmp_path / "data"
    folder.mkdir()
    crowded = [20.0 + step / 100 for step in range(45)]
    # A fine step keeps the reflection column well under PHASE_MAX_FILL, so
    # the parser still reads it as a phase and not as a data column. The
    # peak they all sit on puts every name at the top of the trace, which is
    # where there is no room left to climb into.
    angles = np.arange(10.0, 50.0, 0.1)
    obs = 100.0 + 900.0 * np.exp(-((angles - 20.2) ** 2) / 0.5)
    lines = [HEADER]
    for index, (angle, value) in enumerate(zip(angles, obs)):
        reflection = f"{crowded[index]:.2f}" if index < len(crowded) else ""
        lines.append(f"{angle:.2f};{value:.4f};{value:.4f};100.0;0.1;"
                     f"{reflection}")
    (folder / "crowded.csv").write_text("\n".join(lines) + "\n",
                                        encoding="utf-8")
    traces, phases, _colors = ms.load_series(
            no_labels("crowded.csv"), folder, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=crowded,
                         guide_labels=["(111)"] * len(crowded))
    assert ("more room than the figure can grow to hold"
            in capsys.readouterr().out)
    # Grown as far as it is allowed to and no further.
    assert fig.get_size_inches()[1] <= 2 * xp.FIGURE_HEIGHT
    plt_close(fig)


def test_the_name_height_is_measured_from_the_peak(series):
    """The one value is the clearance over the peak, not a height above it."""
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    on_it = ms.plot_series(traces, phases, guide_lines=[20.2],
                           guide_labels=["(110)"], guide_label_height=0.0)
    above = ms.plot_series(traces, phases, guide_lines=[20.2],
                           guide_labels=["(110)"], guide_label_height=0.40)
    peak = 1.0
    assert guide_names(on_it.axes[0])[0].get_position()[1] == pytest.approx(
        peak)
    assert guide_names(above.axes[0])[0].get_position()[1] == pytest.approx(
        peak + 0.40)
    plt_close(on_it)
    plt_close(above)


def test_the_name_height_falls_back_to_its_setting(series, monkeypatch):
    monkeypatch.setattr(ms, "GUIDE_LABEL_HEIGHT", 0.50)
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, label_height=0.90,
                         guide_lines=[20.2], guide_labels=["(110)"])
    assert guide_names(fig.axes[0])[0].get_position()[1] == pytest.approx(1.50)
    plt_close(fig)


def test_the_names_fall_back_to_the_settings(series, monkeypatch):
    monkeypatch.setattr(ms, "GUIDE_LINES", [20.2])
    monkeypatch.setattr(ms, "GUIDE_LABELS", ["(110)"])
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases)
    assert [text.get_text()
            for text in guide_names(fig.axes[0])] == ["(110)"]
    plt_close(fig)


def test_the_name_weight_is_its_own_control(series, monkeypatch):
    """The names carry their own weight, not the one the trace labels have.

    A name sits alone above the stack and a trace label sits over a pattern,
    so the reader picks what separates each of them best.
    """
    monkeypatch.setattr(ms, "GUIDE_LABEL_WEIGHT", "normal")
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.2],
                         guide_labels=["(110)"], label_weight="bold",
                         guide_label_weight="semibold")
    ax = fig.axes[0]
    assert ax.texts[0].get_fontweight() == "bold"       # the trace label
    assert guide_names(ax)[0].get_fontweight() == "semibold"
    plt_close(fig)


def test_the_name_weight_falls_back_to_its_setting(series, monkeypatch):
    monkeypatch.setattr(ms, "GUIDE_LABEL_WEIGHT", "normal")
    traces, phases, _colors = ms.load_series(
            no_labels("s1.csv"), series, pd.DataFrame())
    fig = ms.plot_series(traces, phases, guide_lines=[20.2],
                         guide_labels=["(110)"])
    assert guide_names(fig.axes[0])[0].get_fontweight() == "normal"
    plt_close(fig)


def test_a_typed_guide_carries_its_own_name():
    positions, labels, unreadable = ms.parse_guides("30.95=(111), 37.04")
    assert positions == [30.95, 37.04]
    assert labels == ["(111)", ""]
    assert unreadable == []


def test_a_bare_list_of_positions_still_reads_as_guides():
    # The syntax section 5 had before the names: every position keeps
    # working and gets an empty name.
    positions, labels, _unreadable = ms.parse_guides("20.1; 30.9")
    assert positions == [20.1, 30.9]
    assert labels == ["", ""]


def test_a_decimal_comma_is_read_as_two_guides():
    # The tooltip has always warned about this: the comma separates, so
    # '28,4' is two guides and not one at 28.4.
    positions, _labels, _unreadable = ms.parse_guides("28,4")
    assert positions == [28.0, 4.0]


def test_an_unreadable_entry_is_named_and_left_out():
    positions, labels, unreadable = ms.parse_guides("here=(110), 30.95")
    assert positions == [30.95]
    assert labels == [""]
    assert unreadable == ["here=(110)"]


def test_a_name_is_taken_exactly_as_it_was_typed():
    # Mathtext is how an overbar is written, so nothing may be stripped from
    # inside a name or wrapped around it.
    _positions, labels, _unreadable = ms.parse_guides(r"30.95=$(\bar{1}11)$")
    assert labels == [r"$(\bar{1}11)$"]


def test_the_height_of_a_reflection_is_the_top_of_its_peak():
    theta = np.array([19.8, 19.9, 20.0, 20.1, 20.2])
    values = np.array([0.1, 0.4, 0.9, 0.3, 0.1])
    assert ms.peak_height(theta, values, 20.0, 0.3) == pytest.approx(0.9)


def test_a_peak_that_moved_is_still_found_inside_the_window():
    """The ticks come from the first sample; the top trace has moved on.

    A guide sits where the reflection started, so the height it is written
    at has to be read over a window and not at that one angle, or a name
    would sit on the shoulder of the peak it belongs to.
    """
    theta = np.array([19.8, 19.9, 20.0, 20.1, 20.2])
    values = np.array([0.1, 0.2, 0.3, 0.8, 0.2])
    assert ms.peak_height(theta, values, 20.0, 0.3) == pytest.approx(0.8)
    # Outside the window the peak is not this reflection's own.
    assert ms.peak_height(theta, values, 20.0, 0.05) == pytest.approx(0.3)


def test_a_reflection_with_no_pattern_near_it_has_no_height():
    theta = np.array([19.8, 20.0, 20.2])
    values = np.array([0.1, 0.9, 0.2])
    assert ms.peak_height(theta, values, 40.0, 0.3) == 0.0


def test_two_names_that_do_not_touch_are_not_lifted():
    assert ms.stagger([(0.0, 0.0, 10.0, 5.0), (20.0, 0.0, 30.0, 5.0)],
                      [], step=6.0) == [0.0, 0.0]


def test_a_name_over_another_is_lifted_clear_of_it():
    assert ms.stagger([(0.0, 0.0, 10.0, 5.0), (5.0, 0.0, 15.0, 5.0)],
                      [], step=6.0) == [0.0, 6.0]


def test_a_name_beside_a_taller_one_stays_where_its_peak_is():
    """Two names at different heights do not have to be on the same row.

    This is the whole point of writing a name over its own peak: a name on a
    small reflection sits low, and one on a tall reflection beside it passes
    over it without either being moved.
    """
    assert ms.stagger([(0.0, 0.0, 10.0, 5.0), (5.0, 20.0, 15.0, 25.0)],
                      [], step=6.0) == [0.0, 0.0]


def test_a_name_is_lifted_clear_of_an_obstacle_it_did_not_place():
    # The trace labels are on the axes before any name is, and a name has to
    # clear them as it clears another name.
    assert ms.stagger([(0.0, 0.0, 10.0, 5.0)], [(0.0, 0.0, 10.0, 5.0)],
                      step=6.0) == [6.0]


def test_the_lowest_name_settles_first():
    """Order of lifting is by height, so the answer does not depend on typing.

    The lower of two overlapping names keeps its place and the higher one
    moves, whichever order they were given in.
    """
    high = (5.0, 3.0, 15.0, 8.0)
    low = (0.0, 0.0, 10.0, 5.0)
    assert ms.stagger([high, low], [], step=6.0) == [6.0, 0.0]
    assert ms.stagger([low, high], [], step=6.0) == [0.0, 6.0]


def shipped_setting(name):
    """One module-level constant as the repository's own file assigns it.

    Read from the source, not from the imported module: the autouse fixture
    above pins every setting to its shipped value, so an assertion about
    ms.NAME would pass on a working copy that carries a tuned one, which is
    exactly the copy a commit is made from.
    """
    source = Path(ms.__file__).read_text(encoding="utf-8")
    for node in ast.parse(source).body:
        if isinstance(node, ast.AnnAssign):
            targets = [node.target]
        elif isinstance(node, ast.Assign):
            targets = node.targets
        else:
            continue
        for target in targets:
            if isinstance(target, ast.Name) and target.id == name:
                return ast.literal_eval(node.value)
    raise AssertionError(f"{name} is no longer assigned in make_series.py")


def test_no_setting_in_degrees_is_about_to_be_committed():
    """The four settings in degrees 2theta must ship unset.

    A position in degrees is a measurement, not a preference: it gives a
    lattice spacing through Bragg's law and a handful of them identify the
    phase. This is the gate that keeps a tuned one out of a commit, and it
    reads the file rather than the module so the fixture cannot hide it.
    See docs/privacy.md.
    """
    unset = {"GUIDE_LINES": [], "PLOT_X_MIN": None, "PLOT_X_MAX": None,
             "LABEL_X": None}
    tuned = {name: shipped_setting(name) for name, blank in unset.items()
             if shipped_setting(name) != blank}
    assert not tuned, (
        f"{tuned} carries a measured 2theta position. Empty it before "
        "committing, or set it in section 5 of the notebook, where it "
        "stays in the widget. See docs/privacy.md")


def test_no_reflection_name_is_about_to_be_committed():
    """The names pair with GUIDE_LINES, so they ship empty for its reason.

    A Miller index is not a position and gives no lattice spacing on its
    own. It still names a reflection picked out of a pattern nobody else
    has, and it is written beside the position it belongs to, so the gate
    that keeps one out of a commit keeps the other out too.
    """
    assert shipped_setting("GUIDE_LABELS") == [], (
        "GUIDE_LABELS names a reflection from your own pattern. Empty it "
        "before committing, or set it in section 5 of the notebook, where "
        "it stays in the widget. See docs/privacy.md")
