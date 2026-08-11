"""Tests for the stacked series figure script."""
import matplotlib

matplotlib.use("Agg")  # no display in CI, same as test_xrd_plotter.py

import numpy as np
import pytest
from matplotlib.pyplot import close as plt_close

import make_series as ms


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


HEADER = "2theta;Obs;Calc;Bkg;diff/sigma;Phase 1"


def write_export(path, peak_at, height):
    """A synthetic GSAS-II publication export: one peak on a flat background.

    The reflection column holds one position and is padded, so the parser
    reads it as a phase rather than as a data column.
    """
    angles = np.arange(10.0, 50.0, 0.5)
    obs = 100.0 + height * np.exp(-((angles - peak_at) ** 2) / 0.5)
    lines = [HEADER]
    for i, (angle, value) in enumerate(zip(angles, obs)):
        reflection = f"{peak_at:.2f}" if i == 0 else ""
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
    traces, phases = ms.load_series(["s1.csv", "s2.csv", "s3.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    assert len(traces) == 3
    assert phases, "the reflection column should be read as a phase"


def test_a_missing_file_stops_the_run_and_names_itself(series, tmp_path):
    with pytest.raises(SystemExit, match="s9.csv"):
        ms.load_series(["s9.csv"], series, tmp_path / "absent_metadata.csv")


def test_the_trace_label_falls_back_to_the_file_stem(series, tmp_path):
    traces, _ = ms.load_series(["s2.csv"], series,
                               tmp_path / "absent_metadata.csv")
    assert traces[0][2] == "s2"


def test_the_figure_carries_one_line_and_one_label_per_sample(series,
                                                              tmp_path):
    traces, phases = ms.load_series(["s1.csv", "s2.csv", "s3.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases)
    ax = fig.axes[0]
    assert len(ax.lines) == 3
    assert sorted(t.get_text() for t in ax.texts) == ["s1", "s2", "s3"]
    plt_close(fig)


def test_the_strongest_sample_does_not_overrun_the_trace_above(series,
                                                              tmp_path):
    """s3 scatters many times harder than s1; normalisation must hide that."""
    traces, phases = ms.load_series(["s1.csv", "s2.csv", "s3.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases, offset=1.5)
    tops = [float(line.get_ydata().max()) for line in fig.axes[0].lines]
    assert tops[0] < 1.5 and tops[1] < 3.0
    plt_close(fig)


def test_the_intensity_axis_carries_no_numbers(series, tmp_path):
    traces, phases = ms.load_series(["s1.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases)
    ax = fig.axes[0]
    assert ax.get_yticklabels() == [] or not any(
        label.get_visible() for label in ax.get_yticklabels())
    plt_close(fig)


def write_metadata(path, filename, formula):
    """A one-row private metadata file, the format load_metadata expects."""
    path.write_text(f"filename;formula\n{filename};{formula}\n",
                    encoding="utf-8")


def test_a_written_label_wins_over_the_metadata_formula(series, tmp_path):
    meta = tmp_path / "meta.csv"
    write_metadata(meta, "s1.csv", "Formula From Metadata")
    traces, _ = ms.load_series(["s1.csv"], series, meta,
                               labels={"s1.csv": "Label I Wrote"})
    assert traces[0][2] == "Label I Wrote"


def test_without_a_written_label_the_metadata_formula_is_used(series,
                                                              tmp_path):
    meta = tmp_path / "meta.csv"
    write_metadata(meta, "s1.csv", "Formula From Metadata")
    traces, _ = ms.load_series(["s1.csv"], series, meta, labels={})
    assert traces[0][2] == "Formula From Metadata"


def test_a_label_written_for_another_file_does_not_leak(series, tmp_path):
    traces, _ = ms.load_series(["s1.csv"], series,
                               tmp_path / "absent_metadata.csv",
                               labels={"s2.csv": "Label I Wrote"})
    assert traces[0][2] == "s1"


def test_the_module_ships_no_labels_of_its_own():
    assert ms.SERIES_LABELS == {}


def test_each_label_sits_at_the_top_of_its_own_trace(series, tmp_path):
    traces, phases = ms.load_series(["s1.csv", "s2.csv", "s3.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases, offset=1.35)
    placed = {text.get_text(): text.get_position()[1]
              for text in fig.axes[0].texts}
    # stack() normalises every pattern to a span of 1.0, so the top of the
    # trace at index i sits at i * offset + 1.0.
    assert placed == {"s1": 1.0, "s2": 2.35, "s3": 3.7}
    plt_close(fig)


def test_the_labels_are_right_aligned_inside_the_frame(series, tmp_path):
    traces, phases = ms.load_series(["s1.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases)
    label = fig.axes[0].texts[0]
    assert label.get_ha() == "right"
    assert label.get_position()[0] <= 1.0, "a label at x > 1 is outside"
    plt_close(fig)


def test_the_topmost_label_is_not_cut_off_by_the_frame(series, tmp_path):
    traces, phases = ms.load_series(["s1.csv", "s2.csv", "s3.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases, offset=1.35)
    ax = fig.axes[0]
    highest_label = max(text.get_position()[1] for text in ax.texts)
    assert ax.get_ylim()[1] > highest_label + 0.2, (
        "no room above the top label for the text itself")
    plt_close(fig)


def test_the_phase_legend_has_no_frame(series, tmp_path):
    traces, phases = ms.load_series(["s1.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases)
    assert fig.axes[0].get_legend().get_frame_on() is False
    plt_close(fig)


def test_snap_takes_the_nearest_reflection_inside_the_tolerance():
    positions = np.array([20.0, 30.95, 37.04])
    assert ms.snap_to_reflection(31.0, positions, tolerance=0.3) == 30.95


def test_snap_returns_none_when_nothing_is_near_enough():
    positions = np.array([20.0, 37.04])
    assert ms.snap_to_reflection(31.0, positions, tolerance=0.3) is None


def test_snap_on_an_empty_reflection_list_returns_none():
    assert ms.snap_to_reflection(31.0, np.array([]), tolerance=0.3) is None


def test_a_guide_is_drawn_at_the_reflection_not_at_the_typed_value(
        series, tmp_path, monkeypatch):
    # The first file's only reflection sits at 20.0; 20.2 is inside the
    # default 0.3 tolerance, so the guide must land on 20.0.
    monkeypatch.setattr(ms, "GUIDE_LINES", [20.2])
    traces, phases = ms.load_series(["s1.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases)
    guides = [line for line in fig.axes[0].lines
              if line.get_linestyle() in (":", "dotted")]
    assert len(guides) == 1
    assert float(guides[0].get_xdata()[0]) == 20.0
    plt_close(fig)


def test_a_guide_with_no_reflection_near_it_is_reported_and_not_drawn(
        series, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(ms, "GUIDE_LINES", [24.0])
    traces, phases = ms.load_series(["s1.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases)
    guides = [line for line in fig.axes[0].lines
              if line.get_linestyle() in (":", "dotted")]
    assert guides == []
    assert "24" in capsys.readouterr().out
    plt_close(fig)


def test_no_guides_are_drawn_by_default(series, tmp_path):
    traces, phases = ms.load_series(["s1.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases)
    assert [line for line in fig.axes[0].lines
            if line.get_linestyle() in (":", "dotted")] == []
    plt_close(fig)


def test_the_ticks_are_taller_than_the_traces_are_apart_is_not_the_case(
        series, tmp_path):
    """A tick row must fit under the bottom trace without reaching it."""
    traces, phases = ms.load_series(["s1.csv"], series,
                                    tmp_path / "absent_metadata.csv")
    fig = ms.plot_series(traces, phases, offset=1.35)
    ax = fig.axes[0]
    tick_collection = ax.collections[0]
    tops = [segment[:, 1].max() for segment in
            [np.array(s) for s in tick_collection.get_segments()]]
    assert max(tops) < 0.0, "the ticks reach into the bottom trace"
    assert ax.get_ylim()[0] < min(
        np.array(s)[:, 1].min()
        for s in tick_collection.get_segments()), "the ticks are clipped"
    plt_close(fig)
