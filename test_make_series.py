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

    The reflection column holds two positions and is padded, so the parser
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
    ax = ms.plot_series(traces, phases).axes[0]
    assert ax.get_yticklabels() == [] or not any(
        label.get_visible() for label in ax.get_yticklabels())
