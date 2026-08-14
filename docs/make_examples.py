"""Regenerate the example figures in the README from an invented pattern.

No measurement is read. The two phases are called Phase 1 and Phase 2, and
their reflections are computed rather than written down: each one follows
from Bragg's law applied to an invented cubic cell stated below, so a reader
can recompute every position in the figures from three numbers on this page
and confirm that nothing here was read off a diffractometer. The intensities
are made up outright.

The same invented material makes the series of five, which is written out as
CSV exports first and then read back through make_series, so the example
figure comes off the path a real run takes rather than off a shortcut.

Run from the repository root:

    python docs/make_examples.py
"""
import math
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import make_series as ms
import xrd_plotter as xp

N = 2000
WAVELENGTH = 1.5406  # Cu K-alpha 1, in angstrom
# Phase 1 is a face-centred cell, so only the unmixed hkl survive; Phase 2 is
# primitive, so all of them do. Two cells close enough that their strongest
# reflections nearly coincide, which is what gives the series figure a pair
# of guides worth drawing side by side.
CELL_1, HKL_1 = 5.00, ((1, 1, 1), (2, 0, 0), (2, 2, 0), (3, 1, 1),
                       (2, 2, 2), (4, 0, 0), (3, 3, 1))
CELL_2, HKL_2 = 4.20, ((1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 0, 0),
                       (2, 1, 0), (2, 1, 1))
HEIGHTS_1 = (9000, 2600, 3000, 2400, 700, 900, 400)
HEIGHTS_2 = (1500, 500, 700, 600, 300, 250)


def two_theta(cell: float, hkl: tuple[int, int, int]) -> float:
    """Bragg angle of one reflection of a cubic cell, in degrees 2theta.

    >>> round(two_theta(5.00, (1, 1, 1)), 2)
    30.95
    >>> two_theta(5.00, (2, 0, 0)) > two_theta(5.00, (1, 1, 1))
    True
    """
    spacing = cell / math.sqrt(sum(index ** 2 for index in hkl))
    return 2.0 * math.degrees(math.asin(WAVELENGTH / (2.0 * spacing)))


PEAKS_1 = [(two_theta(CELL_1, hkl), height)
           for hkl, height in zip(HKL_1, HEIGHTS_1)]
PEAKS_2 = [(two_theta(CELL_2, hkl), height)
           for hkl, height in zip(HKL_2, HEIGHTS_2)]

# The invented story the series tells, one entry per member: how far the
# Phase 1 reflections have moved, how much Phase 2 there is, and the overall
# count rate. The count rate varies on purpose, so the figure shows the
# per-trace rescaling doing its work.
SERIES = (("x = 0.00", 0.00, 0.05, 0.4),
          ("x = 0.10", 0.07, 0.25, 1.0),
          ("x = 0.20", 0.14, 0.55, 0.7),
          ("x = 0.30", 0.21, 0.80, 2.5),
          ("x = 0.40", 0.28, 1.00, 1.3))


def pattern(shift: float = 0.0, second: float = 1.0, scale: float = 1.0
            ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """One synthetic two-phase pattern with counting noise on the observed.

    'shift' moves the Phase 1 reflections, as a solid solution would;
    'second' scales Phase 2, so it can be made to appear across a series;
    'scale' is the count rate, which differs between real measurements.

    >>> x, obs, calc, bkg = pattern()
    >>> len(x) == len(obs) == 2000
    True
    """
    rng = np.random.default_rng(0)
    x = np.linspace(13.0, 85.0, N)
    bkg = np.full(N, 60.0) + 40.0 * np.exp(-(x - 13.0) / 40.0)
    calc = bkg.copy()
    for centre, height in PEAKS_1:
        calc += height * np.exp(-((x - centre - shift) ** 2)
                                / (2 * 0.09 ** 2))
    for centre, height in PEAKS_2:
        calc += second * height * np.exp(-((x - centre) ** 2)
                                         / (2 * 0.09 ** 2))
    calc *= scale
    obs = calc + rng.normal(0.0, 1.0, N) * np.sqrt(calc)
    return x, obs, calc, bkg


def write_export(path: Path, x: np.ndarray, obs: np.ndarray,
                 calc: np.ndarray, bkg: np.ndarray, shift: float) -> None:
    """Write one pattern in the format the GSAS-II publication plot saves.

    The reflection columns hold one position per row and are blank below,
    which is what marks them as phases rather than as data columns.
    """
    reflections_1 = [c + shift for c, _h in PEAKS_1]
    reflections_2 = [c for c, _h in PEAKS_2]
    lines = ["2theta;Obs;Calc;Bkg;diff/sigma;Phase 1;Phase 2"]
    for i in range(len(x)):
        resid = (obs[i] - calc[i]) / np.sqrt(calc[i])
        first = f"{reflections_1[i]:.4f}" if i < len(reflections_1) else ""
        second = f"{reflections_2[i]:.4f}" if i < len(reflections_2) else ""
        lines.append(f"{x[i]:.4f};{obs[i]:.4f};{calc[i]:.4f};{bkg[i]:.4f};"
                     f"{resid:.4f};{first};{second}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_series(data_folder: Path, metadata_path: Path) -> None:
    """Write the five invented exports and the metadata that orders them.

    The metadata is what puts the series in order and names each trace, the
    same two columns a real run uses, so the folder this writes is a
    complete stand-in for a private one.

    >>> import tempfile
    >>> with tempfile.TemporaryDirectory() as d:
    ...     folder = Path(d) / "data"
    ...     folder.mkdir()
    ...     write_series(folder, Path(d) / "Samples_metadata.csv")
    ...     len(list(folder.glob("*.csv")))
    5
    """
    data_folder.mkdir(parents=True, exist_ok=True)
    rows = ["filename;formula;series_order;series_label"]
    for index, (label, shift, second, scale) in enumerate(SERIES, start=1):
        name = f"series_{index}.csv"
        write_export(data_folder / name, *pattern(shift, second, scale),
                     shift)
        rows.append(f"{name};Synthetic example;{index};{label}")
    metadata_path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def hkl_label(hkl: tuple[int, int, int]) -> str:
    """A Miller index as it is written on a figure: (1, 1, 1) -> '(111)'.

    Built from the same tuple the position is computed from, so the name on
    a guide and the guide under it can never drift apart.

    >>> hkl_label((1, 1, 1))
    '(111)'
    >>> hkl_label((2, 0, 0))
    '(200)'
    """
    return "(" + "".join(str(index) for index in hkl) + ")"


def main():
    x, obs, calc, bkg = pattern()
    phases = {"Phase 1 hkl": np.array([c for c, _ in PEAKS_1]),
              "Phase 2 hkl": np.array([c for c, _ in PEAKS_2])}
    fractions = {"phase 1": 62.0, "phase 2": 38.0}
    out = Path(__file__).resolve().parent

    for name, use_sqrt, weighted in (("example_sqrt", True, True),
                                     ("example_linear", False, False)):
        resid = ((obs - calc) / np.sqrt(calc)) if weighted else (obs - calc)
        data = {"x": x, "obs": obs, "calc": calc, "bkg": bkg, "resid": resid,
                **phases}
        args = xp.prepare_data(data, list(phases), use_sqrt=use_sqrt)
        fig = xp.create_plot(*args, "Synthetic example", fractions,
                             use_sqrt=use_sqrt, weighted=weighted)
        path = out / f"{name}.png"
        fig.savefig(path, dpi=110, bbox_inches="tight", facecolor="white")
        print(f"{path.name}: {path.stat().st_size // 1024} kB")

    # A temporary folder, not data/: the exports exist only to be read back
    # through the same calls a real run makes, and the repository keeps no
    # CSV of its own.
    with tempfile.TemporaryDirectory() as workdir:
        folder = Path(workdir) / "data"
        metadata = Path(workdir) / "Samples_metadata.csv"
        write_series(folder, metadata)
        meta = xp.load_metadata(metadata)
        traces, phases, colors = ms.load_series(ms.series_from_metadata(meta),
                                                folder, meta, use_sqrt=True)
        # Every setting given here, none read off the module: the copy of
        # make_series.py an owner runs carries their own tuning, and an
        # example figure drawn through it would publish that tuning. Two
        # guides, one owned by each phase, so the example also shows a guide
        # taking the colour of the tick row it came from. Both are named
        # from the same hkl their positions are computed from, so the figure
        # shows the labels without quoting anything measured.
        fig = ms.plot_series(traces, phases, colors=colors,
                             offset=1.35, label_height=0.90,
                             label_x=float("nan"), label_weight="bold",
                             tick_height=0.10, show_ticks=True,
                             guide_lines=[two_theta(CELL_1, (1, 1, 1)),
                                          two_theta(CELL_2, (1, 1, 0))],
                             guide_labels=[hkl_label((1, 1, 1)),
                                           hkl_label((1, 1, 0))],
                             guide_label_weight="bold",
                             guide_label_height=0.20,
                             guide_label_rotation=90,
                             guide_snap=0.3,
                             guide_style=":", guide_width=1.2,
                             linewidth=0.9, use_sqrt=True,
                             window=(None, None))
    path = out / "example_series.png"
    fig.savefig(path, dpi=110, bbox_inches="tight", facecolor="white")
    print(f"{path.name}: {path.stat().st_size // 1024} kB")


if __name__ == "__main__":
    main()
