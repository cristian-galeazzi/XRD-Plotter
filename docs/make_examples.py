"""Regenerate the example figures in the README from an invented pattern.

No measurement is read: the pattern below is analytic and the phases are
called Phase 1 and Phase 2, so the images show the layout and nothing else.

Run from the repository root:

    python docs/make_examples.py
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import xrd_plotter as xp

N = 2000
PEAKS_1 = [(30.95, 900), (35.89, 9000), (51.67, 2600), (61.46, 3000),
           (64.51, 2400), (76.08, 700), (84.37, 800), (87.09, 400)]
PEAKS_2 = [(21.14, 1500), (30.07, 500), (37.04, 700), (43.04, 600)]


def pattern():
    """One synthetic two-phase pattern with counting noise on the observed."""
    rng = np.random.default_rng(0)
    x = np.linspace(13.0, 85.0, N)
    bkg = np.full(N, 60.0) + 40.0 * np.exp(-(x - 13.0) / 40.0)
    calc = bkg.copy()
    for centre, height in PEAKS_1 + PEAKS_2:
        calc += height * np.exp(-((x - centre) ** 2) / (2 * 0.09 ** 2))
    obs = calc + rng.normal(0.0, 1.0, N) * np.sqrt(calc)
    return x, obs, calc, bkg


def main():
    x, obs, calc, bkg = pattern()
    phases = {"Phase 1 hkl": np.array([c for c, _ in PEAKS_1]),
              "Phase 2 hkl": np.array([c for c, _ in PEAKS_2])}
    fractions = {"phase 1": 62.0, "phase 2": 38.0}
    out = Path(__file__).resolve().parent

    for name, use_sqrt, weighted in (("example_sqrt", True, True),
                                     ("example_counts", False, False)):
        resid = ((obs - calc) / np.sqrt(calc)) if weighted else (obs - calc)
        data = {"x": x, "obs": obs, "calc": calc, "bkg": bkg, "resid": resid,
                **phases}
        args = xp.prepare_data(data, list(phases), use_sqrt=use_sqrt)
        fig = xp.create_plot(*args, "Synthetic example", fractions,
                             use_sqrt=use_sqrt, weighted=weighted)
        path = out / f"{name}.png"
        fig.savefig(path, dpi=110, bbox_inches="tight", facecolor="white")
        print(f"{path.name}: {path.stat().st_size // 1024} kB")


if __name__ == "__main__":
    main()
