# Privacy: nothing private gets published

This repository holds no experimental or personal data, and none must ever be
committed or published. Three layers enforce this.

## 1. Ignored paths

`data/`, `output/` and `Samples_metadata.csv` are excluded in
[`.gitignore`](../.gitignore). The validation suite uses synthetic patterns
only, and so do the example figures in the README.

## 2. Stripped notebook outputs

Executed cells embed their results, figures and sample names included, inside
the `.ipynb` file. The panels of sections 4 and 5 also save their figures,
their metadata line and their trace names into the notebook's widget state,
which a plain output clear does not touch, so strip both before every commit:

```bash
pip install nbstripout                          # once
python -m nbstripout XRD_Rietveld_Plotter.ipynb
```

## 3. Strip on commit instead (recommended)

With [nbstripout](https://github.com/kynan/nbstripout) installed as above,
register it once per clone. Git then strips outputs at commit time, and a
forgotten manual strip leaks nothing:

```bash
nbstripout --install        # run inside the git repository
```

## Running elsewhere

Uploading a pattern to a hosted notebook service sends it to a third party.
An embargo, a group policy or a collaboration agreement might forbid that,
and a diffraction pattern of an unpublished sample is exactly the kind of
file you do not own alone. Run locally instead. It costs one `pip install`
and your files never leave your machine.

## make_series.py

Unlike `Samples_metadata.csv`, [`make_series.py`](../make_series.py) is
tracked: it is code, not data, so it ships in the repository. Which samples
the stacked series holds, in which order and under which names, is not in
it: that comes from the `series_order` and `series_label` columns of your
own metadata file, which is excluded. What remains is the appearance of the
figure, the offsets, heights and line widths, and none of those says
anything about a sample.

Four of those settings are the exception, and they are easy to miss because
they sit among the appearance ones. **Anything in degrees 2θ is a position
you read off your own patterns, not a preference.** A reflection position is
a measurement: it gives a lattice spacing through Bragg's law, and a handful
of them identify the phase. The four:

| Setting | What a filled value says |
|---|---|
| `GUIDE_LINES` | The sharpest of the four. A guide sits exactly on a reflection, so a filled list is a measured peak table |
| `PLOT_X_MIN`, `PLOT_X_MAX` | A window framing the reflection the series is about says where that reflection is |
| `LABEL_X` | Names an empty stretch of your figure, so it says where the reflections are not. Weaker, but still a position |

All four ship unset and should stay that way. Set them in section 5 of the
notebook, where they stay in the widget and never reach a file, or, for the
window, per sample through the `x_min`/`x_max` columns of your metadata. CI
refuses a commit that carries any of them tuned, so a forgotten one is
caught rather than published.

Everything else in that file is appearance: offsets, heights, line widths,
styles, weights. Those can stay as you tuned them. Keep it that way, and
apply the rule that decides it: a setting that would name one of your files,
or carry a number you measured, belongs in the metadata or in the notebook,
not here.

## Phase names

Phase colours belong in `Samples_metadata.csv`, through the `<phase>_color`
columns, rather than in `PHASE_COLORS` in the notebook. The names printed in
the legend come from the phase column headers of your own uncommitted CSV
files, and the sample name from the `formula` column of the metadata file, so
they stay out of the repository as well. Leave `PHASE_LABELS` empty for the
same reason. See [the metadata reference](metadata.md).
