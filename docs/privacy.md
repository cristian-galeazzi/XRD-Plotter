# Privacy: nothing private gets published

This repository holds no experimental or personal data, and none must ever be
committed or published. Three layers enforce this.

## 1. Ignored paths

`data/`, `output/` and `Samples_metadata.csv` are excluded in
[`.gitignore`](../.gitignore). The validation suite uses synthetic patterns
only, and so do the example figures in the README.

## 2. Stripped notebook outputs

Executed cells embed their results, figures and sample names included, inside
the `.ipynb` file. The section 4 panel also saves its figure and metadata
line into the notebook's widget state, which a plain output clear does not
touch, so strip both before every commit:

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
tracked: it is code, not data, so it ships in the repository. But its
`SERIES` and `SERIES_LABELS` settings get filled in with real file names and
real sample names while you use it, and those are private in the same way a
metadata row is. Put both back to their placeholder values before
committing, and check with `git diff --cached` that the diff you are about
to commit carries no file name or sample name of your own.

## Phase names

Phase colours belong in `Samples_metadata.csv`, through the `<phase>_color`
columns, rather than in `PHASE_COLORS` in the notebook. The names printed in
the legend come from the phase column headers of your own uncommitted CSV
files, and the sample name from the `formula` column of the metadata file, so
they stay out of the repository as well. Leave `PHASE_LABELS` empty for the
same reason. See [the metadata reference](metadata.md).
