# Sample metadata

Your file names are usually internal sample codes, and the figure needs a
real name and the phase fractions. Both come from an optional
`Samples_metadata.csv` placed **beside the notebook**, one row per data file,
matched on the exact file name.

Keeping this file outside `data/` is deliberate. It is the one place where
sample identities live, so it stays a single file you control, and
[`.gitignore`](../.gitignore) excludes it by name.

## The columns

| Column | Required | Content |
|--------|----------|---------|
| `filename` (or `file`) | Yes | Name of the data file this row describes |
| `formula` | No | Display name for the figure legend |
| `<phase>_pct`, one per phase | No | Fraction (%) appended to that phase's legend entry |
| `<phase>_color`, one per phase | No | Tick and legend colour for that phase, as `#RRGGBB` or a matplotlib colour name |
| `x_min`, `x_max` | No | 2θ window for this sample alone, in degrees |

A file with two phases and a fixed window looks like this:

```
filename;formula;phase_1_pct;phase_1_color;phase_2_pct;phase_2_color;x_min;x_max
sample_3.csv;Sample 3;62;#1f77b4;38;#EE8031;15;80
```

Name a fraction or colour column after its phase. `phase_1_pct` feeds the
`Phase 1` entry, `phase_1_color` paints it. The part before the suffix is
matched inside the legend name, with underscores read as spaces, so a
fragment of the name is enough. When two columns match one phase, the longer
name wins, so `phase_1_color` beats a general `phase_color`. Decimal commas
work here too.

## Colours stay private

Pinning a colour in the metadata keeps the phase palette in your own
uncommitted file rather than in the notebook. Each run prints what it used,
so you read the mapping back without opening the figure:

```
[sample_3.csv] phases detected: Phase 1 hkl, Phase 2 hkl
  colours: Phase 1 #1f77b4, Phase 2 #EE8031
```

A cell holding something matplotlib does not know as a colour is reported and
skipped, and that phase falls back to the colour cycle.

`PHASE_COLORS` in `xrd_plotter.py` does the same for every sample at once,
and the metadata wins over it. `PHASE_LABELS` sets the legend name when the
column header is not what you want printed, and a phase renamed there takes
its colour under the new name. Both dictionaries ship empty.

## What happens when a cell is missing or wrong

- A fraction of zero, or a column you left out, prints the phase with no
  percentage rather than `0%`.
- An empty `formula` cell falls back to the file name.
- Two columns matching one phase print no percentage and say so, rather than
  picking one.
- A column matching no phase in a file is reported for that file.
- Without the metadata file, or for a data file with no row in it, the legend
  falls back to the file name and the phase entries carry no percentage.
- Duplicate rows for one file are dropped after the first.

## Finding a window before you run the batch

Section 4 of the notebook draws one file at a time. Pick it from the
dropdown, type the 2θ and intensity limits, tick or untick the sqrt intensity
and the `diff/sigma` panel, press **Apply**. An empty box leaves that end of
the axis to the setting behind it, the `xp.PLOT_X_MIN` and `xp.PLOT_X_MAX`
constants for 2θ and the data itself for the intensity. Nothing is written to
`output/`.

Under the figure the panel prints the metadata line for the window you landed
on:

```
filename;formula;x_min;x_max
sample_3.csv;Sample 3;15;80
```

Paste it into `Samples_metadata.csv` and section 3 draws that sample this way
on every run. The intensity limits stay in the panel, since they depend on
`USE_SQRT` and describe a look rather than a fact about the sample.
