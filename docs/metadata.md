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
| `series_order` | No | A number puts this sample in the stacked series at that position, the smallest at the bottom. Blank leaves it out |
| `series_label` | No | Name for this sample's trace in the stacked series, when it should differ from `formula` |
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
fragment of the name is enough. When two colour columns match one phase, the
longer name wins, so `phase_1_color` beats a general `phase_color`. Two
`_pct` columns matching one phase print no percentage instead, and say so.
Decimal commas work here too.

## Subscripts in a formula

The `formula` cell reaches the legend as it is written, and matplotlib draws
it, so anything between `$` signs is typeset as maths. Write the whole formula
inside `$\mathrm{...}$` and the letters stay upright, as a chemical formula
should be, while `_` starts a subscript and `^` a superscript. Braces group
anything longer than one character:

```
filename;formula
sample_3.csv;$\mathrm{Fe_{1-\delta}O}$
sample_4.csv;$\mathrm{MgAl_2O_4}$ standard
```

Text outside the `$` signs is drawn in the ordinary figure font, so a formula
and a plain phrase can share one cell. `PHASE_LABELS` in `xrd_plotter.py`
accepts the same notation for the phase entries of the legend.

Two things to watch. The cell holds a backslash, so keep it out of a
spreadsheet's formula bar. And when it is a *phase* name that carries maths,
set through `PHASE_LABELS` or a phase column header, its `<phase>_pct` column
matches on the printed name, so key that column to a fragment surviving the
maths, `mgal` rather than `phase 2`. A `formula` cell never affects a
percentage: it names the sample, not a phase.

## Colours stay private

Pinning a colour here keeps the phase palette in your own uncommitted file
rather than in the notebook. Each run prints what it used, so you read the
mapping back without opening the figure:

```
sample_3.csv
  phases: Phase 1 #1f77b4, Phase 2 #EE8031
  2theta window: 15 to 80 (metadata)
```

A cell holding something matplotlib does not know as a colour is reported and
skipped, and that phase falls back to the colour cycle.

`PHASE_COLORS` in `xrd_plotter.py` does the same for every sample at once,
and this file wins over it. `PHASE_LABELS` sets the legend name when the
column header is not what you want printed, and a phase renamed there takes
its colour under the new name. Both dictionaries ship empty.

## What happens when a cell is missing or wrong

- A fraction of zero, or a column you left out, prints the phase with no
  percentage rather than `0%`.
- An empty `formula` cell falls back to the file name without its extension.
- Two columns matching one phase print no percentage and say so, rather than
  picking one.
- A column matching no phase in a file is reported for that file.
- A file with no `filename` or `file` column, or one with no header row at
  all, is reported and ignored whole. Every sample then falls back to its file
  name, and the batch still runs.
- Without this file, or for a data file with no row in it, the legend falls
  back to the file name and the phase entries carry no percentage.
- Duplicate rows for one file are dropped after the first.

## Finding a window before you run the batch

Section 4 of the notebook draws one file at a time. Pick it from the
dropdown, type the 2θ and intensity limits, tick or untick the sqrt intensity
and the `diff/sigma` panel. Picking a file fills the 2θ boxes from its row
here, so it opens on the window section 3 draws. The figure redraws on every
change, a text box when you press Enter or leave it, a checkbox and the
dropdown at once.

Clear a 2θ box and that end of the axis falls back to the `xp.PLOT_X_MIN` and
`xp.PLOT_X_MAX` constants, which is the full measured range while they are
`None`. An empty intensity box leaves that end to the data.

**Save to output** writes the window on screen to `output/`, at full
resolution, under the name section 3 uses for that file. The two checkboxes
are part of that name: untick one and the figure is written as `_linear` or
`_unweighted`, so a preview never overwrites a batch figure. Every control
other than that button previews only.

Above the figure the panel prints the row for the window you landed on:

```
filename;formula;x_min;x_max
sample_3.csv;Sample 3;15;80
```

Paste it into `Samples_metadata.csv` and section 3 draws that sample this way
on every run. The intensity limits stay in the panel, since they depend on
`USE_SQRT` and describe a look rather than a fact about the sample.

## The stacked series

`series_order` is what decides the stacked figure of section 5: a sample
joins it by carrying a number, and the numbers only have to sort, so 10, 20,
30 leaves room to insert one later without renumbering the rest. A sample
with no number keeps its own figure and stays out of the series.

```
filename;formula;series_order;series_label
sample_1.csv;$\mathrm{Mg_2SiO_4}$;1;x = 0.00
sample_2.csv;$\mathrm{Mg_2SiO_4}$;2;x = 0.10
sample_9.csv;$\mathrm{SiO_2}$;;
```

![Five diffraction patterns stacked one above another on a square-root intensity axis, each labelled at its left with a composition, one row of orange and one of blue reflection ticks below the lowest pattern, and two dotted vertical guides rising through the stack in the colour of the phase that owns the reflection they mark.](example_series.png)

Five rows of a metadata file like the one above, drawn. `series_order` put
them bottom to top, `series_label` wrote `x = 0.00` to `x = 0.40` inside the
frame, and the two phases keep the names their column headers carry, here
`Phase 1` and `Phase 2`. The pattern is invented, so the figure shows the
layout and no measurement.

Everything else about that figure, the settings, the controls of section 5
and the guides: [the stacked series reference](series.md).

`series_label` exists because the two figures want different names: the
single figure has room for a formula, a stacked trace usually does not. Left
blank it falls back to `formula`, and then to the file name. It accepts the
same maths notation as `formula`.

A cell holding something that is not a number is reported by name and that
sample is left out, rather than being dropped in silence: a typo in one cell
would otherwise remove a sample from the figure without saying so. A blank
cell is a choice and is passed over quietly.

Nothing about your samples reaches [`make_series.py`](../make_series.py),
which holds the appearance of the figure and nothing else.
