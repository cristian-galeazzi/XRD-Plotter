# The stacked series figure

One figure for a whole sample series: the observed patterns, stacked, with
the reflection ticks underneath. No fit, no background, no residual, since
this figure answers a different question from the per-sample ones. Which
reflection appears, moves or goes away across the series is hard to read
from a pile of separate figures and easy to read from one stack.

![Five diffraction patterns stacked one above another on a square-root intensity axis, each labelled at its left with a composition, one row of orange and one of blue reflection ticks below the lowest pattern, and two dotted vertical guides rising through the stack in the colour of the phase that owns the reflection they mark.](example_series.png)

Drawn from an invented pattern with two phases named `Phase 1` and `Phase 2`
by [`make_examples.py`](make_examples.py), so it shows the layout and no
measurement. Every reflection in it follows from Bragg's law on a cubic cell
stated in that script, so the positions can be recomputed from the source
rather than taken on trust.

There are two ways to draw it, and they share every setting:

- **Section 5 of the notebook**, which redraws in place on every change and
  writes the two files when you press *Save PDF and PNG*. Each box there
  shows its default greyed inside it and the useful range beside it, so an
  empty box is never a mystery. Use this while you are still deciding what
  the figure should look like.
- **`python make_series.py`** from the repository root, which draws it once
  from the settings written at the top of that file. Use this to reproduce a
  figure you have already settled.

## Which samples, in which order, under which names

Not in the code. Three columns of `Samples_metadata.csv`, the private file
beside the notebook that `.gitignore` excludes:

| Column | Effect |
|---|---|
| `series_order` | A number puts the sample in the series at that position, the smallest at the bottom. No number leaves it out |
| `series_label` | The name written on that trace |
| `formula` | Used for the name when `series_label` is blank |

The numbers only have to sort, so `10, 20, 30` leaves room to insert one
later without renumbering the rest. A sample left out of the series still
gets its own figure from section 3.

A cell holding something that is not a number is reported by name and that
sample is left out, rather than being dropped in silence: a typo in one cell
would otherwise remove a sample from the figure without saying so. A blank
cell is a choice and passes quietly.

Full reference for the file: [the metadata reference](metadata.md).

## The settings

Every one of them, bar `GUIDE_SNAP`, is a control in section 5, a constant
at the top of `make_series.py`, and a keyword argument of `plot_series`. The
three are the same setting reached three ways: a control left blank falls
back to the constant, and so does an argument left out.

Two units run through the table:

- **degrees** are positions on the 2θ axis;
- **trace heights** are the vertical unit of this figure. Every pattern is
  rescaled to span exactly 1.0 from its own baseline to its tallest
  reflection, so 0.5 is half a pattern tall wherever that pattern sits.

Widths are in points, matplotlib's own unit for a line.

| Control | Constant | Default | What it does |
|---|---|---|---|
| 2θ from, 2θ to | `PLOT_X_MIN`, `PLOT_X_MAX` | measured range | The window every trace is drawn in. It also sets what each pattern is rescaled by, since a pattern is normalised over the part of it the window shows. An end left unset falls back to `xp.PLOT_X_MIN`/`xp.PLOT_X_MAX`, and only then to the measured range, so a window set in section 3 reaches this figure too |
| √ intensity | `USE_SQRT` | `True` | Square root of the intensity on the drawn copy, so weak reflections survive beside a strong one. The axis label and the output file name follow |
| reflection ticks | `SHOW_TICKS` | `True` | One row of ticks per phase below the lowest trace, and the legend that names them. Off hides those two and nothing else: the guides are still drawn and the list of positions still printed |
| trace spacing | `OFFSET` | `1.35` | Baseline to baseline, in trace heights. At 1.0 a pattern touches the one above; higher spreads them out |
| trace width | `LINEWIDTH_TRACE` | `0.9` | Width of a pattern line. A whole series much thicker reads as a solid block |
| label height | `LABEL_HEIGHT` | `0.90` | How high a trace's name sits above its own baseline, in trace heights, so 1.0 is the top of the pattern. Keep it under the trace spacing, or the name lands on the trace above |
| label 2θ | `LABEL_X` | left border | Where a name starts, measured from its left edge. Unset pins every name just inside the left border, whatever window is drawn |
| label weight | `LABEL_WEIGHT` | `normal` | `normal`, `medium`, `semibold` or `bold`. A heavier name separates itself from the pattern it sits over without a box behind it |
| tick height | `TICK_HEIGHT` | `0.10` | Height of one reflection tick, in trace heights, so it does not move when the trace spacing does. The rows sit `TICK_HEIGHT + 0.02` apart, so raising it moves them apart instead of overlapping them |
| guides 2θ | `GUIDE_LINES` | none | Reflections to follow up through the stack. See below |
| guide style | `GUIDE_STYLE` | `:` dotted | `:` dotted, `--` dashed, `-.` dash-dot, `-` solid |
| guide width | `GUIDE_WIDTH` | `1.2` | Width of a guide. Past the width of the traces it crosses, a guide hides the peak it points at |
| (not a control) | `GUIDE_SNAP` | `0.3` | How far a typed guide may sit from a reflection and still snap to it, in degrees |

`plot_series` takes each of them as a keyword, so a script can draw a variant
without writing to the module:

```python
import make_series as ms
import xrd_plotter as xp

meta = xp.load_metadata(ms.METADATA_FILE)
traces, phases, colors = ms.load_series(ms.series_from_metadata(meta),
                                        ms.DATA_FOLDER, meta, use_sqrt=False)
fig = ms.plot_series(traces, phases, colors=colors, use_sqrt=False,
                     window=(13, 85), guide_style="--", guide_width=1.8)
```

`use_sqrt` appears twice on purpose. `load_series` applies the square root,
because it transforms the data; `plot_series` only names the axis. Passing
them apart would label a transform nobody took.

## The guides

A guide is a vertical line rising through the whole stack, in front of the
traces, so it is not lost in the patterns it is meant to be followed
through. Each value you give snaps to the nearest real reflection, so a
position read off a printed figure by eye still lands exactly on its tick,
and each guide is drawn in the colour of the phase that owns the reflection
it landed on. A value with nothing within `GUIDE_SNAP` degrees is reported
and not drawn, rather than being drawn where no reflection is. Turning the
ticks off does not turn the guides off: a guide carries a reflection up
through the stack whether or not a tick marks it at the bottom.

Do not read the positions off the picture. Every run prints the reflections
in the window, one line per phase, and that list is what the guides should
be picked from:

```
reflections between 13 and 85 deg, to pick the guides from:
  Phase 1: 30.95, 35.89, 51.67, 61.46, 64.51, 76.08, 84.37
  Phase 2: 21.14, 30.07, 37.04, 43.04, 48.42, 53.39
```

The list is printed whether or not the ticks are drawn, since it is what a
guide value is picked from.

A phase whose export carries two reflection columns has both listed together
under one name, so a position can appear twice in its line.

## Three things about the figure that are deliberate

- **No height may be compared with another.** Each pattern is rescaled to
  its own span, which is what keeps a weakly scattering sample visible
  beside a strong one. It also destroys any comparison of heights between
  traces, and the caption of the figure has to say so.
- **The ticks and the guides come from the first member of the series.** One
  row of ticks per phase is drawn for the whole stack rather than one per
  sample. In a series drawn precisely because a reflection moves, they mark
  where it started, which is why a guide on a moving reflection sits beside
  the peak in the upper traces instead of on it. Put a well refined pattern
  first: a first file whose reflection column holds nothing but padding wins
  that choice all the same and leaves the series with no ticks at all.
- **A member that cannot be read stops the run.** A series figure silently
  missing one sample is worse than no figure, so the run names the file and
  stops instead of leaving a gap.

## Where the figure is written

`python make_series.py` writes `output/pdf/series_stacked_XRD_analysis_sqrt.pdf`
and the matching PNG, and *Save PDF and PNG* in section 5 writes the same
pair. `_linear` replaces `_sqrt` when the square root is off, so the two
versions never overwrite each other. `OUTPUT_BASENAME` sets the first word.

## Privacy

Nothing in `make_series.py` names a sample: membership, order and names all
come from the metadata file, which is excluded from the repository.

Four settings are the exception, easy to miss because they sit among the
appearance ones. **Anything in degrees 2θ is a position read off your own
patterns**: `GUIDE_LINES`, `PLOT_X_MIN`, `PLOT_X_MAX` and `LABEL_X`. A
reflection position is a measurement, it gives a lattice spacing through
Bragg's law, and a handful of them identify the phase. Set them in section 5,
where they stay in the widget, and leave the constants as shipped. The test
suite refuses a tuned one, so a forgotten setting fails before it is pushed.
Everything else in that file is appearance and can be edited freely.

Procedure: [the privacy reference](privacy.md).
