# Input format

## Which GSAS-II export to use

> [!IMPORTANT]
> GSAS-II has two CSV exports. One of them works here, after you check its
> header by hand. Take the one behind the plot, not the one in the menu bar.

- **Use this.** On the refined powder pattern, open the publication-plot
  dialog and *Save* it with `csv` as the format. Internally this is
  `CopyRietveld2csv`. It writes one row per data point under a single header
  row. Check the header as below, then the file is read as it is.

  > [!WARNING]
  > Save it from a plot on a linear intensity axis. `obs` is copied from the
  > data, but `calc` and `bkg` are copied from the drawn lines, so a
  > square-root plot exports the square roots of those two beside an untouched
  > `obs`. Nothing in the file marks it, and the symptom is a calculated curve
  > drawn far below the observed one. (`diff` is spoiled the same way, and is
  > the one column this notebook never reads.)

- **Not this.** *Export → Powder data as → histogram CSV file* writes a block
  of quoted `"Histogram"`, `"Instparm: …"` and `"Samparm: …"` lines above the
  data, and names its columns `x, y_obs, weight, y_calc, y_bkg, Q`. The
  preamble alone makes it unreadable here. It also carries no residuals and
  no reflection positions, so editing recovers nothing.

## Check the header first

> [!WARNING]
> The usable export usually writes one header name more than it writes data
> fields, so every name sits one column to the left of the data it describes.
> The angles arrive under `used` and the intensities arrive under
> `x, 2theta (deg)`. A file left this way draws a figure that looks convincing
> and is wrong.

The header it writes:

```
used;x, 2theta (deg);obs;calc;bkg;diff;<Phase A>;<Phase B>;tick-pos;diff/sigma;Axis-limits
```

**Usually, not always.** In `CopyRietveld2csv`, the exporter behind that
dialog, the name `used` is written and the column itself is built inside a
bare `try`, from the mask of the x array:

```python
lblList.append('used')
try:
    valueList.append([0 if i else 1 for i in savedX.mask])
except:
    pass
```

A histogram with no excluded points carries no mask array, only the scalar
`False` numpy uses to mean nothing is masked, and iterating that raises. The
bare `except` swallows it, the name stays, the column never arrives, and every
following name sits one place left of its data. Exclude a region in the
refinement and the same export comes out aligned. So the file in front of you
may be shifted or not, and reading its first row is the only way to know.

**Look at the first cell of the first data row.** It tells the two apart on
its own, because `used` holds a flag and the angle column holds your starting
angle:

```
used;x, 2theta (deg);obs;calc;bkg;diff;...
1;10.0;4200.0;4185.0;4180.0;-585.0;...     aligned, change nothing
10.0;4200.0;4185.0;4180.0;-585.0;...       shifted, fix the header row
```

| The first cell reads | The file is | Do this |
|---|---|---|
| `1` or `0` | aligned | Nothing. Save it into `data/` as it is |
| your lowest 2θ, near 10 or 20 | shifted | Open it in a spreadsheet, delete the `used` cell of the **header row only**, shift the rest of that row one place left, save as CSV |

Do not delete the `used` cell of an aligned file. That creates the shift this
section exists to undo.

**Then read one whole data row across and check every name against the values
under it.** The angle column climbs through your scan range, `obs` and `calc`
hold intensities of the same size, `bkg` is the small smooth one, and a phase
column is nearly empty. This takes a moment and catches an export that shifted
a different way from the one described here.

Deleting the columns the figure never uses is safe and keeps the file
readable, `diff`, `tick-pos` and `Axis-limits` among them. The 2θ column has
to survive, and with it either `diff/sigma` or the `obs` and `calc` pair,
depending on which panel you draw.

The notebook refuses a file it still finds shifted, rather than drawing it:

```
FAILED: the 'x, 2theta (deg)' column does not run from one end of the scan
to the other, so it holds something other than angles and the header names
do not line up with the columns they sit above. See docs/input-format.md
```

Two properties a 2θ axis always has and a column of intensities never has are
checked. It runs one way from the first row to the last, measured by distance
travelled, so a handful of repeated or out-of-order points is tolerated while
a column that doubles back is not. And it lies inside 0 to 180 degrees, which
catches a shift that lands a smooth column such as the background under the
2θ header:

```
FAILED: the 'x, 2theta (deg)' column runs from 40 to 400, outside the 0 to
180 degrees a 2theta axis occupies, so it holds something other than angles.
See docs/input-format.md
```

Nothing downstream recovers a shifted file, so the run stops on it and
carries on with the others.

Once it lines up the export holds, in this order: `x, 2theta (deg)`, `obs`,
`calc`, `bkg`, `diff`, then **one column per phase, labelled with the phase
name from your project**, then `tick-pos`, `diff/sigma` and `Axis-limits`.
A file that arrived aligned carries `used` ahead of all of them. The notebook
reads what it needs and ignores the rest.

## Where GSAS-II gets each column

`CopyRietveld2csv`, behind *Save* in the publication-plot dialog, writes some
columns from the refinement's arrays and others from the lines drawn on the
plot. The difference decides which ones can be trusted as numbers, because a
drawn line carries whatever the plot did to it.

| Column | Written from | Holds |
|---|---|---|
| `used` | the mask of the x array | 1 where the point was refined, 0 where it was excluded |
| `x, 2theta (deg)` | the x array | the angle, or TOF, Q, d-spacing or energy on another instrument |
| `obs` | the data array | the observed intensity, untouched |
| `calc` | the drawn line | the calculated intensity as plotted |
| `bkg` | the drawn line | the background as plotted |
| `diff` | the drawn line | `Iobs - Icalc` shifted down by the plot's offset |
| one column per phase | the drawn tick row | the reflection positions of that phase |
| `tick-pos` | the tick rows | the height each row was drawn at, a y value and not an angle |
| `diff/sigma` | recomputed from the arrays | `(Iobs - Icalc) * sqrt(wtFactor * w)` |
| `excluded` | the mask | present only when points were excluded |
| `Axis-limits` | the axes | the four plot limits |
| `mag`, `initial-mag` | the plot annotations | present only when a magnification region was used |

Only six are read here: the angle, `obs`, `calc`, `bkg`, `diff/sigma` and the
phase columns. `GOF`, `Rw` and the rest of the fit statistics are not written
by this exporter at all. They appear in a file because someone pasted them in,
which is why they sit on the blocklist rather than in this table.

`w` in the `diff/sigma` row is the weight of the point, `1/σ²`, so that column
is the residual in units of the uncertainty of that point. On a pattern of raw
counts `σ` is `√Iobs`, and `wtFactor` is 1 unless the histogram was reweighted
by hand.

## The columns

| Column | Header must | Required | Content |
|---|---|---|---|
| 2θ | name the angle: `2theta`, `2 theta`, `2θ`, `two-theta`, `x` or `x, 2theta (deg)` | Yes | Diffraction angle in degrees |
| Observed | be `obs` | No, though the figure is empty without it, and yes when `WEIGHTED_RESIDUALS` is `False` | Measured pattern |
| Calculated | be `calc` | No, and yes when `WEIGHTED_RESIDUALS` is `False` | Refined pattern |
| Background | be `bkg` | No | Refined background |
| Residuals | be `diff/sigma` | Only when `WEIGHTED_RESIDUALS` is `True`, where it is the one pattern column that is not filled with zeros when absent | Drawn in the lower panel |
| Reflection positions, one column per phase | carry the phase name from your project | No | 2θ positions for that phase's tick row |

Case never matters and surrounding spaces are trimmed, so `Obs`, `obs` and
` OBS ` are one column. Beyond that the angle is the tolerant one and the
rest are exact:

| Header | Also accepted | Not accepted |
|---|---|---|
| the angle | `2theta`, `2 theta`, `2-theta`, `two-theta`, `2θ`, `2theta_deg`, `x`, `x, deg`, `x, 2theta (deg)` | `theta`, `angle` |
| `obs` | nothing else | `observed`, `y_obs`, `Iobs` |
| `calc` | nothing else | `calculated`, `y_calc` |
| `bkg` | nothing else | `background` |
| `diff/sigma` | nothing else | `diff / sigma`, `difference` |

**So keep the names of the pattern columns.** A renamed `obs`, `calc` or
`bkg` is filled with zeros and drawn flat, reported on screen but drawn all
the same, and a renamed `diff/sigma` stops the file. Only the angle header
is free, which is what lets an export from another program run unedited.

**Rename the phase columns.** A phase column arrives with whatever GSAS-II
had for it, the name you gave the phase in the project or the file it came
from. Put the real phase there and it goes straight into the legend, with no
setting to change: `Phase 1` becomes `Rutile` by editing that header alone.
A trailing `hkl` is dropped, so `Rutile hkl` prints the same.

A missing `obs`, `calc` or `bkg` is filled with zeros and reported on screen,
so a partial file still draws. The residual is the exception. Zeros there draw
a flat lower panel, which reads as a perfect fit, so a file that cannot supply
one is skipped with `residual column 'diff/sigma' not found`.

**The raw residual is not read from the export at all.** With
`WEIGHTED_RESIDUALS = False` the lower panel is `obs - calc`, subtracted here,
and the `diff` column is ignored wherever it appears. GSAS-II copies that
column off the plotted difference curve, which the plot holds below the
pattern so the two do not overlap: 2% of the largest observed intensity by
default, any value at all once you drag the curve, and the difference of the
square roots in a square-root plot. The same shift sits on every point of the
column. `diff/sigma` carries none of it, because GSAS-II recomputes that one
from the data rather than reading it off the plot. Drawing `diff` would put
the trace off the zero the panel holds at its middle. That mode needs `obs`
and `calc` instead, and says so when one of them is absent:

```
FAILED: the raw residual is computed as obs - calc, and the 'calc' column
was not found
```

## How the phases are found

**By elimination, not from a list of names.** The engine knows the headers
GSAS-II writes itself, `used`, `diff`, `tick-pos`, `tick pos`, `Axis-limits`,
`axis limits`, `excluded` and the pattern columns above, the per-point
uncertainties, `weight`, `weights`, `sig` and `sigma`, and the fit statistics
people keep beside the pattern, `GOF`, `Rw`, `Rwp`, `Rexp` and `chi2`. All of
them are in `NON_PHASE_COLUMNS`, and `Rp` deliberately is not, since a header
that short is as likely to abbreviate a phase. Any other column counts as a
phase when it fills at most half the rows, the share in `PHASE_MAX_FILL`. A
column above that share is named on screen and left undrawn rather than
dropped in silence, so
a reflection list dense enough to cross the ceiling tells you why it is
missing. A header repeated in the export, which pandas renames `tick-pos.1`,
is matched against the blocklist without its suffix.

Blocklist entries are compared folded, ignoring case, spaces and the
punctuation people add, so the single entry `rw` covers `Rw`, `Rw%` and
`Rw / %`. A phase whose name merely begins the same way, `Rw phase`, is
still drawn.

> [!WARNING]
> Every blocklist entry is a name a phase itself could carry, `sigma` for the
> intermetallic among them. A phase column headed exactly `sigma`, `weight`,
> `GOF` or any other entry is read as a bookkeeping column and left undrawn.
> Rename that header, `Sigma phase` is enough, since only an exact folded
> match is blocked. The `phases:` line each run prints is where you catch it.

To keep a column of your own out of the legend, add its header:

```python
xp.NON_PHASE_COLUMNS = xp.NON_PHASE_COLUMNS | {"scan note"}
```

Your phase names need no configuration. Name them as you like in GSAS-II and
they arrive in the legend. Every run prints the phases it found, per file, so
a column read the wrong way shows up at once.

The phase columns are shorter than the pattern. They hold their few
reflection positions, the remaining cells stay blank, and each column is
trimmed on its own.

## What the parser tolerates

| Feature | Behaviour |
|---------|-----------|
| Field separator | `;` and `,` are both detected |
| Decimal mark | `.` and `,` are both accepted, so a locale export needs no editing |
| Text encoding | UTF-8, with Latin-1 fallback |
| Unreadable or malformed file | Reported with its reason and skipped, and the rest of the batch still runs |
| Rows with more fields than the header | Skipped and counted on screen, and the rest of the file keeps its precision |
| Non-numeric cells inside a valid column | Become NaN, counted, reported, and masked out of the plot |
| A 2θ or `obs` column with no numeric cell at all | The file is isolated with that reason, since nothing is left to draw |
| A header row left one place out of step | The file is isolated, because the 2θ column turns instead of running one way |

Nothing in the parser is specific to GSAS-II. An export from another
refinement program is read as soon as its columns carry the headers above, in
any order.

## Phase fractions come from you

No export carries them. The weight fractions the refinement produced are read
off the GSAS-II phase table and typed into the `<phase>_pct` columns of
`Samples_metadata.csv`, one row per data file. See
[the metadata reference](metadata.md). A phase with no fraction is drawn and
named in the legend without a percentage, so the figure is complete either
way.

## Fit statistics are not drawn

No CSV export carries the goodness of fit or a weighted profile R factor.
Reading one back from a column you added by hand would attach a number to the
figure with nothing in this notebook to check it against. Put the numbers in
the sample name instead. The `formula` field of the metadata file is free
text and reaches the legend verbatim, so `Sample 3, wR = 4.2%` renders as
written.
