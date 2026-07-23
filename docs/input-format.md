# Input format

## Which GSAS-II export to use

> [!IMPORTANT]
> GSAS-II has two CSV exports and one of them works here. Take the one behind
> the plot, not the one in the menu bar.

- **Use this.** On the refined powder pattern, open the publication-plot
  dialog and *Save* it with `csv` as the format. Internally this is
  `CopyRietveld2csv`. It writes one row per data point under a single header
  row. Files from it are read unmodified. You edit nothing.
- **Not this.** *Export → Powder data as → histogram CSV file* writes a block
  of quoted `"Histogram"`, `"Instparm: …"` and `"Samparm: …"` lines above the
  data, and names its columns `x, y_obs, weight, y_calc, y_bkg, Q`. The
  preamble alone makes it unreadable here. It also carries no residuals and
  no reflection positions, so editing recovers nothing.

The usable export contains, in this order: `used`, `x, 2theta (deg)`, `obs`,
`calc`, `bkg`, `diff`, then **one column per phase, labelled with the phase
name from your project**, then `tick-pos`, `diff/sigma` and `Axis-limits`.
The notebook reads what it needs and ignores the rest.

## The columns

| Column | Header must | Required | Content |
|---|---|---|---|
| 2θ | contain `2theta` or `x,`, as `x, 2theta (deg)` does | Yes | Diffraction angle in degrees |
| Observed | be `obs` | No, though the figure is empty without it | Measured pattern |
| Calculated | be `calc` | No | Refined pattern |
| Background | be `bkg` | No | Refined background |
| Residuals | be `diff/sigma`, or `diff` with `WEIGHTED_RESIDUALS = False` | Yes | Drawn in the lower panel |
| Reflection positions, one column per phase | carry the phase name from your project | No | 2θ positions for that phase's tick row |

Matching ignores case and tolerates extra text, so `obs` and `Obs` name the
same column, and both `Phase 1` and `Phase 1 hkl` are recognised.

Two columns are required, the 2θ and the residuals. A missing `obs`, `calc`
or `bkg` is filled with zeros and reported on screen, so a partial file still
draws. The residual column is the exception. Zeros there draw a flat lower
panel, which reads as a perfect fit, so the file is skipped with `residual
column 'diff/sigma' not found`.

## How the phases are found

**By elimination, not from a list of names.** The engine knows the headers
GSAS-II writes itself, `used`, `diff`, `tick-pos`, `Axis-limits`, `excluded`
and the pattern columns above, listed in `NON_PHASE_COLUMNS`. Any other
column counts as a phase when it fills at most half the rows, the share in
`PHASE_MAX_FILL`. A column above that share is named on screen and left
undrawn rather than dropped in silence, so a reflection list dense enough to
cross the ceiling tells you why it is missing. A header repeated in the
export, which pandas renames `tick-pos.1`, is matched against the blocklist
without its suffix.

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

Nothing in the parser is specific to GSAS-II. An export from another
refinement program is read as soon as its columns carry the headers above, in
any order.

## Fit statistics are not drawn

No CSV export carries the goodness of fit or a weighted profile R factor.
Reading one back from a column you added by hand would attach a number to the
figure with nothing in this notebook to check it against. Put the numbers in
the sample name instead. The `formula` field of the metadata file is free
text and reaches the legend verbatim, so `Sample 3, wR = 4.2%` renders as
written.
