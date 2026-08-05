# How it works and how it was checked

## Precision

Parsing is the one place where numbers take damage, so it avoids the single
lossy step available to it. Every numeric column is read as text and
converted with Python's built-in `float()` rather than pandas' fast C
converter, which is not correctly rounded and lands up to one unit in the
last place away from the nearest double. Both decimal marks go through the
same conversion, so an export with decimal commas parses to the same doubles
as the same export with decimal points. Everything downstream is IEEE-754
double precision with no intermediate rounding. Numbers are rounded only when
printed on the figure, and √I is applied to the drawn copy alone.

## Running the suite

[`test_xrd_plotter.py`](../test_xrd_plotter.py) builds synthetic exports at
run time from an analytic pattern with a fixed seed (`default_rng(0)`). It
reads nothing from `data/`, so it runs on any machine, including a fresh
clone with an empty data folder.

```bash
pytest -q
```

Section 2 of the notebook calls the same suite, so running the notebook is a
test run as well:

```bash
python -m nbconvert --to notebook --execute --output executed.ipynb XRD_Rietveld_Plotter.ipynb
```

CI runs both on every push, on the Python floor the README claims and on the
version this is developed with.

## What the suite asserts

| Case | Purpose |
|------|---------|
| Semicolon separator, decimal commas | Bit-exact parsing of a locale export, asserted with `array_equal`, at zero tolerance |
| Comma separator, decimal points, `x, deg` header | The same doubles from the other header and locale variant |
| Both variants | Both phase columns detected and trimmed independently |
| The native export, 11 header names above 10 data fields | Refused, because the 2θ header sits above the intensities, so the figure would look right and be wrong |
| A 2θ axis running from 90° down to 10° | Still accepted: the guard rejects an axis that turns, not one that runs backwards |
| Eleven angle-header spellings, and headers that must not match | Only the angle header is flexible; a renamed `obs` is drawn flat and a phase header holding a comma is never read as the axis |
| An axis with one swapped pair, and one with repeated points | Both accepted: neither travels backwards far enough to be anything but an axis |
| A column of alternating intensities, two scans pasted into one file, a column of one repeated value | All three refused, since the distance travelled backwards is not a rounding detail |
| A smooth monotone column running 400 to 40 | Refused, outside the 0 to 180 degrees a 2θ axis occupies |
| A file of arbitrary bytes | Isolated with a reason instead of raising |
| A CSV with no 2θ column | Isolated as `2theta column not found` |
| A CSV with only 2θ and `Obs` | Isolated as `residual column 'diff/sigma' not found` |
| A CSV whose `Obs` column is all text | Isolated as `no valid data in the obs column` |
| A CSV whose `diff/sigma` column is all text | Isolated as `no valid data in the 'diff/sigma' column` |
| A CSV with 2θ, `Obs` and `diff/sigma` | `calc` and `bkg` filled with zeros, still plotted |
| A CSV with one over-long row | The row skipped and counted, every other row still bit-exact |
| A complete export, phases named `Alpha` and `Beta` | Only the two phases detected among eleven columns, `tick-pos` and `Axis-limits` left alone |
| A header repeated in the export | The `tick-pos.1` pandas invents is still blocklisted, not drawn as a phase |
| Two tick columns of different length, tail-padded with a space as the export pads them | All 8 and all 6 positions reach the drawn marks, bit-exact and in order, so a missing tick is never made here |
| `Rw`, `Rw%`, `Rw / %`, `Rwp / %`, `chi2`, `GOF` and `Tick-Pos`, beside `Rp` | The folded blocklist catches every spelling, `Rp` stays undrawn as a phase, and a name added to `NON_PHASE_COLUMNS` takes effect on the next call, not at import |
| A column too full to be a reflection list | Reported by name and left undrawn |
| A synthetic metadata row | Name and both phase fractions arrive in the legend text |
| A third phase, and two columns of one phase | Legend in alphabetical order, one tick row each, the cycle in that order, a repeated label taking one entry and one colour |
| A `<phase>_color` cell, and one holding text that is not a colour | The colour follows its phase alone and beside another, the invalid cell falls back to the cycle |
| A general and a specific `_color` key, and a column named only `_pct` | The longer key wins, the nameless column matches no phase |
| `PHASE_LABELS` and `PHASE_COLORS` set together | The renamed phase prints its new name and keeps the colour keyed to it |
| Metadata with a decimal comma, an empty `formula`, an unmatched `_pct` column, two columns matching one phase | 60,5 read as 60.5, the file name used as the legend name, both mismatches reported, no percentage guessed |
| The window unset, then the constants, then `x_min`/`x_max` in the metadata | Each overrides the one before, with the intensity axis rescaled to the window |
| A calc peak taller than the obs scatter it overshoots | The top margin follows the obs max and still clears the taller calc peak |
| Both intensity labels and the 2θ label | `Intensity / a.u.`, `√Intensity / a.u.` and `2θ / °`, the IUPAC `quantity / unit` form with no parentheses |
| The ticks of the finished figure | Neither y axis numbered nor marked, the 2θ numbers and marks drawn once, under the lower panel |
| The residual panel of both modes | One line drawn on it, the trace |
| A well-fitted pattern, one with a residual past the floor, and a raw residual | Limits symmetric about zero in all three, held at `RESIDUAL_SPAN` for the first, widened rather than clipped for the second, unfloored for the third |
| The same export read with `weighted=False` | The residual subtracted as `obs - calc` bit-exact, the panel relabelled, `WEIGHTED_RESIDUALS` left alone |
| An export whose `diff` column carries the offset GSAS-II puts on the plotted curve, 2% of the largest observed intensity | The offset column never read, the panel drawn from `obs - calc`, and a file with `obs` but no `calc` isolated by name |
| The three ways the unweighted panel can be refused: `calc` absent, both `obs` and `calc` absent, `calc` present but holding no number | Each names what is actually wrong, both absent columns reported at once, and no message sends the reader to `diff` |
| The function behind the section 4 panel | Limits applied to both axes, the metadata line returned, no setting mutated, an unreadable file raising |
| The section 4 save reusing the batch name and writer | `_sqrt`, `_linear` and the `_unweighted` suffix chosen from the toggles, both files written under that name |
| The window line the batch prints | The drawn 2theta range, marked `metadata` when a row set it and `auto` otherwise |
| The window section 4 prefills from a metadata row | The `x_min`/`x_max` pair read back, blank cells and an absent file giving the widen-to-full sentinel |
| A folder whose middle file is unusable | It is reported with its reason and the files after it are still drawn |
| A batch of a good, a broken and a good file | Each file prints its own block, and the failure is repeated in the summary at the end |

## What validation does not cover

The plotting is checked. The crystallography is yours. The weighted panel is
copied from the refinement and the unweighted one is `obs - calc`, so neither
is checked against the structure and a refinement converged on the wrong one
still produces a clean-looking figure.
