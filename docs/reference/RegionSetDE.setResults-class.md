# RegionSetDE.setResults class

S4 class holding the outcome of a set level test. The table in the
`results` slot has one row per region set and carries the effect size
with its confidence interval next to the p-values, since with tens of
thousands of regions in a set the p-value stops discriminating long
before the effect size does.

## Slots

- `results`:

  Data.frame with one row per region set, or per pair of sets for a set
  contrast.

- `regionStats`:

  Data.frame with the per-region statistics the set level test was
  computed on.

- `contrast`:

  String describing the contrast that was tested.

- `contrast.groups`:

  List with the `column` of the `colData` the contrast separates and the
  two `groups` it compares, the first one being the level the fold
  change is positive for. Empty when the contrast is not a difference
  between two levels of one variable.

- `test`:

  String with the kind of test, either `"set"` or `"setContrast"`.

- `methods`:

  Character vector with the tests that were run.

- `universe`:

  `RegionSetDE.universe` object with the rows every set was compared
  against. Empty when the test was self-contained only.

- `engine`:

  String with the engine that produced the per-region statistics.

- `thresholds`:

  List with the FDR cut-off and the correction method.

- `counts`:

  `RegionSetDE.counts` object the test was computed on, so that the
  signal behind a set level claim can be drawn without carrying a second
  object around. Empty when `carryCounts` was `FALSE`.

## Author

Sebastian Gregoricchio
