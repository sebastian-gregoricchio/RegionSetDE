# RegionSetDE.fit class

S4 class storing a model fitted on a `RegionSetDE.counts` object. The
object returned by the engine is kept as it is inside the `fit` slot, so
that any function of `edgeR`, `limma`, `variancePartition` or `DESeq2`
remains usable on it, while the design, the offsets and the samples that
produced it travel with it. The counts are carried along as well: the
per-set tests read the fit and never the counts directly, which is what
guarantees that the two levels share the same dispersion and the same
normalisation.

## Slots

- `fit`:

  List with a single element, `object`, containing the fit returned by
  the engine, together with the auxiliary objects needed to test it.

- `engine`:

  String with the engine used, one of `"edgeR"`, `"voom"`, `"dream"` and
  `"deseq2"`.

- `design`:

  Matrix of the fixed effects, with one row per sample used.

- `design.formula`:

  List holding the formula when the design was declared as such, empty
  otherwise.

- `blocking`:

  List with the blocking variable, the random effect formula and the
  consensus correlation, when they apply.

- `dispersion`:

  List with the dispersion estimates of the engine, summarised.

- `counts`:

  `RegionSetDE.counts` object on which the model has been fitted,
  restricted to the samples used.

- `universe`:

  `RegionSetDE.universe` object with the rows every set will be compared
  against at the set level. It depends on the rows and not on the
  contrast, so it is built once here and reused by every test run on
  this fit. Empty when none could be built.

- `samples`:

  Character vector with the names of the samples used.

- `counting.level`:

  String indicating whether the rows are regions (`"region"`) or tiles
  (`"tile"`).

## Author

Sebastian Gregoricchio
