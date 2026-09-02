# testRegionSets

Asks whether a region set responds to a contrast as a whole. Two
questions are answered side by side: whether the regions of the set move
away from zero, which is a self-contained claim, and whether they move
more than the regions they are compared against, which is a competitive
one. Both are computed from the per-region statistics of the same fit,
so they never disagree with
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
on the design, the offsets or the dispersion.

## Usage

``` r
testRegionSets(
  fit,
  contrast,
  method = c("camera", "fry"),
  universe = NULL,
  matchOn = c("width", "abundance"),
  universeRatio = 5,
  interRegionCor = NULL,
  useRanks = FALSE,
  FDR = 0.05,
  adjustMethod = "BH",
  regionSets = NULL,
  carryCounts = TRUE,
  verbose = TRUE
)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- contrast:

  Contrast to test, in the syntax accepted by
  [`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
  or a named list of contrasts to run on the same fit.

- method:

  Character vector with the tests to run, among `"camera"` (competitive)
  and `"fry"` (self-contained). Default: `c("camera", "fry")`.

- universe:

  What each set is compared against in the competitive test. Default:
  `NULL`, the universe carried by the fit. A `RegionSetDE.universe`
  object, or the strings `"matched"` and `"all"`, override it and are
  built here.

- matchOn:

  Character vector with the covariates the comparison rows are matched
  on, when one has to be built here. Default: `c("width", "abundance")`.

- universeRatio:

  Numeric value with the number of comparison rows drawn per region of
  the set, when one has to be built here. Default: `5`.

- interRegionCor:

  Numeric value with the correlation between regions, used to inflate
  the variance of both the set and the rows it is compared against.
  Default: `NULL`, estimated separately for each of the two from the
  residuals of the fit, or held at 0.01 when the design leaves no
  residual to estimate it from.

- useRanks:

  Logical value to indicate whether `camera` must work on the ranks
  rather than on the statistics, which is more robust and less powerful.
  Default: `FALSE`.

- FDR:

  Numeric value with the adjusted p-value cut-off reported in the
  output. Default: `0.05`.

- adjustMethod:

  String with the multiple testing correction across the sets. Default:
  `"BH"`.

- regionSets:

  Character vector with the names of the sets to test. Default: `NULL`,
  all of them.

- carryCounts:

  Logical value to indicate whether the counts must travel inside the
  result, so that
  [`plotSetSignal`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md)
  can draw the signal without being handed the counts object again.
  Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.setResults` object, or a `RegionSetDE.setResultsList`
when `contrast` is a named list.

## Details

The effect size, not the p-value, is the primary output here. A set of
30,000 promoters tested as if its regions were independent returns a
p-value below anything a computer will print for a mean shift of 0.05
log2, which says nothing about whether the shift matters. The regions of
a set are not independent either, since neighbouring elements inside the
same domain move together, so the variance of the mean log2 fold change
is inflated by the factor `1 + (n - 1) * rho`, with `rho` estimated from
the residuals of the fit through
[`limma::interGeneCorrelation`](https://rdrr.io/pkg/limma/man/camera.html).
The confidence interval in the output carries that inflation; read it
before reading the p-value. The two tests answer different questions and
the pattern between them is informative. `camera` is competitive: it
asks whether the regions of the set moved more than the regions they are
compared against, and it is invariant to a scaling error affecting every
region equally. `fry` is self-contained: it asks whether they moved away
from zero at all, which a global shift in the mark, or a residual
normalisation error, will satisfy for every set at once. When camera
separates the sets and fry does not, the sets redistributed the signal
between them; when fry is significant everywhere and camera nowhere,
everything moved together and the normalisation deserves a second look
before the biology does. The comparison universe comes from the fit,
which built it once, and travels on into the result, so
[`plotUniverseMatching`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotUniverseMatching.md)
can check the matching afterwards without anything being kept on the
side. Passing a `RegionSetDE.universe` object, or one of the two
keywords, overrides it for this test alone. The interval on
`delta.log2FC` carries the correlation of both sides. The regions of the
comparison are no less correlated than the regions of the set, so
treating their mean as if it were known would leave the interval
narrower than the data supports, by around a factor of the square root
of two when the two sides are of similar size.

A fit with no replicates loses the self-contained test. `fry` builds a
linear model inside each set and needs a residual to measure it against,
which a design with one sample per level does not have, so it is dropped
with a message and only the competitive test runs. The correlation
between regions goes the same way: it is estimated from the residuals of
the fit, and without them it falls back to 0.01, the value `limma` uses
when nothing better is available. That number sets how much the
confidence interval is widened, so on such a fit the interval is as
assumed as the dispersion is, and `interRegionCor` is worth setting by
hand from a replicated experiment on the same assay when one exists. The
competitive test runs through
[`limma::cameraPR`](https://rdrr.io/pkg/limma/man/camera.html) on the
per-region statistics, which is what makes it work identically for the
four engines. The self-contained test needs the values themselves and is
computed on the log-CPM matrix of the fit; for `edgeR` and `DESeq2` that
matrix is a transformation of the counts rather than the quantity the
model was fitted on, so the two are close but not identical, and the
competitive test is the one to lead with.

## See also

[`testSetContrast`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md),
[`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md),
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`plotSetEffect`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetEffect.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fitRegions(counts, design = ~ replicate + condition, engine = "edgeR")

# The universe comes from the fit and travels into the result
setRes <- testRegionSets(fit, contrast = "conditionCOMBO")

plotUniverseMatching(setRes)
plotSetEffect(setRes)

# Overriding it for one test
setRes <- testRegionSets(fit, contrast = "conditionCOMBO", universe = "all")
} # }
```
