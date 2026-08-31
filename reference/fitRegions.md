# fitRegions

Fits a linear model on the counts of a `RegionSetDE.counts` object, one
model per row. The row is the region, or the tile when the counts were
tiled, and the model is the same one that the set level tests will read
later, so the two levels never disagree on the design, on the offsets or
on the dispersion.

## Usage

``` r
fitRegions(
  counts,
  design,
  engine = "edgeR",
  samples = NULL,
  block = NULL,
  random = NULL,
  assay = "counts",
  useOffsets = TRUE,
  dispersion = NULL,
  nullSource = "background",
  nullRegionSets = NULL,
  robust = TRUE,
  universe = "matched",
  matchOn = c("width", "abundance"),
  universeRatio = 5,
  BPPARAM = NULL,
  verbose = TRUE
)
```

## Arguments

- counts:

  `RegionSetDE.counts` object.

- design:

  Formula evaluated on the `colData`, e.g. `~ replicate + condition`,
  the same formula written as a string, e.g.
  `"~ replicate + condition"`, or a design matrix with one row per
  sample. For the `"dream"` engine the formula must be given as a
  formula or a string and may contain random terms, e.g.
  `~ condition + (1|donor)`.

- engine:

  String with the model to fit, one of `"edgeR"` (quasi-likelihood
  negative binomial), `"voom"` (limma on log-CPM with precision
  weights), `"dream"` (limma with random effects) and `"deseq2"`
  (negative binomial Wald test). Default: `"edgeR"`.

- samples:

  Character vector with the names of the samples to keep, or conditions
  already applied with
  [`selectSamples`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/selectSamples.md).
  Default: `NULL`, all the samples.

- block:

  String with the name of a `colData` column holding a blocking variable
  whose effect is estimated as a correlation rather than as a
  coefficient. Only for `engine = "voom"`. Default: `NULL`.

- random:

  Formula with the random terms, e.g. `~ (1|donor)`, appended to
  `design`. Only for `engine = "dream"`. Default: `NULL`.

- assay:

  String with the name of the assay holding the values to model.
  Default: `"counts"`.

- useOffsets:

  Logical value to indicate whether the normalisation stored in the
  object must enter the model as offsets. Default: `TRUE`.

- dispersion:

  Dispersion to fit with instead of taking it from the residual
  variation, which is what a design with no replicates needs. A numeric
  value, the list returned by
  [`estimateNullDispersion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md),
  or one of the strings `"background"` and `"regionSet"`, which run that
  function here. Only for `engine = "edgeR"`. Default: `NULL`, from the
  residual variation when there is any and from `nullSource` when there
  is none.

- nullSource:

  String with where the null rows come from when a dispersion has to be
  estimated, either `"background"` or `"regionSet"`. Default:
  `"background"`.

- nullRegionSets:

  Character vector with the names of the sets used as null rows. Only
  for `nullSource = "regionSet"`. Default: `NULL`.

- robust:

  Logical value to indicate whether the dispersion, or the prior
  variance, must be estimated robustly against outlier regions. Default:
  `TRUE`.

- universe:

  What every region set will be compared against at the set level.
  Either the string `"matched"`, which builds a universe matched on
  `matchOn` through
  [`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md),
  the string `"all"`, which takes every other row as it is, `NULL` to
  build none, or a `RegionSetDE.universe` object. Default: `"matched"`.

- matchOn:

  Character vector with the covariates the comparison rows are matched
  on, among `"width"` and `"abundance"`. Default:
  `c("width", "abundance")`.

- universeRatio:

  Numeric value with the number of comparison rows drawn per region of
  the set. Default: `5`.

- BPPARAM:

  `BiocParallelParam` object passed to the `"dream"` engine. Default:
  `NULL`, sequential.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.fit` object.

## Details

Only one assay type belongs in a model. Marks, antibodies and assays
differ in dispersion, in dynamic range and in what their scaling factors
mean, so a fit that pools them borrows information between rows that
describe unrelated experiments. Split the object with
[`splitSamples`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitSamples.md)
and normalise each piece on its own before coming here. Offsets are read
from the object rather than recomputed. When
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
produced a matrix of log offsets, in the `offset` assay, that matrix is
used; otherwise the per-sample `scaling.factor` is expanded into one. In
both cases
[`edgeR::scaleOffset`](https://rdrr.io/pkg/edgeR/man/scaleOffset.html)
puts the offsets back on the scale of the library sizes, which leaves
the coefficients readable as log2 fold changes and keeps the fitted
values on the scale of the raw counts. Running without offsets is
possible and is almost never what you want: the library sizes alone
assume that the depth is the only difference between the samples. A
design with one sample per condition has no residual degree of freedom,
so nothing in the data says how much two libraries differ for reasons
unrelated to the treatment. In that case the number is taken from rows
assumed not to respond, through
[`estimateNullDispersion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md),
and the fit switches from the quasi-likelihood F test to a likelihood
ratio test with the dispersion held fixed. It happens here rather than
being asked for, since the alternative is a fit that cannot be produced
at all, but it is announced when it does: the result is conditional on
that number, and
[`checkNullCalibration`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md)
is what turns the assumption into something checkable. `edgeR`'s own
guidance, when no null rows exist either, is to pick a BCV by
experience, around 0.4 for human samples, 0.1 for genetically identical
model organisms and 0.01 for technical replicates, and to read the
output as descriptive. The engines answer slightly different questions.
`"edgeR"` is the default and the safest with few replicates, since the
quasi-likelihood F test carries the uncertainty of the dispersion
estimate into the p-value. `"voom"` is faster on large objects and more
flexible on the design, and it is the only one of the four that can
absorb a repeated-measures structure through `block` without spending a
coefficient on it. `"dream"` extends the same machinery to explicit
random effects, which is what a design with several samples per donor
asks for. `"deseq2"` is included for comparison and for the shrunken
fold changes; on a few thousand regions it agrees with edgeR almost
everywhere. The comparison universe of the set level tests is built here
rather than there, because it depends on the rows and not on the
contrast: one fit, one universe, however many contrasts are run on it
afterwards. An object holding a single region set has nothing to compare
that set against, and in that case the universe comes out empty with a
message; the per-region analysis is unaffected. Low-count rows are not
removed here. Filter them with
[`filterRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md)
first, or the dispersion trend is fitted on rows that carry no
information.

## See also

[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md),
[`makeSetUniverse`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md),
[`selectSamples`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/selectSamples.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes
counts <- filterRegions(counts, verbose = FALSE)

fit <- fitRegions(counts, design = ~ condition, engine = "edgeR", verbose = FALSE)
fit
#> An object of class 'RegionSetDE.fit'
#>   engine          : edgeR 
#>   rows            : 1895 (region level) 
#>   samples         : 4 (lv-H3K4me3-BN-female-bio1-tech1, lv-H3K4me3-BN-male-bio2-tech1, lv-H3K4me3-SHR-male-bio2-tech1, lv-H3K4me3-SHR-male-bio3-tech1) 
#>   coefficients    : (Intercept), conditionSHR 
#>   set universe    : otherSets (matched on width and abundance) 
#>   common disp.    : 0.139 (BCV 0.373)  

# \donttest{
# limma-voom on the same design
voomFit <- fitRegions(counts, design = ~ condition, engine = "voom", verbose = FALSE)
# }

if (FALSE) { # \dontrun{
# Random effects need the dream engine and the formula given as such
mixedFit <- fitRegions(counts, design = ~ condition + (1|donor), engine = "dream")
} # }
```
