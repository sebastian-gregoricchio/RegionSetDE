# testRegions

Tests a contrast on a `RegionSetDE.fit` object and returns one row per
region. When the counts were tiled, every tile is tested on its own and
the p-values are then combined back to the region, so that the region
stays the unit of inference even though the signal was measured at a
finer scale.

## Usage

``` r
testRegions(
  fit,
  contrast,
  combine = TRUE,
  combineMethod = "simes",
  lfcThreshold = 0,
  FDR = 0.05,
  log2FC = 0,
  adjustMethod = "BH",
  regionSets = NULL,
  extraColumns = TRUE,
  carryCounts = TRUE,
  verbose = TRUE
)
```

## Arguments

- fit:

  `RegionSetDE.fit` object.

- contrast:

  Contrast to test, given in one of four ways. A character vector of
  length three, `c("column", "groupA", "groupB")`, naming a column of
  the `colData` and two of its levels, which is the form to reach for
  when the design uses a reference level. A string with the name of a
  design column, e.g. `"conditionCOMBO"`. A string written as an
  expression over the design columns, e.g.
  `"conditionCOMBO - conditionEPZ"`. Or a numeric vector with one
  coefficient per column of the design. A named list of any of these
  runs every contrast on the same fit and returns a
  `RegionSetDE.resultsList`.

- combine:

  Logical value to indicate whether the tile level p-values must be
  combined into one value per region. Ignored when the counts were not
  tiled. Default: `TRUE`.

- combineMethod:

  String with the method used to combine the tiles, among those accepted
  by
  [`csaw::combineTests`](https://rdrr.io/pkg/csaw/man/combineTests.html):
  `"simes"`, `"holm-min"`, `"wilcoxon"` and `"stouffer"`. Default:
  `"simes"`.

- lfcThreshold:

  Numeric value with the log2 fold change against which the null
  hypothesis is tested. A value above zero moves the threshold inside
  the test, through
  [`edgeR::glmTreat`](https://rdrr.io/pkg/edgeR/man/glmTreat.html),
  [`limma::treat`](https://rdrr.io/pkg/limma/man/ebayes.html) or the
  `lfcThreshold` of
  [`DESeq2::results`](https://rdrr.io/pkg/DESeq2/man/results.html),
  which is stricter and better calibrated than filtering the output
  afterwards. Default: `0`.

- FDR:

  Numeric value with the adjusted p-value cut-off used to fill the
  `diff.status` column. Default: `0.05`.

- log2FC:

  Numeric value with the absolute log2 fold change cut-off used to fill
  the `diff.status` column. Default: `0`.

- adjustMethod:

  String with the multiple testing correction, passed to
  [`stats::p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default:
  `"BH"`.

- regionSets:

  Character vector with the names of the region sets to keep in the
  output. Default: `NULL`, all of them.

- extraColumns:

  Annotation carried by the regions that must be appended to the result,
  at the end of the table. Either `TRUE` for every column of the
  `rowData` beyond the ones the package writes itself, `FALSE` for none,
  or a character vector naming the ones wanted. Default: `TRUE`.

- carryCounts:

  Logical value to indicate whether the counts must travel inside the
  result, so that
  [`plotRegion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md)
  and
  [`plotTopHeatmap`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotTopHeatmap.md)
  can draw the values without being handed the counts object again.
  Several contrasts run on one fit share the same copy in memory.
  Default: `TRUE`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `RegionSetDE.results` object.

## Details

The multiple testing correction is applied over all the rows of the
object, across the region sets, and `regionSets` subsets the output
afterwards. Correcting inside each set separately would make the FDR of
a set depend on how many other sets were loaded, which is not a property
anyone wants in a result.

Two things follow from the combination step. The p-value of a tiled
region is a Simes combination, so it answers "does any part of this
region change" rather than "does the whole region change", and a long
domain that moves over one tile out of forty will come out with a small
p-value and a small overall fold change. The `log2FC` reported for a
combined region is the fold change of the most significant tile, not an
average, which is the quantity that matches the p-value. The tile level
table stays available in the `tiles` slot, and
[`plotRegion`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md)
draws it.

A design written as `~ condition` spends one coefficient per level
except the first, so a level can be a coefficient in the design or the
reference the others are measured against, depending on how the factor
was ordered. Naming a coefficient that turns out to be the reference is
the usual source of confusion, and it is what
`c("column", "groupA", "groupB")` avoids: that form averages the design
rows of each group and takes the difference, which gives the same
contrast whatever the reference is and whether the design was written as
`~ condition` or `~ 0 + condition`. With other covariates in the design
the averaging picks up their imbalance between the two groups, so it
describes what it says only when the design is reasonably balanced.

Whatever the regions were loaded with travels through to the result. A
gene name, a peak score or any other column attached to the `rowData`
comes out at the end of the table, which is what makes
[`topRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md)
readable and lets `plotVolcano(labelColumn = )` label the points with
something other than an identifier. On a tiled object the value is read
off the tile the combination reported, the same one the fold change
comes from, so a row describes one place rather than an average over
several.

The `diff.status` column is a labelling convenience, not a claim. It is
filled from `FDR` and `log2FC` and used by the plotting functions; the
thresholds are stored in the object so that a figure can state them.

## See also

[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md),
[`topRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md),
[`plotVolcano`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotVolcano.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- fitRegions(counts, design = ~ replicate + condition, engine = "edgeR")

res <- testRegions(fit, contrast = "conditionCOMBO")

# Two levels of a column, whichever of them the design took as reference
res <- testRegions(fit, contrast = c("condition", "COMBO", "DMSO"))

# Difference between two coefficients of the design
res <- testRegions(fit, contrast = "conditionCOMBO - conditionEPZ")

# Several contrasts on the same fit
resList <- testRegions(fit, contrast = list(combo = c("condition", "COMBO", "DMSO"),
                                            epz = c("condition", "EPZ", "DMSO")))
resList
topRegions(resList, contrast = "combo")

# Threshold inside the test rather than on the output
resStrict <- testRegions(fit, contrast = "conditionCOMBO", lfcThreshold = 1)
} # }
```
