# plotRegion

Draws the signal of a single region across the samples. On a tiled
object the values are drawn along the coordinates, one line per sample,
which shows whether the whole region moved or only part of it. On a
non-tiled object the region carries one value per sample and those are
drawn as points, with brackets between the groups when `pairwiseTest` is
set.

## Usage

``` r
plotRegion(
  object,
  region,
  counts = NULL,
  contrast = NULL,
  groupBy = NULL,
  assay = NULL,
  log2Scale = TRUE,
  summarise = FALSE,
  pairwiseTest = NULL,
  pairBy = NULL,
  comparisons = NULL,
  pAdjustMethod = "none",
  pLabel = "p.value",
  pDecimals = 2,
  colours = NULL,
  pointSize = 3,
  rotateX = TRUE,
  title = NULL,
  subtitle = NULL,
  legendPosition = "right",
  baseSize = 12
)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit`, `RegionSetDE.results` or
  `RegionSetDE.resultsList` object. A result carries both the values and
  the statistics, so nothing else has to be passed.

- region:

  String identifying the region, written as `"set|id"` or as the region
  identifier alone when it is unique across the sets. A `GRanges` of
  length one is accepted as well, in which case the overlapping rows are
  drawn.

- counts:

  `RegionSetDE.counts` object holding the values, when `object` carries
  none. Default: `NULL`.

- contrast:

  String with the name of the contrast to annotate with, or its
  position, when `object` holds several of them. With
  `pairwiseTest = "model"` it keeps the bracket of that contrast only.
  Default: `NULL`, which with `pairwiseTest = "model"` draws one bracket
  per contrast of the object.

- groupBy:

  String with the name of a `colData` column driving the colour, e.g.
  `"condition"`. Its levels are also the groups compared by
  `pairwiseTest`. Default: `NULL`, one colour per sample.

- assay:

  String with the name of the assay to draw. Default: `NULL`, the
  normalised assay when present, the raw counts otherwise.

- log2Scale:

  Logical value to indicate whether the values must be drawn on a log2
  scale. Default: `TRUE`.

- summarise:

  Logical value to indicate whether the replicates of a group must be
  summarised rather than drawn one by one: a mean line with a ribbon
  along a tiled region, a mean with its spread next to the individual
  points on a region counted as a single row. Requires `groupBy`.
  Default: `FALSE`.

- pairwiseTest:

  String with what the brackets between the groups of `groupBy` report,
  one among `"t.test"`, `"wilcox.test"` and `"model"`. The first two
  test the plotted values, with the Welch t-test or the Wilcoxon
  rank-sum test, or with their paired versions when `pairBy` is given.
  `"model"` writes instead the log2 fold change and the FDR of the
  fitted contrasts, and needs a results object. Only for a region
  counted as a single row. Default: `NULL`, no brackets.

- pairBy:

  String with the name of a `colData` column matching the samples of two
  groups, e.g. `"replicate"`, for a paired test. A value can appear once
  per group, and a sample whose value has no partner in the other group
  is left out of that comparison. Default: `NULL`, unpaired tests.

- comparisons:

  List of character vectors of length two, naming the pairs of groups
  joined by a bracket. Default: `NULL`, every pair of groups, or with
  `pairwiseTest = "model"` every pair tested by a contrast.

- pAdjustMethod:

  String with the correction applied to the p-values of the brackets,
  one of
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html). It
  covers the comparisons drawn in the plot, not the regions tested, and
  it is not used with `pairwiseTest = "model"`. Default: `"none"`.

- pLabel:

  String indicating how the values are written on the brackets:
  `"p.value"` for the number, `"stars"` for the significance symbols
  (`ns` above 0.05, then `*`, `**`, `***` and `****` up to 0.05, 0.01,
  0.001 and 0.0001). Default: `"p.value"`.

- pDecimals:

  Numeric value with the number of decimals of the values written on the
  brackets. Values below 0.1 are written in scientific notation, with
  the exponent as a superscript. Default: `2`.

- colours:

  Named character vector with the colours. Default: `NULL`.

- pointSize:

  Numeric value with the size of the points, on a non-tiled region.
  Default: `3`.

- rotateX:

  Logical value to indicate whether the labels of the x axis must be
  angled. Default: `TRUE`.

- title:

  String with the title of the plot, rendered as markdown. Default:
  `NULL`, the region identifier.

- subtitle:

  String with the subtitle of the plot, rendered as markdown. Default:
  `NULL`, the statistics of the region when they are available.

- legendPosition:

  String with the position of the legend. Default: `"right"`.

- baseSize:

  Numeric value with the base font size. Default: `12`.

## Value

A `ggplot` object.

## Details

The values come from the object and nothing is re-read from the BAM or
bigWig files, so the resolution of the plot is the resolution of the
counting. A region counted as a single row gives a single point per
sample, which is the honest picture of what the model saw.

The brackets drawn by `"t.test"` and `"wilcox.test"` answer a smaller
question than the model does. The test runs on the values as they are
drawn, on the log2 scale when `log2Scale = TRUE`, and it knows nothing
of the offsets, of the moderated dispersion, of the covariates in the
design or of the thousands of regions tested next to this one. Its
p-value and the FDR in the subtitle will often disagree. The region is
also usually picked from
[`topRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md)
on the same samples, so the test describes the region and does not
confirm the result. Where it helps is on counts and fit objects, which
carry no statistics of their own, and between groups that no contrast
compared. On raw counts the function warns, since a difference in
sequencing depth between the groups is read as a difference in signal.

Few replicates limit these tests more than they limit the model. With
three samples per group the exact Wilcoxon test cannot go below p = 0.1,
and with three pairs the signed-rank test cannot go below 0.25; a
message says so when the groups are that small. With `pairBy` the
samples are matched on the value of that column, whatever their order in
the object. A test paired by replicate asks what
`design = ~ replicate + condition` asks in
[`fitRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md),
and the model is the place for that pairing when the numbers are meant
to be reported.

`"model"` writes on each bracket the log2 fold change and the FDR that
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
gave to the region, one bracket per contrast comparing two levels of
`groupBy`. A contrast qualifies when it was written as
`c("column", "groupA", "groupB")`, or when it reduces to that
difference, which the result records in its `contrast.groups` slot. The
caption says which test, or which fit, the brackets come from.

## See also

[`topRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md),
[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`plotTopHeatmap`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotTopHeatmap.md),
[`plotSetSignal`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

topRegion <- topRegions(results, n = 1, FDR = 1)$region.id

plotRegion(results, region = topRegion, groupBy = "condition")


# Summarised to one point per group rather than one per sample
plotRegion(results, region = topRegion, groupBy = "condition", summarise = TRUE)


# Fold change and FDR of the fitted contrast on the bracket
plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "model")


# Welch t-test on the plotted values
plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "t.test")


# Paired t-test; the example has no matched replicates, so the pairs are made up to show the call
pairedCounts <- resultCounts(results)
pairedCounts$pair <- c("p1", "p2", "p1", "p2")

plotRegion(pairedCounts, region = topRegion, groupBy = "condition",
           pairwiseTest = "t.test", pairBy = "pair")

```
