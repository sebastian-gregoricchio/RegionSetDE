# plotSetEffect

Draws the effect size of every region set with its confidence interval,
which is the figure the set level conclusion should rest on. The
interval carries the inflation for the correlation between the regions,
so a set of thirty thousand promoters does not come out looking thirty
thousand times more certain than a set of two hundred enhancers.

## Usage

``` r
plotSetEffect(
  setResults,
  value = "delta.log2FC",
  contrast = NULL,
  colourBy = "FDR",
  FDR = NULL,
  showN = TRUE,
  orderBy = "effect",
  colours = NULL,
  pointSize = 2.5,
  title = NULL,
  subtitle = NULL,
  legendPosition = "right",
  baseSize = 12
)
```

## Arguments

- setResults:

  `RegionSetDE.setResults` or `RegionSetDE.setResultsList` object.

- value:

  String with the quantity on the axis, either `"delta.log2FC"` (the
  difference with the background) or `"mean.log2FC"` (the shift away
  from zero). Default: `"delta.log2FC"`.

- contrast:

  String with the name of the contrast to draw, or its position, when
  `setResults` holds several of them. Default: `NULL`.

- colourBy:

  String with the variable driving the colour, either `"FDR"`,
  `"direction"` or `"none"`. Default: `"FDR"`.

- FDR:

  Numeric value with the adjusted p-value cut-off used by the colouring.
  Default: `NULL`, the threshold stored in the object.

- showN:

  Logical value to indicate whether the number of regions must be
  written next to each set. Default: `TRUE`.

- orderBy:

  String with the ordering of the sets, either `"effect"` or `"name"`.
  Default: `"effect"`.

- colours:

  Character vector of length two with the colours of the significant and
  non-significant sets. Default: `NULL`.

- pointSize:

  Numeric value with the size of the points. Default: `2.5`.

- title:

  String with the title of the plot, rendered as markdown. Default:
  `NULL`, the contrast.

- subtitle:

  String with the subtitle of the plot, rendered as markdown. Default:
  `NULL`.

- legendPosition:

  String with the position of the legend. Default: `"right"`.

- baseSize:

  Numeric value with the base font size. Default: `12`.

## Value

A `ggplot` object.

## Details

`"mean.log2FC"` is a self-contained quantity and moves with any global
shift of the mark, including one left behind by an imperfect
normalisation. `"delta.log2FC"` is the difference between the set and
what it was compared against, and a scaling error common to both cancels
out of it. When the two tell different stories, the second is the one
that survives a reviewer.

Both levels of the colour scale are kept in the legend even when only
one of them occurs, so that a figure in which nothing reaches the
cut-off still says what the cut-off was.

## See also

[`testRegionSets`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md),
[`plotSetDistribution`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetDistribution.md),
[`plotSetSignal`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
setResults <- testRegionSets(fit, contrast = c("condition", "SHR", "BN"),
                             verbose = FALSE)

# One point per set, with its confidence interval
plotSetEffect(setResults)


plotSetEffect(setResults, value = "mean.log2FC", orderBy = "name")

```
