# plotVolcano

Draws the log2 fold change of a contrast against the significance, one
panel per region set, with the number of changing regions written in the
top corners of each panel.

## Usage

``` r
plotVolcano(
  results,
  set = NULL,
  contrast = NULL,
  colourBy = "diff.status",
  facetBySet = TRUE,
  facetScales = "fixed",
  yValue = "FDR",
  FDR = NULL,
  log2FC = NULL,
  showCounts = TRUE,
  labelTop = 0,
  labelColumn = "region.id",
  colours = NULL,
  pointSize = 0.8,
  maxPoints = 20000,
  title = NULL,
  subtitle = NULL,
  legendPosition = "right",
  baseSize = 12
)
```

## Arguments

- results:

  `RegionSetDE.results` or `RegionSetDE.resultsList` object.

- set:

  Character vector with the names of the region sets to draw. Default:
  `NULL`, all of them.

- contrast:

  String with the name of the contrast to draw, or its position, when
  `results` holds several of them. Default: `NULL`.

- colourBy:

  String with the variable driving the colour, either `"diff.status"` or
  `"region.set"`. Default: `"diff.status"`.

- facetBySet:

  Logical value to indicate whether each region set must get its own
  panel. Default: `TRUE`.

- facetScales:

  String with the scales of the panels, one among `"fixed"`, `"free"`,
  `"free_x"` and `"free_y"`. Default: `"fixed"`.

- yValue:

  String with the quantity on the y axis, either `"FDR"` or `"p.value"`.
  Default: `"FDR"`.

- FDR:

  Numeric value with the adjusted p-value cut-off drawn as a line.
  Default: `NULL`, the threshold stored in the object.

- log2FC:

  Numeric value with the absolute log2 fold change cut-off drawn as a
  line. Default: `NULL`, the threshold stored in the object.

- showCounts:

  Logical value to indicate whether the number of changing regions must
  be written in the top corners of each panel. Default: `TRUE`.

- labelTop:

  Numeric value with the number of top regions to label, per panel.
  Default: `0`.

- labelColumn:

  String with the column holding the labels. Default: `"region.id"`.

- colours:

  Named character vector with the colours. Default: `NULL`, a grey, blue
  and red palette for `"diff.status"`.

- pointSize:

  Numeric value with the size of the points. Default: `0.8`.

- maxPoints:

  Numeric value with the number of non-changing points drawn per panel.
  Default: `20000`.

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

The y axis is the adjusted p-value by default. Showing the raw p-value
while drawing the cut-off line at the FDR puts two different quantities
on the same figure, which is where most misread volcano plots come from.

Only the points labelled `"null"` are thinned by `maxPoints`, and the
thinning happens inside each panel so that a small set keeps all of its
points. Everything passing the thresholds is drawn. The thinning is
deterministic, so the figure does not change between calls.

The counts in the corners come from the full table, before any thinning,
and they are the number of regions labelled `"down"` on the left and
`"up"` on the right.

## See also

[`testRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md),
[`plotResultsMA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotResultsMA.md),
[`topRegions`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
fit <- loadExampleData("fit", verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

plotVolcano(results)


# One set, with the strongest regions labelled
plotVolcano(results, set = "promoterCpG", labelTop = 5)

```
