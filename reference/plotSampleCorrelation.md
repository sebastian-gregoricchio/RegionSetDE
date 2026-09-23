# plotSampleCorrelation

Draws the correlation between the samples as a heatmap, clustered and
annotated with any column of the sample table, such as the condition,
the treatment or the replicate. The correlation can be computed on the
normalised or on the raw signal, on every region or on the most variable
ones.

## Usage

``` r
plotSampleCorrelation(
  object,
  set = NULL,
  contrast = NULL,
  method = "spearman",
  groupBy = NULL,
  annotationColumns = NULL,
  annotationColours = NULL,
  useOffsets = TRUE,
  compareOffsets = FALSE,
  facetBySet = FALSE,
  cluster = TRUE,
  clusteringMethod = "complete",
  showDendrogram = TRUE,
  topRegions = NULL,
  excludeDiagonal = FALSE,
  palette = NULL,
  limits = NULL,
  showValues = TRUE,
  valuesColour = NULL,
  digits = 2,
  title = NULL,
  fontSize = 9,
  verbose = TRUE
)
```

## Arguments

- object:

  `RegionSetDE.counts`, `RegionSetDE.fit` or any result object of the
  package, or the list returned by
  [`computeSampleCorrelation`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeSampleCorrelation.md).

- set:

  Character vector with the names of the region sets used. Default:
  `NULL`, all of them.

- contrast:

  String with the name of a contrast, or its position, when `object`
  holds several of them. Default: `NULL`.

- method:

  String with the correlation, one of `"spearman"`, `"pearson"` and
  `"kendall"`. Default: `"spearman"`.

- groupBy:

  String with the name of a `colData` column defining the groups whose
  within and between correlations are summarised in the title. Default:
  `NULL`.

- annotationColumns:

  Character vector with the `colData` columns drawn as annotation bars
  above and beside the heatmap, for instance
  `c("condition", "replicate")`. Default: `NULL`, the `groupBy` column
  when there is one.

- annotationColours:

  Named list with the colours of the annotation columns, as
  `ComplexHeatmap` takes them: a named vector for a categorical column,
  a
  [`circlize::colorRamp2`](https://rdrr.io/pkg/circlize/man/colorRamp2.html)
  function for a numeric one. The columns left out take the palette of
  the package. Default: `NULL`.

- useOffsets:

  Logical value to indicate whether the normalisation stored in the
  object must be applied, `FALSE` scaling the samples by their library
  sizes alone. Default: `TRUE`.

- compareOffsets:

  Logical value to indicate whether the same heatmap must be drawn twice
  side by side, once with the normalisation and once on the library
  sizes alone. Default: `FALSE`.

- facetBySet:

  Logical value to indicate whether each region set must get its own
  heatmap, side by side. Default: `FALSE`.

- cluster:

  Logical value to indicate whether the samples must be ordered by
  hierarchical clustering rather than kept in the order of the object.
  Default: `TRUE`.

- clusteringMethod:

  String with the agglomeration passed to
  [`stats::hclust`](https://rdrr.io/r/stats/hclust.html). Default:
  `"complete"`.

- showDendrogram:

  Logical value to indicate whether the dendrogram of the clustering
  must be drawn. Default: `TRUE`.

- topRegions:

  Numeric value with the number of most variable regions the correlation
  is computed on. Default: `NULL`, all of them.

- excludeDiagonal:

  Logical value to indicate whether the diagonal must be left empty.
  Default: `FALSE`.

- palette:

  Character vector with the colours of the scale. Default: `NULL`,
  `viridisLite::mako(100, direction = -1)`.

- limits:

  Numeric vector of length two with the range of the colour scale,
  either value possibly `NA` to take that end from the data. Values
  outside are drawn at the nearest end rather than dropped, and how many
  were is reported. Default: `NULL`, the range of the values off the
  diagonal.

- showValues:

  Logical value to indicate whether the correlations must be written in
  the cells. Default: `TRUE`.

- valuesColour:

  String with the colour of the written values. Default: `NULL`, black
  or white on each cell depending on how dark it is.

- digits:

  Numeric value with the number of decimals written. Default: `2`.

- title:

  String with the title, written above the first heatmap. Default:
  `NULL`.

- fontSize:

  Numeric value with the font size of the names, the values being
  written slightly smaller. Default: `9`.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `Heatmap` built by `ComplexHeatmap`, or a `HeatmapList` when
`compareOffsets` or `facetBySet` draw several of them, drawn when
printed.
[`ComplexHeatmap::draw`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html)
opens the layout, for instance `draw(x, heatmap_legend_side = "bottom")`
or `draw(x, ht_gap = grid::unit(5, "mm"))` to push several panels apart.
The matrices themselves come out of
[`computeSampleCorrelation`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeSampleCorrelation.md).

## Details

The scale runs over the values off the diagonal rather than from zero to
one, because every sample correlates with itself perfectly and every
pair of libraries from the same assay correlates highly. A scale
anchored at zero turns the whole matrix one shade and hides the
differences that matter. `limits` takes that decision back, and either
end can be left as `NA` to be read from the data: `c(NA, 1)` fixes the
top at one and lets the bottom follow the values. Values outside
`limits` are drawn at the nearest end, which hides how far past it they
went, so the number of cells concerned is reported.

The palette is sequential, since a correlation has a low end and a high
end and nothing meaningful in the middle. A diverging scale with white
at the centre reads that midpoint as an absence, which on a matrix where
everything sits between 0.9 and 1 is exactly wrong.

With `groupBy`, the mean correlation within a group and between groups
is written in the title. Within above between is what a usable
experiment looks like; the two being equal says the condition effect is
small next to the replicate noise, and that is the answer about whether
to block, regardless of what an ordination suggests.

The clustering, on one minus the correlation, comes from the first
heatmap and is reused by the others, so a comparison across
`compareOffsets` or `facetBySet` shows the values changing rather than
the samples moving. The dendrogram on the side is drawn once for the
same reason, and so are the sample names, the columns of every panel
following the order of the rows. Numeric annotation columns get a grey
gradient, and turning one into a factor colours it by level instead,
which suits a replicate number.

`compareOffsets` answers less here than it does on an ordination. A
correlation does not see a single factor per sample, so the two panels
come out identical unless the normalisation holds one offset per region,
as `method = "loess"` and offsets supplied from outside do. Whether the
scaling factors are driving a grouping is a question for
[`plotRegionPCA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegionPCA.md).

## See also

[`computeSampleCorrelation`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeSampleCorrelation.md),
[`plotRegionPCA`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegionPCA.md),
[`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
counts <- loadExampleData("counts", verbose = FALSE)
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#> calcNormFactors has been renamed to normLibSizes

plotSampleCorrelation(counts, groupBy = "condition", annotationColumns = c("condition", "sex"))


# The same samples before and after the normalisation, on the CpG island promoters only
plotSampleCorrelation(counts, set = "promoterCpG", method = "pearson", compareOffsets = TRUE)

```
