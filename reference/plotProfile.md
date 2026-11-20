# plotProfile

Draws the signal around the centre of a group of regions, by default the
regions going up and those going down in a contrast, as DiffBind does
with `dba.plotProfile`. The `"heatmap"` style draws one heatmap per
sample or group, one row per region and the mean profile on top; the
`"lines"` style draws the mean profiles alone, with their standard
error.

## Usage

``` r
plotProfile(
  object,
  style = "heatmap",
  colours = NULL,
  groupColours = NULL,
  limits = NULL,
  title = NULL,
  fontSize = 9,
  baseSize = 12,
  legendPosition = "right",
  regions = NULL,
  contrast = NULL,
  set = NULL,
  samples = NULL,
  direction = "both",
  FDR = NULL,
  log2FC = NULL,
  maxRegions = 1000,
  blacklist = NULL,
  whitelist = NULL,
  groupBy = NULL,
  signal = "auto",
  signalFiles = NULL,
  distance = 1500,
  binWidth = 50,
  centre = "summit",
  useOffsets = TRUE,
  BPPARAM = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  Any object accepted by
  [`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md),
  or the list it returns.

- style:

  String with the drawing, either `"heatmap"` or `"lines"`. Default:
  `"heatmap"`.

- colours:

  Character vector with the colours of the heatmap scale, from low to
  high. Default: `NULL`, `viridisLite::mako(100, direction = -1)`.

- groupColours:

  Character vector with the colours of the groups of rows in the
  heatmap, or of the samples or groups in the `"lines"` style, named
  after them or given in their order. Default: `NULL`, the palette of
  the package, with the regions going up in red and those going down in
  blue as in
  [`plotVolcano`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotVolcano.md).

- limits:

  Numeric vector of length two with the range of the colour scale of the
  heatmap. Default: `NULL`, from zero to the 99th percentile of the
  values.

- title:

  String with the title of the `"lines"` style. A `HeatmapList` takes
  its title when it is drawn, as
  `ComplexHeatmap::draw(x, column_title = "title")`. Default: `NULL`,
  the contrast when there is one.

- fontSize:

  Numeric value with the font size of the heatmap. Default: `9`.

- baseSize:

  Numeric value with the base font size of the `"lines"` style. Default:
  `12`.

- legendPosition:

  String with the position of the legend of the `"lines"` style.
  Default: `"right"`.

- regions:

  Regions to profile instead of those picked from `object`: a `GRanges`,
  a named list of `GRanges` or a `GRangesList`, one group of rows per
  element, or a `RegionSetDE` object, one group per set. Default:
  `NULL`.

- contrast:

  String with the name of a contrast, or its position, when `object`
  holds several of them. Default: `NULL`.

- set:

  Character vector with the names of the region sets used. Default:
  `NULL`, all of them.

- samples:

  Samples profiled: a character vector with their names, a numeric
  vector with their positions, or a logical vector with one value per
  sample. The scaling factors of the whole analysis are kept. Default:
  `NULL`, all of them.

- direction:

  String with the differential regions profiled, one among `"both"`,
  which draws the regions going up and those going down as two groups of
  rows, `"up"` and `"down"`. Only for results. Default: `"both"`.

- FDR:

  Numeric value with the adjusted p-value cut-off defining a
  differential region. Default: `NULL`, the one used by the test.

- log2FC:

  Numeric value with the absolute log2 fold change cut-off defining a
  differential region. Default: `NULL`, the one used by the test.

- maxRegions:

  Numeric value with the maximum number of regions per group of rows:
  the most significant ones for results, the strongest ones for the
  region sets of a counts object, and regions spread evenly along the
  list for regions given through `regions`. It is applied after
  `blacklist` and `whitelist`. Default: `1000`.

- blacklist:

  Regions whose signal must stay out of the figure: a row is left out
  when its drawn window, `distance` on each side of the centre, overlaps
  them. A `GRanges`, a path to a BED-like file, a data.frame, a list of
  them, or `TRUE` for the blacklist stored in the object by
  [`applyBlacklist`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)
  or
  [`loadConsensusPeaks`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md).
  Default: `NULL`.

- whitelist:

  Regions the rows are restricted to: a row is kept when its region
  overlaps them. Accepts the same forms as `blacklist`, except `TRUE`.
  Default: `NULL`.

- groupBy:

  String with the name of a `colData` column whose groups are averaged
  into one profile each, for instance `"condition"`. Default: `NULL`,
  one profile per sample.

- signal:

  String with the files the signal is read from: `"bam"`, the BAM files
  of
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
  `"bigwig"`, the bigWig files of the sample sheet or of
  [`countBigwig`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md),
  or `"auto"`, the BAM files when they are known and the bigWig files
  otherwise. Default: `"auto"`.

- signalFiles:

  Character vector with one BAM or bigWig file per sample, in the order
  of the samples, overriding the ones recorded in the object. Default:
  `NULL`.

- distance:

  Numeric value with the number of base pairs drawn on each side of the
  centre. Default: `1500`.

- binWidth:

  Numeric value with the width of the bins the signal is averaged over,
  in base pairs. Default: `50`.

- centre:

  String with the point the regions are aligned on: `"summit"`, the
  `summit` column written by
  [`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
  with `summits`, falling back to the midpoint where there is none, or
  `"midpoint"`. Default: `"summit"`.

- useOffsets:

  Logical value to indicate whether the signal read from BAM files must
  be scaled by the scaling factors of
  [`normalizeCounts`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md),
  `FALSE` scaling it by the library sizes alone. Ignored for bigWig
  files, whose signal is taken as it is. Default: `TRUE`.

- BPPARAM:

  A `BiocParallelParam` object, used as it is, or a number of cores as a
  shorthand for one. One file goes to each worker. Default: `NULL`, a
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html),
  so that nothing runs in parallel unasked.

- verbose:

  Logical value to indicate whether the messages must be printed.
  Default: `TRUE`.

## Value

A `HeatmapList` built by `ComplexHeatmap`, drawn when printed, for the
`"heatmap"` style, or a `ggplot` object for the `"lines"` style.
[`ComplexHeatmap::draw`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html)
opens the layout of the heatmaps, for instance
`draw(x, column_title = "R1881 against DMSO", heatmap_legend_side = "bottom")`.

## Details

Given a results, counts or fit object, `plotProfile()` calls
[`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md)
and draws what it returns. The arguments shared with `computeProfiles`,
from `regions` to `verbose`, are passed on to it, and they are ignored
when `object` is already the list returned by `computeProfiles`. The
reading of the files is the slow part, so profiles drawn in both styles,
or with different colours, are better computed once with
[`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md)
and handed over as they are.

The rows of each group are sorted by their mean signal over all the
samples, strongest on top, so that the same row holds the same region in
every heatmap, and the eye can run across the panels to see a region
gain or lose signal. The colour scale is shared by the panels for the
same reason. The profile on top of each heatmap is the mean of the rows
of every group, drawn on a common axis.

A change found by the test should show here as a difference in height
between the conditions, centred on the summit. A difference spread
evenly across the whole window points to background rather than binding,
and a profile whose peak sits away from the centre says the regions were
not aligned on the signal, which recentring with
`countReads(summits = )` fixes.

## See also

[`computeProfiles`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md),
[`plotTopHeatmap`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotTopHeatmap.md),
[`countReads`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)

## Author

Sebastian Gregoricchio

## Examples

``` r
sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)

counts <- countReads(peakRegions, sampleSheet = sampleSheet, summits = 200, countInput = FALSE, verbose = FALSE)
counts <- normalizeCounts(counts, method = "TMM", verbose = FALSE)
fit <- fitRegions(counts, design = ~condition, verbose = FALSE)
results <- testRegions(fit, contrast = c("condition", "R1881_24h", "DMSO"), verbose = FALSE)

# The regions going up and down, one heatmap per condition
plotProfile(results, groupBy = "condition", distance = 1000, verbose = FALSE)


# Computed once, drawn twice: the BAM files are read a single time
profileData <- computeProfiles(results, groupBy = "condition", distance = 1000, verbose = FALSE)
plotProfile(profileData)

plotProfile(profileData, style = "lines")

```
