# Changelog

## RegionSetDE 0.99.0

First version.

### Regions and counting

- [`loadRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)
  reads region sets from BED, narrowPeak, broadPeak or GRanges and keeps
  them as named sets of arbitrary width inside a `RegionSetDE` object.
- [`applyBlacklist()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md)
  and
  [`applyWhitelist()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyWhitelist.md)
  restrict the regions, and record what they removed in the
  `filtering.log` slot.
- [`countReads()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
  and
  [`countBigwig()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)
  count over the regions, either one row per region or one row per tile,
  and return a `RegionSetDE.counts` object.
- [`countBackground()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)
  counts genome-wide bins alongside the regions, for the normalisation
  and for the null estimates.

### Normalisation and filtering

- [`normalizeCounts()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
  estimates scaling factors from the background bins, from the regions
  themselves, or takes them from outside, for instance from a spike-in
  or a greenlist.
- [`plotSetMA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetMA.md)
  and
  [`plotNormComparison()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotNormComparison.md)
  show what a normalisation did before anything is fitted on it.
- [`filterRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md)
  removes the rows that carry too little signal to say anything, on
  average abundance alone so that the choice is independent of the
  contrast tested afterwards. Width-adjusted by default, since a
  threshold in reads otherwise keeps every broad region and drops every
  narrow one.

### Sample selection

- [`selectSamples()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/selectSamples.md)
  filters the samples with
  [`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html)
  syntax on the `colData`, and
  [`splitSamples()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitSamples.md)
  splits an object into one piece per mark or assay. Both drop the
  stored normalisation by default, since factors estimated across marks
  describe a library composition that no longer exists once the object
  is subset.

### Fitting and testing

- [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
  fits one model per region with `edgeR`, `limma-voom`,
  [`variancePartition::dream`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
  or `DESeq2`, reading the normalisation out of the object as offsets
  rather than recomputing it.
- [`testRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
  tests a contrast, or a named list of them, and combines tiled regions
  back to one row per region through
  [`csaw::combineTests`](https://rdrr.io/pkg/csaw/man/combineTests.html).
- Contrasts can be given as a coefficient name, an expression over the
  design columns, a numeric vector, or as
  `c("column", "groupA", "groupB")`, which works whatever the reference
  level is.

### Region sets

- [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  asks whether a set responds as a whole, competitively through
  [`limma::cameraPR`](https://rdrr.io/pkg/limma/man/camera.html) and
  self-contained through
  [`limma::fry`](https://rdrr.io/pkg/limma/man/roast.html), with the
  variance inflated for the correlation between regions so that a large
  set does not come out certain by virtue of being large.
- [`testSetContrast()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md)
  asks whether a contrast affects one set differently from another.
- [`makeSetUniverse()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md)
  builds the comparison universe, matched on width and abundance;
  [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
  builds one automatically and keeps it in the fit.

### Designs without replicates

- [`estimateNullDispersion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md)
  reads the between-sample variation off rows assumed not to respond, so
  that a design with one sample per condition has a dispersion to be
  tested against.
  [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
  calls it on its own when the design leaves no residual.
- [`checkNullCalibration()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md)
  runs the same contrast on rows that should not respond and reports how
  many come out significant anyway, broken down by abundance, with a
  suggested dispersion when the current one is off.

### Plots

- Per region:
  [`plotVolcano()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotVolcano.md),
  [`plotResultsMA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotResultsMA.md),
  [`plotRegion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md),
  [`plotTopHeatmap()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotTopHeatmap.md).
- Per set:
  [`plotSetEffect()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetEffect.md),
  [`plotSetDistribution()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetDistribution.md),
  [`plotSetSignal()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md),
  [`plotUniverseMatching()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotUniverseMatching.md).
- Samples:
  [`plotRegionPCA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegionPCA.md)
  and
  [`plotSampleCorrelation()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSampleCorrelation.md),
  both able to draw the same figure with and without the normalisation
  so that a grouping caused by the scaling factors can be told apart
  from one in the data.

### Export

- [`asDGEList()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDGEList.md),
  `as(x, "DGEList")` and
  [`asDESeqDataSet()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDESeqDataSet.md)
  hand the counts to `edgeR` or `DESeq2` with the offsets attached the
  right way round.
- [`exportResults()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/exportResults.md)
  writes the table, a BED coloured by direction, and every parameter the
  analysis was run with.
