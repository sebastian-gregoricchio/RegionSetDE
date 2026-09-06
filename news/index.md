# Changelog

## RegionSetDE 0.99.1

Statistical and interpretive corrections, from a review of the set-level
inference. Three of these change what the output means rather than how
it is computed, so a result produced with 0.99.0 should be re-read
against them.

### Set-level effect sizes

- [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  and
  [`testSetContrast()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md)
  gain `effectMethod`. With the default `"sample"` the confidence
  interval is built from one set score per library, the mean signal over
  the set minus the mean signal over its comparison, run through the
  design of the experiment. The replication behind the interval is now
  the biological samples.
- The interval of 0.99.0 is still reported, under
  `heterogeneity.CI.lower` and `heterogeneity.CI.upper`, and named for
  what it measures: how far the effect varies between the loci of a set,
  conditional on the libraries at hand. It was documented as a
  confidence interval on the set effect, which it is not, since its
  sampling units are genomic loci.
- `sample.delta.log2FC`, its standard error, degrees of freedom, p-value
  and adjusted p-value are reported alongside.

### Interpretation of camera and fry

- The reading of a significant competitive test with a non-significant
  self-contained one as redistribution has been removed from the manual,
  the vignette and the README. Failing to reject a self-contained null
  is not evidence that the absolute change is zero, and unequal power
  between the two tests produces the same pattern. The four outcomes are
  now described as evidence rather than as mechanism.
- A redistribution claim belongs to
  [`testSetContrast()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md),
  which compares two sets directly. The documentation points there.
- Also stated: a centring normalisation removes a genuinely global shift
  before `fry` ever sees the data, so the self-contained test is not
  where to look for one.

### What normalisation can and cannot decide

- The claim that
  [`plotNormComparison()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotNormComparison.md)
  and
  [`plotSetMA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetMA.md)
  say whether a global shift is technical or biological has been
  dropped. The two are not separable from endogenous data alone. Those
  plots answer whether the conclusions are sensitive to the
  normalisation assumption, which is what they are now documented as
  answering.

### No-replicate analysis

- The “no replicates needed” heading is gone.
  [`estimateNullDispersion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md)
  measures how two libraries differ over rows assumed not to respond,
  which is not the biological variability that was never sampled, and
  the documentation says so.
- [`normalizeCounts()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
  gains `backgroundHoldout`, which keeps a fraction of the background
  bins out of the estimation of the scaling factors.
  [`estimateNullDispersion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md)
  picks up that split rather than drawing its own, so the calibration
  rows sit outside the whole preprocessing chain.
- The returned list gains `holdout.type`, distinguishing rows held out
  of the dispersion alone from rows held out of both steps.

### Counts and coverage

- [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
  gains the `"limma"` engine, limma-trend on the log2 signal, for values
  that are not counts.
- [`countBigwig()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)
  gains `countLike` and no longer rounds by default. Rounding coverage
  to integers does not make it a fragment count, and the negative
  binomial and voom engines now refuse an object built from bigWig files
  unless `countLike` was declared or `assumeCountLike` overrides it in
  [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md).
- [`loadCounts()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md)
  gains the same `countLike`, for external matrices holding coverage
  rather than counts.

### Tiles

- [`countReads()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md),
  [`countBigwig()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)
  and
  [`loadCounts()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md)
  mark a tiled object as tiled. They previously stored
  `counting.level = "region"` whatever `tileWidth` was, so
  [`testRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
  never recombined the tiles of a region and the tile-level output path
  could not be reached.
- [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  gains `tileHandling`, collapsing the tiles of a region into one row
  before the set is assembled. Without it a 40 kb region counted at 1 kb
  weighs forty times a 1 kb one, which makes the set effect an average
  over base pairs rather than over regions.

### Overlapping sets

- [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  gains `overlapPolicy`, and
  [`testSetContrast()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md)
  now detects shared regions through
  [`IRanges::findOverlaps`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  rather than through region identifiers. Two sets can cover the same
  chromatin without sharing an identifier, and the shared reads pull the
  difference between them towards zero.
- The number of overlapping comparison rows is reported in
  `n.comparison.overlapping`.

### The comparison universe

- [`makeSetUniverse()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md),
  [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
  and
  [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  gain `universeSets`, naming the sets the comparison rows are drawn
  from.
- `RegionSetDE.universe` gains the `comparison.sets` slot and prints it.
  A competitive p-value is relative to the sets that happen to be
  loaded, and that is now recorded with the result instead of having to
  be reconstructed.

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
