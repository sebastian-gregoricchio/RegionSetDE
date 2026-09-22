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
  and return a `RegionSetDE.counts` object. A tiled object is marked as
  tiled, so
  [`testRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
  recombines the tiles of a region instead of treating each of them as a
  region of its own.
- [`countBigwig()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)
  takes `countLike` and does not round. Rounding coverage to integers
  does not make it a fragment count, and the negative binomial and voom
  engines refuse an object built from bigWig files unless `countLike`
  was declared or `assumeCountLike` overrides it in
  [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md).
  [`loadCounts()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md)
  takes the same argument, for external matrices holding coverage rather
  than counts.
- [`countBackground()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)
  counts genome-wide bins alongside the regions, for the normalisation
  and for the null estimates.
- [`countReads()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
  and
  [`countBackground()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)
  rebuild each paired-end fragment from the first mate of the pair and
  its template length, and read the BAM files in pieces shared among the
  threads. `excludeChromosomes` keeps chromosomes such as chrM out of
  the library sizes and of the background bins, while
  `fullLibrarySize = FALSE` reads only the chromosomes carrying regions,
  which is faster but gives library sizes that are not meant for
  normalisation.

### Normalisation and filtering

- [`normalizeCounts()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeCounts.md)
  estimates scaling factors from the background bins, from the regions
  themselves, or takes them from outside, for instance from a spike-in
  or a greenlist. `backgroundHoldout` keeps a fraction of the bins out
  of that estimation, so that the rows used to check the calibration sit
  outside the whole preprocessing chain.
- [`plotSetMA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetMA.md)
  and
  [`plotNormComparison()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotNormComparison.md)
  show what a normalisation did before anything is fitted on it. What
  they answer is whether the conclusions are sensitive to the
  normalisation assumption. Whether a global shift is technical or
  biological is not separable from endogenous data alone, and neither
  plot decides it.
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
  fits one model per region with `edgeR`, `limma-voom`, `limma-trend`,
  [`variancePartition::dream`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
  or `DESeq2`, reading the normalisation out of the object as offsets
  rather than recomputing it. The `"limma"` engine runs limma-trend on
  the log2 signal, for values that are not counts.
- [`testRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
  tests a contrast, or a named list of them, and combines tiled regions
  back to one row per region through
  [`csaw::combineTests`](https://rdrr.io/pkg/csaw/man/combineTests.html).
- Contrasts can be given as a coefficient name, an expression over the
  design columns, a numeric vector, or as
  `c("column", "groupA", "groupB")`, which works whatever the reference
  level is.
- `RegionSetDE.results` records in `contrast.groups` the column and the
  two levels a contrast compares, as `RegionSetDE.setResults` already
  did. A contrast written as a coefficient name, an expression or a
  vector is assigned to the design variable it belongs to, never to a
  column naming each sample on its own, such as the sample names.

### Region sets

- [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  asks whether a set responds as a whole, competitively through
  [`limma::cameraPR`](https://rdrr.io/pkg/limma/man/camera.html) and
  self-contained through
  [`limma::fry`](https://rdrr.io/pkg/limma/man/roast.html), with the
  variance inflated for the correlation between regions so that a large
  set does not come out certain by virtue of being large. The four
  combinations of the two outcomes are documented as evidence rather
  than as mechanism: failing to reject a self-contained null is not
  evidence that the absolute change is zero, unequal power between the
  two tests produces the same pattern, and a centring normalisation
  removes a genuinely global shift before `fry` ever sees the data.
- `effectMethod` decides what the confidence interval describes. With
  the default `"sample"` it is built from one set score per library, the
  mean signal over the set minus the mean signal over its comparison,
  run through the design of the experiment, so the replication behind it
  is the biological samples. `sample.delta.log2FC`, its standard error,
  degrees of freedom, p-value and adjusted p-value are reported
  alongside. `heterogeneity.CI.lower` and `heterogeneity.CI.upper`
  report the other quantity, how far the effect varies between the loci
  of a set, conditional on the libraries at hand.
- [`testSetContrast()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md)
  asks whether a contrast affects one set differently from another,
  which is where a redistribution claim belongs.
- `tileHandling` collapses the tiles of a region into one row before the
  set is assembled. Without it a 40 kb region counted at 1 kb weighs
  forty times a 1 kb one, which makes the set effect an average over
  base pairs rather than over regions.
- `overlapPolicy` handles sets covering the same chromatin, found
  through
  [`IRanges::findOverlaps`](https://rdrr.io/pkg/IRanges/man/findOverlaps-methods.html)
  rather than through region identifiers, since two sets can cover the
  same chromatin without sharing an identifier and the shared reads pull
  the difference between them towards zero. `n.comparison.overlapping`
  reports how many comparison rows overlap.
- [`makeSetUniverse()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md)
  builds the comparison universe, matched on width and abundance;
  [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
  builds one automatically and keeps it in the fit. `universeSets` names
  the sets the comparison rows are drawn from, and
  `RegionSetDE.universe` records them in `comparison.sets` and prints
  them, since a competitive p-value is relative to the sets that happen
  to be loaded.

### Region set scores

- [`scoreRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md)
  compares region sets to each other in an experiment holding a single
  condition, where there is no contrast for
  [`testSetContrast()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md)
  to run on. It computes one score per set per library, the summarised
  signal over the set divided by a reference measured in the same
  library, and compares the sets library by library.
- The replication is the libraries. The regions of a set collapse to one
  number per library rather than being tested across, so the interval
  rests on the biological samples and the p-value is not the vanishing
  one a per-region test over the same data returns.
- The score is a ratio taken inside one library, so the sequencing depth
  and the scaling factors cancel before it exists. With
  `reference = "background"` the regions are read from the raw `counts`
  assay, since the bins hold raw counts, and the normalisation never
  enters. The reference cancels a second time in the difference between
  two sets, which therefore depends on neither.
- What the comparison does not separate is the factor from the
  composition of the sets. Width, mappability, GC content and
  accessibility produce coverage in a library where nothing is bound,
  reproducibly across replicates. `perBasepair` removes the width term,
  and the documentation states plainly that the rest remain, so that the
  step from more signal to more factor is made as an argument rather
  than read off the output.
- `RegionSetDE.setScores` holds the per-library scores, the comparisons
  and the `colData` of the counts object. `show` prints the sets, the
  libraries, the reference and the comparisons;
  [`resultsTable()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/resultsTable.md)
  returns the comparisons,
  [`scoreTable()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreTable.md)
  the scores behind them, and
  [`regionSetNames()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/regionSetNames.md)
  works on the object.

### Designs without replicates

- [`estimateNullDispersion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md)
  reads the between-sample variation off rows assumed not to respond, so
  that a design with one sample per condition has a dispersion to be
  tested against.
  [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
  calls it on its own when the design leaves no residual. What it
  measures is how two libraries differ over those rows, which is not the
  biological variability that was never sampled, and the documentation
  says so.
- The returned list carries `holdout.type`, distinguishing rows held out
  of the dispersion alone from rows held out of both the dispersion and
  the normalisation.
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
- [`plotRegion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md)
  draws brackets between the groups of `groupBy` through `pairwiseTest`.
  With `"model"` each bracket carries the fold change and the FDR of a
  fitted contrast comparing two of those groups. `"t.test"` and
  `"wilcox.test"` test the plotted values instead, paired through a
  `colData` column with `pairBy`, and `pAdjustMethod` corrects across
  the brackets of the plot. The caption names the test, the function
  warns when it runs on raw counts, and a message says when the groups
  are too small for the exact Wilcoxon test to reach 0.05.
- Per set:
  [`plotSetEffect()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetEffect.md),
  [`plotSetDistribution()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetDistribution.md),
  [`plotSetSignal()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md),
  [`plotUniverseMatching()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotUniverseMatching.md).
  [`plotSetSignal()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md)
  also draws a `RegionSetDE.setScores` object, one point per library per
  set with the points of a library joined across the sets, since the
  comparison behind the bracket is paired.
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
