# Changelog

## RegionSetDE 0.99.0

First version, submitted to Bioconductor.

### Regions

- [`loadRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadRegions.md)
  reads region sets from BED, narrowPeak, broadPeak, `GRanges` or
  data.frames into named sets of arbitrary width, kept in a
  `RegionSetDE` object, and
  [`splitLoadRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitLoadRegions.md)
  splits one file into sets by one of its columns. A single `GRanges` is
  accepted wherever a region set is expected.
- [`loadSampleSheet()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadSampleSheet.md)
  reads a table with one row per sample: the BAM or bigWig file, the
  peaks called on it, the input it was sequenced against and any
  annotation column. A DiffBind sample sheet is read as it is, relative
  paths are resolved against the folder of the sheet, and missing or
  unindexed files are reported together before anything is counted.
- [`loadConsensusPeaks()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadConsensusPeaks.md)
  builds the regions from the peaks of the samples: a consensus within
  each group through *consensusRegions*, pooled into a total consensus.
  `technicalBy` joins the libraries of a biological sample first, so
  that one sample does not count twice. Every region carries the groups
  and the samples that called a peak on it (`peak.<group>`,
  `peak.groups`, `peak.samples`), which follow it into the counts and
  the results.
  [`consensusData()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusData.md)
  and
  [`consensusGroupList()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/consensusGroupList.md)
  read the consensus back, and
  [`plotPeakUpset()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotPeakUpset.md)
  draws the occupancy.
- [`loadBlacklist()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadBlacklist.md),
  [`loadGreenlist()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadGreenlist.md)
  and
  [`availableRegionLists()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/availableRegionLists.md)
  give the region lists installed with the package: the ENCODE blacklist
  version 2 for hg19, hg38, mm10, dm3, dm6, ce10 and ce11, the T2T-CHM13
  exclusion set of *excluderanges*, and the CUT&RUN and CUT&Tag
  greenlists and high signal regions of de Mello et al. for hg38 and
  mm39. Nothing is downloaded.
- [`makeGreylist()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeGreylist.md)
  builds a greylist from the input libraries following *GreyListChIP*.
  The negative binomial is fitted to every window rather than to
  bootstrap samples, so the result does not depend on a random seed.
- [`applyBlacklist()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyBlacklist.md),
  [`applyWhitelist()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyWhitelist.md)
  and
  [`applyGreylist()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/applyGreylist.md)
  restrict the regions and record what they removed in the
  `filtering.log` slot. A list built for another assembly than the one
  of the regions is refused.

### Counting

- [`countReads()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countReads.md)
  counts BAM files over the regions, once per region or per fixed-width
  tile through `tileWidth`, and returns a `RegionSetDE.counts` object,
  which is a `RangedSummarizedExperiment`. Paired-end fragments are
  rebuilt from the first mate and its template length, and
  `pairedEnd = "auto"` reads the layout from the first records of each
  file.
- `summits` recentres the regions on the summit of the pooled pileup, in
  windows of `2 * summits + 1` bp as `DiffBind::dba.count()` does, and
  `summitSource = "peaks"` takes the summits written in the narrowPeak
  files instead.
- `fragmentLength` takes one value per sample, a column of the sheet or
  `"auto"`, in which case
  [`estimateFragmentLength()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateFragmentLength.md)
  estimates it from the strand cross-correlation of the reads.
- The inputs of the samples are counted into an `input` assay, and
  nothing subtracts them from the counts.
- `blacklist`, `greylist` and `discardRegions` remove regions or ignore
  their reads at the counting, and `excludeChromosomes` keeps
  chromosomes such as chrM out of the library sizes.
- BAM files, bigWig files, regions and lists naming the chromosomes of
  one assembly in different styles, `chr1` in some and `1` in others,
  are matched name by name. CSI indices are found as well as BAI ones.
- [`countBigwig()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBigwig.md)
  counts coverage from bigWig files without rounding it, and the
  negative binomial and voom engines refuse such an object unless
  `countLike` declares the coverage count-like.
  [`loadCounts()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadCounts.md)
  takes a count matrix computed elsewhere.
- [`countBackground()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countBackground.md)
  counts genome-wide bins for the normalisation and the null estimates,
  and
  [`countGreenlist()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countGreenlist.md)
  counts the libraries over the greenlist.
- [`libInfo()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/libInfo.md)
  summarises the libraries in one table: the records of every file, the
  fragments kept, the fragments in the regions, the FRiP, the fragment
  length and the input of each sample.
- [`countTable()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/countTable.md)
  returns the values of any object as a table, raw or normalised, per
  region or per tile, wide or long with the `colData` attached.
- `progressBar` draws a bar while the files are read, by default only in
  interactive sessions.

### Normalisation and filtering

- [`normalizeRegionCounts()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/normalizeRegionCounts.md)
  estimates the scaling factors from background bins, from the regions
  themselves (TMM, RLE or the reads in the regions), from the greenlist,
  or takes them from outside, a spike-in for instance.
  `backgroundHoldout` keeps part of the bins out of the estimate so that
  the calibration can be checked on them later.
- [`plotNormComparison()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotNormComparison.md)
  and
  [`plotSetMA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetMA.md)
  show how far the choice of normalisation moves the samples and the
  sets before anything is fitted on it.
- [`filterRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/filterRegions.md)
  removes the rows carrying too little signal, on their average
  abundance alone, so that the choice is independent of the contrast
  tested afterwards.
- [`selectSamples()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/selectSamples.md)
  filters the samples with
  [`dplyr::filter`](https://dplyr.tidyverse.org/reference/filter.html)
  syntax on the `colData`, and
  [`splitSamples()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/splitSamples.md)
  splits an object by mark or assay.

### Fitting and testing

- [`fitRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/fitRegions.md)
  fits one model per region with *edgeR*, *limma*-voom, *limma*-trend,
  [`variancePartition::dream`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
  or *DESeq2*, reading the normalisation out of the object as offsets.
- [`testRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
  tests a contrast or a named list of them. A contrast can be a
  coefficient name, an expression over the design columns, a numeric
  vector or `c("column", "groupA", "groupB")`, and
  [`pairwiseContrasts()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/pairwiseContrasts.md)
  writes out every pair of levels of a column. Tiled regions are
  combined back to one row per region through
  [`csaw::combineTests()`](https://rdrr.io/pkg/csaw/man/combineTests.html)
  or a representative tile.
- The results carry the distribution of `stat` and its degrees of
  freedom, for a power calculation for instance, and
  [`contrastInfo()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/contrastInfo.md)
  summarises every contrast of an object.
- [`updateThresholds()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/updateThresholds.md)
  labels the regions again with other `FDR` and `log2FC` cut-offs,
  without testing again.
- [`peakOccupancyTable()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/peakOccupancyTable.md)
  and
  [`plotPeakOccupancy()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotPeakOccupancy.md)
  cross the groups that called a region with what the test said about
  it.

### Region sets

- [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  tests whole sets, competitively through
  [`limma::cameraPR`](https://rdrr.io/pkg/limma/man/camera.html) and
  self-contained through
  [`limma::fry`](https://rdrr.io/pkg/limma/man/roast.html). `fry.p` is
  the directional test and `fry.mixed.p` the one for a set changing in
  either direction, and `n.up` and `n.down` count the regions changing
  inside each set. The effect of a set comes with a confidence interval
  built on the libraries and a second interval for how much it varies
  between regions.
- [`makeSetUniverse()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/makeSetUniverse.md)
  builds the rows a set is compared against, matched on width and
  abundance.
- [`testSetContrast()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md)
  asks whether a contrast affects two sets differently, which is where a
  redistribution claim belongs.
- `tileHandling`, `overlapPolicy` and `overlapWithinSet` decide how
  tiles, sets covering the same chromatin and overlapping regions inside
  one set are handled.
- [`scoreRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/scoreRegionSets.md)
  compares the sets to each other in an experiment holding a single
  condition, through one score per set per library.

### Region to gene links

- [`annotateRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/annotateRegions.md)
  stores a many-to-many region-to-gene table in the object, built from a
  `TxDb`, an `EnsDb`, gene ranges or a table of contacts of your own,
  Hi-C loops or ABC scores for instance, with optional weights.
  `maxDistance`, `flankDistance` and `tssRegion` control the links, and
  [`regionLinks()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/regionLinks.md)
  reads and replaces them.
- [`geneSetsToRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/geneSetsToRegionSets.md)
  turns gene sets into region sets for
  [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md),
  so that the region stays the unit of the test.
- [`aggregateToGenes()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/aggregateToGenes.md)
  writes one row per gene when a gene table is needed. `direction`
  reports a gene whose regions disagree as `mixed` rather than deciding
  for it, next to combined p-values and a weighted score for ranked
  enrichment.

### Designs without replicates

- [`estimateNullDispersion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md)
  estimates the dispersion from rows assumed not to respond and holds
  half of them back, so that
  [`checkNullCalibration()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md)
  and
  [`plotNullCalibration()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotNullCalibration.md)
  can check the assumption on rows the estimate never saw.

### Plots

- Regions:
  [`plotVolcano()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotVolcano.md),
  [`plotResultsMA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotResultsMA.md),
  [`plotRegion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md)
  and
  [`plotTopHeatmap()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotTopHeatmap.md).
- Sets:
  [`plotSetEffect()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetEffect.md),
  [`plotSetDistribution()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetDistribution.md),
  [`plotSetSignal()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetSignal.md)
  and
  [`plotUniverseMatching()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotUniverseMatching.md).
- Samples:
  [`plotRegionPCA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegionPCA.md)
  and
  [`plotSampleCorrelation()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSampleCorrelation.md),
  on the signal or on the Jaccard index of the peak calls, with
  [`computeSamplePCA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeSamplePCA.md)
  and
  [`computeSampleCorrelation()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeSampleCorrelation.md)
  returning the numbers behind them.
- [`computeProfiles()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/computeProfiles.md)
  and
  [`plotProfile()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotProfile.md)
  draw the signal around the regions, as `DiffBind::dba.plotProfile()`
  does.

### Objects, export and example data

- Every class has a validity method, a `show` method, accessors for its
  slots and an `updateObject()` method, and the objects answer the
  BiocGenerics generics `counts()`, `design()`, `plotMA()` and
  `plotPCA()`.
- [`asDGEList()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDGEList.md)
  and
  [`asDESeqDataSet()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/asDESeqDataSet.md)
  hand the counts to *edgeR* or *DESeq2* with the offsets attached.
- [`exportResults()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/exportResults.md)
  writes the tables, a BED coloured by direction and every parameter the
  analysis was run with.
- `BPPARAM` takes a number of cores or any `BiocParallelParam` object in
  every function reading files.
- [`loadExampleData()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/loadExampleData.md)
  reads the two example datasets, H3K4me3 ChIP-seq of two rat strains
  from the EURATRANS project, as distributed by *chromstaRData*, and a
  slice of androgen receptor ChIP-seq in LNCaP cells (GSE284522). The
  scripts that made them are in `inst/scripts`.
