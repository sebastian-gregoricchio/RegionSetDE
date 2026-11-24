# RegionSetDE 0.99.0

First version, submitted to Bioconductor.

## Regions
* `loadRegions()` reads region sets from BED, narrowPeak, broadPeak, `GRanges` or data.frames into named sets of arbitrary width, kept in a `RegionSetDE` object, and `splitLoadRegions()` splits one file into sets by one of its columns. A single `GRanges` is accepted wherever a region set is expected.
* `loadSampleSheet()` reads a table with one row per sample: the BAM or bigWig file, the peaks called on it, the input it was sequenced against and any annotation column. A DiffBind sample sheet is read as it is, relative paths are resolved against the folder of the sheet, and missing or unindexed files are reported together before anything is counted.
* `loadConsensusPeaks()` builds the regions from the peaks of the samples: a consensus within each group through *consensusRegions*, pooled into a total consensus. `technicalBy` joins the libraries of a biological sample first, so that one sample does not count twice. Every region carries the groups and the samples that called a peak on it (`peak.<group>`, `peak.groups`, `peak.samples`), which follow it into the counts and the results. `consensusData()` and `consensusGroupList()` read the consensus back, and `plotPeakUpset()` draws the occupancy.
* `loadBlacklist()`, `loadGreenlist()` and `availableRegionLists()` give the region lists installed with the package: the ENCODE blacklist version 2 for hg19, hg38, mm10, dm3, dm6, ce10 and ce11, the T2T-CHM13 exclusion set of *excluderanges*, and the CUT&RUN and CUT&Tag greenlists and high signal regions of de Mello et al. for hg38 and mm39. Nothing is downloaded.
* `makeGreylist()` builds a greylist from the input libraries following *GreyListChIP*. The negative binomial is fitted to every window rather than to bootstrap samples, so the result does not depend on a random seed.
* `applyBlacklist()`, `applyWhitelist()` and `applyGreylist()` restrict the regions and record what they removed in the `filtering.log` slot. A list built for another assembly than the one of the regions is refused.

## Counting
* `countReads()` counts BAM files over the regions, once per region or per tile through `tileWidth`. By default each region is split into equal tiles of about `tileWidth` bp, so a region narrower than one tile stays whole and no region ends with a sliver; `tileMethod = "fixed"` cuts exact widths from the start of the region, and returns a `RegionSetDE.counts` object, which is a `RangedSummarizedExperiment`. Paired-end fragments run from the 5' end of the forward mate to the 5' end of the reverse mate, as in csaw: the two mates are matched by read name rather than trusting the template length, so that pairs whose mates run past each other are placed right, and the mapping quality of the second mate is read from its own record, with no need for the `MQ` tag. `pairedEnd = "auto"` reads the layout from the first records of each file.
* `summits` recentres the regions on the summit of the pooled pileup, in windows of `2 * summits + 1` bp as `DiffBind::dba.count()` does, and `summitSource = "peaks"` takes the summits written in the narrowPeak files instead.
* `fragmentLength` takes one value per sample, a column of the sheet or `"auto"`, in which case `estimateFragmentLength()` estimates it from the strand cross-correlation of the reads.
* The inputs of the samples are counted into an `input` assay, and nothing subtracts them from the counts.
* `blacklist`, `greylist` and `discardRegions` remove regions or ignore their reads at the counting, and `excludeChromosomes` keeps chromosomes such as chrM out of the library sizes.
* BAM files, bigWig files, regions and lists naming the chromosomes of one assembly in different styles, `chr1` in some and `1` in others, are matched name by name. CSI indices are found as well as BAI ones.
* CRAM files are read wherever BAM files are, in `countReads()`, `countBackground()`, `countGreenlist()`, `makeGreylist()`, `estimateFragmentLength()`, `computeProfiles()`, `libInfo()` and the sample sheet, and both formats can be mixed. BAM and CRAM files are read by the same C code, built on Rhtslib, so they give the same counts; CRAM files are read without the reference genome, decoding only the fields CRAM stores apart from the bases. For a CRAM file `libInfo()` reads the whole file once, since a CRAM index holds no read counts.
* `countBigwig()` counts coverage from bigWig files without rounding it, and the negative binomial and voom engines refuse such an object unless `countLike` declares the coverage count-like. `loadCounts()` takes a count matrix computed elsewhere.
* `countBackground()` counts genome-wide bins for the normalisation and the null estimates, and `countGreenlist()` counts the libraries over the greenlist.
* `libInfo()` summarises the libraries in one table: the records of every file, the fragments kept, the fragments in the regions, the FRiP, the fragment length and the input of each sample.
* `countTable()` returns the values of any object as a table, raw or normalised, per region or per tile, wide or long with the `colData` attached.
* `progressBar` draws a bar while the files are read, by default only in interactive sessions.

## Normalisation and filtering
* `normalizeRegionCounts()` estimates the scaling factors from background bins, from the regions themselves (TMM, RLE or the reads in the regions), from the greenlist, or takes them from outside, a spike-in for instance. `backgroundHoldout` keeps part of the bins out of the estimate so that the calibration can be checked on them later.
* `plotNormComparison()` and `plotSetMA()` show how far the choice of normalisation moves the samples and the sets before anything is fitted on it.
* `filterRegions()` removes the rows carrying too little signal, on their average abundance alone, so that the choice is independent of the contrast tested afterwards.
* `selectSamples()` filters the samples with `dplyr::filter` syntax on the `colData`, and `splitSamples()` splits an object by mark or assay.

## Fitting and testing
* `fitRegions()` fits one model per region with *edgeR*, *limma*-voom, *limma*-trend, `variancePartition::dream` or *DESeq2*, reading the normalisation out of the object as offsets. voom and dream receive it as effective library sizes, since the offsets of a `DGEList` are read differently, or not at all, by each version of `limma::voom`, and the `"loess"` offsets of every region go through voom with one library size per observation.
* `testRegions()` tests a contrast or a named list of them. A contrast can be a coefficient name, an expression over the design columns, a numeric vector or `c("column", "groupA", "groupB")`, and `pairwiseContrasts()` writes out every pair of levels of a column. Tiled regions are combined back to one row per region through `csaw::combineTests()` or a representative tile.
* The results carry the distribution of `stat` and its degrees of freedom, for a power calculation for instance, and `contrastInfo()` summarises every contrast of an object.
* `updateThresholds()` labels the regions again with other `FDR` and `log2FC` cut-offs, without testing again.
* `peakOccupancyTable()` and `plotPeakOccupancy()` cross the groups that called a region with what the test said about it.

## Region sets
* `testRegionSets()` tests whole sets, competitively through `limma::cameraPR` and self-contained through `limma::fry`. `fry.p` is the directional test and `fry.mixed.p` the one for a set changing in either direction. `camera.mixed.p` is the competitive version of the latter, `limma::cameraPR` on the ranks of the absolute statistics, for a set changing more than its comparison in either direction, and `n.up` and `n.down` count the regions changing inside each set. The effect of a set comes with a confidence interval built on the libraries and a second interval for how much it varies between regions.
* `makeSetUniverse()` builds the rows a set is compared against, matched on width and abundance.
* `testSetContrast()` asks whether a contrast affects two sets differently, which is where a redistribution claim belongs.
* `tileHandling`, `overlapPolicy` and `overlapWithinSet` decide how tiles, sets covering the same chromatin and overlapping regions inside one set are handled.
* `scoreRegionSets()` compares the sets to each other in an experiment holding a single condition, through one score per set per library.

## Region to gene links
* `annotateRegions()` stores a many-to-many region-to-gene table in the object, built from a `TxDb`, an `EnsDb`, gene ranges or a table of contacts of your own, Hi-C loops or ABC scores for instance, with optional weights. `maxDistance`, `flankDistance` and `tssRegion` control the links, and `regionLinks()` reads and replaces them.
* `geneSetsToRegionSets()` turns gene sets into region sets for `testRegionSets()`, so that the region stays the unit of the test.
* `aggregateToGenes()` writes one row per gene when a gene table is needed. `direction` reports a gene whose regions disagree as `mixed` rather than deciding for it, next to combined p-values and a weighted score for ranked enrichment.

## Designs without replicates
* `estimateNullDispersion()` estimates the dispersion from rows assumed not to respond and holds half of them back, so that `checkNullCalibration()` and `plotNullCalibration()` can check the assumption on rows the estimate never saw. `checkNullCalibration()` also fits the regions the same way and compares the median of their fold changes with the one of the null rows, and returns `"regions shifted"` when the two do not sit on the same baseline, the case of libraries differing in their share of fragments in the regions under background normalisation.

## Plots
* Regions: `plotVolcano()`, `plotResultsMA()`, `plotRegion()` and `plotTopHeatmap()`.
* Sets: `plotSetEffect()`, `plotSetDistribution()`, `plotSetSignal()` and `plotUniverseMatching()`.
* Samples: `plotRegionPCA()` and `plotSampleCorrelation()`, on the signal or on the Jaccard index of the peak calls, with `computeSamplePCA()` and `computeSampleCorrelation()` returning the numbers behind them.
* `computeProfiles()` and `plotProfile()` draw the signal around the regions, as `DiffBind::dba.plotProfile()` does.

## Objects, export and example data
* Every class has a validity method, a `show` method, accessors for its slots and an `updateObject()` method, and the objects answer the BiocGenerics generics `counts()`, `design()`, `plotMA()` and `plotPCA()`.
* `asDGEList()` and `asDESeqDataSet()` hand the counts to *edgeR* or *DESeq2* with the offsets attached.
* `exportResults()` writes the tables, a BED coloured by direction and every parameter the analysis was run with.
* `BPPARAM` takes a number of cores or any `BiocParallelParam` object in every function reading files.
* `loadExampleData()` reads the two example datasets, H3K4me3 ChIP-seq of two rat strains from the EURATRANS project, as distributed by *chromstaRData*, and a slice of androgen receptor ChIP-seq in LNCaP cells (GSE284522). The scripts that made them are in `inst/scripts`.
