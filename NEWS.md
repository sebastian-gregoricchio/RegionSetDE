# RegionSetDE 0.99.0
First version.

## Regions and counting
* `loadRegions()` reads region sets from BED, narrowPeak, broadPeak or GRanges and keeps them as named sets of arbitrary width inside a `RegionSetDE` object.
* `applyBlacklist()` and `applyWhitelist()` restrict the regions, and record what they removed in the `filtering.log` slot.
* `countReads()` and `countBigwig()` count over the regions, either one row per region or one row per tile, and return a `RegionSetDE.counts` object.
* `countBackground()` counts genome-wide bins alongside the regions, for the normalisation and for the null estimates.


## Normalisation and filtering
* `normalizeCounts()` estimates scaling factors from the background bins, from the regions themselves, or takes them from outside, for instance from a spike-in or a greenlist.
* `plotSetMA()` and `plotNormComparison()` show what a normalisation did before anything is fitted on it.
* `filterRegions()` removes the rows that carry too little signal to say anything, on average abundance alone so that the choice is independent of the contrast tested afterwards. Width-adjusted by default, since a threshold in reads otherwise keeps every broad region and drops every narrow one.

## Sample selection
* `selectSamples()` filters the samples with `dplyr::filter` syntax on the `colData`, and `splitSamples()` splits an object into one piece per mark or assay. Both drop the stored normalisation by default, since factors estimated across marks describe a library composition that no longer exists once the object is subset.


## Fitting and testing
* `fitRegions()` fits one model per region with `edgeR`, `limma-voom`, `variancePartition::dream` or `DESeq2`, reading the normalisation out of the object as offsets rather than recomputing it.
* `testRegions()` tests a contrast, or a named list of them, and combines tiled regions back to one row per region through `csaw::combineTests`.
* Contrasts can be given as a coefficient name, an expression over the design columns, a numeric vector, or as `c("column", "groupA", "groupB")`, which works whatever the reference level is.


## Region sets
* `testRegionSets()` asks whether a set responds as a whole, competitively through `limma::cameraPR` and self-contained through `limma::fry`, with the variance inflated for the correlation between regions so that a large set does not come out certain by virtue of being large.
* `testSetContrast()` asks whether a contrast affects one set differently from another.
* `makeSetUniverse()` builds the comparison universe, matched on width and abundance; `fitRegions()` builds one automatically and keeps it in the fit.


## Designs without replicates
* `estimateNullDispersion()` reads the between-sample variation off rows assumed not to respond, so that a design with one sample per condition has a dispersion to be tested against. `fitRegions()` calls it on its own when the design leaves no residual.
* `checkNullCalibration()` runs the same contrast on rows that should not respond and reports how many come out significant anyway, broken down by abundance, with a suggested dispersion when the current one is off.


## Plots
* Per region: `plotVolcano()`, `plotResultsMA()`, `plotRegion()`, `plotTopHeatmap()`.
* Per set: `plotSetEffect()`, `plotSetDistribution()`, `plotSetSignal()`, `plotUniverseMatching()`.
* Samples: `plotRegionPCA()` and `plotSampleCorrelation()`, both able to draw the same figure with and without the normalisation so that a grouping caused by the scaling factors can be told apart from one in the data.


## Export
* `asDGEList()`, `as(x, "DGEList")` and `asDESeqDataSet()` hand the counts to `edgeR` or `DESeq2` with the offsets attached the right way round.
* `exportResults()` writes the table, a BED coloured by direction, and every parameter the analysis was run with.
