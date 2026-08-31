# RegionSetDE

Differential analysis of chromatin signal over user-defined genomic
region sets, without peak calling

## Introduction

The concept behind `RegionSetDE` (*Region Set Differential Enrichment*)
is to ask whether chromatin signal changes over regions that you define,
rather than over regions that a peak caller happens to find. Regions are
grouped into named sets of arbitrary width and number, such as promoter
classes, enhancer catalogues, chromatin states or motif-derived binding
sites, and a single model fit is then tested at two levels: one region
at a time, and one whole set at a time.

The second level is the reason the package exists. A great deal of
chromatin biology is phrased as a question about a class of elements
rather than about individual loci. *Does H3K4me3 redistribute away from
CpG island promoters in this mutant?* is a question about a set, and
answering it by counting how many members of the set crossed an FDR
threshold conflates the size of the effect with the power to detect it.
Here the set is the unit of the test, and the effect size comes with a
confidence interval that accounts for regions inside a set not being
independent of one another.

Because no peak calling happens anywhere in the pipeline, the region
definitions come from you and stay fixed across conditions. Nothing is
redefined when the signal moves, which is what makes a redistribution
visible instead of being absorbed into a new set of peak boundaries.

  

### Main features

- **No peak calling.** Regions are supplied as BED, narrowPeak,
  broadPeak, `GRanges` or data.frames, and are the same in every
  condition.
- **Two levels of testing from one fit.**
  [`testRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md)
  for individual regions,
  [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md)
  for whole sets, and
  [`testSetContrast()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md)
  for the difference between two sets.
- **Effect sizes before p-values.** Set-level results carry a confidence
  interval inflated for the correlation between neighbouring regions, so
  a tiny shift over 30,000 promoters is not mistaken for a finding.
- **Normalisation you can defend.** Scaling factors from background
  bins, from the regions themselves, or supplied from a spike-in or
  greenlist, with
  [`plotNormComparison()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotNormComparison.md)
  to show how far the choice moves the samples before anything is fitted
  on it.
- **Four engines.** `edgeR`, `limma-voom`,
  [`variancePartition::dream`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
  for random effects, and `DESeq2`.
- **Calibration you can check.**
  [`checkNullCalibration()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md)
  runs the contrast on rows known to be null and shows whether the
  p-values come out uniform.
- **Tiled counting.** Regions can be split into fixed-width tiles, so a
  change confined to part of a wide domain is not diluted across the
  whole of it.
- **No replicates needed.** One library per condition is testable: the
  dispersion is estimated from rows assumed not to respond, half of them
  held back so
  [`checkNullCalibration()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md)
  can check the assumption rather than confirm it.

  

### Citation

If you use this package, please cite:

*RegionSetDE*: differential analysis of chromatin signal over
user-defined genomic region sets.  
Gregoricchio S.  
*Manuscript in preparation*  

  

## Installation

### Bioconductor

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("RegionSetDE")
```

### Developmental version

``` r

## Install remotes from CRAN (if not already installed)
if (!require("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install the RegionSetDE package
remotes::install_github("sebastian-gregoricchio/RegionSetDE",
                        build_manual = TRUE,
                        build_vignettes = TRUE)
```

### Possible installation issues

Most of the dependencies live on Bioconductor rather than on the CRAN,
so `BiocManager` has to be installed and its repositories active before
`remotes::install_github()` is able to resolve them:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

setRepositories(ind = 1:8)
```

Building the vignette additionally requires
[`BiocStyle`](https://bioconductor.org/packages/BiocStyle/), and a few
of the plotting functions require packages listed under `Suggests`:

``` r

BiocManager::install(c("BiocStyle", "ComplexHeatmap", "circlize", "ggrepel"))

# Only for the random-effects engine
BiocManager::install(c("variancePartition", "lme4"))

# Only for the DESeq2 engine
BiocManager::install("DESeq2")
```

  

## Quick start

``` r

library(RegionSetDE)

# Region sets: a named list of BED files, GRanges or data.frames
regions <- loadRegions(list(promoters = "annotation/promoters.bed",
                            enhancers = "annotation/enhancers.bed",
                            CTCF      = "peaks/CTCF.narrowPeak"),
                       genomeAssembly = "hg38")

regions <- applyBlacklist(regions, blacklist = encodeBlacklist)

# Counting, plus the background bins used for the normalisation and the null
counts <- countReads(regions,
                     bamFiles = bamPaths,
                     sampleMetadata = sampleTable,
                     pairedEnd = TRUE,
                     nThreads = 4)

counts <- countBackground(counts, binSize = 10000, nThreads = 4)

# Normalisation and filtering
counts <- normalizeCounts(counts, method = "background")
counts <- filterRegions(counts, method = "background", foldChange = 2)

# One fit, two levels of testing
fit <- fitRegions(counts, design = ~ condition, engine = "edgeR")

results    <- testRegions(fit,     contrast = c("condition", "treated", "control"))
setResults <- testRegionSets(fit,  contrast = c("condition", "treated", "control"))

topRegions(results, n = 20)
resultsTable(setResults)
```

Every object in the package ships with a worked example, so the whole
pipeline can be run without any data of your own:

``` r

counts <- loadExampleData("counts")
fit    <- loadExampleData("fit")

results <- testRegions(fit, contrast = c("condition", "SHR", "BN"))
plotVolcano(results)
```

  

## Which analysis for which question

| Question | Function |
|:---|:---|
| Which individual regions changed? | [`testRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegions.md), then [`topRegions()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/topRegions.md) |
| Did this class of regions respond as a class? | [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md) |
| Did the effect differ between two classes? | [`testSetContrast()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testSetContrast.md) |
| Did the mark redistribute rather than change globally? | [`testRegionSets()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/testRegionSets.md), competitive test significant and self-contained not |
| Is a global shift technical or biological? | [`plotNormComparison()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotNormComparison.md), [`plotSetMA()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotSetMA.md) |
| Are my p-values trustworthy? | [`checkNullCalibration()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/checkNullCalibration.md) |
| Can I test without replicates? | [`estimateNullDispersion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/estimateNullDispersion.md), then `fitRegions(dispersion = )` |
| Is the competitive comparison fair? | [`plotUniverseMatching()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotUniverseMatching.md) |
| Where inside the region did the change happen? | `countReads(tileWidth = )`, then [`plotRegion()`](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/plotRegion.md) |

  

## Documentation

With the package there are available:

- [web-manual](https://sebastian-gregoricchio.github.io/RegionSetDE/reference/index.html):
  technical manual of each function;
- [overview
  vignette](https://sebastian-gregoricchio.github.io/RegionSetDE/articles/RegionSetDE.vignette.html):
  presents the whole workflow, what each object contains, how to extract
  its parts, and how to read the results.

Note that the vignette can be inspected on R as well by typing
`browseVignettes("RegionSetDE")`.

  

## Package history and releases

A list of all releases and the respective description of the changes
applied can be found
[here](https://sebastian-gregoricchio.github.io/RegionSetDE/news/index.html).

  

------------------------------------------------------------------------

## Contact

For any suggestion, bug fixing or commentary please report it in the
[issues](https://github.com/sebastian-gregoricchio/RegionSetDE/issues)/[request](https://github.com/sebastian-gregoricchio/RegionSetDE/pulls)
tab of this repository.

## License

This package is under a GNU General Public License (version 3).

  

#### Contributors

![contributors](https://badges.pufler.dev/contributors/sebastian-gregoricchio/RegionSetDE?size=50&padding=5&bots=true)

contributors
