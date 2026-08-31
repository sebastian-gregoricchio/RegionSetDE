## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = ">",
                      warning = FALSE, message = FALSE, fig.align = "center",
                      dev = "png", dpi = 96, fig.width = 7, fig.height = 4.5)

library(RegionSetDE)
library(dplyr)

## ----installation, eval = FALSE-----------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install("RegionSetDE")
# 
# # Development version
# # devtools::install_github("sebastian-gregoricchio/RegionSetDE")

## ----load_example-------------------------------------------------------------
counts <- loadExampleData("counts")

## ----make_data_path, eval = FALSE---------------------------------------------
# file.edit(system.file("scripts", "make-data.R", package = "RegionSetDE"))

## ----load_regions_generic, eval = FALSE---------------------------------------
# regions <- loadRegions(list(promoters = "annotation/promoters.bed",
#                             enhancers = enhancerGRanges,
#                             CTCF      = "peaks/CTCF.narrowPeak"),
#                        genomeAssembly = "hg38")

## ----split_load_regions-------------------------------------------------------
regionTable <- loadExampleData("regions", verbose = FALSE)

regionRanges <- GenomicRanges::makeGRangesFromDataFrame(regionTable,
                                                        keep.extra.columns = TRUE)

regions <- splitLoadRegions(regionRanges,
                            splitBy = "setName",
                            genomeAssembly = "rn4",
                            verbose = FALSE)

regions

## ----regions_accessors--------------------------------------------------------
regionSetNames(regions)

lengths(regions@regions)

head(regions@regions$promoterCpG, 3)

## ----apply_blacklist----------------------------------------------------------
exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)

regions <- applyBlacklist(regions,
                          blacklist = exclusionRegions,
                          verbose = FALSE)

regions@filtering.log

## ----count_reads, eval = FALSE------------------------------------------------
# counts <- countReads(regions,
#                      bamFiles = bamPaths,
#                      sampleNames = c("BN_1", "BN_2", "SHR_1", "SHR_2"),
#                      sampleMetadata = data.frame(sample = c("BN_1", "BN_2", "SHR_1", "SHR_2"),
#                                                  condition = c("BN", "BN", "SHR", "SHR")),
#                      pairedEnd = FALSE,
#                      fragmentLength = 180,
#                      minMapq = 10,
#                      nThreads = 4)

## ----counts_structure---------------------------------------------------------
counts <- loadExampleData("counts", verbose = FALSE)

dim(counts)

head(SummarizedExperiment::assay(counts, "counts"), 3)

SummarizedExperiment::colData(counts)

## ----counts_rowdata-----------------------------------------------------------
head(SummarizedExperiment::rowData(counts), 3)

table(SummarizedExperiment::rowData(counts)$region.set)

## ----counts_provenance--------------------------------------------------------
counts@counting.level

counts@genome.assembly

## ----count_background, eval = FALSE-------------------------------------------
# counts <- countBackground(counts, binSize = 10000, nThreads = 4)

## ----background_access--------------------------------------------------------
backgroundBins <- S4Vectors::metadata(counts)$background

dim(backgroundBins)

head(backgroundBins, 3)

## ----normalize----------------------------------------------------------------
counts <- normalizeCounts(counts, method = "background", verbose = FALSE)

SummarizedExperiment::colData(counts)

## ----norm_assays--------------------------------------------------------------
SummarizedExperiment::assayNames(counts)

## ----plot_norm_comparison, fig.height = 4-------------------------------------
plotNormComparison(loadExampleData("counts", verbose = FALSE))

## ----plot_set_ma, fig.height = 4.5--------------------------------------------
plotSetMA(counts, groupBy = "condition")

## ----filter_regions-----------------------------------------------------------
nrow(counts)

counts <- filterRegions(counts, method = "background", foldChange = 2, verbose = FALSE)

nrow(counts)

table(SummarizedExperiment::rowData(counts)$region.set)

## ----plot_qc, fig.height = 4--------------------------------------------------
plotRegionPCA(counts, colourBy = "condition", shapeBy = "sex")

## ----plot_correlation, fig.height = 4.5---------------------------------------
plotSampleCorrelation(counts, groupBy = "condition")

## ----fit_regions--------------------------------------------------------------
fit <- fitRegions(counts,
                  design = ~ condition,
                  engine = "edgeR",
                  verbose = FALSE)

fit

## ----fit_accessors------------------------------------------------------------
fitCounts(fit)

class(fitObject(fit))

## ----plot_universe, fig.height = 4--------------------------------------------
setResults <- testRegionSets(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)

plotUniverseMatching(setResults, set = "promoterCpG", covariate = "abundance")

## ----test_regions-------------------------------------------------------------
results <- testRegions(fit,
                       contrast = c("condition", "SHR", "BN"),
                       verbose = FALSE)

results

## ----results_table------------------------------------------------------------
resultTable <- resultsTable(results)

head(resultTable, 3)

## ----results_accessors--------------------------------------------------------
contrastName(results)

head(resultRanges(results), 2)

resultCounts(results)

## ----top_regions--------------------------------------------------------------
topRegions(results, n = 5, FDR = 1)

## ----top_regions_set----------------------------------------------------------
topRegions(results, n = 3, set = "promoterCpG", FDR = 1, sortBy = "log2FC")

## ----tile_table, eval = FALSE-------------------------------------------------
# head(tileTable(tiledResults))

## ----plot_volcano, fig.height = 5---------------------------------------------
plotVolcano(results, labelTop = 3)

## ----plot_ma, fig.height = 5--------------------------------------------------
plotResultsMA(results)

## ----plot_region, fig.height = 4----------------------------------------------
topRegion <- topRegions(results, n = 1, FDR = 1)$region.id

plotRegion(results, region = topRegion, groupBy = "condition")

## ----plot_heatmap, fig.height = 5---------------------------------------------
plotTopHeatmap(results, n = 20, FDR = 1)

## ----test_region_sets---------------------------------------------------------
setResults <- testRegionSets(fit,
                             contrast = c("condition", "SHR", "BN"),
                             method = c("camera", "fry"),
                             verbose = FALSE)

setResults

## ----set_results_table--------------------------------------------------------
setTable <- resultsTable(setResults)

setTable

## ----plot_set_effect, fig.height = 4------------------------------------------
plotSetEffect(setResults)

## ----plot_set_distribution, fig.height = 4------------------------------------
plotSetDistribution(setResults)

## ----plot_set_signal, fig.height = 4.5----------------------------------------
plotSetSignal(setResults, groupBy = "condition")

## ----test_set_contrast--------------------------------------------------------
setContrast <- testSetContrast(fit,
                               contrast = c("condition", "SHR", "BN"),
                               set1 = "promoterCpG",
                               set2 = "intergenic",
                               verbose = FALSE)

setContrast

## ----null_dispersion----------------------------------------------------------
nullDispersion <- estimateNullDispersion(loadExampleData("counts", verbose = FALSE) |>
                                           normalizeCounts(method = "background", verbose = FALSE),
                                         source = "background",
                                         verbose = FALSE)

nullDispersion

## ----null_calibration, fig.height = 4-----------------------------------------
calibration <- checkNullCalibration(fit,
                                    contrast = c("condition", "SHR", "BN"),
                                    source = "background",
                                    verbose = FALSE)

plotNullCalibration(calibration)

## ----null_from_intergenic, error = TRUE---------------------------------------
try({
estimateNullDispersion(counts, source = "regionSet",
                       regionSets = "intergenic", verbose = FALSE)
})

## ----export_results-----------------------------------------------------------
outputDirectory <- file.path(tempdir(), "regionSetDE-vignette")
dir.create(outputDirectory, showWarnings = FALSE)

exportResults(results,
              path = outputDirectory,
              prefix = "SHR_vs_BN",
              verbose = FALSE)

list.files(outputDirectory)

## ----read_provenance----------------------------------------------------------
parameterFile <- list.files(outputDirectory, pattern = "parameters", full.names = TRUE)

head(read.delim(parameterFile), 10)

## ----cleanup, include = FALSE-------------------------------------------------
unlink(outputDirectory, recursive = TRUE)

## ----session_info-------------------------------------------------------------
sessionInfo()

