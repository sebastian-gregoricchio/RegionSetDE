## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = ">",
                      warning = FALSE, message = FALSE, fig.align = "center",
                      dev = "png", dpi = 96, fig.width = 7, fig.height = 4.5)

library(RegionSetDE)
library(dplyr)
library(ggplot2)

## ----load_sheet---------------------------------------------------------------
sampleSheet <- loadExampleData("peakSheet")

sampleSheet %>%
  dplyr::select(sample, condition, treatment, time, replicate)

## ----sheet_columns------------------------------------------------------------
colnames(sampleSheet)

## ----available_lists----------------------------------------------------------
availableRegionLists(genome = "hg38")

## ----blacklist----------------------------------------------------------------
blacklist <- loadBlacklist("hg38", seqlevelsStyle = "Ensembl")
blacklist

## ----greylist_default---------------------------------------------------------
greylist <- makeGreylist(sampleSheet)

S4Vectors::metadata(greylist)$thresholds

## ----greylist-----------------------------------------------------------------
greylist <- makeGreylist(sampleSheet, binSize = 10000)

S4Vectors::metadata(greylist)$thresholds %>%
  dplyr::select(input, fragments, mean, threshold, regions, greylisted.bp)

## ----artefacts----------------------------------------------------------------
artefactRegions <- c(GenomicRanges::granges(blacklist), GenomicRanges::granges(greylist))

## ----consensus----------------------------------------------------------------
consensus <- loadConsensusPeaks(sampleSheet,
                                groupBy = "condition",
                                excludeRegions = artefactRegions,
                                seqlevelsStyle = "Ensembl",
                                genomeAssembly = "hg38")

consensus

## ----min_replicates-----------------------------------------------------------
strictConsensus <- loadConsensusPeaks(sampleSheet,
                                      groupBy = "condition",
                                      excludeRegions = artefactRegions,
                                      seqlevelsStyle = "Ensembl",
                                      minReplicates = 3,
                                      verbose = FALSE)

data.frame(group = names(consensusData(consensus)$groups),
           twoOfThree = lengths(consensusData(consensus)$groups),
           threeOfThree = lengths(consensusData(strictConsensus)$groups))

## ----consensus_data-----------------------------------------------------------
consensusInfo <- consensusData(consensus)

consensusInfo$samples %>%
  dplyr::select(sample, group, n.peaks, n.excluded)

## ----consensus_parts----------------------------------------------------------
lengths(consensusInfo$groups)           # the consensus of every group
length(consensusInfo$total)             # the pooled region set that will be counted

## ----occupancy_columns--------------------------------------------------------
as.data.frame(S4Vectors::mcols(consensus@regions$consensus)) %>%
  dplyr::select(dplyr::starts_with("peak.")) %>%
  head(4)

## ----upset, fig.height = 4----------------------------------------------------
plotPeakUpset(consensus)

## ----counting-----------------------------------------------------------------
counts <- countReads(consensus,
                     sampleSheet = sampleSheet,
                     nThreads = 2)

counts

## ----count_table_wide---------------------------------------------------------
countTable(counts) %>%
  head(3)

## ----count_table_matrix-------------------------------------------------------
countTable(counts, format = "matrix")[1:3, 1:4]

## ----count_table_long---------------------------------------------------------
countTable(counts, format = "long") %>%
  dplyr::select(region.id, sample, condition, counts) %>%
  head(4)

## ----lib_info-----------------------------------------------------------------
libInfo(counts, annotationColumns = "condition")

## ----background---------------------------------------------------------------
counts <- countBackground(counts, binSize = 10000, nThreads = 2)

## ----norm_comparison, fig.height = 5------------------------------------------
plotNormComparison(counts,
                   methods = c("librarySize", "background", "TMM", "readsInRegions"))

## ----norm_effect--------------------------------------------------------------
normalisationEffect <-
  lapply(c("librarySize", "background", "TMM", "RLE", "readsInRegions"),
         function(method) {
           normalised <- normalizeCounts(counts, method = method, verbose = FALSE)
           fitted <- fitRegions(normalised, design = ~ condition, engine = "edgeR", verbose = FALSE)
           tested <- resultsTable(testRegions(fitted,
                                              contrast = c("condition", "R1881_24h", "DMSO"),
                                              verbose = FALSE))

           data.frame(method = method,
                      up = sum(tested$diff.status == "up"),
                      down = sum(tested$diff.status == "down"),
                      median.log2FC = round(median(tested$log2FC), 2))
         })

dplyr::bind_rows(normalisationEffect)

## ----greenlist, eval = FALSE--------------------------------------------------
# greenlist <- loadGreenlist("hg38", assay = "cuttag")
# 
# counts <- countGreenlist(counts, greenlist = greenlist)
# counts <- normalizeCounts(counts, method = "greenlist")

## ----greenlist_check, eval = FALSE--------------------------------------------
# greenlistCounts <- SummarizedExperiment::assay(S4Vectors::metadata(counts)$greenlist, "counts")
# colSums(greenlistCounts)

## ----manual_factors-----------------------------------------------------------
backgroundCounts <- normalizeCounts(counts, method = "background", verbose = FALSE)

backgroundFactors <- setNames(sampleInfo(backgroundCounts)$scaling.factor,
                              colnames(backgroundCounts))
round(backgroundFactors, 2)

manualCounts <- normalizeCounts(counts,
                                method = "manual",
                                scalingFactors = backgroundFactors,
                                factorType = "division",
                                verbose = FALSE)

all.equal(countTable(manualCounts, normalized = TRUE, format = "matrix"),
          countTable(backgroundCounts, normalized = TRUE, format = "matrix"))

## ----normalise----------------------------------------------------------------
counts <- normalizeCounts(counts, method = "background")

## ----scaling_factors----------------------------------------------------------
sampleInfo(counts, columns = c("sample", "condition", "library.size", "scaling.factor"))

## ----normalised_counts--------------------------------------------------------
countTable(counts, normalized = TRUE) %>%
  head(3)

## ----pca, fig.height = 4.5----------------------------------------------------
plotRegionPCA(counts, colourBy = "condition", shapeBy = "time")

## ----correlation, fig.height = 5.5--------------------------------------------
plotSampleCorrelation(counts,
                      annotationColumns = c("condition", "replicate"),
                      method = "spearman")

## ----fit----------------------------------------------------------------------
fit <- fitRegions(counts, design = ~ condition, engine = "edgeR")

## ----single_contrast----------------------------------------------------------
resultsEarly <- testRegions(fit,
                            contrast = c("condition", "R1881_4h", "DMSO"),
                            FDR = 0.05,
                            log2FC = 1)

resultsEarly

## ----contrasts----------------------------------------------------------------
contrasts <- pairwiseContrasts(fit, column = "condition")
contrasts

## ----test---------------------------------------------------------------------
results <- testRegions(fit,
                       contrast = contrasts,
                       FDR = 0.05,
                       log2FC = 1,
                       carryCounts = TRUE)

results

## ----results_table------------------------------------------------------------
resultsLate <- results$R1881_24h_vs_DMSO

resultsTable(resultsLate) %>%
  head(3)

## ----results_all--------------------------------------------------------------
resultsTable(results) %>%
  dplyr::filter(diff.status != "null") %>%
  dplyr::count(contrast, diff.status)

## ----top_regions--------------------------------------------------------------
topRegions(results, contrast = "R1881_24h_vs_R1881_4h", n = 5) %>%
  dplyr::select(region.id, log2FC, FDR, peak.R1881_4h, peak.R1881_24h)

## ----contrast_info------------------------------------------------------------
contrastInfo(results)

## ----volcano------------------------------------------------------------------
plotVolcano(results, contrast = "R1881_24h_vs_DMSO", labelTop = 5)

## ----volcano_late-------------------------------------------------------------
plotVolcano(results, contrast = "R1881_24h_vs_R1881_4h", labelTop = 5)

## ----ma-----------------------------------------------------------------------
plotResultsMA(results, contrast = "R1881_24h_vs_DMSO")

## ----heatmap, fig.height = 6--------------------------------------------------
plotTopHeatmap(results, contrast = "R1881_24h_vs_DMSO", n = 30, annotationColumns = "condition")

## ----occupancy_table----------------------------------------------------------
peakOccupancyTable(results, contrast = "R1881_24h_vs_DMSO")

## ----occupancy_plot, fig.height = 4.5-----------------------------------------
plotPeakOccupancy(results, contrast = "R1881_24h_vs_DMSO")

## ----export, eval = FALSE-----------------------------------------------------
# exportResults(results,
#               path = "results",
#               onlyChanging = TRUE,
#               splitByDirection = TRUE)

## ----tiles_count, eval = FALSE------------------------------------------------
# tiledCounts <- countReads(consensus,
#                           sampleSheet = sampleSheet,
#                           tileWidth = 200)

## ----tiles_test, eval = FALSE-------------------------------------------------
# tiledCounts <- countBackground(tiledCounts, binSize = 10000)
# tiledCounts <- normalizeCounts(tiledCounts, method = "background")
# 
# tiledFit <- fitRegions(tiledCounts, design = ~ condition, engine = "edgeR")
# 
# tiledResults <- testRegions(tiledFit,
#                             contrast = c("condition", "R1881_24h", "DMSO"),
#                             combineMethod = "simes")
# 
# resultsTable(tiledResults)     # one row per region
# tileTable(tiledResults)        # one row per tile

## ----recentre, eval = FALSE---------------------------------------------------
# recentredConsensus <- loadConsensusPeaks(sampleSheet,
#                                          groupBy = "condition",
#                                          excludeRegions = artefactRegions,
#                                          seqlevelsStyle = "Ensembl",
#                                          recentre = TRUE,
#                                          width = 400)

## ----session_info-------------------------------------------------------------
sessionInfo()

