## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = ">",
                      warning = FALSE, message = FALSE, fig.align = "center",
                      dev = "png", dpi = 96, fig.width = 7, fig.height = 4.5)

library(RegionSetDE)
library(dplyr)
library(ggplot2)

## ----load_sheet---------------------------------------------------------------
sampleSheet <- loadExampleData("peakSheet")

# Every column, with the file columns shown as names alone: the paths point at
# wherever the package was installed and say nothing useful here
sampleSheet %>%
  dplyr::mutate(dplyr::across(dplyr::any_of(c("bam", "peaks", "input")),
                              basename))

## ----sheet_columns------------------------------------------------------------
colnames(sampleSheet)

## ----fragment_column, eval = FALSE--------------------------------------------
# # A 'fragmentLength' column in the sheet, one value per library
# counts <- countReads(consensus,
#                      sampleSheet = sampleSheet,
#                      fragmentLength = "fragmentLength")

## ----fragment_estimate, fig.height = 5----------------------------------------
fragmentEstimate <- estimateFragmentLength(sampleSheet = sampleSheet,
                                           pairedEnd = FALSE,
                                           verbose = FALSE)

insertSizes <- estimateFragmentLength(sampleSheet = sampleSheet,
                                      verbose = FALSE)

data.frame(sample = fragmentEstimate$table$sample,
           crossCorrelation = fragmentEstimate$table$fragment.length,
           insertSize = insertSizes$table$fragment.length)

fragmentEstimate$plot

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

S4Vectors::metadata(greylist)$thresholds

## ----consensus----------------------------------------------------------------
consensus <- loadConsensusPeaks(sampleSheet,
                                groupBy = "condition",
                                blacklist = blacklist,
                                greylist = greylist,
                                seqlevelsStyle = "Ensembl",
                                nThreads = 1)

consensus

## ----filtering_log------------------------------------------------------------
filteringLog(consensus) %>% head(3)

## ----min_replicates-----------------------------------------------------------
strictConsensus <- loadConsensusPeaks(sampleSheet,
                                      groupBy = "condition",
                                      blacklist = blacklist,
                                      greylist = greylist,
                                      seqlevelsStyle = "Ensembl",
                                      minReplicates = 3,
                                      nThreads = 1,
                                      verbose = FALSE)

data.frame(group = names(consensusData(consensus)$groups),
           twoOfThree = lengths(consensusData(consensus)$groups),
           threeOfThree = lengths(consensusData(strictConsensus)$groups))

## ----consensus_data-----------------------------------------------------------
consensusInfo <- consensusData(consensus)

consensusInfo$samples

## ----consensus_parts----------------------------------------------------------
lengths(consensusInfo$groups)           # the consensus of every group
length(consensusInfo$total)             # the pooled set that is counted

## ----group_list, eval = FALSE-------------------------------------------------
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# groupList <- consensusGroupList(consensus, seqlevelsStyle = "UCSC")
# 
# annotationList <- lapply(groupList,
#                          ChIPseeker::annotatePeak,
#                          TxDb = txdb)
# 
# ChIPseeker::plotAnnoBar(annotationList)

## ----occupancy_columns--------------------------------------------------------
as.data.frame(S4Vectors::mcols(regionRanges(consensus)$consensus)) %>%
  head(4)

## ----upset--------------------------------------------------------------------
plotPeakUpset(consensus)

## ----jaccard, fig.height = 5.5------------------------------------------------
plotSampleCorrelation(consensus,
                      method = "jaccard",
                      groupBy = "condition",
                      annotationColumns = c("condition", "replicate"))

## ----counting-----------------------------------------------------------------
counts <- countReads(consensus,
                     sampleSheet = sampleSheet)

counts

## ----input_assay--------------------------------------------------------------
countTable(counts, input = TRUE, format = "matrix")[1:3, 1:3]

## ----count_table_wide---------------------------------------------------------
countTable(counts) %>%
  head(3)

## ----count_table_matrix-------------------------------------------------------
countTable(counts, format = "matrix")[1:3, 1:4]

## ----count_table_long---------------------------------------------------------
countTable(counts, format = "long") %>%
  head(4)

## ----lib_info-----------------------------------------------------------------
libInfo(counts, annotationColumns = "condition")

## ----background---------------------------------------------------------------
counts <- countBackground(counts, binSize = 10000)

## ----norm_comparison, fig.height = 5------------------------------------------
plotNormComparison(counts,
                   methods = c("librarySize", "background",
                               "TMM", "readsInRegions"))

## ----norm_effect--------------------------------------------------------------
normalisationEffect <-
  lapply(c("librarySize", "background", "TMM", "RLE", "readsInRegions"),
         function(method) {
           normalised <- normalizeCounts(counts,
                                         method = method,
                                         verbose = FALSE)

           fitted <- fitRegions(normalised,
                                design = ~ condition,
                                engine = "edgeR",
                                verbose = FALSE)

           tested <- testRegions(fitted,
                                 contrast = c("condition", "R1881_24h", "DMSO"),
                                 verbose = FALSE)
           tested <- resultsTable(tested)

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
# SummarizedExperiment::colData(S4Vectors::metadata(counts)$greenlist)

## ----manual_factors-----------------------------------------------------------
backgroundCounts <- normalizeCounts(counts,
                                    method = "background",
                                    verbose = FALSE)

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
sampleInfo(counts)

## ----normalised_counts--------------------------------------------------------
countTable(counts, normalized = TRUE) %>%
  head(3)

## ----pca, fig.height = 4.5----------------------------------------------------
plotRegionPCA(counts, colourBy = "condition", shapeBy = "time")

## ----pca_treated, fig.height = 4.5--------------------------------------------
plotRegionPCA(counts,
              samples = sampleInfo(counts)$treatment == "R1881",
              colourBy = "condition")

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
topRegions(results, contrast = "R1881_24h_vs_R1881_4h", n = 5)

## ----contrast_info------------------------------------------------------------
contrastInfo(results)

## ----volcano------------------------------------------------------------------
plotVolcano(results, contrast = "R1881_24h_vs_DMSO", labelTop = 5)

## ----volcano_late-------------------------------------------------------------
plotVolcano(results, contrast = "R1881_24h_vs_R1881_4h", labelTop = 5)

## ----ma-----------------------------------------------------------------------
plotResultsMA(results, contrast = "R1881_24h_vs_DMSO")

## ----heatmap, fig.height = 6--------------------------------------------------
plotTopHeatmap(results,
               contrast = "R1881_24h_vs_DMSO",
               n = 30,
               annotationColumns = "condition")

## ----profile, fig.height = 6, fig.width = 8-----------------------------------
plotProfile(results,
            contrast = "R1881_24h_vs_DMSO",
            groupBy = "condition")

## ----profile_compute----------------------------------------------------------
profileData <- computeProfiles(results,
                               contrast = "R1881_24h_vs_DMSO",
                               groupBy = "condition")

names(profileData$profiles)
dim(profileData$profiles$DMSO)

## ----profile_lines, fig.height = 3.5------------------------------------------
plotProfile(profileData, style = "lines")

## ----profile_regions, fig.height = 3.5----------------------------------------
consensusRanges <- regionRanges(consensus)$consensus

plotProfile(counts,
            regions = list(
              boundInDMSO = consensusRanges[consensusRanges$peak.DMSO],
              newWithR1881 = consensusRanges[!consensusRanges$peak.DMSO]),
            groupBy = "condition",
            style = "lines")

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

## ----summits------------------------------------------------------------------
summitCounts <- countReads(consensus,
                           sampleSheet = sampleSheet,
                           summits = 200,
                           verbose = FALSE)

head(SummarizedExperiment::rowRanges(summitCounts), 3)

table(BiocGenerics::width(SummarizedExperiment::rowRanges(summitCounts)))

## ----annotate_txdb, eval = FALSE----------------------------------------------
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# annotatedFit <- annotateRegions(fit,
#                                 txdb = txdb,
#                                 orgdb = org.Hs.eg.db::org.Hs.eg.db,
#                                 maxDistance = 1e5)

## ----annotate_contacts, eval = FALSE------------------------------------------
# # One row per link. link.weight is optional and carries no assumed scale.
# annotatedFit <- annotateRegions(fit,
#                                 linkTable = abcLinks,
#                                 weightSource = "ABC score")

## ----annotate_chipseeker, eval = FALSE----------------------------------------
# peakAnnotation <- ChIPseeker::annotatePeak(regionRanges(consensus)$consensus,
#                                            TxDb = txdb)
# 
# peakLinks <- GenomicRanges::granges(peakAnnotation, use.mcols = TRUE)
# 
# annotatedFit <- annotateRegions(fit, linkTable = peakLinks)

## ----annotate_real------------------------------------------------------------
geneSites <- readRDS(system.file("extdata", "lncapAR", "chr19_genes.rds",
                                 package = "RegionSetDE"))

annotatedFit <- annotateRegions(fit,
                                genes = geneSites,
                                maxDistance = 1e5,
                                verbose = FALSE)

head(regionLinks(annotatedFit), 4)

## ----gene_sets_peaks----------------------------------------------------------
geneSets <- readRDS(system.file("extdata", "lncapAR", "chr19_geneSets.rds",
                                package = "RegionSetDE"))

regionSets <- geneSetsToRegionSets(annotatedFit,
                                   geneSets = geneSets,
                                   verbose = FALSE)

head(attr(regionSets, "geneSetSummary"))

## ----test_gene_sets_peaks-----------------------------------------------------
geneSetResults <- testRegionSets(annotatedFit,
                                 contrast = c("condition", "R1881_24h", "DMSO"),
                                 regionSets = regionSets,
                                 verbose = FALSE)

resultsTable(geneSetResults)

## ----aggregate_peaks----------------------------------------------------------
annotatedLate <- testRegions(annotatedFit,
                             contrast = c("condition", "R1881_24h", "DMSO"),
                             FDR = 0.05, log2FC = 1, verbose = FALSE)

geneResults <- aggregateToGenes(annotatedLate,
                                method = c("simes", "potential"),
                                verbose = FALSE)
geneResults

## ----gene_directions_peaks----------------------------------------------------
table(resultsTable(geneResults)$direction)

## ----gene_ora_peaks, eval = FALSE---------------------------------------------
# geneTable <- resultsTable(geneResults)
# 
# upGenes <- dplyr::filter(geneTable, direction == "up", FDR <= 0.05)$gene.id
# 
# clusterProfiler::enrichGO(upGenes,
#                           universe = geneTable$gene.id,
#                           OrgDb = org.Hs.eg.db::org.Hs.eg.db,
#                           keyType = "ENTREZID",
#                           ont = "BP")

## ----gsea_peaks, eval = FALSE-------------------------------------------------
# geneRanks <- stats::setNames(geneTable$score, geneTable$gene.id)
# 
# fgsea::fgseaMultilevel(pathways = geneSets,
#                        stats = geneRanks,
#                        minSize = 5)

## ----gsea_signed, eval = FALSE------------------------------------------------
# geneRanks <- stats::setNames(sign(geneTable$score) * -log10(geneTable$p.value),
#                              geneTable$gene.id)

## ----annotate_flank-----------------------------------------------------------
annotatedFlank <- annotateRegions(fit,
                                  genes = geneSites,
                                  flankDistance = 5e4,
                                  tssRegion = c(-2000, 500),
                                  verbose = FALSE)

table(regionLinks(annotatedFlank)$link.type, useNA = "ifany")

## ----link_type_split----------------------------------------------------------
promoterSets <- geneSetsToRegionSets(annotatedFlank,
                                     geneSets = geneSets,
                                     linkType = "promoter",
                                     verbose = FALSE)

distalSets <- geneSetsToRegionSets(annotatedFlank,
                                   geneSets = geneSets,
                                   linkType = "distal",
                                   verbose = FALSE)

c(promoter = length(promoterSets), distal = length(distalSets))

## ----session_info-------------------------------------------------------------
sessionInfo()

