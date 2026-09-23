# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

testthat::skip_if_not_installed("consensusRegions")

sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)
counts <- countReads(consensus, sampleSheet = sampleSheet, verbose = FALSE)


test_that("libInfo puts the reads of the files next to the ones of the regions", {

  infoTable <- libInfo(counts, annotationColumns = "condition")

  expect_identical(colnames(infoTable),
                   c("sample", "condition", "paired.end", "bam.reads", "bam.mapped", "library.size", "reads.in.regions", "FRiP"))
  expect_identical(infoTable$sample, colnames(counts))
  expect_identical(infoTable$condition, sampleSheet$condition)

  # The index gives every record of the file, which is two per fragment on paired-end data
  expect_identical(infoTable$bam.reads[1],
                   as.numeric(Rsamtools::countBam(RegionSetDE:::.bamWithIndex(sampleSheet$bam[1]))$records))
  expect_true(all(infoTable$bam.mapped > infoTable$library.size))

  expect_equal(infoTable$FRiP, round(infoTable$reads.in.regions / infoTable$library.size, 4))
  expect_lt(max(infoTable$FRiP[infoTable$condition == "DMSO"]), min(infoTable$FRiP[infoTable$condition == "R1881_24h"]))
})


test_that("libInfo counts a region shared by two sets once", {

  regionRanges <- consensus@regions$consensus
  twoSets <- GenomicRanges::GRangesList(first = regionRanges, second = regionRanges)

  doubledCounts <- countReads(twoSets, sampleSheet = sampleSheet, verbose = FALSE)

  expect_identical(libInfo(doubledCounts)$reads.in.regions, libInfo(counts)$reads.in.regions)
})


test_that("libInfo refuses what it cannot summarise", {

  expect_error(libInfo(counts, annotationColumns = "tissue"), "not in the colData")
  expect_error(libInfo(counts, bamFiles = sampleSheet$bam[1:2]), "does not match")
  expect_error(libInfo(data.frame()), "must be a RegionSetDE.counts")
})


test_that("pairwiseContrasts writes every pair, or every level against a reference", {

  allPairs <- pairwiseContrasts(counts, column = "condition")

  expect_length(allPairs, 3L)
  expect_identical(names(allPairs), c("R1881_4h_vs_DMSO", "R1881_24h_vs_DMSO", "R1881_24h_vs_R1881_4h"))
  expect_identical(allPairs$R1881_24h_vs_DMSO, c("condition", "R1881_24h", "DMSO"))

  againstVehicle <- pairwiseContrasts(sampleSheet, column = "condition", reference = "DMSO")
  expect_identical(names(againstVehicle), c("R1881_4h_vs_DMSO", "R1881_24h_vs_DMSO"))

  # The order of the levels sets the sign of the fold change, so it can be turned around
  reversed <- pairwiseContrasts(counts, column = "condition", levels = c("R1881_24h", "DMSO"))
  expect_identical(reversed, list(DMSO_vs_R1881_24h = c("condition", "DMSO", "R1881_24h")))

  expect_error(pairwiseContrasts(counts, column = "tissue"), "not in the sample table")
  expect_error(pairwiseContrasts(counts, column = "condition", reference = "EtOH"), "not among the levels")
  expect_error(pairwiseContrasts(counts, column = "condition", levels = "DMSO"), "two levels")
})


test_that("the contrasts written out run through testRegions as a list", {

  fit <- fitRegions(normalizeCounts(counts, method = "librarySize", verbose = FALSE),
                    design = ~ condition, engine = "edgeR", verbose = FALSE)
  results <- testRegions(fit, contrast = pairwiseContrasts(fit, column = "condition"), verbose = FALSE)

  expect_s4_class(results, "RegionSetDE.resultsList")
  expect_identical(names(results), names(pairwiseContrasts(fit, column = "condition")))
  expect_s4_class(results$R1881_24h_vs_DMSO, "RegionSetDE.results")

  # The stacked table is filtered with the same string every other function takes
  stackedTable <- resultsTable(results)
  expect_identical(colnames(stackedTable)[1:2], c("contrast", "contrast.description"))
  expect_setequal(unique(stackedTable$contrast), names(results))
  expect_identical(nrow(dplyr::filter(stackedTable, .data$contrast == "R1881_24h_vs_DMSO")),
                   nrow(resultsTable(results$R1881_24h_vs_DMSO)))
  expect_identical(unique(stackedTable$contrast.description[stackedTable$contrast == "R1881_24h_vs_DMSO"]),
                   contrastName(results$R1881_24h_vs_DMSO))
})


test_that("single-end and paired-end libraries are counted together on one scale", {

  # The first mates of a paired-end library, written back as single-end records
  pairedFile <- sampleSheet$bam[sampleSheet$sample == "AR_R1881_24h_r1"]
  readFields <- Rsamtools::scanBam(RegionSetDE:::.bamWithIndex(pairedFile),
                                   param = Rsamtools::ScanBamParam(what = c("qname", "flag", "rname", "strand", "pos", "mapq", "cigar"),
                                                                   flag = Rsamtools::scanBamFlag(isFirstMateRead = TRUE)))[[1]]

  singleFlag <- ifelse(readFields$strand == "-", 16L, 0L) + ifelse(bitwAnd(readFields$flag, 1024L) > 0, 1024L, 0L)
  chromosomeLengths <- Rsamtools::scanBamHeader(pairedFile)[[1]]$targets

  samFile <- tempfile(fileext = ".sam")
  writeLines(c("@HD\tVN:1.6\tSO:coordinate",
               paste0("@SQ\tSN:", names(chromosomeLengths), "\tLN:", chromosomeLengths),
               paste(readFields$qname, singleFlag, readFields$rname, readFields$pos, readFields$mapq, readFields$cigar,
                     "*", 0, 0, "*", "*", sep = "\t")),
             samFile)
  # asBam prints the number of records it wrote
  utils::capture.output(singleFile <- Rsamtools::asBam(samFile, sub("\\.sam$", "", samFile), overwrite = TRUE, indexDestination = TRUE))

  # The layout is read from the files, one per sample
  mixedCounts <- countReads(consensus, bamFiles = c(pairedFile, singleFile), sampleNames = c("paired", "single"), verbose = FALSE)
  expect_identical(SummarizedExperiment::colData(mixedCounts)$paired.end, c(TRUE, FALSE))

  # One fragment, one count, whichever way it was sequenced
  mixedInfo <- libInfo(mixedCounts)
  expect_identical(mixedInfo$library.size[1], mixedInfo$library.size[2])
  expect_lt(abs(mixedInfo$reads.in.regions[2] / mixedInfo$reads.in.regions[1] - 1), 0.05)

  countMatrix <- countTable(mixedCounts, format = "matrix")
  expect_gt(stats::cor(countMatrix[, 1], countMatrix[, 2]), 0.95)
})


test_that("sampleInfo returns the sample table of counts, fits and results alike", {

  sampleTable <- sampleInfo(counts)
  expectedTable <- as.data.frame(SummarizedExperiment::colData(counts), optional = TRUE)
  rownames(expectedTable) <- NULL

  expect_s3_class(sampleTable, "data.frame")
  expect_identical(sampleTable$sample, colnames(counts))
  expect_identical(sampleTable, expectedTable)

  expect_identical(colnames(sampleInfo(counts, columns = c("condition", "sample"))), c("condition", "sample"))
  expect_error(sampleInfo(counts, columns = "tissue"), "not in the sample table")

  # The normalisation adds its factors, and the fit and the results reach the same table
  normalisedCounts <- normalizeCounts(counts, method = "librarySize", verbose = FALSE)
  expect_true("scaling.factor" %in% colnames(sampleInfo(normalisedCounts)))

  fit <- fitRegions(normalisedCounts, design = ~ condition, engine = "edgeR", verbose = FALSE)
  results <- testRegions(fit, contrast = c("condition", "R1881_24h", "DMSO"), verbose = FALSE)

  expect_identical(sampleInfo(fit), sampleInfo(normalisedCounts))
  expect_identical(sampleInfo(results), sampleInfo(normalisedCounts))

  expect_error(sampleInfo(testRegions(fit, contrast = c("condition", "R1881_24h", "DMSO"), carryCounts = FALSE, verbose = FALSE)),
               "carries no counts")
  expect_error(sampleInfo(data.frame()), "must be a RegionSetDE.counts")
})
