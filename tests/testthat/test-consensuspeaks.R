# A narrowPeak file with strong peaks at the given starts, shifted a little per replicate as real calls would be
writePeaks <- function(fileName,
                       starts,
                       chromosome = "chr1",
                       width = 500L,
                       shift = 0L) {
  peakStarts <- as.integer(starts + shift)
  peakTable <- data.frame(chromosome, peakStarts - 1L, peakStarts - 1L + width, paste0("peak_", seq_along(starts)),
                          1000L, ".", 10, 30, 29, width %/% 2L)

  peakPath <- file.path(tempdir(), fileName)
  utils::write.table(peakTable, peakPath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  return(peakPath)
}


# Three replicates of A and two of B share ten peaks, and each group has two peaks of its own
consensusSheet <- function(withSingleGroup = FALSE) {
  sharedStarts <- seq(10000L, 100000L, by = 10000L)
  aStarts <- c(sharedStarts, 200000L, 210000L)
  bStarts <- c(sharedStarts, 300000L, 310000L)

  sheetTable <- data.frame(sample = c("A_1", "A_2", "A_3", "B_1", "B_2"),
                           bam = paste0(c("A_1", "A_2", "A_3", "B_1", "B_2"), ".bam"),
                           peaks = c(writePeaks("A_1.narrowPeak", aStarts, shift = 0L),
                                     writePeaks("A_2.narrowPeak", aStarts, shift = 20L),
                                     writePeaks("A_3.narrowPeak", aStarts, shift = 40L),
                                     writePeaks("B_1.narrowPeak", bStarts, shift = 10L),
                                     writePeaks("B_2.narrowPeak", bStarts, shift = 30L)),
                           condition = c("A", "A", "A", "B", "B"))

  if (isTRUE(withSingleGroup)) {
    sheetTable <- rbind(sheetTable, data.frame(sample = "C_1", bam = "C_1.bam",
                                               peaks = writePeaks("C_1.narrowPeak", c(sharedStarts, 400000L)),
                                               condition = "C"))
  }

  return(loadSampleSheet(sheetTable, checkFiles = FALSE, verbose = FALSE))
}


test_that("the consensus is built per group and pooled into the total one", {

  skip_if_not_installed("consensusRegions")
  regions <- loadConsensusPeaks(consensusSheet(), groupBy = "condition", verbose = FALSE)

  expect_s4_class(regions, "RegionSetDE")
  expect_identical(names(regions@regions), "consensus")

  consensusRanges <- regions@regions$consensus
  expect_length(consensusRanges, 14)
  expect_true(all(c("peak.A", "peak.B", "peak.groups", "peak.samples") %in% colnames(S4Vectors::mcols(consensusRanges))))

  # The regions of one group only are kept by the union, and flagged as such
  specificA <- IRanges::overlapsAny(consensusRanges, GenomicRanges::GRanges("chr1", IRanges::IRanges(c(200000, 210000), width = 500)))
  expect_true(all(consensusRanges$peak.A[specificA]))
  expect_false(any(consensusRanges$peak.B[specificA]))
  expect_true(all(consensusRanges$peak.samples[specificA] == 3))
  expect_identical(sum(consensusRanges$peak.groups == 2), 10L)

  consensusList <- consensusData(regions)
  expect_identical(names(consensusList$groups), c("A", "B"))
  expect_identical(lengths(consensusList$groups), c(A = 12L, B = 12L))
  expect_s4_class(consensusList$objects$A, "ConsensusRegions")
  expect_identical(consensusList$mode, "consensus")
  expect_identical(consensusList$samples$group, c("A", "A", "A", "B", "B"))
})


test_that("without groups all the samples form a single one", {

  skip_if_not_installed("consensusRegions")
  regions <- loadConsensusPeaks(consensusSheet(), verbose = FALSE)

  expect_identical(names(consensusData(regions)$groups), "all")
  expect_true("peak.all" %in% colnames(S4Vectors::mcols(regions@regions$consensus)))
})


test_that("a group with a single sample keeps its peaks as they are", {

  skip_if_not_installed("consensusRegions")
  regions <- loadConsensusPeaks(consensusSheet(withSingleGroup = TRUE), groupBy = "condition", verbose = FALSE)

  consensusList <- consensusData(regions)
  expect_identical(length(consensusList$groups$C), 11L)
  expect_null(consensusList$objects$C)
  expect_length(regions@regions$consensus, 15)
})


test_that("the excluded regions take their peaks away before the consensus", {

  skip_if_not_installed("consensusRegions")
  blacklist <- GenomicRanges::GRanges("chr1", IRanges::IRanges(300000, width = 1000))
  greylist <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10000, width = 1000))

  regions <- loadConsensusPeaks(consensusSheet(), groupBy = "condition",
                                excludeRegions = list(blacklist, greylist), verbose = FALSE)

  consensusRanges <- regions@regions$consensus
  expect_length(consensusRanges, 12)
  expect_false(any(IRanges::overlapsAny(consensusRanges, c(blacklist, greylist))))

  sampleTable <- consensusData(regions)$samples
  expect_identical(sampleTable$n.excluded, c(1L, 1L, 1L, 2L, 2L))
  expect_identical(regions@parameters$loadConsensusPeaks$n.excluded.regions, 2L)
})


test_that("the regions of the user split the consensus, the first set taking the shared regions", {

  skip_if_not_installed("consensusRegions")
  firstStretch <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 50000))
  secondStretch <- GenomicRanges::GRanges("chr1", IRanges::IRanges(40000, 250000))

  regions <- loadConsensusPeaks(consensusSheet(), groupBy = "condition",
                                regionSets = list(first = firstStretch, second = secondStretch), verbose = FALSE)

  expect_identical(names(regions@regions), c("first", "second", "other"))
  expect_identical(as.integer(lengths(regions@regions)), c(5L, 7L, 2L))
  expect_identical(consensusData(regions)$mode, "split")

  droppedRegions <- loadConsensusPeaks(consensusSheet(), groupBy = "condition", unassignedSet = NULL,
                                       regionSets = list(first = firstStretch), verbose = FALSE)
  expect_identical(names(droppedRegions@regions), "first")

  expect_error(loadConsensusPeaks(consensusSheet(), groupBy = "condition", unassignedSet = "first",
                                  regionSets = list(first = firstStretch), verbose = FALSE), "already used")
})


test_that("the regions of the user can replace the consensus and be annotated by it", {

  skip_if_not_installed("consensusRegions")
  userRegions <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10000, 200000, 500000), width = 1000))

  regions <- loadConsensusPeaks(consensusSheet(), groupBy = "condition", regionMode = "replace",
                                regionSets = list(targets = userRegions), verbose = FALSE)

  targetRanges <- regions@regions$targets
  expect_length(targetRanges, 3)
  expect_identical(targetRanges$peak.groups, c(2L, 1L, 0L))
  expect_identical(targetRanges$peak.samples, c(5L, 3L, 0L))
  expect_identical(consensusData(regions)$mode, "replace")
  expect_length(consensusData(regions)$total, 14)
})


test_that("the arguments are checked", {

  skip_if_not_installed("consensusRegions")
  sampleSheet <- consensusSheet()

  expect_error(loadConsensusPeaks(sampleSheet[, c("sample", "bam")], verbose = FALSE), "'peaks' columns")
  expect_error(loadConsensusPeaks(sampleSheet, groupBy = "absent", verbose = FALSE), "groupBy")
  expect_error(loadConsensusPeaks(sampleSheet, regionMode = "merge", verbose = FALSE), "regionMode")
  expect_error(loadConsensusPeaks(sampleSheet, sampleNames = "a", verbose = FALSE), "set from the sample sheet")
  expect_error(loadConsensusPeaks(sampleSheet[1, ], verbose = FALSE), "At least two samples")
})


test_that("consensusData refuses regions that were not built from peaks", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
                              splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

  expect_error(consensusData(regions), "not built from peaks")
  expect_identical(regions@consensus, list())
})


test_that("the UpSet plot is drawn by group and by sample", {

  skip_if_not_installed("consensusRegions")
  skip_if_not_installed("ComplexHeatmap")
  regions <- loadConsensusPeaks(consensusSheet(), groupBy = "condition", verbose = FALSE)

  groupPlot <- plotPeakUpset(regions, by = "group")
  samplePlot <- plotPeakUpset(regions, by = "sample", groupColours = c(A = "darkgreen", B = "orange"))

  expect_s4_class(groupPlot, "Heatmap")
  expect_s4_class(samplePlot, "Heatmap")
  expect_identical(nrow(groupPlot@matrix), 2L)
  expect_identical(nrow(samplePlot@matrix), 5L)
  expect_identical(ncol(groupPlot@matrix), 3L)

  expect_error(plotPeakUpset(regions, by = "condition"), "'group' or 'sample'")
  expect_error(plotPeakUpset(regions, set = "absent"), "absent from the object")
  expect_error(plotPeakUpset(regions, groupColours = c(A = "red")), "no colour")
})


test_that("the UpSet plot reports the regions of the user without any peak", {

  skip_if_not_installed("consensusRegions")
  skip_if_not_installed("ComplexHeatmap")
  userRegions <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10000, 200000, 500000), width = 1000))

  regions <- loadConsensusPeaks(consensusSheet(), groupBy = "condition", regionMode = "replace",
                                regionSets = list(targets = userRegions), verbose = FALSE)

  upsetPlot <- plotPeakUpset(regions, by = "group")
  expect_match(upsetPlot@column_title, "1 without any peak")
})


test_that("the occupancy travels with the regions into the counts", {

  skip_if_not_installed("consensusRegions")
  regions <- loadConsensusPeaks(consensusSheet(), groupBy = "condition", verbose = FALSE)

  # consensusRegions keeps the standard chromosomes only, so the library is simulated on chr1
  bamFile <- syntheticInput("consensusLibrary", chromosome = "chr1", contigLength = 400000L)
  counts <- countReads(regions, bamFiles = bamFile, sampleNames = "library", pairedEnd = FALSE, verbose = FALSE)
  rowTable <- as.data.frame(SummarizedExperiment::rowData(counts))

  expect_identical(nrow(rowTable), 14L)
  expect_true(all(c("peak.A", "peak.B", "peak.groups", "peak.samples") %in% colnames(rowTable)))
  expect_identical(sum(rowTable$peak.groups == 1), 4L)
})
