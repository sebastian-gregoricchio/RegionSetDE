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
consensusSheet <- function(withSingleGroup = FALSE,
                           bamFiles = NULL) {
  sharedStarts <- seq(10000L, 100000L, by = 10000L)
  aStarts <- c(sharedStarts, 200000L, 210000L)
  bStarts <- c(sharedStarts, 300000L, 310000L)

  sheetTable <- data.frame(sample = c("A_1", "A_2", "A_3", "B_1", "B_2"),
                           bam = if (is.null(bamFiles)) {paste0(c("A_1", "A_2", "A_3", "B_1", "B_2"), ".bam")} else {bamFiles},
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

  return(loadSampleSheet(sheetTable, checkFiles = !is.null(bamFiles), verbose = FALSE))
}


# One library per sample on the contig the peaks sit on, the A samples piling reads over the first shared peak
consensusLibraries <- function() {
  sampleNames <- c("A_1", "A_2", "A_3", "B_1", "B_2")

  return(vapply(seq_along(sampleNames),
                function(i) {
                  syntheticInput(paste0("consensus_", sampleNames[i]),
                                 pileStart = if (i <= 3) {10000L} else {NULL},
                                 seed = i,
                                 chromosome = "chr1",
                                 contigLength = 400000L,
                                 readCount = 8000L)
                },
                character(1)))
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


test_that("the counting takes the files and the annotation from the sheet of the consensus", {

  skip_if_not_installed("consensusRegions")
  sampleSheet <- consensusSheet(bamFiles = consensusLibraries())
  regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)

  expect_message(counts <- countReads(regions, pairedEnd = FALSE),
                 "sample sheet the consensus was built from")

  expect_identical(dim(counts), c(14L, 5L))
  expect_identical(colnames(counts), sampleSheet$sample)

  sampleTable <- SummarizedExperiment::colData(counts)
  expect_identical(as.character(sampleTable$condition), sampleSheet$condition)
  expect_identical(as.character(sampleTable$peaks), sampleSheet$peaks)
  expect_true(all(sampleTable$library.size > 0))

  # The pile sits on the first shared peak, which the A libraries therefore cover far better
  countMatrix <- SummarizedExperiment::assay(counts, "counts")
  pileRow <- which(GenomicRanges::start(SummarizedExperiment::rowRanges(counts)) < 11000)
  expect_gt(min(countMatrix[pileRow, 1:3]), max(countMatrix[pileRow, 4:5]))
})


test_that("the peaks can be counted per tile as well as per region", {

  skip_if_not_installed("consensusRegions")
  sampleSheet <- consensusSheet(bamFiles = consensusLibraries())
  regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)

  tiledCounts <- countReads(regions, sampleSheet = sampleSheet, tileWidth = 100,
                            pairedEnd = FALSE, verbose = FALSE)

  expect_identical(tiledCounts@counting.level, "tile")
  expect_gt(nrow(tiledCounts), 14L)

  # The occupancy of the consensus follows the tiles, and countTable puts the tiles back together
  tileTable <- countTable(tiledCounts, level = "tile", verbose = FALSE)
  regionTable <- countTable(tiledCounts, level = "region", verbose = FALSE)

  expect_true(all(c("peak.A", "peak.B", "peak.groups", "peak.samples") %in% colnames(tileTable)))
  expect_identical(nrow(regionTable), 14L)
  expect_identical(sum(regionTable$peak.groups == 2), 10L)
})


test_that("a tiled object is tested tile by tile and combined back to its regions", {

  skip_if_not_installed("consensusRegions")
  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)

  tiledCounts <- countReads(regions, sampleSheet = sampleSheet, tileWidth = 100, verbose = FALSE)
  tiledFit <- fitRegions(normalizeCounts(tiledCounts, method = "librarySize", verbose = FALSE),
                         design = ~ condition, engine = "edgeR", verbose = FALSE)

  # Both procedures csaw offers, one function each
  for (combineMethod in c("simes", "holm-min")) {
    tiledResults <- testRegions(tiledFit, contrast = c("condition", "R1881_24h", "DMSO"),
                                combineMethod = combineMethod, verbose = FALSE)

    regionTable <- resultsTable(tiledResults)
    expect_identical(nrow(regionTable), length(regions@regions$consensus))
    expect_true(all(c("n.tiles", "n.tiles.up", "n.tiles.down") %in% colnames(regionTable)))
    expect_identical(sum(regionTable$n.tiles), nrow(tiledCounts))
    expect_identical(nrow(tileTable(tiledResults)), nrow(tiledCounts))
  }

  # A region carried by a strong gain has tiles counted as moving up. csaw counts them against a
  # false discovery rate within the region, and handing it the fold change threshold of the test
  # instead, zero by default, left every count at zero
  simesTable <- resultsTable(testRegions(tiledFit, contrast = c("condition", "R1881_24h", "DMSO"), verbose = FALSE))
  strongGains <- simesTable$FDR < 1e-5 & simesTable$log2FC > 0
  expect_gt(sum(strongGains), 0L)
  expect_true(all(simesTable$n.tiles.up[strongGains] > 0))
  expect_true(all(simesTable$n.tiles.down[strongGains] == 0))

  # holm-min asks for several tiles to agree, so it can only be as permissive as Simes or stricter
  simesTable <- resultsTable(testRegions(tiledFit, contrast = c("condition", "R1881_24h", "DMSO"), combineMethod = "simes", verbose = FALSE))
  holmTable <- resultsTable(testRegions(tiledFit, contrast = c("condition", "R1881_24h", "DMSO"), combineMethod = "holm-min", verbose = FALSE))
  expect_true(all(holmTable$p.value >= simesTable$p.value - 1e-12))

  expect_error(testRegions(tiledFit, contrast = c("condition", "R1881_24h", "DMSO"), combineMethod = "stouffer", verbose = FALSE),
               "either 'simes' or 'holm-min'")
})


test_that("the peak counts are normalised like any other region set", {

  skip_if_not_installed("consensusRegions")
  sampleSheet <- consensusSheet(bamFiles = consensusLibraries())
  regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)

  counts <- countReads(regions, sampleSheet = sampleSheet, pairedEnd = FALSE, verbose = FALSE)
  counts <- countBackground(counts, binSize = 20000, verbose = FALSE)

  backgroundNormalised <- normalizeCounts(counts, method = "background", verbose = FALSE)

  # A handful of peaks is a thin basis for TMM, and the function says so
  expect_warning(peakNormalised <- normalizeCounts(counts, method = "TMM", verbose = FALSE),
                 "consider 'background'")

  for (normalisedCounts in list(backgroundNormalised, peakNormalised)) {
    scalingFactors <- SummarizedExperiment::colData(normalisedCounts)$scaling.factor
    expect_length(scalingFactors, 5L)
    expect_true(all(is.finite(scalingFactors) & scalingFactors > 0))
  }

  # The normalised table divides the counts by the factors of the samples, occupancy columns included
  normalisedTable <- countTable(backgroundNormalised, normalized = TRUE, verbose = FALSE)
  rawMatrix <- countTable(backgroundNormalised, format = "matrix", verbose = FALSE)
  normalisedMatrix <- countTable(backgroundNormalised, normalized = TRUE, format = "matrix", verbose = FALSE)

  expect_true(all(c("peak.A", "peak.groups") %in% colnames(normalisedTable)))
  expect_equal(normalisedMatrix,
               sweep(rawMatrix, MARGIN = 2, STATS = backgroundNormalised$scaling.factor, FUN = "/"))

  # Five libraries in two groups are enough for the sample level plots
  samplePCA <- computeSamplePCA(backgroundNormalised, verbose = FALSE)
  expect_identical(samplePCA$scores$sample, sampleSheet$sample)
  expect_identical(as.character(samplePCA$scores$condition), sampleSheet$condition)

  sampleCorrelation <- computeSampleCorrelation(backgroundNormalised, verbose = FALSE)
  expect_identical(dim(sampleCorrelation$correlation), c(5L, 5L))
})


test_that("the counting refuses two sources of files at once", {

  skip_if_not_installed("consensusRegions")
  bamFiles <- consensusLibraries()
  sampleSheet <- consensusSheet(bamFiles = bamFiles)
  regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)

  expect_error(countReads(regions, bamFiles = bamFiles, sampleSheet = sampleSheet, verbose = FALSE),
               "not both")
  expect_error(countReads(regions, sampleSheet = sampleSheet, sampleNames = sampleSheet$sample, verbose = FALSE),
               "leave 'sampleNames'")
  expect_error(countReads(regions, sampleSheet = data.frame(sample = c("A_1", "A_2"),
                                                            bam = c(bamFiles[1], NA)), verbose = FALSE),
               "have no bam file")
  expect_error(countBigwig(regions, sampleSheet = sampleSheet, verbose = FALSE),
               "'sample' and 'bigwig' columns")

  # Regions that were not built from peaks have no sheet to fall back on
  expect_error(countReads(toyRegionSet(), verbose = FALSE), "No bam file was given")
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


# The whole peak workflow, built once and reused by the tests below
peakAnalysisCache <- new.env(parent = emptyenv())

peakAnalysis <- function() {
  if (is.null(peakAnalysisCache$results)) {
    sampleSheet <- consensusSheet(bamFiles = consensusLibraries())
    regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)

    counts <- countReads(regions, sampleSheet = sampleSheet, pairedEnd = FALSE, verbose = FALSE)
    counts <- countBackground(counts, binSize = 20000, verbose = FALSE)
    counts <- normalizeCounts(counts, method = "background", verbose = FALSE)

    peakFit <- fitRegions(counts, design = ~ condition, engine = "edgeR", verbose = FALSE)
    peakAnalysisCache$results <- testRegions(peakFit, contrast = c("condition", "A", "B"),
                                             FDR = 0.05, verbose = FALSE)
  }

  return(peakAnalysisCache$results)
}


test_that("the peaks go through the fit and the test carrying their occupancy", {

  skip_if_not_installed("consensusRegions")
  results <- peakAnalysis()

  resultTable <- resultsTable(results)
  expect_identical(nrow(resultTable), 14L)
  expect_true(all(c("peak.A", "peak.B", "peak.groups", "peak.samples") %in% colnames(resultTable)))
  expect_identical(results@contrast.groups$groups, c("A", "B"))

  # The pile sits on the first shared peak in the A libraries, and it is the region that comes out
  changedRegions <- dplyr::filter(resultTable, .data$diff.status != "null")
  expect_identical(nrow(changedRegions), 1L)
  expect_lt(changedRegions$start[1], 11000)
  expect_gt(changedRegions$log2FC[1], 2)

  # The occupancy columns follow the results out to the files
  outputDirectory <- file.path(tempdir(), "peakExport")
  exportResults(results, path = outputDirectory, prefix = "peaks", verbose = FALSE)

  writtenTable <- utils::read.delim(file.path(outputDirectory, "peaks_regions.tsv.gz"), stringsAsFactors = FALSE)
  expect_true(all(c("peak.A", "peak.B") %in% colnames(writtenTable)))
})


test_that("the occupancy table crosses the consensus with the direction", {

  skip_if_not_installed("consensusRegions")
  occupancyCounts <- peakOccupancyTable(peakAnalysis())

  expect_s3_class(occupancyCounts, "data.frame")
  expect_identical(as.character(occupancyCounts$occupancy), c("shared", "A only", "B only"))
  expect_identical(occupancyCounts$n.regions, c(10L, 2L, 2L))
  expect_identical(occupancyCounts$up + occupancyCounts$down + occupancyCounts$null, occupancyCounts$n.regions)
  expect_identical(sum(occupancyCounts$up), 1L)
  expect_equal(occupancyCounts$percent.changed[1], 10)

  # The number of samples carrying a peak is the other reading of the same regions
  sampleCounts <- peakOccupancyTable(peakAnalysis(), by = "samples")
  expect_identical(sum(sampleCounts$n.regions), 14L)
  expect_true(all(as.integer(as.character(sampleCounts$occupancy)) %in% 1:5))

  # A single group leaves the regions it covers and the ones it does not
  singleGroup <- peakOccupancyTable(peakAnalysis(), groups = "A")
  expect_identical(as.character(singleGroup$occupancy), c("A", "none"))

  expect_error(peakOccupancyTable(peakAnalysis(), by = "replicate"), "'group' or 'samples'")
  expect_error(peakOccupancyTable(peakAnalysis(), groups = "C"), "no occupancy column")
  expect_error(peakOccupancyTable(exampleResults()), "not built by loadConsensusPeaks")
})


test_that("the occupancy plot is drawn from that table", {

  skip_if_not_installed("consensusRegions")
  occupancyPlot <- plotPeakOccupancy(peakAnalysis())

  expect_s3_class(occupancyPlot, "ggplot")
  expect_s3_class(plotPeakOccupancy(peakAnalysis(), proportion = TRUE, showCounts = FALSE), "ggplot")
  expect_identical(plotPeakOccupancy(peakAnalysis(), returnData = TRUE), peakOccupancyTable(peakAnalysis()))

  # A threshold given here relabels the regions without touching the p-values
  expect_gt(sum(peakOccupancyTable(peakAnalysis(), FDR = 1)$up + peakOccupancyTable(peakAnalysis(), FDR = 1)$down),
            sum(peakOccupancyTable(peakAnalysis())$up + peakOccupancyTable(peakAnalysis())$down))
})


test_that("the androgen receptor example builds a consensus of its own", {

  skip_if_not_installed("consensusRegions")
  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)

  expect_identical(nrow(sampleSheet), 9L)
  expect_true(all(file.exists(sampleSheet$peaks)))
  expect_true(all(file.exists(sampleSheet$bam)))
  expect_identical(sort(unique(sampleSheet$condition)), c("DMSO", "R1881_24h", "R1881_4h"))

  regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)
  consensusRanges <- regions@regions$consensus

  expect_gt(length(consensusRanges), 80L)
  expect_identical(unique(as.character(GenomeInfoDb::seqnames(consensusRanges))), "19")
  expect_true(all(c("peak.DMSO", "peak.R1881_4h", "peak.R1881_24h") %in% colnames(S4Vectors::mcols(consensusRanges))))

  # The androgen brings the receptor to chromatin, so the three groups differ in how much they hold
  expect_lt(sum(consensusRanges$peak.DMSO), sum(consensusRanges$peak.R1881_24h))
  expect_gt(sum(consensusRanges$peak.groups == 1), 5L)

  # Every region of a consensus comes from the peaks of at least one sample, a zero here means the
  # regions and the peaks ended up in two different naming styles
  expect_true(all(consensusRanges$peak.samples > 0))
  expect_lte(max(consensusRanges$peak.samples), nrow(sampleSheet))
})


test_that("the naming style asked for reaches the peaks and the consensus alike", {

  skip_if_not_installed("consensusRegions")
  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)

  for (style in list("Ensembl", "UCSC")) {
    regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = style, verbose = FALSE)
    consensusInfo <- consensusData(regions)

    expect_identical(GenomeInfoDb::seqlevels(consensusInfo$total),
                     GenomeInfoDb::seqlevels(consensusInfo$peaks[[1]]))
  }

  # NULL leaves the files alone, and the files are written the Ensembl way
  untouched <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = NULL, verbose = FALSE)
  expect_identical(GenomeInfoDb::seqlevels(consensusData(untouched)$total), "19")
})


test_that("the example alignments are counted through their CSI index", {

  skip_if_not_installed("consensusRegions")
  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)

  # The shipped files carry a CSI and no BAI, which Rsamtools does not find on its own
  expect_true(all(grepl("\\.csi$", vapply(sampleSheet$bam, RegionSetDE:::.bamIndexPath, character(1), USE.NAMES = FALSE))))
  expect_true(all(RegionSetDE:::.hasBamIndex(sampleSheet$bam)))

  regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)
  counts <- countReads(regions, sampleSheet = sampleSheet, pairedEnd = TRUE, verbose = FALSE)
  countMatrix <- countTable(counts, format = "matrix")

  expect_identical(ncol(countMatrix), 9L)
  expect_true(all(colSums(countMatrix) > 0))

  # The share of the library sitting in the regions has to follow the treatment, or the reads and
  # the peaks are not describing the same experiment
  frip <- colSums(countMatrix) / SummarizedExperiment::colData(counts)$library.size
  expect_lt(max(frip[sampleSheet$condition == "DMSO"]), min(frip[sampleSheet$condition == "R1881_24h"]))
})
