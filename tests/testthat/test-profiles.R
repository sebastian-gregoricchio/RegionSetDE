# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

testthat::skip_if_not_installed("consensusRegions")

sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)
counts <- normalizeCounts(countReads(consensus, summits = 200, countInput = FALSE, verbose = FALSE), method = "TMM", verbose = FALSE)
results <- testRegions(fitRegions(counts, design = ~ condition, verbose = FALSE),
                       contrast = c("condition", "R1881_24h", "DMSO"), verbose = FALSE)


test_that("computeProfiles bins the signal around the differential regions", {

  profileData <- computeProfiles(results, groupBy = "condition", distance = 1000, binWidth = 50, verbose = FALSE)

  expect_named(profileData, c("profiles", "regions", "bins", "samples", "parameters"))
  expect_identical(names(profileData$profiles), c("DMSO", "R1881_4h", "R1881_24h"))
  expect_equal(profileData$bins, seq(-975, 975, by = 50))

  statusTable <- table(resultsTable(results)$diff.status)
  expect_identical(as.integer(table(profileData$regions$row.group)[c("up", "down")]), as.integer(statusTable[c("up", "down")]))
  expect_identical(dim(profileData$profiles$DMSO), c(length(profileData$regions), 40L))
  expect_true(all(BiocGenerics::width(profileData$regions) == 2000L))

  # The regions going up are aligned on their summit, where the treated samples pile up
  upRows <- profileData$regions$row.group == "up"
  treatedProfile <- colMeans(profileData$profiles$R1881_24h[upRows, , drop = FALSE])
  expect_gt(mean(treatedProfile[20:21]), 5 * mean(treatedProfile[c(1:3, 38:40)]))
  expect_gt(mean(treatedProfile[20:21]), mean(colMeans(profileData$profiles$DMSO[upRows, , drop = FALSE])[20:21]))
})


test_that("the BAM coverage is scaled by the normalisation, a bigWig is read as it is", {

  bamProfiles <- computeProfiles(counts, distance = 500, useOffsets = FALSE, verbose = FALSE)
  normalisedProfiles <- computeProfiles(counts, distance = 500, useOffsets = TRUE, verbose = FALSE)

  # One factor per sample, so the two differ by a constant for every sample
  profileRatio <- normalisedProfiles$profiles[[1]] / bamProfiles$profiles[[1]]
  expect_lt(stats::sd(profileRatio[is.finite(profileRatio)]), 1e-8)

  # A bigWig written from the same fragments gives the same shape
  chromosomeLengths <- Rsamtools::scanBamHeader(sampleSheet$bam[1])[[1]]$targets
  bigwigFiles <- vapply(seq_len(nrow(sampleSheet)), function(i) {
    sampleCoverage <- RegionSetDE:::.windowCoverage(bamFile = sampleSheet$bam[i],
                                                    windows = GenomicRanges::GRanges(names(chromosomeLengths)[1], IRanges::IRanges(1, chromosomeLengths[1])),
                                                    chromosomeLengths = chromosomeLengths, isPairedEnd = TRUE, fragmentLength = NA,
                                                    maxFragmentLength = 1000, minMapq = 20, removeDuplicates = TRUE)$coverage
    bigwigFile <- file.path(tempdir(), paste0(sampleSheet$sample[i], ".bw"))
    rtracklayer::export.bw(sampleCoverage, bigwigFile)
    return(bigwigFile)
  }, character(1))

  bigwigProfiles <- computeProfiles(counts, distance = 500, signalFiles = bigwigFiles, verbose = FALSE)
  expect_identical(bigwigProfiles$parameters$signal, "bigwig")
  expect_gt(stats::cor(as.numeric(bigwigProfiles$profiles[[9]]), as.numeric(bamProfiles$profiles[[9]])), 0.999)
})


test_that("regions on the minus strand are read from right to left", {

  centreRanges <- SummarizedExperiment::rowRanges(counts)[1:5]
  plusRanges <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(centreRanges), IRanges::ranges(centreRanges), strand = "+")
  minusRanges <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(centreRanges), IRanges::ranges(centreRanges), strand = "-")

  plusProfiles <- computeProfiles(counts, regions = list(plus = plusRanges), distance = 500, verbose = FALSE)
  minusProfiles <- computeProfiles(counts, regions = list(minus = minusRanges), distance = 500, verbose = FALSE)

  expect_equal(minusProfiles$profiles[[1]], plusProfiles$profiles[[1]][, rev(seq_len(20))], ignore_attr = TRUE)
})


test_that("plotProfile draws heatmaps or lines", {

  heatmapPlot <- plotProfile(results, groupBy = "condition", distance = 500, verbose = FALSE)
  expect_s4_class(heatmapPlot, "HeatmapList")
  expect_length(heatmapPlot@ht_list, 3)

  profileData <- computeProfiles(results, distance = 500, verbose = FALSE)
  linePlot <- plotProfile(profileData, style = "lines")
  expect_s3_class(linePlot, "ggplot")
  expect_setequal(unique(as.character(linePlot$data$row.group)), c("up", "down"))

  expect_error(plotProfile(profileData, style = "bars"), "'heatmap' or 'lines'")
  expect_error(computeProfiles(results, direction = "sideways", verbose = FALSE), "'both', 'up', 'down'")
  expect_error(computeProfiles(results, FDR = 1e-30, verbose = FALSE), "No region passes")
})


test_that("computeProfiles draws the samples asked for", {

  profileData <- computeProfiles(results, samples = c("AR_DMSO_r1", "AR_R1881_24h_r1"), distance = 500, verbose = FALSE)
  expect_identical(names(profileData$profiles), c("AR_DMSO_r1", "AR_R1881_24h_r1"))

  # The same sample gives the same profile whether the others are there or not
  allProfiles <- computeProfiles(results, distance = 500, verbose = FALSE)
  expect_equal(profileData$profiles$AR_R1881_24h_r1, allProfiles$profiles$AR_R1881_24h_r1)
})


test_that("a blacklist takes out the rows whose window it touches, a whitelist keeps the regions it covers", {

  allProfiles <- computeProfiles(results, distance = 500, verbose = FALSE)
  profileRegions <- allProfiles$regions

  # A blacklist sitting just outside the first region, but inside its window
  firstWindow <- profileRegions[1]
  nearRegion <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(firstWindow), IRanges::IRanges(BiocGenerics::end(firstWindow) - 10, width = 5))

  blacklistedProfiles <- computeProfiles(results, distance = 500, blacklist = nearRegion, verbose = FALSE)
  expect_false(profileRegions$region.key[1] %in% blacklistedProfiles$regions$region.key)
  expect_identical(length(blacklistedProfiles$regions), length(profileRegions) - 1L)

  # The object carries no blacklist here
  expect_error(computeProfiles(results, distance = 500, blacklist = TRUE, verbose = FALSE), "there is none")

  # A whitelist on the first two regions keeps those two
  whitelist <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(profileRegions[1:2]),
                                      IRanges::IRanges(start = profileRegions$centre.position[1:2], width = 1))
  whitelistedProfiles <- computeProfiles(results, distance = 500, whitelist = whitelist, verbose = FALSE)
  expect_setequal(whitelistedProfiles$regions$region.key, profileRegions$region.key[1:2])

  expect_error(computeProfiles(results, distance = 500, whitelist = GenomicRanges::GRanges("19", IRanges::IRanges(1, 10)), verbose = FALSE),
               "No region is left")
})
