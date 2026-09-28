# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

test_that("countReads recentres the regions on the summit of the reads", {

  bamFile <- syntheticFragments("summit_sites", fragmentLength = 200L)
  sitePositions <- attr(bamFile, "sites")
  siteRegions <- syntheticSiteRegions(sitePositions, shift = 700L)

  counts <- countReads(siteRegions, bamFiles = bamFile, fragmentLength = 200, summits = 150, verbose = FALSE)
  rowRanges <- SummarizedExperiment::rowRanges(counts)

  expect_true(all(BiocGenerics::width(rowRanges) == 301L))
  expect_true(all(abs(rowRanges$summit - sitePositions) < 30))
  expect_equal(BiocGenerics::start(rowRanges) + 150L, rowRanges$summit)

  # The identifiers still point at the regions the windows came from
  expect_identical(rowRanges$region.id,
                   paste0("chrT:", sitePositions + 700L - 2000L, "-", sitePositions + 700L - 2000L + 3999L))
  expect_identical(counts@parameters$countReads$summits, 150)
  expect_identical(counts@parameters$countReads$summitSource, "reads")
})


test_that("summits = 0 locates the summits without moving the regions", {

  bamFile <- syntheticFragments("summit_zero", fragmentLength = 200L)
  siteRegions <- syntheticSiteRegions(attr(bamFile, "sites"), shift = 700L)

  counts <- countReads(siteRegions, bamFiles = bamFile, fragmentLength = 200, summits = 0, verbose = FALSE)
  plainCounts <- countReads(siteRegions, bamFiles = bamFile, fragmentLength = 200, verbose = FALSE)

  expect_identical(BiocGenerics::start(SummarizedExperiment::rowRanges(counts)), BiocGenerics::start(SummarizedExperiment::rowRanges(plainCounts)))
  expect_equal(SummarizedExperiment::assay(counts, "counts"), SummarizedExperiment::assay(plainCounts, "counts"))
  expect_true("summit" %in% colnames(SummarizedExperiment::rowData(counts)))
  expect_false("summit" %in% colnames(SummarizedExperiment::rowData(plainCounts)))
})


test_that("the summits of the peak calls and of the reads agree on a consensus", {

  testthat::skip_if_not_installed("consensusRegions")

  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)
  consensusRanges <- regionRanges(consensus)$consensus

  readSummits <- countReads(consensus, summits = 200, summitSource = "reads", countInput = FALSE, verbose = FALSE)
  peakSummits <- countReads(consensus, summits = 200, summitSource = "peaks", countInput = FALSE, verbose = FALSE)

  readPosition <- SummarizedExperiment::rowData(readSummits)$summit
  peakPosition <- SummarizedExperiment::rowData(peakSummits)$summit

  # Every summit lies within the region it was found for
  expect_true(all(peakPosition >= BiocGenerics::start(consensusRanges) & peakPosition <= BiocGenerics::end(consensusRanges)))
  expect_true(all(readPosition >= BiocGenerics::start(consensusRanges) & readPosition <= BiocGenerics::end(consensusRanges)))
  expect_lt(stats::median(abs(readPosition - peakPosition)), 100)

  # The inputs are counted over the recentred windows
  withInput <- countReads(consensus, summits = 200, verbose = FALSE)
  expect_true(all(BiocGenerics::width(SummarizedExperiment::rowRanges(withInput)) == 401L))
  expect_true("input" %in% SummarizedExperiment::assayNames(withInput))
})


test_that("the summit options are checked", {

  bamFile <- syntheticFragments("summit_checks", fragmentLength = 200L)
  siteRegions <- syntheticSiteRegions(attr(bamFile, "sites"))

  expect_error(countReads(siteRegions, bamFiles = bamFile, summits = 100, tileWidth = 50, verbose = FALSE), "cannot be used together")
  expect_error(countReads(siteRegions, bamFiles = bamFile, summits = 100, summitSource = "peaks", verbose = FALSE), "loadConsensusPeaks")
  expect_error(countReads(siteRegions, bamFiles = bamFile, summits = -1, verbose = FALSE), "positive number")
  expect_error(countReads(siteRegions, bamFiles = bamFile, summits = 100, summitSource = "caller", verbose = FALSE), "'reads' or 'peaks'")
})
