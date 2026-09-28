# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

test_that("estimateFragmentLength recovers the fragment length of single-end libraries", {

  shortFragments <- syntheticFragments("fragments_200", fragmentLength = 200L, readLength = 50L)
  longFragments <- syntheticFragments("fragments_300", fragmentLength = 300L, readLength = 36L, seed = 2L)

  fragmentEstimate <- estimateFragmentLength(c(shortFragments, longFragments), sampleNames = c("short", "long"), verbose = FALSE)

  expect_named(fragmentEstimate, c("table", "profile", "plot"))
  expect_identical(fragmentEstimate$table$sample, c("short", "long"))
  expect_identical(fragmentEstimate$table$method, rep("cross-correlation", 2))
  expect_equal(fragmentEstimate$table$read.length, c(50, 36))
  expect_lt(abs(fragmentEstimate$table$fragment.length[1] - 200), 10)
  expect_lt(abs(fragmentEstimate$table$fragment.length[2] - 300), 10)

  expect_setequal(unique(as.character(fragmentEstimate$profile$sample)), c("short", "long"))
  expect_s3_class(fragmentEstimate$plot, "ggplot")
})


test_that("estimateFragmentLength looks past the read length and within the regions", {

  bamFile <- syntheticFragments("fragments_regions", fragmentLength = 180L, readLength = 50L)
  siteRegions <- syntheticSiteRegions(attr(bamFile, "sites"), shift = 0L)

  withRegions <- estimateFragmentLength(bamFile, regions = siteRegions, verbose = FALSE)
  expect_lt(abs(withRegions$table$fragment.length - 180), 10)

  # The search starts past the read length, where the mappability peak sits
  expect_gt(withRegions$table$fragment.length, 50 + 15)
  expect_error(estimateFragmentLength(bamFile, maxDistance = 10, verbose = FALSE), "at least 50")
})


test_that("paired-end libraries report their insert size", {

  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  insertSizes <- estimateFragmentLength(sampleSheet = sampleSheet[1:2, ], verbose = FALSE)

  expect_identical(insertSizes$table$method, rep("insert size", 2))
  expect_true(all(insertSizes$table$paired.end))
  expect_true(all(insertSizes$table$fragment.length > 100 & insertSizes$table$fragment.length < 400))
  expect_null(insertSizes$plot)
})


test_that("countReads takes the fragment length from a column of the sample sheet", {

  firstFile <- syntheticFragments("column_a", fragmentLength = 200L)
  secondFile <- syntheticFragments("column_b", fragmentLength = 200L, seed = 3L)
  siteRegions <- syntheticSiteRegions(attr(firstFile, "sites"))

  sampleTable <- data.frame(sample = c("a", "b"), phantom = c(120, 280))

  counts <- countReads(siteRegions, bamFiles = c(firstFile, secondFile), sampleNames = c("a", "b"),
                       sampleMetadata = sampleTable, fragmentLength = "phantom", verbose = FALSE)

  expect_equal(SummarizedExperiment::colData(counts)$fragment.length, c(120, 280))
  expect_equal(counts@parameters$countReads$fragmentLength, c(120, 280))
  expect_identical(counts@parameters$countReads$fragmentLengthSource, "column phantom")

  # A named vector is matched to the samples, whatever its order
  namedCounts <- countReads(siteRegions, bamFiles = c(firstFile, secondFile), sampleNames = c("a", "b"),
                            fragmentLength = c(b = 280, a = 120), verbose = FALSE)
  expect_equal(SummarizedExperiment::assay(namedCounts, "counts"), SummarizedExperiment::assay(counts, "counts"))

  expect_error(countReads(siteRegions, bamFiles = firstFile, fragmentLength = "absentColumn", verbose = FALSE),
               "names no column")
  expect_error(countReads(siteRegions, bamFiles = c(firstFile, secondFile), sampleNames = c("a", "b"),
                          sampleMetadata = data.frame(sample = c("a", "b"), phantom = c(150, NA)),
                          fragmentLength = "phantom", verbose = FALSE),
               "No valid fragment length")
})


test_that("the length the reads are extended to changes what the regions collect", {

  bamFile <- syntheticFragments("extension", fragmentLength = 200L)

  # A narrow window right of a site catches the reads extended from its left only when they are long enough
  sitePosition <- attr(bamFile, "sites")[1]
  rightWindow <- loadRegions(list(right = GenomicRanges::GRanges("chrT", IRanges::IRanges(sitePosition + 150L, width = 50L))),
                             seqlevelsStyle = NULL, verbose = FALSE)

  shortCounts <- countReads(rightWindow, bamFiles = bamFile, fragmentLength = 60, verbose = FALSE)
  longCounts <- countReads(rightWindow, bamFiles = bamFile, fragmentLength = 300, verbose = FALSE)

  expect_gt(SummarizedExperiment::assay(longCounts, "counts")[1, 1], SummarizedExperiment::assay(shortCounts, "counts")[1, 1])
})


test_that("fragmentLength = 'auto' estimates every single-end sample and leaves the paired-end ones alone", {

  bamFile <- syntheticFragments("auto_single", fragmentLength = 220L)
  siteRegions <- syntheticSiteRegions(attr(bamFile, "sites"), shift = 0L)

  counts <- countReads(siteRegions, bamFiles = bamFile, fragmentLength = "auto", verbose = FALSE)
  expect_lt(abs(SummarizedExperiment::colData(counts)$fragment.length - 220), 10)
  expect_identical(counts@parameters$countReads$fragmentLengthSource, "cross-correlation")

  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), verbose = FALSE)
  pairedCounts <- countReads(peakRegions, sampleSheet = sampleSheet[1:2, ], fragmentLength = "auto", countInput = FALSE, verbose = FALSE)

  expect_true(all(is.na(SummarizedExperiment::colData(pairedCounts)$fragment.length)))
})
