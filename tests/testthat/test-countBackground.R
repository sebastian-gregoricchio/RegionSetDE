## countBackground needs an object counted from BAM files, since it recovers the
## files and the read filters from the parameters of countReads. The Rsamtools
## example alignment is the only one available at check time, and its contigs are
## about 1.5 kb, so the bins have to be small.

backgroundBinWidth <- 100L


test_that("countBackground stores the bins in the metadata", {

  counts <- toyCounts()
  counts <- countBackground(counts, binSize = backgroundBinWidth, verbose = FALSE)

  backgroundBins <- S4Vectors::metadata(counts)$background

  expect_false(is.null(backgroundBins))
  expect_gt(nrow(backgroundBins), 0)
  expect_equal(ncol(backgroundBins), ncol(counts))
  expect_identical(colnames(backgroundBins), colnames(counts))
})


test_that("countBackground returns the same object rather than a new one", {

  counts <- toyCounts()
  withBackground <- countBackground(counts, binSize = backgroundBinWidth, verbose = FALSE)

  expect_s4_class(withBackground, "RegionSetDE.counts")
  expect_equal(nrow(withBackground), nrow(counts))
  expect_equal(ncol(withBackground), ncol(counts))

  # The region counts are untouched by the background step
  expect_equal(SummarizedExperiment::assay(withBackground, "counts"),
               SummarizedExperiment::assay(counts, "counts"))
})


test_that("excludeRegions drops the bins overlapping the analysed regions", {

  counts <- toyCounts()

  withExclusion <- countBackground(counts, binSize = backgroundBinWidth,
                                   excludeRegions = TRUE, verbose = FALSE)
  withoutExclusion <- countBackground(counts, binSize = backgroundBinWidth,
                                      excludeRegions = FALSE, verbose = FALSE)

  expect_lte(nrow(S4Vectors::metadata(withExclusion)$background),
             nrow(S4Vectors::metadata(withoutExclusion)$background))
})


test_that("a wider bin gives fewer of them", {

  counts <- toyCounts()

  narrowBins <- countBackground(counts, binSize = 100L, excludeRegions = FALSE,
                                verbose = FALSE)
  wideBins <- countBackground(counts, binSize = 400L, excludeRegions = FALSE,
                              verbose = FALSE)

  expect_gt(nrow(S4Vectors::metadata(narrowBins)$background),
            nrow(S4Vectors::metadata(wideBins)$background))
})


test_that("minCount removes the empty bins", {

  counts <- toyCounts()

  keepingEmpty <- countBackground(counts, binSize = backgroundBinWidth,
                                  minCount = 0, excludeRegions = FALSE,
                                  verbose = FALSE)
  droppingEmpty <- countBackground(counts, binSize = backgroundBinWidth,
                                   minCount = 5, excludeRegions = FALSE,
                                   verbose = FALSE)

  expect_gte(nrow(S4Vectors::metadata(keepingEmpty)$background),
             nrow(S4Vectors::metadata(droppingEmpty)$background))
})


test_that("restrictChromosomes narrows the bins to one contig", {

  counts <- toyCounts()

  counts <- countBackground(counts, binSize = backgroundBinWidth,
                            restrictChromosomes = "seq1",
                            excludeRegions = FALSE, verbose = FALSE)

  backgroundBins <- S4Vectors::metadata(counts)$background

  expect_gt(nrow(backgroundBins), 0)
})


test_that("the bins carry counts of the expected shape", {

  counts <- toyCounts()
  counts <- countBackground(counts, binSize = backgroundBinWidth,
                            excludeRegions = FALSE, verbose = FALSE)

  backgroundBins <- S4Vectors::metadata(counts)$background

  # The bins may be a SummarizedExperiment or a rectangular S4 table
  countMatrix <- if (methods::is(backgroundBins, "SummarizedExperiment")) {
    SummarizedExperiment::assay(backgroundBins)
  } else {
    as.matrix(as.data.frame(backgroundBins))
  }

  expect_true(is.numeric(countMatrix))
  expect_true(all(countMatrix >= 0))
  expect_equal(ncol(countMatrix), ncol(counts))
})


test_that("bamFiles can be given rather than recovered", {

  counts <- toyCounts()

  counts <- countBackground(counts,
                            bamFiles = toyBamFile(),
                            binSize = backgroundBinWidth,
                            excludeRegions = FALSE,
                            verbose = FALSE)

  expect_gt(nrow(S4Vectors::metadata(counts)$background), 0)
})


test_that("countBackground refuses a file that is not there", {

  expect_error(
    countBackground(toyCounts(), bamFiles = "there/is/no/such.bam", verbose = FALSE)
  )
})


test_that("the packaged counts carry bins matching their samples", {

  counts <- exampleCounts()
  backgroundBins <- S4Vectors::metadata(counts)$background

  expect_equal(ncol(backgroundBins), ncol(counts))
  expect_identical(colnames(backgroundBins), colnames(counts))
  expect_gt(nrow(backgroundBins), 100)
})
