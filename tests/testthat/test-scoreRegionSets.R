## scoreRegionSets covers the design that has no contrast: one condition, several
## libraries, and a question about which set carries more signal. The example
## object supplies four libraries, four disjoint sets of 1 kb windows and the
## background bins, which is everything the function reads.


test_that("scoreRegionSets returns one score per set per library", {

  counts <- exampleCounts()
  setScores <- scoreRegionSets(counts, verbose = FALSE)

  expect_s4_class(setScores, "RegionSetDE.setScores")

  scores <- scoreTable(setScores)

  expect_equal(nrow(scores), length(regionSetNames(setScores)) * ncol(counts))
  expect_true(all(c("sample", "region.set", "n.regions", "set.signal",
                    "reference.signal", "score") %in% colnames(scores)))
  expect_true(all(is.finite(scores$score)))
  expect_setequal(unique(scores$sample), colnames(counts))
  expect_setequal(regionSetNames(setScores), regionSetNames(counts))
})


test_that("every pair of sets is compared by default", {

  counts <- exampleCounts()
  setScores <- scoreRegionSets(counts, verbose = FALSE)
  comparisons <- resultsTable(setScores)

  expect_equal(nrow(comparisons), choose(length(regionSetNames(setScores)), 2))
  expect_true(all(comparisons$n.libraries == ncol(counts)))
  expect_true(all(comparisons$df == ncol(counts) - 1))
  expect_true(all(comparisons$p.value >= 0 & comparisons$p.value <= 1))
  expect_true(all(comparisons$FDR >= 0 & comparisons$FDR <= 1))
  expect_true(all(comparisons$FDR >= comparisons$p.value))
})


test_that("the reported difference is the mean of the paired differences", {

  setScores <- scoreRegionSets(exampleCounts(), verbose = FALSE)

  scores <- scoreTable(setScores)
  comparisonRow <- resultsTable(setScores)[1, ]

  firstScores <- dplyr::filter(scores, region.set == comparisonRow$set.1)
  secondScores <- dplyr::filter(scores, region.set == comparisonRow$set.2)
  pairedDifference <- firstScores$score[match(secondScores$sample, firstScores$sample)] -
    secondScores$score

  expect_equal(comparisonRow$mean.delta.score, mean(pairedDifference))
  expect_lte(comparisonRow$CI.lower, comparisonRow$mean.delta.score)
  expect_gte(comparisonRow$CI.upper, comparisonRow$mean.delta.score)
})


test_that("the comparison between two sets does not depend on the reference", {

  counts <- exampleCounts()

  orderedComparisons <- function(...) {
    comparisons <- resultsTable(scoreRegionSets(counts, verbose = FALSE, ...))
    return(dplyr::arrange(comparisons, set.1, set.2))
  }

  # Both scores of a library were divided by the same number, so it leaves the difference
  againstBackground <- orderedComparisons(reference = "background")
  againstRegions <- orderedComparisons(reference = "regions")
  withoutReference <- orderedComparisons(reference = "none")

  expect_equal(againstBackground$mean.delta.score, againstRegions$mean.delta.score)
  expect_equal(againstBackground$mean.delta.score, withoutReference$mean.delta.score)
  expect_equal(againstBackground$p.value, withoutReference$p.value)

  # The score itself is a different number under each reference, only the difference is shared
  expect_false(isTRUE(all.equal(scoreTable(scoreRegionSets(counts, reference = "background", verbose = FALSE))$score,
                                scoreTable(scoreRegionSets(counts, reference = "none", verbose = FALSE))$score)))
})


test_that("perBasepair leaves the comparison alone when the widths are equal", {

  counts <- exampleCounts()
  regionWidths <- BiocGenerics::width(SummarizedExperiment::rowRanges(counts))

  skip_if_not(length(unique(regionWidths)) == 1,
              "the example regions are not all of the same width")

  # Dividing by a width shared by every region rescales both sets by the same number
  perRegion <- resultsTable(scoreRegionSets(counts, perBasepair = FALSE, verbose = FALSE))
  perBasepair <- resultsTable(scoreRegionSets(counts, perBasepair = TRUE, verbose = FALSE))

  expect_equal(dplyr::arrange(perRegion, set.1, set.2)$mean.delta.score,
               dplyr::arrange(perBasepair, set.1, set.2)$mean.delta.score)
})


test_that("comparisons restricts the pairs that are tested", {

  counts <- exampleCounts()
  requestedPair <- regionSetNames(counts)[1:2]

  setScores <- scoreRegionSets(counts, comparisons = list(requestedPair), verbose = FALSE)
  comparisons <- resultsTable(setScores)

  expect_equal(nrow(comparisons), 1)
  expect_setequal(c(comparisons$set.1, comparisons$set.2), requestedPair)

  # The scores still cover every set, only the comparisons were restricted
  expect_setequal(regionSetNames(setScores), regionSetNames(counts))
})


test_that("regionSets restricts the sets that are scored", {

  counts <- exampleCounts()
  requestedSets <- regionSetNames(counts)[1:3]

  setScores <- scoreRegionSets(counts, regionSets = requestedSets, verbose = FALSE)

  expect_setequal(regionSetNames(setScores), requestedSets)
  expect_equal(nrow(resultsTable(setScores)), choose(length(requestedSets), 2))
})


test_that("the object carries the provenance and the samples of the counts", {

  counts <- exampleCounts()
  setScores <- scoreRegionSets(counts, summary = "median", verbose = FALSE)

  expect_identical(setScores@genome.assembly, counts@genome.assembly)
  expect_identical(setScores@seqlevels.style, counts@seqlevels.style)
  expect_true("scoreRegionSets" %in% names(setScores@parameters))
  expect_identical(setScores@summary, "median")
  expect_identical(setScores@reference, "background")
  expect_identical(setScores@assay, "counts")

  expect_equal(nrow(setScores@sample.metadata), ncol(counts))
  expect_true("condition" %in% colnames(setScores@sample.metadata))
})


test_that("scoreRegionSets refuses the arguments it cannot honour", {

  counts <- exampleCounts()

  expect_error(scoreRegionSets(exampleFit(), verbose = FALSE), "RegionSetDE.counts")
  expect_error(scoreRegionSets(counts, summary = "mode", verbose = FALSE), "'mean' or 'median'")
  expect_error(scoreRegionSets(counts, reference = "input", verbose = FALSE), "'background', 'regions' or 'none'")
  expect_error(scoreRegionSets(counts, confLevel = 1, verbose = FALSE), "between 0 and 1")
  expect_error(scoreRegionSets(counts, regionSets = "absentSet", verbose = FALSE), "absent from the object")
  expect_error(scoreRegionSets(counts, comparisons = list("promoterCpG"), verbose = FALSE), "length two")
  expect_error(scoreRegionSets(counts, comparisons = list(c("promoterCpG", "absentSet")), verbose = FALSE),
               "were not scored")

  # The bins hold raw counts, so a ratio against them cannot read a normalised assay
  expect_error(scoreRegionSets(counts, assay = "offset", reference = "background", verbose = FALSE),
               "cannot be set")
})


test_that("a single library leaves nothing to compare across", {

  expect_error(scoreRegionSets(toyCounts(), verbose = FALSE), "two libraries")
})


test_that("reference = 'background' needs the bins to exist", {

  counts <- countReads(toyRegionSet(),
                       bamFiles = rep(toyBamFile(), 2),
                       sampleNames = c("first", "second"),
                       verbose = FALSE)

  expect_error(scoreRegionSets(counts, minRegions = 3, verbose = FALSE), "countBackground")
})


test_that("two identical libraries leave no spread to test", {

  counts <- countReads(toyRegionSet(),
                       bamFiles = rep(toyBamFile(), 2),
                       sampleNames = c("first", "second"),
                       verbose = FALSE)

  setScores <- scoreRegionSets(counts, reference = "regions", minRegions = 3, verbose = FALSE)
  comparisons <- resultsTable(setScores)

  # The same file counted twice gives the same difference twice, which is degenerate and not decisive
  expect_true(all(is.na(comparisons$t.statistic)))
  expect_true(all(is.na(comparisons$p.value)))
  expect_true(all(is.finite(comparisons$mean.delta.score)))
})


test_that("a set below minRegions is dropped, and two of them must survive", {

  regionRanges <- GenomicRanges::GRanges(
    seqnames = c(rep("seq1", 6), rep("seq2", 3)),
    ranges = IRanges::IRanges(start = c(1L, 200L, 400L, 600L, 800L, 1000L, 1L, 300L, 600L),
                              width = 150L)
  )

  regionRanges$setName <- c(rep("bigSet", 4), rep("otherSet", 3), rep("tinySet", 2))

  counts <- countReads(splitLoadRegions(regionRanges, splitBy = "setName",
                                        seqlevelsStyle = NULL, verbose = FALSE),
                       bamFiles = rep(toyBamFile(), 2),
                       sampleNames = c("first", "second"),
                       verbose = FALSE)

  expect_warning(setScores <- scoreRegionSets(counts, reference = "regions",
                                              minRegions = 3, verbose = FALSE),
                 "tinySet")

  expect_setequal(regionSetNames(setScores), c("bigSet", "otherSet"))

  expect_error(suppressWarnings(scoreRegionSets(counts, reference = "regions",
                                                minRegions = 5, verbose = FALSE)),
               "must survive")
})


test_that("scoreRegionSets reports what it is doing when asked", {

  expect_message(scoreRegionSets(exampleCounts()), "Comparing")
})


test_that("the accessors and the print method return what they say", {

  setScores <- scoreRegionSets(exampleCounts(), verbose = FALSE)

  expect_identical(resultsTable(setScores), setScores@comparisons)
  expect_identical(scoreTable(setScores), setScores@scores)
  expect_identical(regionSetNames(setScores), sort(unique(setScores@scores$region.set)))

  expect_output(show(setScores), "RegionSetDE.setScores")
  expect_output(show(setScores), "region sets")
})


test_that("plotSetSignal draws a set scores object", {

  counts <- exampleCounts()
  setScores <- scoreRegionSets(counts, verbose = FALSE)
  twoSets <- regionSetNames(setScores)[1:2]

  expect_s3_class(plotSetSignal(setScores), "ggplot")
  expect_s3_class(plotSetSignal(setScores, groupBy = "condition"), "ggplot")
  expect_s3_class(plotSetSignal(setScores, groupBy = "condition", legendPosition = "right"), "ggplot")
  expect_s3_class(plotSetSignal(setScores, annotate = FALSE), "ggplot")
  expect_s3_class(plotSetSignal(setScores, set = twoSets), "ggplot")
  expect_s3_class(plotSetSignal(setScores, comparisons = list(twoSets)), "ggplot")
  expect_s3_class(plotSetSignal(setScores, title = "A title", subtitle = "A subtitle"), "ggplot")

  expect_error(plotSetSignal(setScores, set = "absentSet"), "not scored")
  expect_error(plotSetSignal(setScores, groupBy = "absentColumn"), "absent from the metadata")
})
