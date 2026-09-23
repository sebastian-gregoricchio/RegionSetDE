test_that("the negative binomial fit recovers known parameters", {

  set.seed(3)
  simulatedCounts <- stats::rnbinom(2e5, size = 2, mu = 10)
  fittedDistribution <- RegionSetDE:::.fitNegativeBinomial(simulatedCounts)

  expect_equal(fittedDistribution$mu, 10, tolerance = 0.01)
  expect_equal(fittedDistribution$size, 2, tolerance = 0.03)

  # Counts without overdispersion fall back on the Poisson
  expect_identical(RegionSetDE:::.fitNegativeBinomial(rep(5L, 100))$size, Inf)
})


test_that("the windows overlap by half and stop at the end of the chromosome", {

  windowRanges <- RegionSetDE:::.greylistWindows(chromosomeLengths = c(chrA = 2500), binSize = 1000L)

  expect_identical(BiocGenerics::start(windowRanges), c(1L, 501L, 1001L, 1501L, 2001L))
  expect_identical(BiocGenerics::end(windowRanges), c(1000L, 1500L, 2000L, 2500L, 2500L))
})


test_that("a pile of reads in the input is greylisted and nothing else", {

  inputFile <- syntheticInput("inputPile", pileStart = 100000L)
  greylist <- makeGreylist(inputFile, quantile = 0.9999, verbose = FALSE)

  expect_s4_class(greylist, "GRanges")
  expect_length(greylist, 1)
  expect_true(IRanges::overlapsAny(GenomicRanges::GRanges("chrT:100000-100900"), greylist))
  expect_lt(sum(BiocGenerics::width(greylist)), 5000)

  thresholdTable <- S4Vectors::metadata(greylist)$thresholds
  expect_identical(thresholdTable$input, "inputPile")
  expect_identical(thresholdTable$fragments, 23000)
  expect_identical(S4Vectors::metadata(greylist)$parameters$binSize, 1024L)
})


test_that("the inputs are pooled, and minInputs keeps what several of them flag", {

  firstInput <- syntheticInput("inputFirst", pileStart = 100000L, seed = 1L)
  secondInput <- syntheticInput("inputSecond", pileStart = 150000L, seed = 2L)
  sharedInput <- syntheticInput("inputShared", pileStart = 100000L, seed = 3L)

  pooledList <- makeGreylist(c(firstInput, secondInput), quantile = 0.9999, pairedEnd = FALSE, verbose = FALSE)
  expect_length(pooledList, 2)
  expect_identical(pooledList$inputs, c("inputFirst", "inputSecond"))
  expect_true(all(pooledList$n.inputs == 1))

  expect_length(makeGreylist(c(firstInput, secondInput), quantile = 0.9999, minInputs = 2, pairedEnd = FALSE, verbose = FALSE), 0)

  sharedList <- makeGreylist(c(firstInput, secondInput, sharedInput), quantile = 0.9999, minInputs = 2, pairedEnd = FALSE, verbose = FALSE)
  expect_length(sharedList, 1)
  expect_identical(sharedList$n.inputs, 2L)
  expect_true(IRanges::overlapsAny(GenomicRanges::GRanges("chrT:100000-100900"), sharedList))
})


test_that("a sample sheet brings each distinct input once, under its name", {

  firstInput <- syntheticInput("inputFirst", pileStart = 100000L, seed = 1L)
  secondInput <- syntheticInput("inputSecond", pileStart = 150000L, seed = 2L)

  sampleSheet <- loadSampleSheet(data.frame(sample = c("a", "b", "c", "d"),
                                            bam = c("a.bam", "b.bam", "c.bam", "d.bam"),
                                            input = c(firstInput, firstInput, secondInput, NA)),
                                 checkFiles = FALSE, verbose = FALSE)

  greylist <- makeGreylist(sampleSheet, quantile = 0.9999, pairedEnd = FALSE, verbose = FALSE)

  expect_identical(S4Vectors::metadata(greylist)$thresholds$input, c("inputFirst", "inputSecond"))
  expect_error(makeGreylist(data.frame(sample = "a", bam = "a.bam"), pairedEnd = FALSE, verbose = FALSE), "no 'input' column")
})


test_that("an empty chromosome widens the fit unless it is excluded", {

  inputFile <- syntheticInput("inputEmptyContig", pileStart = 100000L, emptyContig = TRUE)

  withEmpty <- S4Vectors::metadata(makeGreylist(inputFile, pairedEnd = FALSE, verbose = FALSE))$thresholds
  withoutEmpty <- S4Vectors::metadata(makeGreylist(inputFile, excludeChromosomes = "chrU", pairedEnd = FALSE, verbose = FALSE))$thresholds

  expect_lt(withEmpty$mean, withoutEmpty$mean)
  expect_lt(withEmpty$size, withoutEmpty$size)
  expect_lt(withoutEmpty$windows, withEmpty$windows)
})


test_that("the arguments are checked", {

  inputFile <- syntheticInput("inputChecks")

  expect_error(makeGreylist(file.path(tempdir(), "absent.bam"), pairedEnd = FALSE, verbose = FALSE), "do not exist")
  expect_error(makeGreylist(c(inputFile, inputFile), pairedEnd = FALSE, verbose = FALSE), "more than once")
  expect_error(makeGreylist(inputFile, quantile = 1, pairedEnd = FALSE, verbose = FALSE), "quantile")
  expect_error(makeGreylist(inputFile, minInputs = 2, pairedEnd = FALSE, verbose = FALSE), "minInputs")
  expect_error(makeGreylist(inputFile, binSize = 1, pairedEnd = FALSE, verbose = FALSE), "binSize")
  expect_error(makeGreylist(inputFile, excludeChromosomes = "chrT", pairedEnd = FALSE, verbose = FALSE), "no window is left")
})


test_that("applyGreylist removes the regions, logs the step and keeps the blacklist", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)

  regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
                              splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

  # The first half of the exclusion list acts as a blacklist, the second half as a greylist
  blacklist <- exclusionRegions[seq_len(length(exclusionRegions) %/% 2)]
  greylist <- exclusionRegions[-seq_len(length(exclusionRegions) %/% 2)]

  blacklisted <- applyBlacklist(regions, blacklist = blacklist, verbose = FALSE)
  greylisted <- applyGreylist(blacklisted, greylist = greylist, verbose = FALSE)

  expect_lte(sum(lengths(greylisted@regions)), sum(lengths(blacklisted@regions)))
  expect_true("greylist" %in% greylisted@filtering.log$step)
  expect_identical(length(greylisted@blacklist), length(blacklisted@blacklist))
  expect_identical(greylisted@parameters$greylist$n.regions, length(greylist))

  trimmed <- applyGreylist(regions, greylist = greylist, trimRegions = TRUE, verbose = FALSE)
  expect_gte(sum(lengths(trimmed@regions)), sum(lengths(applyGreylist(regions, greylist = greylist, verbose = FALSE)@regions)))

  expect_error(applyGreylist(exampleCounts(), greylist = greylist, verbose = FALSE), "greylist must be applied before the counting")
})
