## The packaged data has two biological replicates per strain, and both strains
## have a bio2 library, so subsetting to it gives a genuine one-versus-one. Both
## are male, which drops the sex difference out of the comparison as well.
singleSampleCounts <- function() {

  counts <- RegionSetDE::selectSamples(exampleCounts(),
                                       biologicalReplicate == "bio2",
                                       verbose = FALSE)

  counts <- RegionSetDE::normalizeCounts(counts, method = "background", verbose = FALSE)

  RegionSetDE::filterRegions(counts, verbose = FALSE)
}


test_that("the example data reduces to one library per strain", {

  counts <- singleSampleCounts()

  expect_equal(ncol(counts), 2)
  expect_setequal(as.character(SummarizedExperiment::colData(counts)$condition),
                  c("BN", "SHR"))
})


test_that("estimateNullDispersion returns a usable dispersion and a holdout", {

  counts <- singleSampleCounts()

  nullDispersion <- estimateNullDispersion(counts, source = "background",
                                           holdout = 0.5, verbose = FALSE)

  expect_type(nullDispersion, "list")
  expect_true(all(c("dispersion", "bcv", "source", "holdout.index") %in%
                    names(nullDispersion)))

  expect_true(is.finite(nullDispersion$dispersion))
  expect_gt(nullDispersion$dispersion, 0)

  # The bcv is the dispersion on the fold change scale
  expect_equal(nullDispersion$bcv, sqrt(nullDispersion$dispersion), tolerance = 1e-6)

  expect_identical(nullDispersion$source, "background")
  expect_gt(length(nullDispersion$holdout.index), 0)
})


test_that("the holdout keeps rows out of the estimate", {

  counts <- singleSampleCounts()

  withHoldout <- estimateNullDispersion(counts, source = "background",
                                        holdout = 0.5, verbose = FALSE)
  withoutHoldout <- estimateNullDispersion(counts, source = "background",
                                           holdout = 0, verbose = FALSE)

  expect_gt(length(withHoldout$holdout.index), length(withoutHoldout$holdout.index))
})


test_that("fitRegions accepts a supplied dispersion with no replicates", {

  counts <- singleSampleCounts()

  nullDispersion <- estimateNullDispersion(counts, source = "background",
                                           verbose = FALSE)

  fit <- fitRegions(counts, design = ~ condition, engine = "edgeR",
                    dispersion = nullDispersion, verbose = FALSE)

  expect_s4_class(fit, "RegionSetDE.fit")
  expect_equal(ncol(fitCounts(fit)), 2)
})


test_that("the dispersion keyword runs the estimation inside fitRegions", {

  counts <- singleSampleCounts()

  fit <- fitRegions(counts, design = ~ condition, engine = "edgeR",
                    dispersion = "background", verbose = FALSE)

  expect_s4_class(fit, "RegionSetDE.fit")
})


test_that("a design with no replicates falls back to the background dispersion", {

  counts <- singleSampleCounts()

  fit <- fitRegions(counts, design = ~ condition, engine = "edgeR", verbose = FALSE)

  expect_s4_class(fit, "RegionSetDE.fit")

  # The fallback is recorded rather than left to be guessed at
  expect_true(fit@dispersion$no.replicates)
  expect_true(fit@dispersion$fixed)
  expect_identical(fit@dispersion$source, "background")
  expect_gt(fit@dispersion$common, 0)

  # It is the same estimate estimateNullDispersion returns on the same object
  standalone <- estimateNullDispersion(counts, source = "background", verbose = FALSE)
  expect_equal(fit@dispersion$common, standalone$dispersion, tolerance = 1e-8)

  # And the held-out rows travel with the fit, so the calibration check has them
  expect_gt(length(fit@dispersion$holdout.index), 0)
})


test_that("both levels of testing work on a single pair of libraries", {

  counts <- singleSampleCounts()

  fit <- fitRegions(counts, design = ~ condition, engine = "edgeR",
                    dispersion = "background", verbose = FALSE)

  results <- testRegions(fit, contrast = exampleContrast(), verbose = FALSE)
  setResults <- testRegionSets(fit, contrast = exampleContrast(), verbose = FALSE)

  resultTable <- resultsTable(results)
  setTable <- resultsTable(setResults)

  expect_s4_class(results, "RegionSetDE.results")
  expect_true(all(resultTable$p.value >= 0 & resultTable$p.value <= 1))

  expect_s4_class(setResults, "RegionSetDE.setResults")
  expect_equal(nrow(setTable), length(regionSetNames(fit)))
})


test_that("the fold changes agree with the replicated analysis", {

  counts <- singleSampleCounts()

  fit <- fitRegions(counts, design = ~ condition, engine = "edgeR",
                    dispersion = "background", verbose = FALSE)

  singleTable <- resultsTable(testRegions(fit, contrast = exampleContrast(),
                                          verbose = FALSE))
  replicatedTable <- resultsTable(exampleResults())

  sharedRegions <- intersect(singleTable$region.id, replicatedTable$region.id)
  expect_gt(length(sharedRegions), 100)

  singleValues <- singleTable$log2FC[match(sharedRegions, singleTable$region.id)]
  replicatedValues <- replicatedTable$log2FC[match(sharedRegions, replicatedTable$region.id)]

  # The dispersion does not enter the point estimate, so the effect sizes track
  # each other, though not tightly: the two fits filter independently and a single
  # library is noisiest exactly where the fold changes are largest
  expect_gt(stats::cor(singleValues, replicatedValues, method = "spearman"), 0.5)
})


test_that("the calibration check runs on the rows held out of the estimate", {

  counts <- singleSampleCounts()

  fit <- fitRegions(counts, design = ~ condition, engine = "edgeR", verbose = FALSE)

  calibration <- checkNullCalibration(fit,
                                      contrast = exampleContrast(),
                                      source = "background",
                                      index = fit@dispersion$holdout.index,
                                      verbose = FALSE)

  expect_false(is.null(calibration))
  expect_s3_class(plotNullCalibration(calibration), "ggplot")
})


test_that("the dispersion is recorded in the provenance", {

  counts <- singleSampleCounts()

  fit <- fitRegions(counts, design = ~ condition, engine = "edgeR",
                    dispersion = "background", verbose = FALSE)

  expect_false(is.null(fit@dispersion))
})
