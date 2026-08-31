test_that("testRegions returns a table with the expected columns", {

  results <- exampleResults()
  resultTable <- resultsTable(results)

  expect_s4_class(results, "RegionSetDE.results")
  expect_s3_class(resultTable, "data.frame")

  expectedColumns <- c("region.set", "region.id", "seqnames", "start", "end",
                       "log2FC", "p.value", "FDR", "diff.status")
  expect_true(all(expectedColumns %in% colnames(resultTable)))
})


test_that("the p-values and the FDR stay inside the unit interval", {

  resultTable <- resultsTable(exampleResults())

  expect_true(all(resultTable$p.value >= 0 & resultTable$p.value <= 1))
  expect_true(all(resultTable$FDR >= 0 & resultTable$FDR <= 1))
  expect_true(all(resultTable$FDR >= resultTable$p.value - 1e-8))
})


test_that("reversing the contrast flips the sign of the fold change", {

  fit <- exampleFit()

  forward <- resultsTable(testRegions(fit, contrast = c("condition", "SHR", "BN"),
                                      verbose = FALSE))
  reverse <- resultsTable(testRegions(fit, contrast = c("condition", "BN", "SHR"),
                                      verbose = FALSE))

  expect_equal(forward$log2FC, -reverse$log2FC, tolerance = 1e-6)
  expect_equal(forward$p.value, reverse$p.value, tolerance = 1e-6)
})


test_that("a fold change threshold cannot increase the number of hits", {

  fit <- exampleFit()

  plain <- resultsTable(testRegions(fit, contrast = exampleContrast(),
                                    verbose = FALSE))
  thresholded <- resultsTable(testRegions(fit, contrast = exampleContrast(),
                                          lfcThreshold = 1, verbose = FALSE))

  countSignificant <- function(resultTable) {
    sum(dplyr::pull(dplyr::filter(resultTable, FDR < 0.05), FDR) < 0.05)
  }

  expect_lte(countSignificant(thresholded), countSignificant(plain))
})


test_that("testRegions errors on a level that is not in the design", {

  expect_error(
    testRegions(exampleFit(), contrast = c("condition", "WKY", "BN"), verbose = FALSE)
  )
})


test_that("testRegionSets returns one row per region set", {

  setResults <- exampleSetResults()
  setTable <- resultsTable(setResults)

  expect_s4_class(setResults, "RegionSetDE.setResults")
  expect_s3_class(setTable, "data.frame")
  expect_equal(nrow(setTable), length(regionSetNames(exampleFit())))
  expect_true(all(c("region.set", "n.regions") %in% colnames(setTable)))
})


test_that("the set level confidence intervals bracket the effect", {

  setTable <- resultsTable(exampleSetResults())

  expect_true(all(setTable$CI.lower <= setTable$delta.log2FC + 1e-8))
  expect_true(all(setTable$CI.upper >= setTable$delta.log2FC - 1e-8))
})


test_that("testSetContrast compares two named sets", {

  setContrast <- testSetContrast(exampleFit(),
                                 contrast = exampleContrast(),
                                 set1 = "promoterCpG",
                                 set2 = "intergenic",
                                 verbose = FALSE)

  expect_false(is.null(setContrast))
})


test_that("checkNullCalibration runs against the background", {

  calibration <- checkNullCalibration(exampleFit(),
                                      contrast = exampleContrast(),
                                      source = "background",
                                      verbose = FALSE)

  expect_false(is.null(calibration))
})


test_that("contrastName reports the contrast at both levels", {

  expect_type(contrastName(exampleResults()), "character")
  expect_type(contrastName(exampleSetResults()), "character")
  expect_length(contrastName(exampleResults()), 1)
})


test_that("regionSetNames agrees across the classes", {

  expectedSets <- c("promoterCpG", "promoterNonCpG", "geneBody", "intergenic")

  expect_setequal(regionSetNames(exampleCounts()), expectedSets)
  expect_setequal(regionSetNames(exampleFit()), expectedSets)
  expect_setequal(regionSetNames(exampleResults()), expectedSets)
  expect_setequal(regionSetNames(exampleSetResults()), expectedSets)
})


test_that("resultCounts and resultRanges return the carried objects", {

  results <- exampleResults()

  expect_s4_class(resultCounts(results), "RegionSetDE.counts")
  expect_s4_class(resultRanges(results), "GRanges")
  expect_equal(length(resultRanges(results)), nrow(resultsTable(results)))

  expect_s4_class(resultCounts(exampleSetResults()), "RegionSetDE.counts")
})


test_that("resultRanges has no method for the set level", {
  expect_error(resultRanges(exampleSetResults()))
})


test_that("tileTable refuses an object that was not tiled", {
  expect_error(tileTable(exampleResults()))
})


test_that("topRegions ranks and filters as asked", {

  results <- exampleResults()

  ranked <- topRegions(results, n = 10, FDR = 1)
  expect_s3_class(ranked, "data.frame")
  expect_lte(nrow(ranked), 10)
  expect_false(is.unsorted(ranked$FDR))

  oneSet <- topRegions(results, n = 5, set = "promoterCpG", FDR = 1)
  expect_setequal(as.character(oneSet$region.set), "promoterCpG")

  upOnly <- topRegions(results, n = 5, FDR = 1, direction = "up")
  expect_true(all(upOnly$log2FC > 0))
})
