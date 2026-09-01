## The testing files are around 1600 lines, most of it branching on how the test
## is run. The list classes and their accessors had no coverage at all.


test_that("testRegions honours its correction and threshold arguments", {

  fit <- exampleFit()

  for (correction in c("BH", "BY", "holm", "bonferroni")) {

    results <- testRegions(fit, contrast = exampleContrast(),
                           adjustMethod = correction, verbose = FALSE)

    resultTable <- resultsTable(results)

    expect_true(all(resultTable$FDR >= resultTable$p.value - 1e-8),
                info = paste("adjustment below the raw p-value with", correction))
  }

  # A stricter threshold cannot call more regions
  loose <- resultsTable(testRegions(fit, contrast = exampleContrast(),
                                    FDR = 0.5, verbose = FALSE))
  strict <- resultsTable(testRegions(fit, contrast = exampleContrast(),
                                     FDR = 0.01, verbose = FALSE))

  expect_lte(sum(strict$diff.status != "null"), sum(loose$diff.status != "null"))
})


test_that("testRegions can be restricted to a subset of the sets", {

  results <- testRegions(exampleFit(), contrast = exampleContrast(),
                         regionSets = c("promoterCpG", "intergenic"),
                         verbose = FALSE)

  expect_setequal(unique(as.character(resultsTable(results)$region.set)),
                  c("promoterCpG", "intergenic"))
})


test_that("carryCounts decides whether the counts travel with the results", {

  fit <- exampleFit()

  carried <- testRegions(fit, contrast = exampleContrast(),
                         carryCounts = TRUE, verbose = FALSE)
  notCarried <- testRegions(fit, contrast = exampleContrast(),
                            carryCounts = FALSE, verbose = FALSE)

  expect_s4_class(resultCounts(carried), "RegionSetDE.counts")
  expect_s4_class(notCarried, "RegionSetDE.results")
})


test_that("a named list of contrasts gives a results list", {

  fit <- exampleFit()

  resultsList <- testRegions(fit,
                             contrast = list(strainEffect = c("condition", "SHR", "BN"),
                                             reversed = c("condition", "BN", "SHR")),
                             verbose = FALSE)

  expect_s4_class(resultsList, "RegionSetDE.resultsList")
  expect_length(resultsList, 2)
  expect_setequal(names(resultsList), c("strainEffect", "reversed"))

  # The accessors added for the list classes
  expect_type(contrastName(resultsList), "character")
  expect_length(contrastName(resultsList), 2)

  stackedTable <- resultsTable(resultsList)
  expect_s3_class(stackedTable, "data.frame")
  expect_true("contrast" %in% colnames(stackedTable))
  expect_equal(colnames(stackedTable)[1], "contrast")
  expect_equal(length(unique(stackedTable$contrast)), 2)

  expect_s4_class(resultCounts(resultsList), "RegionSetDE.counts")
  expect_setequal(regionSetNames(resultsList), regionSetNames(fit))

  # Single elements come back as ordinary results
  expect_s4_class(resultsList[[1]], "RegionSetDE.results")
  expect_s4_class(resultsList$strainEffect, "RegionSetDE.results")
})


test_that("a named list of contrasts gives a set results list", {

  setResultsList <- testRegionSets(exampleFit(),
                                   contrast = list(strainEffect = c("condition", "SHR", "BN"),
                                                   reversed = c("condition", "BN", "SHR")),
                                   verbose = FALSE)

  expect_s4_class(setResultsList, "RegionSetDE.setResultsList")
  expect_length(setResultsList, 2)

  stackedTable <- resultsTable(setResultsList)
  expect_true("contrast" %in% colnames(stackedTable))
  expect_equal(nrow(stackedTable), 2 * length(regionSetNames(exampleFit())))

  expect_length(contrastName(setResultsList), 2)
  expect_setequal(regionSetNames(setResultsList), regionSetNames(exampleFit()))
})


test_that("testRegionSets runs either test on its own", {

  fit <- exampleFit()

  cameraOnly <- resultsTable(testRegionSets(fit, contrast = exampleContrast(),
                                            method = "camera", verbose = FALSE))
  fryOnly <- resultsTable(testRegionSets(fit, contrast = exampleContrast(),
                                         method = "fry", verbose = FALSE))

  expect_true(any(grepl("camera", colnames(cameraOnly))))
  expect_true(any(grepl("fry", colnames(fryOnly))))
})


test_that("testRegionSets honours the universe arguments", {

  fit <- exampleFit()

  matchedOnBoth <- testRegionSets(fit, contrast = exampleContrast(),
                                  matchOn = c("width", "abundance"), verbose = FALSE)
  matchedOnWidth <- testRegionSets(fit, contrast = exampleContrast(),
                                   matchOn = "width", verbose = FALSE)
  unmatched <- testRegionSets(fit, contrast = exampleContrast(),
                              matchOn = character(0), verbose = FALSE)

  for (setResults in list(matchedOnBoth, matchedOnWidth, unmatched)) {
    expect_s4_class(setResults, "RegionSetDE.setResults")
    expect_equal(nrow(resultsTable(setResults)), length(regionSetNames(fit)))
  }

  wideUniverse <- testRegionSets(fit, contrast = exampleContrast(),
                                 universeRatio = 2, verbose = FALSE)
  expect_s4_class(wideUniverse, "RegionSetDE.setResults")
})


test_that("testRegionSets can be given the correlation rather than estimating it", {

  setResults <- testRegionSets(exampleFit(), contrast = exampleContrast(),
                               interRegionCor = 0.05, verbose = FALSE)

  expect_s4_class(setResults, "RegionSetDE.setResults")
})


test_that("testSetContrast compares sets both ways round", {

  fit <- exampleFit()

  forward <- testSetContrast(fit, contrast = exampleContrast(),
                             set1 = "promoterCpG", set2 = "intergenic",
                             verbose = FALSE)
  reverse <- testSetContrast(fit, contrast = exampleContrast(),
                             set1 = "intergenic", set2 = "promoterCpG",
                             verbose = FALSE)

  expect_s4_class(forward, "RegionSetDE.setResults")
  expect_s4_class(reverse, "RegionSetDE.setResults")
})


test_that("testSetContrast refuses a set that is not there", {

  expect_error(
    testSetContrast(exampleFit(), contrast = exampleContrast(),
                    set1 = "promoterCpG", set2 = "notASet", verbose = FALSE)
  )
})


test_that("makeSetUniverse honours its matching arguments", {

  fit <- exampleFit()

  defaultUniverse <- makeSetUniverse(fit, verbose = FALSE)
  widthOnly <- makeSetUniverse(fit, match = "width", verbose = FALSE)
  unmatched <- makeSetUniverse(fit, match = character(0), verbose = FALSE)
  narrowRatio <- makeSetUniverse(fit, ratio = 2, strata = 3, verbose = FALSE)

  for (universe in list(defaultUniverse, widthOnly, unmatched, narrowRatio)) {
    expect_s4_class(universe, "RegionSetDE.universe")
  }
})


test_that("makeSetUniverse can be built for one set", {

  universe <- makeSetUniverse(exampleFit(), regionSets = "promoterCpG",
                              verbose = FALSE)

  expect_s4_class(universe, "RegionSetDE.universe")
})


test_that("the voom engine fits and tests", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)
  counts <- filterRegions(counts, verbose = FALSE)

  fit <- fitRegions(counts, design = ~ condition, engine = "voom", verbose = FALSE)

  expect_s4_class(fit, "RegionSetDE.fit")

  results <- testRegions(fit, contrast = exampleContrast(), verbose = FALSE)
  resultTable <- resultsTable(results)

  expect_true(all(resultTable$p.value >= 0 & resultTable$p.value <= 1))

  setResults <- testRegionSets(fit, contrast = exampleContrast(), verbose = FALSE)
  expect_equal(nrow(resultsTable(setResults)), length(regionSetNames(fit)))
})


test_that("the DESeq2 engine fits and tests", {

  skip_if_not_installed("DESeq2")

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)
  counts <- filterRegions(counts, verbose = FALSE)

  fit <- fitRegions(counts, design = ~ condition, engine = "deseq2", verbose = FALSE)

  expect_s4_class(fit, "RegionSetDE.fit")

  results <- testRegions(fit, contrast = exampleContrast(), verbose = FALSE)
  expect_s4_class(results, "RegionSetDE.results")
})


test_that("fitRegions refuses an engine it does not have", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  expect_error(
    fitRegions(counts, design = ~ condition, engine = "notAnEngine", verbose = FALSE)
  )
})


test_that("estimateNullDispersion refuses a source it does not have", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  expect_error(
    estimateNullDispersion(counts, source = "notASource", verbose = FALSE)
  )
})


test_that("checkNullCalibration runs from a region set as well", {

  calibration <- checkNullCalibration(exampleFit(),
                                      contrast = exampleContrast(),
                                      source = "regionSet",
                                      regionSets = "geneBody",
                                      verbose = FALSE)

  expect_false(is.null(calibration))
})
