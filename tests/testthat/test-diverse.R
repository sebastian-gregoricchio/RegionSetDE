test_that("a tiled object says that it is tiled", {

  plainCounts <- toyCounts()
  tiledCounts <- toyCounts(tileWidth = 100L)

  expect_identical(plainCounts@counting.level, "region")
  expect_identical(tiledCounts@counting.level, "tile")
})


test_that("the sample interval rests on the samples and not on the regions", {

  fit <- exampleFit()

  setResults <- testRegionSets(fit, contrast = exampleContrast(),
                               effectMethod = "sample", verbose = FALSE)
  setTable <- resultsTable(setResults)

  expect_true(all(c("sample.delta.log2FC", "sample.delta.SE", "sample.delta.df",
                    "heterogeneity.CI.lower", "CI.type") %in% colnames(setTable)))

  # The degrees of freedom come from the design, so they are the same whatever the set holds
  expect_true(all(setTable$sample.delta.df == (nrow(fit@design) - ncol(fit@design))))
})


test_that("effectMethod switches which interval is reported without dropping the other", {

  regionResults <- testRegionSets(exampleFit(), contrast = exampleContrast(),
                                  effectMethod = "region", verbose = FALSE)
  regionTable <- resultsTable(regionResults)

  expect_true(all(regionTable$CI.type == "region"))
  expect_equal(regionTable$CI.lower, regionTable$heterogeneity.CI.lower)
  expect_false("sample.delta.log2FC" %in% colnames(regionTable))
})


test_that("overlapping comparison rows are found on the coordinates", {

  setResults <- testRegionSets(exampleFit(), contrast = exampleContrast(),
                               overlapPolicy = "drop", verbose = FALSE)

  expect_true("n.comparison.overlapping" %in% colnames(resultsTable(setResults)))
  expect_true(all(resultsTable(setResults)$n.comparison.overlapping >= 0))
})


test_that("the comparison pool can be named and is recorded", {

  setNames <- regionSetNames(exampleFit())
  skip_if(length(setNames) < 3, "at least three sets are needed to narrow the pool")

  universe <- makeSetUniverse(exampleFit(),
                              universeSets = setNames[2:3],
                              verbose = FALSE)

  expect_setequal(universe@comparison.sets, setNames[2:3])
})


test_that("bigWig values are kept out of the count engines unless declared", {

  counts <- exampleCounts()
  S4Vectors::metadata(counts)$signal.type <- "bigwig"
  S4Vectors::metadata(counts)$count.like <- FALSE

  expect_error(fitRegions(counts, design = ~ condition, engine = "edgeR", verbose = FALSE),
               "limma")

  # The assertion is the way through, and it has to be made deliberately
  expect_s4_class(fitRegions(counts, design = ~ condition, engine = "edgeR",
                             assumeCountLike = TRUE, verbose = FALSE),
                  "RegionSetDE.fit")
})


test_that("the limma engine fits a signal that is not counts", {

  counts <- exampleCounts()
  S4Vectors::metadata(counts)$signal.type <- "bigwig"
  S4Vectors::metadata(counts)$count.like <- FALSE

  fit <- fitRegions(counts, design = ~ condition, engine = "limma", verbose = FALSE)
  results <- testRegions(fit, contrast = exampleContrast(), verbose = FALSE)

  expect_identical(fit@engine, "limma")
  expect_equal(nrow(resultsTable(results)), nrow(counts))
})


test_that("the normalisation holdout travels into the dispersion estimate", {

  counts <- normalizeCounts(exampleCounts(), method = "background",
                            backgroundHoldout = 0.5, verbose = FALSE)

  expect_gt(length(S4Vectors::metadata(counts)$background.holdout), 0)

  nullDispersion <- estimateNullDispersion(counts, source = "background", verbose = FALSE)

  expect_identical(nullDispersion$holdout.type, "normalization and dispersion")
  expect_true(all(nullDispersion$holdout.index %in% S4Vectors::metadata(counts)$background.holdout))
})


test_that("without that split the holdout is a dispersion holdout and says so", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)
  nullDispersion <- estimateNullDispersion(counts, source = "background", verbose = FALSE)

  expect_identical(nullDispersion$holdout.type, "dispersion")
})
