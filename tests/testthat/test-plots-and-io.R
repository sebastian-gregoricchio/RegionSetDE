test_that("loadExampleData returns each object with the right class", {

  expect_s4_class(loadExampleData("counts", verbose = FALSE), "RegionSetDE.counts")
  expect_s4_class(loadExampleData("fit", verbose = FALSE), "RegionSetDE.fit")
  expect_s4_class(loadExampleData("exclusionRegions", verbose = FALSE), "GRanges")
  expect_s3_class(loadExampleData("regions", verbose = FALSE), "data.frame")
  expect_s3_class(loadExampleData("sampleSheet", verbose = FALSE), "data.frame")
  expect_type(loadExampleData("buildMetadata", verbose = FALSE), "list")
})


test_that("blacklist is accepted as a synonym", {

  expect_identical(loadExampleData("blacklist", verbose = FALSE),
                   loadExampleData("exclusionRegions", verbose = FALSE))
})


test_that("loadExampleData rejects an unknown dataset", {
  expect_error(loadExampleData("notAnObject", verbose = FALSE))
})


test_that("the packaged region table covers the four sets", {

  regionTable <- loadExampleData("regions", verbose = FALSE)

  expect_setequal(unique(regionTable$setName),
                  c("promoterCpG", "promoterNonCpG", "geneBody", "intergenic"))
  expect_true(all(regionTable$end > regionTable$start))
  expect_false(any(duplicated(regionTable$regionId)))
})


test_that("the build metadata records what the script used", {

  buildMetadata <- loadExampleData("buildMetadata", verbose = FALSE)

  expect_identical(buildMetadata$genomeBuild, "rn4")
  expect_identical(buildMetadata$targetChromosomes, "chr12")
  expect_true(is.numeric(buildMetadata$fragmentLength))
})


test_that("the counts level plots build", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  expect_s3_class(plotNormComparison(exampleCounts()), "ggplot")
  expect_s3_class(plotSetMA(counts, groupBy = "condition"), "ggplot")
  expect_s3_class(plotRegionPCA(counts, colourBy = "condition"), "ggplot")
  expect_s3_class(plotSampleCorrelation(counts, groupBy = "condition"), "ggplot")
})


test_that("plotNormComparison can hand back its data", {

  plotData <- plotNormComparison(exampleCounts(), returnData = TRUE)

  expect_s3_class(plotData, "data.frame")
  expect_gt(nrow(plotData), 0)
})


test_that("the region level plots build", {

  results <- exampleResults()
  topRegion <- topRegions(results, n = 1, FDR = 1)$region.id

  expect_s3_class(plotVolcano(results), "ggplot")
  expect_s3_class(plotResultsMA(results), "ggplot")
  expect_s3_class(plotRegion(results, region = topRegion, groupBy = "condition"),
                  "ggplot")
})


test_that("the set level plots build", {

  setResults <- exampleSetResults()

  expect_s3_class(plotSetEffect(setResults), "ggplot")
  expect_s3_class(plotSetDistribution(setResults), "ggplot")
  expect_s3_class(plotSetSignal(setResults, groupBy = "condition"), "ggplot")
  expect_s3_class(plotUniverseMatching(setResults, set = "promoterCpG"), "ggplot")
})


test_that("plotNullCalibration builds from a calibration object", {

  calibration <- checkNullCalibration(exampleFit(),
                                      contrast = exampleContrast(),
                                      source = "background",
                                      verbose = FALSE)

  expect_s3_class(plotNullCalibration(calibration), "ggplot")
})


test_that("plotTopHeatmap builds", {

  heatmap <- plotTopHeatmap(exampleResults(), n = 10, FDR = 1)
  expect_false(is.null(heatmap))
})


test_that("exportResults writes files and cleans up after itself", {

  outputDirectory <- file.path(tempdir(), "regionSetDE-test-export")
  dir.create(outputDirectory, showWarnings = FALSE)
  on.exit(unlink(outputDirectory, recursive = TRUE), add = TRUE)

  exportResults(exampleResults(), path = outputDirectory,
                prefix = "test", verbose = FALSE)

  writtenFiles <- list.files(outputDirectory)

  expect_gt(length(writtenFiles), 0)
  expect_true(any(grepl("^test", writtenFiles)))
})


test_that("exportResults creates the output directory when it is missing", {

  parentDirectory <- file.path(tempdir(), "regionSetDE-test-nested")
  outputDirectory <- file.path(parentDirectory, "deeper")
  on.exit(unlink(parentDirectory, recursive = TRUE), add = TRUE)

  unlink(parentDirectory, recursive = TRUE)
  expect_false(dir.exists(outputDirectory))

  exportResults(exampleResults(), path = outputDirectory,
                prefix = "nested", verbose = FALSE)

  expect_true(dir.exists(outputDirectory))
  expect_gt(length(list.files(outputDirectory)), 0)
})
