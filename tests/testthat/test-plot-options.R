## The plotting files are about a quarter of the package, and almost all of it is
## branches on the display arguments. One call per branch is what covers them.


test_that("plotSetEffect draws each of its variants", {

  setResults <- exampleSetResults()

  for (plotValue in c("delta.log2FC", "mean.log2FC")) {
    expect_s3_class(plotSetEffect(setResults, value = plotValue), "ggplot")
  }

  for (colourVariable in c("FDR", "direction", "none")) {
    expect_s3_class(plotSetEffect(setResults, colourBy = colourVariable), "ggplot")
  }

  for (setOrder in c("effect", "name")) {
    expect_s3_class(plotSetEffect(setResults, orderBy = setOrder), "ggplot")
  }

  expect_s3_class(plotSetEffect(setResults, showN = FALSE, FDR = 0.5), "ggplot")
})


test_that("plotSetDistribution draws each of its styles", {

  setResults <- exampleSetResults()

  for (plotStyle in c("violin", "boxplot", "ecdf")) {
    expect_s3_class(plotSetDistribution(setResults, style = plotStyle), "ggplot")
  }

  expect_s3_class(plotSetDistribution(setResults, set = "promoterCpG",
                                      annotate = FALSE), "ggplot")
})


test_that("plotSetSignal draws each of its styles", {

  setResults <- exampleSetResults()

  for (plotStyle in c("violin", "boxplot")) {
    expect_s3_class(plotSetSignal(setResults, style = plotStyle,
                                  groupBy = "condition"), "ggplot")
  }

  expect_s3_class(plotSetSignal(setResults, set = "promoterCpG",
                                log2Scale = FALSE, annotate = FALSE), "ggplot")
})


test_that("plotVolcano draws each of its variants", {

  results <- exampleResults()

  for (colourVariable in c("diff.status", "region.set")) {
    expect_s3_class(plotVolcano(results, colourBy = colourVariable), "ggplot")
  }

  for (yVariable in c("FDR", "p.value")) {
    expect_s3_class(plotVolcano(results, yValue = yVariable), "ggplot")
  }

  expect_s3_class(plotVolcano(results, facetBySet = FALSE, showCounts = FALSE), "ggplot")
  expect_s3_class(plotVolcano(results, set = "promoterCpG", labelTop = 3), "ggplot")
  expect_s3_class(plotVolcano(results, FDR = 0.5, log2FC = 0.5), "ggplot")
})


test_that("plotResultsMA draws each of its variants", {

  results <- exampleResults()

  for (colourVariable in c("diff.status", "region.set")) {
    expect_s3_class(plotResultsMA(results, colourBy = colourVariable), "ggplot")
  }

  expect_s3_class(plotResultsMA(results, facetBySet = FALSE, showTrend = FALSE), "ggplot")
  expect_s3_class(plotResultsMA(results, set = "geneBody", showCounts = FALSE), "ggplot")
})


test_that("plotTopHeatmap draws each of its rankings and directions", {

  results <- exampleResults()

  for (ranking in c("log2FC", "FDR", "stat")) {
    expect_false(is.null(plotTopHeatmap(results, n = 8, sortBy = ranking, FDR = 1)))
  }

  for (changeDirection in c("both", "up", "down")) {
    expect_false(is.null(plotTopHeatmap(results, n = 8, direction = changeDirection,
                                        FDR = 1)))
  }

  expect_false(is.null(plotTopHeatmap(results, n = 8, FDR = 1, scaleRows = FALSE,
                                      showRegionNames = TRUE, showLog2FC = FALSE)))
})


test_that("plotSampleCorrelation draws each correlation", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  for (correlationMethod in c("spearman", "pearson", "kendall")) {
    expect_s3_class(plotSampleCorrelation(counts, method = correlationMethod), "ggplot")
  }

  expect_s3_class(plotSampleCorrelation(counts, cluster = FALSE,
                                        showValues = FALSE), "ggplot")
  expect_s3_class(plotSampleCorrelation(counts, set = "promoterCpG",
                                        excludeDiagonal = TRUE), "ggplot")
  expect_s3_class(plotSampleCorrelation(counts, facetBySet = TRUE), "ggplot")
})


test_that("plotRegionPCA draws each of its variants", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  expect_s3_class(plotRegionPCA(counts, colourBy = "condition",
                                shapeBy = "sex"), "ggplot")
  expect_s3_class(plotRegionPCA(counts, labelBy = NULL), "ggplot")
  expect_s3_class(plotRegionPCA(counts, labelBy = "condition"), "ggplot")
  expect_s3_class(plotRegionPCA(counts, dimensions = c(1, 3)), "ggplot")
  expect_s3_class(plotRegionPCA(counts, topRegions = 500), "ggplot")
  expect_s3_class(plotRegionPCA(counts, useOffsets = FALSE), "ggplot")
})


test_that("plotSetMA draws with and without a chosen set", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  expect_s3_class(plotSetMA(counts, groupBy = "condition"), "ggplot")
  expect_s3_class(plotSetMA(counts, set = "geneBody", groupBy = "condition",
                            showTrend = FALSE, showMedian = FALSE), "ggplot")

  plotData <- plotSetMA(counts, groupBy = "condition", returnData = TRUE)
  expect_s3_class(plotData, "data.frame")
  expect_gt(nrow(plotData), 0)
})


test_that("plotRegion draws with and without summarising", {

  results <- exampleResults()
  topRegion <- topRegions(results, n = 1, FDR = 1)$region.id

  expect_s3_class(plotRegion(results, region = topRegion,
                             groupBy = "condition"), "ggplot")
  expect_s3_class(plotRegion(results, region = topRegion, groupBy = "condition",
                             summarise = TRUE), "ggplot")
  expect_s3_class(plotRegion(results, region = topRegion, log2Scale = FALSE,
                             rotateX = FALSE), "ggplot")
})


test_that("plotNullCalibration draws each of its styles", {

  calibration <- checkNullCalibration(exampleFit(),
                                      contrast = exampleContrast(),
                                      source = "background",
                                      verbose = FALSE)

  expect_s3_class(plotNullCalibration(calibration, style = "histogram"), "ggplot")
  expect_s3_class(plotNullCalibration(calibration, bins = 20,
                                      title = "Calibration"), "ggplot")
})


test_that("plotUniverseMatching draws on each covariate", {

  setResults <- exampleSetResults()

  for (matchingCovariate in c("width", "abundance")) {
    expect_s3_class(plotUniverseMatching(setResults, set = "promoterCpG",
                                         covariate = matchingCovariate), "ggplot")
  }
})


test_that("the plots refuse a set that is not there", {

  results <- exampleResults()

  expect_error(plotVolcano(results, set = "notASet"))
  expect_error(plotResultsMA(results, set = "notASet"))
})
