# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

## The results table carries the average signal of every group, and enough about the test statistic,
## its distribution and its degrees of freedom, to take the analysis further, into a power analysis
## for instance. The degrees of freedom are checked the only way that proves them right: handed back to
## the distribution they belong to, with the statistic, they have to give the p-value of the engine.

exampleCounts <- normalizeCounts(loadExampleData("counts", verbose = FALSE), method = "TMM", verbose = FALSE)
exampleContrast <- c("condition", "SHR", "BN")

fitWith <- function(engine, ...) {
  suppressMessages(fitRegions(exampleCounts, design = ~ condition, engine = engine, verbose = FALSE, ...))
}


test_that("the statistic, its distribution and its degrees of freedom give back the p-value of every engine", {

  # edgeR, quasi-likelihood F: numerator and denominator degrees of freedom
  edgeTable <- resultsTable(testRegions(fitWith("edgeR"), contrast = exampleContrast, verbose = FALSE))
  expect_true(all(edgeTable$stat.distribution == "f"))
  expect_true(all(edgeTable$df1 == 1))
  expect_equal(stats::pf(edgeTable$stat, edgeTable$df1, edgeTable$df2, lower.tail = FALSE), edgeTable$p.value, tolerance = 1e-8)

  # edgeR with a fixed dispersion, likelihood ratio: one degree of freedom per coefficient tested
  lrtTable <- resultsTable(testRegions(fitWith("edgeR", dispersion = 0.05), contrast = exampleContrast, verbose = FALSE))
  expect_true(all(lrtTable$stat.distribution == "chisq"))
  expect_true(all(is.na(lrtTable$df2)))
  expect_equal(stats::pchisq(lrtTable$stat, lrtTable$df1, lower.tail = FALSE), lrtTable$p.value, tolerance = 1e-8)

  # voom and limma-trend, moderated t: residual plus prior degrees of freedom
  for (engine in c("voom", "limma")) {
    limmaTable <- resultsTable(testRegions(fitWith(engine), contrast = exampleContrast, verbose = FALSE))
    expect_true(all(limmaTable$stat.distribution == "t"))
    expect_equal(2 * stats::pt(-abs(limmaTable$stat), limmaTable$df1), limmaTable$p.value, tolerance = 1e-8)
  }

  # DESeq2, Wald: a standard normal, no degrees of freedom at all
  deseqTable <- resultsTable(testRegions(fitWith("DESeq2"), contrast = exampleContrast, verbose = FALSE))
  expect_true(all(deseqTable$stat.distribution == "norm"))
  expect_true(all(is.na(deseqTable$df1)) & all(is.na(deseqTable$df2)))
  expect_equal(2 * stats::pnorm(-abs(deseqTable$stat)), deseqTable$p.value, tolerance = 1e-8)
})


test_that("a threshold test says that its statistic follows no standard distribution", {

  thresholdTable <- resultsTable(testRegions(fitWith("voom"), contrast = exampleContrast, lfcThreshold = 1, verbose = FALSE))
  expect_true(all(is.na(thresholdTable$stat.distribution)))
})


test_that("the average signal of each group is the one of the engine, computed on that group alone", {

  for (engine in c("edgeR", "voom", "DESeq2")) {
    fit <- fitWith(engine)
    resultTable <- resultsTable(testRegions(fit, contrast = exampleContrast, verbose = FALSE))
    groupColumns <- c("average.signal.SHR", "average.signal.BN")

    # Right after the overall average, in the order of the levels in the sample table rather than of the
    # contrast, so that every contrast of a list has the same columns in the same place
    expect_setequal(colnames(resultTable)[match("average.signal", colnames(resultTable)) + 1:2], groupColumns)

    shrSamples <- which(SummarizedExperiment::colData(fit@counts)$condition == "SHR")
    expectedSignal <- switch(fit@engine,
                             "edgeR" = edgeR::aveLogCPM(fit@fit$dge[, shrSamples]),
                             "voom" = rowMeans(fit@fit$voom$E[, shrSamples]),
                             "deseq2" = log2(rowMeans(DESeq2::counts(fit@fit$object, normalized = TRUE)[, shrSamples]) + 1))
    expect_equal(resultTable$average.signal.SHR, as.numeric(expectedSignal), tolerance = 1e-10)

    # A region called up has more signal in the first group than in the second
    upRows <- resultTable$diff.status == "up"
    if (any(upRows)) {
      expect_true(all(resultTable$average.signal.SHR[upRows] > resultTable$average.signal.BN[upRows]))
    }
  }
})


test_that("signalBy picks the column of the group averages, or none", {

  fit <- fitWith("edgeR")

  noGroups <- resultsTable(testRegions(fit, contrast = exampleContrast, signalBy = FALSE, verbose = FALSE))
  expect_false(any(grepl("^average\\.signal\\.", colnames(noGroups))))

  expect_error(testRegions(fit, contrast = exampleContrast, signalBy = "tissue", verbose = FALSE), "signalBy")

  # Any column of the sample table works, and a level starting with a digit keeps its name: the column
  # name as a whole is what has to be syntactic, and it already starts with a letter
  timedFit <- fit
  SummarizedExperiment::colData(timedFit@counts)$timepoint <- rep(c("4h", "24h"), length.out = ncol(timedFit@counts))

  timedTable <- resultsTable(testRegions(timedFit, contrast = exampleContrast, signalBy = "timepoint", verbose = FALSE))
  expect_true(all(c("average.signal.4h", "average.signal.24h") %in% colnames(timedTable)))
  expect_false(any(c("average.signal.SHR", "average.signal.BN") %in% colnames(timedTable)))
})


test_that("contrastInfo reports the engine, the groups, their sizes and the distribution of the statistic", {

  fit <- fitWith("edgeR")
  results <- testRegions(fit, contrast = list(shrVsBn = exampleContrast, bnVsShr = c("condition", "BN", "SHR")), verbose = FALSE)

  infoTable <- contrastInfo(results)
  groupSizes <- table(as.character(SummarizedExperiment::colData(fit@counts)$condition))

  expect_identical(infoTable$contrast, c("shrVsBn", "bnVsShr"))
  expect_identical(infoTable$engine, c("edgeR", "edgeR"))
  expect_identical(infoTable$group1, c("SHR", "BN"))
  expect_identical(infoTable$n.group1, as.integer(groupSizes[c("SHR", "BN")]))
  expect_identical(infoTable$n.group2, as.integer(groupSizes[c("BN", "SHR")]))
  expect_identical(infoTable$stat.distribution, c("f", "f"))
  expect_identical(infoTable$df1, c(1, 1))

  # The sizes are kept inside the object as well
  expect_identical(as.integer(results$shrVsBn@contrast.groups$n.samples), as.integer(groupSizes[c("SHR", "BN")]))

  # A numeric contrast over the coefficients still finds its two groups when it is a difference of two levels
  coefficientResults <- testRegions(fit, contrast = c(0, 1), verbose = FALSE)
  expect_false(is.na(contrastInfo(coefficientResults)$n.group1))

  expect_error(contrastInfo(data.frame()), "must be a RegionSetDE.results")
})


test_that("the fold change cut-off of the labels is drawn on the volcano and on the MA plot", {

  fit <- fitWith("edgeR")
  withCutoff <- testRegions(fit, contrast = exampleContrast, log2FC = 1, verbose = FALSE)
  withoutCutoff <- testRegions(fit, contrast = exampleContrast, verbose = FALSE)

  # ggplot2 4 names the layers of a plot, and unlist() would carry those names into the values
  lineIntercepts <- function(plotObject, geomClass, aesthetic) {
    layerData <- lapply(plotObject$layers, function(layer) {if (inherits(layer$geom, geomClass)) {layer$data[[aesthetic]]} else {NULL}})
    interceptValues <- unlist(layerData, use.names = FALSE)
    if (is.null(interceptValues)) {numeric(0)} else {sort(as.numeric(interceptValues))}
  }

  expect_identical(lineIntercepts(plotVolcano(withCutoff), "GeomVline", "xintercept"), c(-1, 1))
  expect_length(lineIntercepts(plotVolcano(withoutCutoff), "GeomVline", "xintercept"), 0)

  # The MA plot draws zero as a solid line in any case, and the cut-off on either side of it
  expect_identical(lineIntercepts(plotResultsMA(withCutoff), "GeomHline", "yintercept"), c(-1, 0, 1))
  expect_identical(lineIntercepts(plotResultsMA(withoutCutoff), "GeomHline", "yintercept"), 0)
})
