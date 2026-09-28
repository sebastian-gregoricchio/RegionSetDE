## The input of fry, the overlaps inside a set and the fold change of a tiled region.


test_that("fry runs on the NB z-scores of an edgeR fit, as edgeR's own fry does", {

  fit <- exampleFit()
  expect_identical(fit@engine, "edgeR")

  setResults <- testRegionSets(fit, contrast = exampleContrast(), verbose = FALSE)
  setTable <- resultsTable(setResults)
  expect_identical(setResults@parameters$testRegionSets$fryInput, "NB z-scores")

  # edgeR refits the null model for every call, the package computes it once: the p-values must agree
  contrastObject <- RegionSetDE:::.resolveContrast(contrast = exampleContrast(),
                                                   design = fit@design,
                                                   colData = SummarizedExperiment::colData(fit@counts))
  setColumn <- as.character(SummarizedExperiment::rowData(fit@counts)$region.set)

  for (setName in setTable$region.set) {
    edgeRtable <- limma::fry(fit@fit$dge,
                             index = list(set = which(setColumn == setName)),
                             design = fit@design,
                             contrast = contrastObject$vector,
                             sort = FALSE)
    expect_equal(setTable$fry.p[setTable$region.set == setName], edgeRtable$PValue[1], tolerance = 1e-8)
  }
})


test_that("fryInput = 'logCPM' gives back the test on the log-CPM matrix", {

  fit <- exampleFit()
  logResults <- testRegionSets(fit, contrast = exampleContrast(), fryInput = "logCPM", verbose = FALSE)
  logTable <- resultsTable(logResults)
  expect_identical(logResults@parameters$testRegionSets$fryInput, "logCPM")

  contrastObject <- RegionSetDE:::.resolveContrast(contrast = exampleContrast(),
                                                   design = fit@design,
                                                   colData = SummarizedExperiment::colData(fit@counts))
  logMatrix <- RegionSetDE:::.expressionMatrix(fit = fit)
  setColumn <- as.character(SummarizedExperiment::rowData(fit@counts)$region.set)

  setName <- logTable$region.set[1]
  limmaTable <- limma::fry(logMatrix,
                           index = list(set = which(setColumn == setName)),
                           design = fit@design,
                           contrast = contrastObject$vector,
                           sort = FALSE)
  expect_equal(logTable$fry.p[1], limmaTable$PValue[1], tolerance = 1e-10)

  # The competitive test and the effect sizes do not depend on what fry reads
  zTable <- resultsTable(exampleSetResults())
  zTable <- zTable[match(logTable$region.set, zTable$region.set), ]
  expect_equal(logTable$camera.p, zTable$camera.p)
  expect_equal(logTable$delta.log2FC, zTable$delta.log2FC)

  expect_error(testRegionSets(fit, contrast = exampleContrast(), fryInput = "counts", verbose = FALSE), "fryInput")
})


test_that("the DESeq2 engine gets z-scores from its own dispersions and factors", {

  skip_if_not_installed("DESeq2")
  counts <- filterRegions(exampleCounts(), method = "abundance", verbose = FALSE)
  deseqFit <- suppressWarnings(fitRegions(counts, design = ~ condition, engine = "deseq2", verbose = FALSE))

  zValues <- RegionSetDE:::.fryValues(fit = deseqFit, contrastVector = c(0, 1), fryInput = "auto")
  expect_identical(zValues$standardize, "none")
  expect_identical(dim(zValues$values), dim(SummarizedExperiment::assay(deseqFit@counts, 1)))
  expect_true(all(is.finite(zValues$values)))

  # Under the null model the z-scores are centred and close to unit variance
  expect_lt(abs(mean(zValues$values)), 0.1)
  expect_gt(stats::sd(as.vector(zValues$values)), 0.5)
  expect_lt(stats::sd(as.vector(zValues$values)), 1.5)

  setResults <- testRegionSets(deseqFit, contrast = exampleContrast(), verbose = FALSE)
  expect_true(all(resultsTable(setResults)$fry.p >= 0 & resultsTable(setResults)$fry.p <= 1))
})


test_that("the linear engines keep their own log values for fry", {

  voomFit <- fitRegions(filterRegions(exampleCounts(), method = "abundance", verbose = FALSE),
                        design = ~ condition, engine = "voom", verbose = FALSE)

  fryList <- RegionSetDE:::.fryValues(fit = voomFit, contrastVector = c(0, 1), fryInput = "auto")
  expect_identical(fryList$standardize, "posterior.sd")
  expect_identical(fryList$values, RegionSetDE:::.expressionMatrix(fit = voomFit))
})


test_that("z-scores of the tiles are pooled so a region keeps the scale of one z-score", {

  zMatrix <- matrix(c(1, 1, 1, 1, 2, 2), ncol = 2, byrow = TRUE)
  tileMap <- c(1L, 1L, 2L)

  pooled <- RegionSetDE:::.collapseTileMatrix(expressionMatrix = zMatrix, tileMap = tileMap, pooling = "stouffer")
  expect_equal(pooled[1, ], c(2 / sqrt(2), 2 / sqrt(2)), ignore_attr = TRUE)
  expect_equal(pooled[2, ], c(2, 2), ignore_attr = TRUE)

  averaged <- RegionSetDE:::.collapseTileMatrix(expressionMatrix = zMatrix, tileMap = tileMap)
  expect_equal(averaged[1, ], c(1, 1), ignore_attr = TRUE)
})


test_that("regions overlapping inside a set are counted, warned about or refused", {

  fit <- exampleFit()
  regionRows <- SummarizedExperiment::rowRanges(fit@counts)
  setColumn <- as.character(SummarizedExperiment::rowData(fit@counts)$region.set)

  # The second promoter is moved onto the first, half of it now covering the same bases
  promoterRows <- which(setColumn == "promoterCpG")[1:2]
  movedRanges <- regionRows
  GenomicRanges::ranges(movedRanges)[promoterRows[2]] <- IRanges::IRanges(start = BiocGenerics::start(regionRows)[promoterRows[1]] + 100L,
                                                                           width = BiocGenerics::width(regionRows)[promoterRows[1]])
  GenomeInfoDb::seqnames(movedRanges)[promoterRows[2]] <- GenomeInfoDb::seqnames(regionRows)[promoterRows[1]]
  SummarizedExperiment::rowRanges(fit@counts) <- movedRanges

  expect_warning(setResults <- testRegionSets(fit, contrast = exampleContrast(), verbose = FALSE),
                 "promoterCpG \\(2\\)")
  setTable <- resultsTable(setResults)
  expect_identical(setTable$n.overlapping.within[setTable$region.set == "promoterCpG"], 2L)
  expect_true(all(setTable$n.overlapping.within[setTable$region.set != "promoterCpG"] == 0L))

  expect_error(testRegionSets(fit, contrast = exampleContrast(), overlapWithinSet = "stop", verbose = FALSE),
               "overlap another region of the same set")

  expect_no_warning(testRegionSets(fit, contrast = exampleContrast(), overlapWithinSet = "allow", verbose = FALSE))
  expect_error(testRegionSets(fit, contrast = exampleContrast(), overlapWithinSet = "merge", verbose = FALSE),
               "overlapWithinSet")

  # The pair test counts each side on its own
  expect_warning(pairResults <- testSetContrast(fit, contrast = exampleContrast(),
                                                set1 = "promoterCpG", set2 = "intergenic", verbose = FALSE),
                 "promoterCpG")
  expect_identical(resultsTable(pairResults)$n.overlapping.within.1, 2L)
  expect_identical(resultsTable(pairResults)$n.overlapping.within.2, 0L)
})


test_that("a set without overlaps reports zero and stays silent", {

  expect_no_warning(setResults <- exampleSetResults())
  expect_true(all(resultsTable(setResults)$n.overlapping.within == 0L))
})


test_that("a combined tiled region carries the mean fold change of its tiles, weighted by width", {

  tileTable <- data.frame(region.set = "a",
                          region.id = rep(c("r1", "r2"), c(3, 2)),
                          region.key = rep(c("a|r1", "a|r2"), c(3, 2)),
                          log2FC = c(2, 0, 0, -1, 1),
                          average.signal = 5,
                          stat = c(4, 0.1, 0.1, -2, 2),
                          p.value = c(1e-4, 0.9, 0.9, 0.001, 0.5),
                          stringsAsFactors = FALSE)
  tileRanges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = c(1, 101, 201, 1001, 1101),
                                                                width = c(100, 100, 50, 100, 300)))

  combinedTable <- suppressMessages(RegionSetDE:::.combineTiles(tileTable = tileTable, tileRanges = tileRanges))$results

  # The reported fold change is the one of the tile carrying the p-value, the mean covers the region
  expect_equal(combinedTable$log2FC, c(2, -1))
  expect_equal(combinedTable$mean.tile.log2FC, c(2 * 100 / 250, (-1 * 100 + 1 * 300) / 400))
  expect_equal(combinedTable$rep.tile.start, c(1, 1001))
})
