# Two libraries from the same BAM file, so that the normalisation has more than one factor to work with
twoSampleCounts <- function(...) {
  RegionSetDE::countReads(toyRegionSet(),
                          bamFiles = rep(toyBamFile(), 2),
                          sampleNames = c("first", "second"),
                          sampleMetadata = data.frame(sample = c("first", "second"),
                                                      condition = c("A", "B")),
                          verbose = FALSE,
                          ...)
}


test_that("the wide table holds the raw counts, the coordinates and the annotation", {

  counts <- exampleCounts()
  wideTable <- countTable(counts)

  expect_s3_class(wideTable, "data.frame")
  expect_identical(nrow(wideTable), nrow(counts))
  expect_identical(rownames(wideTable), rownames(counts))
  expect_true(all(c("region.set", "region.id", "seqnames", "start", "end", "width", "regionId") %in% colnames(wideTable)))
  expect_false(any(c("tile.id", "n.tiles") %in% colnames(wideTable)))

  expect_equal(as.matrix(wideTable[, colnames(counts)]),
               as.matrix(SummarizedExperiment::assay(counts, "counts")))

  expect_false("regionId" %in% colnames(countTable(counts, extraColumns = FALSE)))
})


test_that("the matrix format returns the assay itself", {

  counts <- exampleCounts()

  expect_equal(countTable(counts, format = "matrix"),
               as.matrix(SummarizedExperiment::assay(counts, "counts")))
})


test_that("the long format has one row per region and sample, with the colData attached", {

  counts <- exampleCounts()
  longTable <- countTable(counts, format = "long")

  expect_identical(nrow(longTable), nrow(counts) * ncol(counts))
  expect_true(all(c("sample", "counts", "condition", "library.size") %in% colnames(longTable)))
  expect_identical(unique(longTable$sample), colnames(counts))

  # Taken out of the filter, since 'counts' is also a column of the long table
  firstSampleName <- colnames(counts)[1]
  firstSample <- dplyr::filter(longTable, .data$sample == firstSampleName)
  expect_equal(firstSample$counts, as.numeric(SummarizedExperiment::assay(counts, "counts")[, 1]))
  expect_true(all(firstSample$condition == counts$condition[1]))
})


test_that("the normalised values require a normalisation and match its assay", {

  counts <- exampleCounts()
  expect_error(countTable(counts, normalized = TRUE), "normalizeCounts")

  counts <- normalizeCounts(counts, method = "background", verbose = FALSE)

  expect_equal(countTable(counts, normalized = TRUE, format = "matrix"),
               as.matrix(SummarizedExperiment::assay(counts, "norm.counts")))
  expect_true("norm.counts" %in% colnames(countTable(counts, normalized = TRUE, format = "long")))
})


test_that("the region sets can be restricted and an absent one is refused", {

  counts <- exampleCounts()
  setTable <- countTable(counts, set = "promoterCpG")

  expect_true(all(setTable$region.set == "promoterCpG"))
  expect_identical(nrow(setTable), sum(SummarizedExperiment::rowData(counts)$region.set == "promoterCpG"))
  expect_error(countTable(counts, set = "notASet"), "absent")
})


test_that("the arguments are checked", {

  counts <- exampleCounts()

  expect_error(countTable(counts, level = "tile"), "not tiled")
  expect_error(countTable(counts, level = "gene"), "level")
  expect_error(countTable(counts, format = "tidy"), "format")
  expect_error(countTable(counts, normalized = NA), "normalized")
  expect_error(countTable(toyCounts(tileWidth = 100L), tileSummary = "median", verbose = FALSE), "tileSummary")
})


test_that("the tile level returns the tiles as they were counted", {

  tiledCounts <- toyCounts(tileWidth = 100L)
  tileTableOut <- countTable(tiledCounts, level = "tile")

  expect_identical(nrow(tileTableOut), nrow(tiledCounts))
  expect_true("tile.id" %in% colnames(tileTableOut))
  expect_equal(tileTableOut$example, as.numeric(SummarizedExperiment::assay(tiledCounts, "counts")[, 1]))
})


test_that("the tiles are summed into their region, and exceed the count of the region taken whole", {

  plainCounts <- toyCounts()
  tiledCounts <- toyCounts(tileWidth = 100L)

  expect_message(regionTable <- countTable(tiledCounts), "counted in both")
  expect_no_message(countTable(tiledCounts, verbose = FALSE))

  rowTable <- as.data.frame(SummarizedExperiment::rowData(tiledCounts))
  expectedSums <- rowsum(SummarizedExperiment::assay(tiledCounts, "counts"),
                         group = paste(rowTable$region.set, rowTable$region.id, sep = "|"),
                         reorder = FALSE)

  expect_identical(rownames(regionTable), rownames(plainCounts))
  expect_equal(regionTable$example, as.numeric(expectedSums[rownames(regionTable), 1]))
  expect_true(all(regionTable$n.tiles == 3))
  expect_false("tile.id" %in% colnames(regionTable))

  # The span of the tiles is the region, since 300 bp cut into 100 bp tiles leaves no remainder
  plainRanges <- SummarizedExperiment::rowRanges(plainCounts)
  expect_identical(regionTable$start, BiocGenerics::start(plainRanges))
  expect_identical(regionTable$end, BiocGenerics::end(plainRanges))

  # A fragment crossing a border between tiles is counted in both of them
  expect_true(all(regionTable$example >= SummarizedExperiment::assay(plainCounts, "counts")[, 1]))
})


test_that("the normalised region values equal the summed counts divided by the scaling factor", {

  tiledCounts <- twoSampleCounts(tileWidth = 100L)
  tiledCounts <- normalizeCounts(tiledCounts, method = "manual", scalingFactors = c(1, 3), verbose = FALSE)

  rawMatrix <- countTable(tiledCounts, format = "matrix", verbose = FALSE)
  normalizedMatrix <- countTable(tiledCounts, normalized = TRUE, format = "matrix", verbose = FALSE)
  scalingFactors <- tiledCounts$scaling.factor

  expect_equal(normalizedMatrix, sweep(rawMatrix, MARGIN = 2, STATS = scalingFactors, FUN = "/"))
})


test_that("the tiles stay with their region whatever the row order", {

  tiledCounts <- toyCounts(tileWidth = 100L)

  set.seed(1)
  shuffledCounts <- tiledCounts[sample(nrow(tiledCounts)), ]

  orderedTable <- countTable(tiledCounts, format = "matrix", verbose = FALSE)
  shuffledTable <- countTable(shuffledCounts, format = "matrix", verbose = FALSE)

  expect_equal(shuffledTable[rownames(orderedTable), , drop = FALSE], orderedTable)
})


test_that("a mean is weighted by the width of the tiles, maxima and minima are kept", {

  # 300 bp cut into 120 bp tiles leaves a trailing tile of 60 bp
  tiledCounts <- toyCounts(tileWidth = 120L)
  tileTableOut <- countTable(tiledCounts, level = "tile")
  regionKey <- paste(tileTableOut$region.set, tileTableOut$region.id, sep = "|")

  expect_true(all(table(regionKey) == 3))

  expectedMean <- tapply(tileTableOut$example * tileTableOut$width, regionKey, sum) / tapply(tileTableOut$width, regionKey, sum)
  expectedMax <- tapply(tileTableOut$example, regionKey, max)
  expectedMin <- tapply(tileTableOut$example, regionKey, min)

  meanTable <- countTable(tiledCounts, tileSummary = "mean", verbose = FALSE)
  maxTable <- countTable(tiledCounts, tileSummary = "max", verbose = FALSE)
  minTable <- countTable(tiledCounts, tileSummary = "min", verbose = FALSE)

  expect_equal(meanTable$example, as.numeric(expectedMean[rownames(meanTable)]))
  expect_equal(maxTable$example, as.numeric(expectedMax[rownames(maxTable)]))
  expect_equal(minTable$example, as.numeric(expectedMin[rownames(minTable)]))
})


test_that("the rule for bigWig signal follows the summary used when counting", {

  tiledCounts <- toyCounts(tileWidth = 100L)
  S4Vectors::metadata(tiledCounts)$signal.type <- "bigwig"
  tiledCounts@parameters$countBigwig <- list(summaryFunction = "max")

  expect_equal(countTable(tiledCounts, format = "matrix", verbose = FALSE),
               countTable(tiledCounts, format = "matrix", tileSummary = "max", verbose = FALSE))

  # Base-by-base signal adds up exactly, so there is nothing to warn about
  tiledCounts@parameters$countBigwig <- list(summaryFunction = "sum")
  expect_no_message(countTable(tiledCounts))
})


test_that("the fit and the results return the counts they carry", {

  fit <- exampleFit()
  results <- exampleResults()

  expect_equal(countTable(fit), countTable(fitCounts(fit)))
  expect_equal(countTable(results, normalized = TRUE), countTable(fitCounts(fit), normalized = TRUE))
  expect_equal(countTable(exampleSetResults(), format = "matrix"), countTable(fit, format = "matrix"))

  bareResults <- testRegions(fit, contrast = exampleContrast(), carryCounts = FALSE, verbose = FALSE)
  expect_error(countTable(bareResults), "carryCounts")
})


test_that("names shared between the rows and the samples do not produce duplicated columns", {

  counts <- exampleCounts()
  counts$width <- 1

  expect_message(longTable <- countTable(counts, format = "long"), "width")
  expect_true(all(c("width", "width.sample") %in% colnames(longTable)))
  expect_false(any(duplicated(colnames(longTable))))

  renamedCounts <- exampleCounts()
  colnames(renamedCounts)[1] <- "start"
  expect_error(countTable(renamedCounts), "format = 'long'")
})
