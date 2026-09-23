test_that("the components, the variance and the annotation come back together", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)
  samplePCA <- computeSamplePCA(counts, topRegions = 500, verbose = FALSE)

  expect_named(samplePCA, c("scores", "variance", "loadings", "parameters"))
  expect_identical(samplePCA$scores$sample, colnames(counts))
  expect_true(all(c("PC1", "PC2", "condition", "sex") %in% colnames(samplePCA$scores)))

  # Centring costs one degree of freedom, so four samples leave three components carrying the whole variance
  expect_identical(nrow(samplePCA$variance), ncol(counts) - 1L)
  expect_equal(sum(samplePCA$variance$variance.percent), 100)
  expect_equal(samplePCA$variance$cumulative.percent[nrow(samplePCA$variance)], 100)

  expect_identical(dim(samplePCA$loadings), c(500L, ncol(counts) - 1L))
  expect_identical(samplePCA$parameters$n.regions, 500L)
  expect_true(samplePCA$parameters$useOffsets)
})


test_that("the ordination moves with the normalisation and the correlation does not", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  normalisedPCA <- computeSamplePCA(counts, topRegions = 500, verbose = FALSE)
  rawPCA <- computeSamplePCA(counts, topRegions = 500, useOffsets = FALSE, verbose = FALSE)

  # The rows are chosen once on the normalised values, so the two ordinations differ by the transformation alone
  expect_identical(rownames(normalisedPCA$loadings), rownames(rawPCA$loadings))
  expect_false(isTRUE(all.equal(normalisedPCA$scores$PC1, rawPCA$scores$PC1)))

  # A single factor per sample shifts the log values by a constant, which no correlation sees
  normalisedCorrelation <- computeSampleCorrelation(counts, method = "pearson", verbose = FALSE)
  rawCorrelation <- computeSampleCorrelation(counts, method = "pearson", useOffsets = FALSE, verbose = FALSE)

  expect_equal(normalisedCorrelation$correlation, rawCorrelation$correlation)
})


test_that("the correlation comes back with the samples it was computed on", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)
  sampleCorrelation <- computeSampleCorrelation(counts, topRegions = 400, verbose = FALSE)

  correlationMatrix <- sampleCorrelation$correlation
  expect_identical(dim(correlationMatrix), rep(ncol(counts), 2))
  expect_identical(colnames(correlationMatrix), colnames(counts))
  expect_equal(correlationMatrix, t(correlationMatrix))
  expect_equal(diag(correlationMatrix), rep(1, ncol(counts)), ignore_attr = TRUE)

  expect_identical(sampleCorrelation$samples$sample, colnames(counts))
  expect_identical(sampleCorrelation$parameters$method, "spearman")
  expect_identical(sampleCorrelation$parameters$n.regions, 400L)

  expect_error(computeSampleCorrelation(counts, method = "cosine", verbose = FALSE), "'spearman'")
  expect_error(computeSampleCorrelation(selectSamples(counts, condition == "BN" & sex == "male",
                                                      verbose = FALSE), verbose = FALSE),
               "At least two samples")
})


test_that("an object without normalisation is scaled by the library sizes and says so", {

  counts <- exampleCounts()

  expect_message(computeSamplePCA(counts, topRegions = 200), "library sizes alone")
  expect_message(computeSampleCorrelation(counts, topRegions = 200), "library sizes alone")
  expect_silent(computeSampleCorrelation(counts, topRegions = 200, useOffsets = FALSE))
})


test_that("the plots take a computation done beforehand", {

  skip_if_not_installed("ComplexHeatmap")
  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  samplePCA <- computeSamplePCA(counts, topRegions = 500, verbose = FALSE)
  pcaPlot <- plotRegionPCA(samplePCA, colourBy = "condition")

  expect_s3_class(pcaPlot, "ggplot")
  expect_equal(attr(pcaPlot, "pca")$x.value, samplePCA$scores$PC1)

  sampleCorrelation <- computeSampleCorrelation(counts, method = "pearson", verbose = FALSE)
  correlationPlot <- plotSampleCorrelation(sampleCorrelation, groupBy = "condition")

  expect_s4_class(correlationPlot, "Heatmap")
  expect_equal(correlationPlot@matrix, sampleCorrelation$correlation)
  expect_match(correlationPlot@column_title, "within")
})


test_that("the annotation bars carry the columns and the colours asked for", {

  skip_if_not_installed("ComplexHeatmap")
  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  correlationPlot <- plotSampleCorrelation(counts,
                                           annotationColumns = c("condition", "sex", "biologicalReplicate"),
                                           annotationColours = list(sex = c(male = "grey40", female = "orange")),
                                           title = "Sample correlation",
                                           verbose = FALSE)

  annotationNames <- names(correlationPlot@top_annotation@anno_list)
  expect_identical(annotationNames, c("condition", "sex", "biologicalReplicate"))
  expect_match(correlationPlot@column_title, "Sample correlation")

  # The colours given are kept, the others come from the palette of the package
  sexColours <- correlationPlot@top_annotation@anno_list$sex@color_mapping@colors
  expect_identical(grDevices::col2rgb(sexColours[["male"]]), grDevices::col2rgb("grey40"))

  expect_error(plotSampleCorrelation(counts, groupBy = "absent", verbose = FALSE), "absent from the colData")
  expect_error(plotSampleCorrelation(counts, annotationColumns = "absent", verbose = FALSE), "absent from the colData")
})


test_that("the diagonal can be left out of the drawing and of the scale", {

  skip_if_not_installed("ComplexHeatmap")
  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  correlationPlot <- plotSampleCorrelation(counts, excludeDiagonal = TRUE, cluster = FALSE, verbose = FALSE)

  expect_true(all(is.na(diag(correlationPlot@matrix))))
  expect_false(any(is.na(correlationPlot@matrix[row(correlationPlot@matrix) != col(correlationPlot@matrix)])))

  expect_message(plotSampleCorrelation(counts, limits = c(0.99, 1), verbose = FALSE),
                 "drawn at its ends")
})
