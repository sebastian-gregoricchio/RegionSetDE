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


test_that("the Jaccard index compares the peak calls of the samples", {

  peakList <- GenomicRanges::GRangesList(
    a = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 1000, 5000), width = 100)),
    b = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(150, 1000, 9000), width = 100)),
    c = GenomicRanges::GRanges("chr1", IRanges::IRanges(20000, width = 100)))

  # Places: 100-249 (a, b), 1000 (a, b), 5000 (a), 9000 (b), 20000 (c)
  regionJaccard <- computeSampleCorrelation(peakList, method = "jaccard")$correlation
  expect_equal(regionJaccard["a", "b"], 2 / 4)
  expect_equal(regionJaccard["a", "c"], 0)
  expect_equal(unname(diag(regionJaccard)), rep(1, 3))

  # Base pairs: 50 + 100 shared, 150 + 100 + 100 + 100 in the union
  basepairJaccard <- computeSampleCorrelation(peakList, method = "jaccard", jaccardLevel = "basepair")$correlation
  expect_equal(basepairJaccard["a", "b"], 150 / 450)
  expect_true(isSymmetric(basepairJaccard))

  expect_error(computeSampleCorrelation(exampleCounts(), method = "jaccard"), "compares peak calls")
  expect_error(computeSampleCorrelation(peakList, method = "jaccard", jaccardLevel = "window"), "'region' or 'basepair'")
})


test_that("the Jaccard heatmap takes the peaks and the annotation of a consensus", {

  testthat::skip_if_not_installed("consensusRegions")

  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)

  peakJaccard <- computeSampleCorrelation(consensus, method = "jaccard")
  expect_identical(colnames(peakJaccard$correlation), sampleSheet$sample)
  expect_identical(peakJaccard$samples$condition, sampleSheet$condition)
  expect_true(all(peakJaccard$correlation >= 0 & peakJaccard$correlation <= 1))

  # Replicates of a condition share more peaks with each other than with the vehicle
  expect_gt(peakJaccard$correlation["AR_R1881_24h_r1", "AR_R1881_24h_r2"], peakJaccard$correlation["AR_R1881_24h_r1", "AR_DMSO_r1"])

  jaccardHeatmap <- plotSampleCorrelation(consensus, method = "jaccard", groupBy = "condition")
  expect_s4_class(jaccardHeatmap, "Heatmap")
  expect_s4_class(plotSampleCorrelation(peakJaccard), "Heatmap")
})


test_that("the ordination and the correlation can be restricted to some samples", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)
  sampleNames <- colnames(counts)

  # The normalisation of the whole analysis is kept, so a subset of the matrix is a subset of the samples
  fullCorrelation <- computeSampleCorrelation(counts, method = "pearson")$correlation
  subsetCorrelation <- computeSampleCorrelation(counts, method = "pearson", samples = sampleNames[c(1, 2, 4)])$correlation
  expect_equal(subsetCorrelation, fullCorrelation[c(1, 2, 4), c(1, 2, 4)])

  # Names, positions and a logical vector select the same samples
  expect_identical(computeSamplePCA(counts, samples = 1:3, topRegions = 500)$scores$sample, sampleNames[1:3])
  expect_identical(computeSamplePCA(counts, samples = c(TRUE, TRUE, TRUE, FALSE), topRegions = 500)$scores$sample, sampleNames[1:3])

  pcaPlot <- plotRegionPCA(counts, samples = sampleNames[2:4], colourBy = "condition", topRegions = 500)
  expect_setequal(attr(pcaPlot, "pca")$sample, sampleNames[2:4])

  correlationHeatmap <- plotSampleCorrelation(counts, samples = 2:4, groupBy = "condition")
  expect_identical(nrow(correlationHeatmap@matrix), 3L)

  expect_error(computeSamplePCA(counts, samples = 1:2), "three samples")
  expect_error(computeSampleCorrelation(counts, samples = "absentSample"), "absent from the object")
  expect_error(computeSampleCorrelation(counts, samples = c(TRUE, FALSE)), "one TRUE or FALSE per sample")
})


test_that("the Jaccard index can be restricted to some samples", {

  testthat::skip_if_not_installed("consensusRegions")

  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)

  treatedSamples <- sampleSheet$sample[sampleSheet$treatment == "R1881"]
  treatedJaccard <- computeSampleCorrelation(consensus, method = "jaccard", samples = treatedSamples)

  expect_identical(rownames(treatedJaccard$correlation), treatedSamples)
  expect_identical(treatedJaccard$samples$sample, treatedSamples)

  # The places are the union of the peaks of the samples kept, so the index is the one of those samples alone
  pairJaccard <- computeSampleCorrelation(consensus, method = "jaccard", samples = treatedSamples[1:2])$correlation
  expect_equal(pairJaccard[1, 2], treatedJaccard$correlation[1, 2])
})


test_that("selectSamples keeps the files, the layouts and the greenlist counts in line with the samples", {

  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
  peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), verbose = FALSE)
  counts <- countReads(peakRegions, sampleSheet = sampleSheet, countInput = FALSE, verbose = FALSE)
  counts <- countGreenlist(counts, greenlist = GenomicRanges::GRanges("19", IRanges::IRanges(seq(46.5e6, 57.5e6, by = 1e6), width = 2e5)),
                           excludeCounted = FALSE, verbose = FALSE)

  selectedCounts <- selectSamples(counts, samples = c(1, 4, 7), verbose = FALSE)

  expect_identical(selectedCounts@parameters$countReads$bamFiles, sampleSheet$bam[c(1, 4, 7)])
  expect_length(selectedCounts@parameters$countReads$fragmentLength, 3)
  expect_identical(colnames(S4Vectors::metadata(selectedCounts)$greenlist), colnames(selectedCounts))

  # Every step reading the files again works on the selection
  expect_s3_class(libInfo(selectedCounts), "data.frame")
  expect_s4_class(normalizeCounts(selectedCounts, method = "greenlist", verbose = FALSE), "RegionSetDE.counts")
  expect_s4_class(countBackground(selectedCounts, binSize = 10000, verbose = FALSE), "RegionSetDE.counts")
})
