# The effective library sizes of a fit, rebuilt from the object as
# edgeR::scaleOffset rebuilds them.
effectiveSizes <- function(counts, offsetMatrix) {
    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, "counts"))
    librarySizes <- SummarizedExperiment::colData(counts)$library.size
    if (is.null(librarySizes) || any(is.na(librarySizes))) {
        librarySizes <- colSums(countMatrix)
    }
    exp(offsetMatrix - rowMeans(offsetMatrix) + mean(log(librarySizes)))
}


test_that("voom reads the scaling factors as effective library sizes", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "TMM", verbose = FALSE)
    counts <- filterRegions(counts, verbose = FALSE)
    fit <- fitRegions(counts, design = ~condition, engine = "voom", verbose = FALSE)

    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, "counts"))
    logFactors <- log(SummarizedExperiment::colData(counts)$scaling.factor)
    offsetMatrix <- matrix(rep(logFactors, each = nrow(countMatrix)), nrow = nrow(countMatrix))
    expectedValues <- log2((countMatrix + 0.5) / (effectiveSizes(counts, offsetMatrix) + 1) * 1e6)

    expect_equal(unname(fitComponents(fit)$voom$E), unname(expectedValues), tolerance = 1e-8)
})


test_that("voom changes with the normalisation", {
    counts <- filterRegions(exampleCounts(), verbose = FALSE)

    tmmFit <- fitRegions(normalizeRegionCounts(counts, method = "TMM", verbose = FALSE),
        design = ~condition, engine = "voom", verbose = FALSE
    )
    sizeFit <- fitRegions(normalizeRegionCounts(counts, method = "librarySize", verbose = FALSE),
        design = ~condition, engine = "voom", verbose = FALSE
    )

    expect_false(isTRUE(all.equal(fitComponents(tmmFit)$voom$E, fitComponents(sizeFit)$voom$E)))
})


test_that("the voom written for one size per observation matches limma::voom", {
    counts <- exampleCounts()
    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, "counts"))
    designMatrix <- stats::model.matrix(~condition, data = as.data.frame(SummarizedExperiment::colData(counts)))
    librarySizes <- colSums(countMatrix) * seq(0.8, 1.2, length.out = ncol(countMatrix))

    dgeList <- edgeR::DGEList(counts = countMatrix, lib.size = librarySizes)
    logLibraries <- RegionSetDE:::.effectiveLogLibraries(librarySizes = librarySizes, nRows = nrow(countMatrix))

    perSample <- RegionSetDE:::.voomOffsets(dgeList, designMatrix, logLibraries, perObservation = FALSE)
    perObservation <- RegionSetDE:::.voomOffsets(dgeList, designMatrix, logLibraries, perObservation = TRUE)

    expect_equal(unname(perObservation$E), unname(perSample$E), tolerance = 1e-8)
    expect_equal(unname(perObservation$weights), unname(perSample$weights), tolerance = 1e-6)
})


test_that("voom carries the loess offsets of every region", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "loess", verbose = FALSE)
    counts <- filterRegions(counts, verbose = FALSE)
    fit <- fitRegions(counts, design = ~condition, engine = "voom", verbose = FALSE)

    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, "counts"))
    offsetMatrix <- as.matrix(SummarizedExperiment::assay(counts, "offset"))
    expectedValues <- log2((countMatrix + 0.5) / (effectiveSizes(counts, offsetMatrix) + 1) * 1e6)

    expect_equal(unname(fitComponents(fit)$voom$E), unname(expectedValues), tolerance = 1e-8)
    expect_true(all(is.finite(fitComponents(fit)$voom$weights)))
    expect_true(all(fitComponents(fit)$voom$weights > 0))
})


test_that("dream refuses the loess offsets", {
    skip_if_not_installed("variancePartition")

    counts <- normalizeRegionCounts(exampleCounts(), method = "loess", verbose = FALSE)
    counts <- filterRegions(counts, verbose = FALSE)
    SummarizedExperiment::colData(counts)$pair <- rep_len(c("a", "b"), ncol(counts))

    expect_error(
        fitRegions(counts, design = ~condition, engine = "dream", random = ~ (1 | pair), verbose = FALSE),
        "loess"
    )
})


test_that("a new normalisation removes the loess offsets", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "loess", verbose = FALSE)
    counts <- normalizeRegionCounts(counts, method = "TMM", verbose = FALSE)

    expect_false("offset" %in% SummarizedExperiment::assayNames(counts))
})
