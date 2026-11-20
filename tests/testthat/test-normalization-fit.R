test_that("normalizeRegionCounts adds a normalized assay and positive factors", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "background", verbose = FALSE)

    expect_true("norm.counts" %in% SummarizedExperiment::assayNames(counts))

    scalingFactors <- SummarizedExperiment::colData(counts)$scaling.factor
    expect_length(scalingFactors, ncol(counts))
    expect_true(all(scalingFactors > 0))
    expect_true(all(is.finite(scalingFactors)))
})


test_that("the normalization methods do not all give the same factors", {
    counts <- exampleCounts()

    backgroundFactors <- SummarizedExperiment::colData(
        normalizeRegionCounts(counts, method = "background", verbose = FALSE)
    )$scaling.factor

    tmmFactors <- SummarizedExperiment::colData(
        normalizeRegionCounts(counts, method = "TMM", verbose = FALSE)
    )$scaling.factor

    expect_length(tmmFactors, ncol(counts))
    expect_false(isTRUE(all.equal(backgroundFactors, tmmFactors)))
})


test_that("normalizeRegionCounts rejects a wrong number of supplied factors", {
    expect_error(
        normalizeRegionCounts(exampleCounts(), scalingFactors = c(1, 1), verbose = FALSE)
    )
})


test_that("filterRegions removes rows without emptying any set", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "background", verbose = FALSE)
    filtered <- filterRegions(counts, verbose = FALSE)

    expect_s4_class(filtered, "RegionSetDE.counts")
    expect_lt(nrow(filtered), nrow(counts))
    expect_equal(ncol(filtered), ncol(counts))

    # Every set must survive, otherwise the downstream set tests break
    expect_setequal(
        unique(as.character(SummarizedExperiment::rowData(filtered)$region.set)),
        unique(as.character(SummarizedExperiment::rowData(counts)$region.set))
    )
})


test_that("a stricter fold change keeps fewer regions", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "background", verbose = FALSE)

    loose <- filterRegions(counts, foldChange = 1.5, verbose = FALSE)
    strict <- filterRegions(counts, foldChange = 4, verbose = FALSE)

    expect_lte(nrow(strict), nrow(loose))
})


test_that("fitRegions returns a fit whose accessors work", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "background", verbose = FALSE)
    counts <- filterRegions(counts, verbose = FALSE)

    fit <- fitRegions(counts, design = ~condition, engine = "edgeR", verbose = FALSE)

    expect_s4_class(fit, "RegionSetDE.fit")
    expect_s4_class(fitCounts(fit), "RegionSetDE.counts")
    expect_false(is.null(fitObject(fit)))
    expect_equal(nrow(fitCounts(fit)), nrow(counts))
})


test_that("fitRegions errors on a design variable that is not in colData", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "background", verbose = FALSE)

    expect_error(
        fitRegions(counts, design = ~notAVariable, engine = "edgeR", verbose = FALSE)
    )
})


test_that("the packaged fit matches what the script produced", {
    fit <- exampleFit()

    expect_s4_class(fit, "RegionSetDE.fit")
    expect_setequal(
        regionSetNames(fit),
        c("promoterCpG", "promoterNonCpG", "geneBody", "intergenic")
    )
    expect_equal(ncol(fitCounts(fit)), 4)
})


test_that("estimateNullDispersion works from the background", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "background", verbose = FALSE)
    nullDispersion <- estimateNullDispersion(counts,
        source = "background",
        verbose = FALSE
    )

    expect_false(is.null(nullDispersion))
})


test_that("estimateNullDispersion refuses a set with too little signal", {
    counts <- normalizeRegionCounts(exampleCounts(), method = "background", verbose = FALSE)

    # The intergenic set is the low-signal control by construction
    expect_error(
        estimateNullDispersion(counts,
            source = "regionSet",
            regionSets = "intergenic", verbose = FALSE
        )
    )
})


test_that("makeSetUniverse returns a universe object", {
    universe <- makeSetUniverse(exampleFit(), verbose = FALSE)
    expect_s4_class(universe, "RegionSetDE.universe")
})
