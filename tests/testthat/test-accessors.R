## Every slot of every class has an accessor, and the tables come back as data.frames.


slotTable <- function(object, slotName) {
    as.data.frame(methods::slot(object, slotName), optional = TRUE)
}


test_that("the provenance accessors return the slots of every object", {
    regions <- toyRegionSet()
    counts <- exampleCounts()
    fit <- exampleFit()
    results <- exampleResults()

    for (object in list(regions, counts, fit, results)) {
        expect_identical(blacklistRanges(object), methods::slot(object, "blacklist"))
        expect_identical(greylistRanges(object), methods::slot(object, "greylist"))
        expect_identical(whitelistRanges(object), methods::slot(object, "whitelist"))
        expect_identical(chromosomeStyle(object), methods::slot(object, "seqlevels.style"))
        expect_identical(parameterLog(object), methods::slot(object, "parameters"))
        expect_identical(filteringLog(object), slotTable(object, "filtering.log"))
    }

    expect_true(all(c("countReads") %in% names(parameterLog(counts))))
})


test_that("the tables are stored as DataFrame and read as data.frame", {
    results <- exampleResults()
    setResults <- exampleSetResults()

    expect_s4_class(methods::slot(results, "results"), "DataFrame")
    expect_s4_class(methods::slot(results, "filtering.log"), "DataFrame")
    expect_s4_class(methods::slot(setResults, "results"), "DataFrame")

    expect_s3_class(resultsTable(results), "data.frame")
    expect_identical(resultsTable(results), slotTable(results, "results"))
    expect_identical(resultsTable(setResults), slotTable(setResults, "results"))
    expect_identical(regionStats(setResults), slotTable(setResults, "regionStats"))

    # dplyr reads what the accessors return
    expect_s3_class(dplyr::filter(resultsTable(results), .data$FDR <= 1), "data.frame")
})


test_that("the fit and contrast accessors return their slots", {
    fit <- exampleFit()
    results <- exampleResults()
    setResults <- exampleSetResults()

    expect_identical(fitComponents(fit), methods::slot(fit, "fit"))
    expect_identical(fitObject(fit), fitComponents(fit)$object)
    expect_identical(fitEngine(fit), methods::slot(fit, "engine"))
    expect_identical(fitEngine(results), methods::slot(results, "engine"))
    expect_identical(fitEngine(setResults), methods::slot(setResults, "engine"))
    expect_identical(designMatrix(fit), methods::slot(fit, "design"))
    expect_identical(blockingInfo(fit), methods::slot(fit, "blocking"))
    expect_identical(fitSamples(fit), methods::slot(fit, "samples"))
    expect_identical(setUniverse(fit), methods::slot(fit, "universe"))
    expect_identical(setUniverse(setResults), methods::slot(setResults, "universe"))
    expect_true(is.null(designFormula(fit)) || methods::is(designFormula(fit), "formula"))

    expect_identical(contrastVector(results), methods::slot(results, "contrast.vector"))
    expect_identical(contrastGroups(results), methods::slot(results, "contrast.groups"))
    expect_identical(contrastGroups(setResults), methods::slot(setResults, "contrast.groups"))
    expect_identical(thresholdInfo(results), methods::slot(results, "thresholds"))
    expect_identical(thresholdInfo(setResults), methods::slot(setResults, "thresholds"))
    expect_identical(tileCombination(results), methods::slot(results, "combination"))
    expect_identical(regionRanges(results), methods::slot(results, "regions"))

    expect_identical(setTestType(setResults), "set")
    expect_identical(setTestMethods(setResults), methods::slot(setResults, "methods"))
})


test_that("the universe accessors return their slots", {
    universe <- setUniverse(exampleFit())

    expect_identical(universeIndex(universe), methods::slot(universe, "index"))
    expect_identical(universeType(universe), methods::slot(universe, "type"))
    expect_identical(universeMatching(universe), methods::slot(universe, "matching"))
    expect_identical(comparisonSets(universe), methods::slot(universe, "comparison.sets"))
    expect_identical(universeDiagnostics(universe), slotTable(universe, "diagnostics"))
    expect_identical(universeRows(universe), methods::slot(universe, "n.rows"))
})


test_that("the score accessors return their slots", {
    setScores <- scoreRegionSets(exampleCounts(), verbose = FALSE)

    expect_identical(scoreComparisons(setScores), slotTable(setScores, "comparisons"))
    expect_identical(scoreSamples(setScores), slotTable(setScores, "sample.metadata"))
    expect_identical(scoreReference(setScores), methods::slot(setScores, "reference"))
    expect_identical(scoreAssay(setScores), methods::slot(setScores, "assay"))
    expect_identical(scoreSummary(setScores), methods::slot(setScores, "summary"))
    expect_identical(scorePerBasepair(setScores), methods::slot(setScores, "per.basepair"))
})


test_that("a list of results gives back all of its contrasts", {
    resultsList <- testRegions(exampleFit(),
        contrast = list(first = exampleContrast(), second = c("condition", "BN", "SHR")),
        verbose = FALSE
    )

    expect_identical(as.list(resultsList), methods::slot(resultsList, "results"))
    expect_identical(names(as.list(resultsList)), c("first", "second"))
    expect_identical(contrastName(resultsList), methods::slot(resultsList, "contrasts"))
})


test_that("the generics of BiocGenerics answer for the objects of the package", {
    fit <- exampleFit()
    results <- exampleResults()

    expect_identical(counts(fit), SummarizedExperiment::assay(fitCounts(fit), "counts"))
    expect_identical(counts(results), SummarizedExperiment::assay(resultCounts(results), "counts"))
    expect_identical(design(fit), designMatrix(fit))

    expect_s3_class(plotMA(results), "ggplot")
    expect_s3_class(plotPCA(fit), "ggplot")

    # A result without its counts says so
    bareResults <- testRegions(fit, contrast = exampleContrast(), carryCounts = FALSE, verbose = FALSE)
    expect_error(counts(bareResults), "carries no counts")
})


test_that("the slots are written with their class checked and the tables converted", {
    results <- exampleResults()

    updatedResults <- RegionSetDE:::.replaceSlots(results, results = resultsTable(results)[1:3, ], regions = regionRanges(results)[1:3])
    expect_s4_class(methods::slot(updatedResults, "results"), "DataFrame")
    expect_identical(nrow(resultsTable(updatedResults)), 3L)
    expect_true(methods::validObject(updatedResults))

    expect_error(RegionSetDE:::.replaceSlots(results, contrast = 1), "contrast")
})


test_that("the list slots are checked by the validity methods", {
    results <- exampleResults()
    fit <- exampleFit()

    brokenThresholds <- RegionSetDE:::.replaceSlots(results, thresholds = list(0.05))
    expect_error(methods::validObject(brokenThresholds), "must be named")

    brokenGroups <- RegionSetDE:::.replaceSlots(results, contrast.groups = list(column = "condition"))
    expect_error(methods::validObject(brokenGroups), "contrast.groups")

    brokenFit <- RegionSetDE:::.replaceSlots(fit, fit = list(model = 1))
    expect_error(methods::validObject(brokenFit), "'object' element")

    brokenParameters <- RegionSetDE:::.replaceSlots(fit, parameters = list(1))
    expect_error(methods::validObject(brokenParameters), "parameters")
})


test_that("an object saved with data.frame tables is brought to DataFrame", {
    results <- exampleResults()

    # readRDS of a file written before the change returns the tables as data.frames
    oldResults <- results
    attr(oldResults, "results") <- resultsTable(results)
    attr(oldResults, "filtering.log") <- filteringLog(results)
    expect_true(is.data.frame(methods::slot(oldResults, "results")))

    expect_message(updatedResults <- BiocGenerics::updateObject(oldResults, verbose = TRUE), "DataFrame")
    expect_s4_class(methods::slot(updatedResults, "results"), "DataFrame")
    expect_identical(resultsTable(updatedResults), resultsTable(results))
    expect_true(methods::validObject(updatedResults))
})
