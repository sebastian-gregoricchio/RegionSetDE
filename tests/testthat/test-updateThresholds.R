# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

test_that("updateThresholds labels the regions again without touching the statistics", {
    results <- exampleResults()
    strictResults <- updateThresholds(results, FDR = 0.01, log2FC = 1, verbose = FALSE)

    originalTable <- resultsTable(results)
    strictTable <- resultsTable(strictResults)

    expect_identical(
        strictTable[, setdiff(colnames(strictTable), "diff.status")],
        originalTable[, setdiff(colnames(originalTable), "diff.status")]
    )

    expectedStatus <- ifelse(strictTable$FDR < 0.01 & strictTable$log2FC > 1, "up",
        ifelse(strictTable$FDR < 0.01 & strictTable$log2FC < -1, "down", "null")
    )
    expect_identical(as.character(strictTable$diff.status), expectedStatus)
    expect_identical(levels(strictTable$diff.status), c("down", "null", "up"))

    # The same cut-offs given to the test give the same labels
    retested <- testRegions(exampleFit(), contrast = exampleContrast(), FDR = 0.01, log2FC = 1, verbose = FALSE)
    expect_identical(strictTable$diff.status, resultsTable(retested)$diff.status)
})


test_that("updateThresholds records the new cut-offs and keeps the previous ones", {
    results <- exampleResults()
    strictResults <- updateThresholds(results, log2FC = 0.5, verbose = FALSE)

    expect_identical(thresholdInfo(strictResults)$log2FC, 0.5)
    expect_identical(thresholdInfo(strictResults)$FDR, thresholdInfo(results)$FDR)
    expect_identical(parameterLog(strictResults)$testRegions$log2FC, 0.5)
    expect_identical(parameterLog(strictResults)$updateThresholds$previous.log2FC, thresholdInfo(results)$log2FC)

    expect_message(updateThresholds(results, FDR = 0.1), "up and")
    expect_error(updateThresholds(results), "at least one new cut-off")
    expect_error(updateThresholds(results, FDR = 2), "not above 1")
    expect_error(updateThresholds(results, log2FC = -1), "at least 0")
    expect_error(updateThresholds(exampleFit(), FDR = 0.1), "returned by testRegions")
})


test_that("updateThresholds works on several contrasts and on the region sets", {
    fit <- exampleFit()
    contrastList <- list(first = exampleContrast(), second = c("condition", "BN", "SHR"))
    resultsList <- testRegions(fit, contrast = contrastList, verbose = FALSE)

    allUpdated <- updateThresholds(resultsList, FDR = 0.01, verbose = FALSE)
    expect_identical(vapply(as.list(allUpdated), function(x) {
        thresholdInfo(x)$FDR
    }, numeric(1)), c(first = 0.01, second = 0.01))

    oneUpdated <- updateThresholds(resultsList, FDR = 0.01, contrast = "second", verbose = FALSE)
    expect_identical(thresholdInfo(oneUpdated$first)$FDR, thresholdInfo(resultsList$first)$FDR)
    expect_identical(thresholdInfo(oneUpdated$second)$FDR, 0.01)
    expect_error(updateThresholds(resultsList, FDR = 0.01, contrast = "third"), "absent from the object")

    setResults <- exampleSetResults()
    updatedSets <- updateThresholds(setResults, FDR = 0.5, verbose = FALSE)
    expect_identical(thresholdInfo(updatedSets)$FDR, 0.5)
    expect_identical(parameterLog(updatedSets)$testRegionSets$FDR, 0.5)
    expect_error(updateThresholds(setResults, log2FC = 1), "only 'FDR' can be updated")
})
