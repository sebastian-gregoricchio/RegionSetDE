## Every method is a separate branch, and the branches are most of the two files.


test_that("each normalisation method returns usable factors", {

  counts <- exampleCounts()

  # The methods needing no extra input
  for (normalizationMethod in c("TMM", "TMMwsp", "RLE", "upperQuartile",
                                "librarySize", "background", "none")) {

    normalized <- normalizeCounts(counts, method = normalizationMethod, verbose = FALSE)

    scalingFactors <- SummarizedExperiment::colData(normalized)$scaling.factor

    expect_s4_class(normalized, "RegionSetDE.counts")
    expect_length(scalingFactors, ncol(counts))
    expect_true(all(is.finite(scalingFactors)),
                info = paste("non finite factors from", normalizationMethod))
    expect_true(all(scalingFactors > 0),
                info = paste("non positive factors from", normalizationMethod))
  }
})


test_that("'none' leaves every sample on the same footing", {

  counts <- normalizeCounts(exampleCounts(), method = "none", verbose = FALSE)

  scalingFactors <- SummarizedExperiment::colData(counts)$scaling.factor

  expect_equal(length(unique(round(scalingFactors, 8))), 1)
})


test_that("supplied factors are applied in both directions", {

  counts <- exampleCounts()
  suppliedFactors <- c(1, 2, 1, 2)

  byDivision <- normalizeCounts(counts, method = "manual",
                                scalingFactors = suppliedFactors,
                                factorType = "division", verbose = FALSE)

  byMultiplication <- normalizeCounts(counts, method = "manual",
                                      scalingFactors = suppliedFactors,
                                      factorType = "multiplication", verbose = FALSE)

  expect_s4_class(byDivision, "RegionSetDE.counts")
  expect_s4_class(byMultiplication, "RegionSetDE.counts")

  # The two readings of the same numbers cannot give the same normalised counts
  expect_false(isTRUE(all.equal(
    SummarizedExperiment::assay(byDivision, "norm.counts"),
    SummarizedExperiment::assay(byMultiplication, "norm.counts"))))
})


test_that("normalizeCounts rejects the wrong number of supplied factors", {

  expect_error(
    normalizeCounts(exampleCounts(), method = "manual",
                    scalingFactors = c(1, 2), verbose = FALSE)
  )
})


test_that("normalizeCounts rejects a method it does not have", {

  expect_error(
    normalizeCounts(exampleCounts(), method = "notAMethod", verbose = FALSE)
  )
})


test_that("useRegionSets restricts where the factors are estimated", {

  counts <- exampleCounts()

  fromAll <- normalizeCounts(counts, method = "TMM", verbose = FALSE)
  fromPromoters <- normalizeCounts(counts, method = "TMM",
                                   useRegionSets = "promoterCpG", verbose = FALSE)

  expect_false(isTRUE(all.equal(
    SummarizedExperiment::colData(fromAll)$scaling.factor,
    SummarizedExperiment::colData(fromPromoters)$scaling.factor)))
})


test_that("a normalized assay can be given its own name", {

  counts <- normalizeCounts(exampleCounts(), method = "TMM",
                            normalizedAssay = "scaled", verbose = FALSE)

  expect_true("scaled" %in% SummarizedExperiment::assayNames(counts))
})


test_that("each filtering method keeps a sensible subset", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  for (filterMethod in c("background", "abundance", "proportion")) {

    filtered <- filterRegions(counts, method = filterMethod, verbose = FALSE)

    expect_s4_class(filtered, "RegionSetDE.counts")
    expect_gt(nrow(filtered), 0)
    expect_lte(nrow(filtered), nrow(counts))
    expect_equal(ncol(filtered), ncol(counts))
  }
})


test_that("the byExpr filter needs the groups to be declared", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  expect_error(filterRegions(counts, method = "byExpr", verbose = FALSE))

  # A vector of values rather than a column name
  expect_error(
    filterRegions(counts, method = "byExpr",
                  group = SummarizedExperiment::colData(counts)$condition,
                  verbose = FALSE),
    "single column name")

  expect_error(
    filterRegions(counts, method = "byExpr", group = "notAColumn", verbose = FALSE),
    "absent from the colData")

  filtered <- filterRegions(counts, method = "byExpr", group = "condition",
                            verbose = FALSE)

  expect_gt(nrow(filtered), 0)
  expect_lte(nrow(filtered), nrow(counts))
})


test_that("the manual filter keeps exactly what it is given", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  keptRows <- seq_len(200)

  filtered <- filterRegions(counts, method = "manual", keep = keptRows, verbose = FALSE)

  expect_equal(nrow(filtered), length(keptRows))
})


test_that("filterRegions rejects a method it does not have", {

  expect_error(
    filterRegions(exampleCounts(), method = "notAMethod", verbose = FALSE)
  )
})


test_that("bySet filters each set on its own terms", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  perSet <- filterRegions(counts, bySet = TRUE, verbose = FALSE)
  pooled <- filterRegions(counts, bySet = FALSE, verbose = FALSE)

  expect_s4_class(perSet, "RegionSetDE.counts")
  expect_s4_class(pooled, "RegionSetDE.counts")
})


test_that("byWidth accounts for the region length", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  withWidth <- filterRegions(counts, byWidth = TRUE, verbose = FALSE)
  withoutWidth <- filterRegions(counts, byWidth = FALSE, verbose = FALSE)

  expect_s4_class(withWidth, "RegionSetDE.counts")
  expect_s4_class(withoutWidth, "RegionSetDE.counts")
})


test_that("plotNormComparison draws both of its layouts", {

  counts <- exampleCounts()

  expect_s3_class(plotNormComparison(counts, plotType = "factors"), "ggplot")
  expect_s3_class(plotNormComparison(counts, plotType = "ma"), "ggplot")

  expect_s3_class(
    plotNormComparison(counts, methods = c("TMM", "RLE"), plotType = "factors"),
    "ggplot")
})
