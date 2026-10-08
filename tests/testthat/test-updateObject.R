# Objects saved by an earlier version of the package lack the slots added since. Taking the attribute away
# leaves an object in that state, which is what readRDS returns for a file written before the slot existed.
withoutGreylist <- function(object) {
  attr(object, "greylist") <- NULL
  return(object)
}


test_that("updateObject adds the greylist slot to an object saved without it", {

  counts <- toyCounts()
  oldCounts <- withoutGreylist(counts)

  expect_false(methods::.hasSlot(oldCounts, "greylist"))
  expect_error(oldCounts@greylist, "no slot")

  expect_message(updatedCounts <- BiocGenerics::updateObject(oldCounts, verbose = TRUE), "gained the 'greylist' slot, left empty")
  expect_true(methods::.hasSlot(updatedCounts, "greylist"))
  expect_null(updatedCounts@greylist)
  expect_true(methods::validObject(updatedCounts))
  expect_identical(updatedCounts, counts)

  # An object already up to date comes back as it was
  expect_identical(BiocGenerics::updateObject(counts), counts)

  # The region sets, without a consensus, get an empty slot as well
  regions <- toyRegionSet()
  expect_identical(BiocGenerics::updateObject(withoutGreylist(regions)), regions)
})


test_that("a fit and its results update the counts they carry", {

  fit <- exampleFit()
  results <- testRegions(fit, contrast = exampleContrast(), verbose = FALSE)

  # The counts inside are the ones a saved fit would carry from before the slot
  oldFit <- fit
  oldCounts <- withoutGreylist(fit@counts)
  methods::slot(oldFit, "counts", check = FALSE) <- oldCounts
  oldFit <- withoutGreylist(oldFit)

  updatedFit <- BiocGenerics::updateObject(oldFit)
  expect_true(methods::.hasSlot(updatedFit, "greylist"))
  expect_true(methods::.hasSlot(updatedFit@counts, "greylist"))
  expect_identical(updatedFit, fit)

  oldResults <- withoutGreylist(results)
  expect_identical(BiocGenerics::updateObject(oldResults), results)

  # A list of results updates every result in it
  resultsList <- testRegions(fit, contrast = list(strainEffect = exampleContrast(), reversed = c("condition", "BN", "SHR")), verbose = FALSE)
  expect_s4_class(resultsList, "RegionSetDE.resultsList")

  oldList <- resultsList
  oldList@results <- lapply(oldList@results, withoutGreylist)
  expect_identical(BiocGenerics::updateObject(oldList), resultsList)
})


test_that("the example objects come back with the current slots", {

  counts <- exampleCounts()

  expect_true(methods::.hasSlot(counts, "greylist"))
  expect_true(methods::validObject(counts))
  expect_true(methods::.hasSlot(exampleFit(), "greylist"))
})
