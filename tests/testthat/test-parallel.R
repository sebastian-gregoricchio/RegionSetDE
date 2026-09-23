# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

## Socket workers, the only kind Windows has, are new R processes. They receive functions that point
## to the RegionSetDE namespace by name and load it from their own libraries, so a session that found
## the package somewhere else, as R CMD build does with its temporary library, would send new code to
## an old namespace, or to none at all.

test_that("socket workers load RegionSetDE from the same library as the session", {

  testthat::skip_on_cran()

  # Loaded from the sources with devtools::load_all() there is no installed package a worker could load
  sessionPath <- getNamespaceInfo(asNamespace("RegionSetDE"), "path")
  testthat::skip_if(file.exists(file.path(sessionPath, "R", "bamUtils.R")), "RegionSetDE is loaded from its sources, not installed")

  socketParam <- BiocParallel::SnowParam(workers = 2)

  workerPaths <- RegionSetDE:::.bplapplySameLibraries(seq_len(2),
                                                      function(i) {getNamespaceInfo(asNamespace("RegionSetDE"), "path")},
                                                      BPPARAM = socketParam)

  expect_identical(unique(normalizePath(unlist(workerPaths))), normalizePath(sessionPath))
})


test_that("the library setting handed to the workers does not outlive the call", {

  testthat::skip_on_cran()

  previousLibraries <- Sys.getenv("R_LIBS", unset = NA)

  RegionSetDE:::.bplapplySameLibraries(seq_len(2), function(i) {i}, BPPARAM = BiocParallel::SnowParam(workers = 2))
  expect_identical(Sys.getenv("R_LIBS", unset = NA), previousLibraries)

  # Serial and forked workers share the session, nothing is set for them
  expect_identical(RegionSetDE:::.bplapplySameLibraries(seq_len(2), function(i) {i * 2}, BPPARAM = BiocParallel::SerialParam()),
                   list(2, 4))
  expect_identical(Sys.getenv("R_LIBS", unset = NA), previousLibraries)
})
