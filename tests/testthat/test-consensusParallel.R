## The groups of loadConsensusPeaks built in parallel.


test_that("the threads are shared between the groups and the calibration", {

  workerPlan <- RegionSetDE:::.consensusWorkers(nThreads = 6, nGroups = 3, calibrate = FALSE)
  expect_identical(workerPlan$outer, 3L)
  expect_identical(workerPlan$inner, 1L)

  workerPlan <- RegionSetDE:::.consensusWorkers(nThreads = 6, nGroups = 3, calibrate = TRUE)
  expect_identical(workerPlan$inner, 2L)

  # One group leaves all the threads to its permutations
  workerPlan <- RegionSetDE:::.consensusWorkers(nThreads = 4, nGroups = 1, calibrate = TRUE)
  expect_identical(workerPlan$outer, 1L)
  expect_identical(workerPlan$inner, 4L)

  # A back end of the user runs the groups one after the other
  workerPlan <- RegionSetDE:::.consensusWorkers(nThreads = 4, nGroups = 3, userBPPARAM = BiocParallel::SerialParam())
  expect_identical(workerPlan$outer, 1L)
  expect_s4_class(workerPlan$inner, "SerialParam")

  expect_error(RegionSetDE:::.consensusWorkers(nThreads = 0, nGroups = 3), "nThreads")
})


test_that("the consensus is the same on one thread or several", {

  skip_if_not_installed("consensusRegions")
  skip_on_os("windows")
  sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)

  serialRegions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl",
                                      nThreads = 1, verbose = FALSE)
  parallelRegions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl",
                                        nThreads = 2, verbose = FALSE)

  expect_identical(regionRanges(serialRegions), regionRanges(parallelRegions))
  expect_identical(names(consensusData(parallelRegions)$groups), names(consensusData(serialRegions)$groups))
  expect_identical(lengths(consensusData(parallelRegions)$groups), lengths(consensusData(serialRegions)$groups))

  # The messages still come one per group, in the order of the groups
  groupMessages <- testthat::capture_messages(loadConsensusPeaks(sampleSheet, groupBy = "condition",
                                                                 seqlevelsStyle = "Ensembl", nThreads = 2))
  groupLines <- grep("^Group ", groupMessages, value = TRUE)
  expect_length(groupLines, length(unique(sampleSheet$condition)))
})
