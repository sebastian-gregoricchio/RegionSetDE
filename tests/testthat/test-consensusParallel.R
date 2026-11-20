## The groups of loadConsensusPeaks built in parallel.


test_that("the cores are shared between the groups and the calibration", {
    workerPlan <- RegionSetDE:::.consensusWorkers(BPPARAM = 6, nGroups = 3, calibrate = FALSE)
    expect_identical(workerPlan$outer, 3L)
    expect_identical(workerPlan$inner, 1L)

    workerPlan <- RegionSetDE:::.consensusWorkers(BPPARAM = 6, nGroups = 3, calibrate = TRUE)
    expect_identical(workerPlan$inner, 2L)

    # One group leaves all the cores to its permutations
    workerPlan <- RegionSetDE:::.consensusWorkers(BPPARAM = 4, nGroups = 1, calibrate = TRUE)
    expect_identical(workerPlan$outer, 1L)
    expect_identical(workerPlan$inner, 4L)

    # NULL is a single core
    workerPlan <- RegionSetDE:::.consensusWorkers(BPPARAM = NULL, nGroups = 3, calibrate = TRUE)
    expect_identical(workerPlan$outer, 1L)
    expect_identical(workerPlan$inner, 1L)

    # A back end of the user goes to the calibration, and the groups run one after the other
    workerPlan <- RegionSetDE:::.consensusWorkers(BPPARAM = BiocParallel::SerialParam(), nGroups = 3, calibrate = TRUE)
    expect_identical(workerPlan$outer, 1L)
    expect_s4_class(workerPlan$inner, "SerialParam")

    # Without calibration the same back end runs the groups
    workerPlan <- RegionSetDE:::.consensusWorkers(BPPARAM = BiocParallel::SerialParam(), nGroups = 3, calibrate = FALSE)
    expect_s4_class(workerPlan$outer, "SerialParam")
    expect_identical(workerPlan$inner, 1L)

    expect_error(RegionSetDE:::.consensusWorkers(BPPARAM = 0, nGroups = 3), "BPPARAM")
})


test_that("BPPARAM takes a back end, a number of cores or nothing", {
    expect_s4_class(RegionSetDE:::.resolveBPPARAM(NULL), "SerialParam")
    expect_s4_class(RegionSetDE:::.resolveBPPARAM(1), "SerialParam")

    # A back end of the user comes back as it was given
    userParam <- BiocParallel::SnowParam(workers = 2, progressbar = FALSE)
    expect_identical(RegionSetDE:::.resolveBPPARAM(userParam, tasks = 10L, progressBar = TRUE), userParam)

    if (.Platform$OS.type != "windows") {
        expect_s4_class(RegionSetDE:::.resolveBPPARAM(2), "MulticoreParam")
        expect_identical(BiocParallel::bpnworkers(RegionSetDE:::.resolveBPPARAM(2)), 2L)
    }

    expect_error(RegionSetDE:::.resolveBPPARAM(0), "'BPPARAM' parameter must be")
    expect_error(RegionSetDE:::.resolveBPPARAM("two"), "'BPPARAM' parameter must be")
})


test_that("the consensus is the same on one core or several", {
    skip_if_not_installed("consensusRegions")
    skip_on_os("windows")
    sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)

    serialRegions <- loadConsensusPeaks(sampleSheet,
        groupBy = "condition", seqlevelsStyle = "Ensembl",
        BPPARAM = 1, verbose = FALSE
    )
    parallelRegions <- loadConsensusPeaks(sampleSheet,
        groupBy = "condition", seqlevelsStyle = "Ensembl",
        BPPARAM = 2, verbose = FALSE
    )

    expect_identical(regionRanges(serialRegions), regionRanges(parallelRegions))
    expect_identical(names(consensusData(parallelRegions)$groups), names(consensusData(serialRegions)$groups))
    expect_identical(lengths(consensusData(parallelRegions)$groups), lengths(consensusData(serialRegions)$groups))

    # The messages still come one per group, in the order of the groups
    groupMessages <- testthat::capture_messages(loadConsensusPeaks(sampleSheet,
        groupBy = "condition",
        seqlevelsStyle = "Ensembl", BPPARAM = 2
    ))
    groupLines <- grep("^Group ", groupMessages, value = TRUE)
    expect_length(groupLines, length(unique(sampleSheet$condition)))
})
