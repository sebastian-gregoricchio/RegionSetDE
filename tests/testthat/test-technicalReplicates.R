# A narrowPeak file with strong peaks at the given starts, the summit in the middle
writeTechnicalPeaks <- function(fileName,
                                starts,
                                width = 500L,
                                shift = 0L) {
    peakStarts <- as.integer(starts + shift)
    peakTable <- data.frame(
        "chr1", peakStarts - 1L, peakStarts - 1L + width, paste0("peak_", seq_along(starts)),
        1000L, ".", 10, 30, 29, width %/% 2L
    )

    peakPath <- file.path(tempdir(), fileName)
    utils::write.table(peakTable, peakPath, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    return(peakPath)
}


# Two conditions, two biological samples each, two libraries per sample, R1 and R2 reused in both conditions.
# 500000 is found in the two libraries of A_R1 only, 600000 in one library of A_R1 and one of A_R2
technicalSheet <- function() {
    sharedStarts <- seq(10000L, 100000L, by = 10000L)
    libraryNames <- c("A_1a", "A_1b", "A_2a", "A_2b", "B_1a", "B_1b", "B_2a", "B_2b")

    libraryStarts <- list(
        A_1a = c(sharedStarts, 500000L, 600000L),
        A_1b = c(sharedStarts, 500000L),
        A_2a = c(sharedStarts, 600000L),
        A_2b = sharedStarts,
        B_1a = sharedStarts,
        B_1b = sharedStarts,
        B_2a = sharedStarts,
        B_2b = sharedStarts
    )

    peakFiles <- vapply(seq_along(libraryNames), function(i) {
        writeTechnicalPeaks(paste0("technical_", libraryNames[i], ".narrowPeak"),
            starts = libraryStarts[[libraryNames[i]]],
            shift = (i - 1L) * 10L
        )
    }, character(1))

    return(loadSampleSheet(
        data.frame(
            sample = libraryNames,
            bam = paste0(libraryNames, ".bam"),
            peaks = peakFiles,
            condition = rep(c("A", "B"), each = 4),
            bioRep = rep(c("R1", "R1", "R2", "R2"), times = 2)
        ),
        checkFiles = FALSE,
        verbose = FALSE
    ))
}


test_that("technical replicates are joined before the consensus of the group", {
    skip_if_not_installed("consensusRegions")

    sheet <- technicalSheet()
    plainRegions <- loadConsensusPeaks(sheet, groupBy = "condition", verbose = FALSE)
    technicalRegions <- loadConsensusPeaks(sheet, groupBy = "condition", technicalBy = "bioRep", verbose = FALSE)

    probe <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(500100L, 600100L), width = 10L))

    # As libraries, the two libraries of A_R1 support each other and so do A_1a and A_2a
    expect_true(all(IRanges::overlapsAny(probe, regionRanges(plainRegions)$consensus)))

    # As biological samples, A_R1 has no other sample to agree with, and 600000 is in one library of each
    expect_false(any(IRanges::overlapsAny(probe, regionRanges(technicalRegions)$consensus)))
    expect_length(regionRanges(technicalRegions)$consensus, 10L)
})


test_that("the biological samples are named, counted and kept", {
    skip_if_not_installed("consensusRegions")

    expect_message(
        technicalRegions <- loadConsensusPeaks(technicalSheet(), groupBy = "condition", technicalBy = "bioRep"),
        "named after both"
    )
    storedConsensus <- consensusData(technicalRegions)

    expect_named(storedConsensus$technical, c("A_R1", "A_R2", "B_R1", "B_R2"))
    expect_true(all(vapply(storedConsensus$technical, methods::is, logical(1), "ConsensusRegions")))
    expect_identical(storedConsensus$technicalBy, "bioRep")
    expect_identical(storedConsensus$samples$biological, rep(c("A_R1", "A_R2", "B_R1", "B_R2"), each = 2))

    consensusRanges <- regionRanges(technicalRegions)$consensus
    expect_true("peak.biological" %in% colnames(S4Vectors::mcols(consensusRanges)))
    expect_true(all(consensusRanges$peak.biological <= 4L))
    expect_true(all(consensusRanges$peak.samples >= consensusRanges$peak.biological))

    loggedParameters <- parameterLog(technicalRegions)$loadConsensusPeaks
    expect_identical(loggedParameters$technicalBy, "bioRep")
    expect_identical(loggedParameters$n.biological, 4L)
    expect_identical(loggedParameters$technicalArguments$replicateType, "technical")
})


test_that("labels unique to a group are kept as they are", {
    skip_if_not_installed("consensusRegions")

    sheet <- technicalSheet()
    sheet$bioRep <- paste(sheet$condition, sheet$bioRep, sep = "-")

    technicalRegions <- loadConsensusPeaks(sheet, groupBy = "condition", technicalBy = "bioRep", verbose = FALSE)
    expect_named(consensusData(technicalRegions)$technical, c("A-R1", "A-R2", "B-R1", "B-R2"))
})


test_that("a biological sample with a single library keeps its peaks", {
    skip_if_not_installed("consensusRegions")

    sheet <- technicalSheet()
    sheet <- sheet[sheet$sample != "B_2b", ]

    technicalRegions <- loadConsensusPeaks(sheet, groupBy = "condition", technicalBy = "bioRep", verbose = FALSE)
    storedConsensus <- consensusData(technicalRegions)

    expect_null(storedConsensus$technical$B_R2)
    expect_s4_class(storedConsensus$technical$B_R1, "ConsensusRegions")
    expect_true(isTRUE(S4Vectors::mcols(regionRanges(technicalRegions)$consensus)$peak.B[1]))
})


test_that("technicalArguments reach the step joining the libraries", {
    skip_if_not_installed("consensusRegions")

    technicalRegions <- loadConsensusPeaks(technicalSheet(),
        groupBy = "condition",
        technicalBy = "bioRep",
        technicalArguments = list(mergeMethod = "iterative"),
        verbose = FALSE
    )

    loggedArguments <- parameterLog(technicalRegions)$loadConsensusPeaks$technicalArguments
    expect_identical(loggedArguments$mergeMethod, "iterative")
    expect_identical(parameterLog(technicalRegions)$loadConsensusPeaks$consensusArguments$mergeMethod, "reduce")
})


test_that("the summits of the libraries reach the recentred consensus", {
    skip_if_not_installed("consensusRegions")

    sheet <- technicalSheet()
    sheet <- sheet[sheet$condition == "A", ]

    recentredRegions <- loadConsensusPeaks(sheet,
        technicalBy = "bioRep",
        recentre = TRUE,
        width = 200L,
        verbose = FALSE
    )

    consensusRanges <- regionRanges(recentredRegions)$consensus
    expect_true(all(BiocGenerics::width(consensusRanges) == 200L))

    # Every library has its summit 250 bp into a 500 bp peak, so a window lies inside the peaks of its region
    expect_true(all(GenomicRanges::start(consensusRanges) >= 10000L))
})


test_that("the technical arguments are checked", {
    skip_if_not_installed("consensusRegions")

    sheet <- technicalSheet()

    expect_error(loadConsensusPeaks(sheet, groupBy = "condition", technicalBy = "donor", verbose = FALSE), "must name a column")
    expect_error(loadConsensusPeaks(sheet, groupBy = "condition", technicalBy = "condition", verbose = FALSE), "other than 'sample' and 'groupBy'")
    expect_error(loadConsensusPeaks(sheet, technicalBy = "bioRep", technicalArguments = list("technical"), verbose = FALSE), "named list")
    expect_error(
        loadConsensusPeaks(sheet, technicalBy = "bioRep", technicalArguments = list(calibrate = TRUE), verbose = FALSE),
        "cannot be given to the step joining"
    )
    expect_error(loadConsensusPeaks(sheet, technicalBy = "bioRep", weightMethod = "frip", verbose = FALSE), "per biological sample")

    sheet$bioRep[1] <- NA
    expect_error(loadConsensusPeaks(sheet, technicalBy = "bioRep", verbose = FALSE), "no value in the 'bioRep' column")
})


test_that("libraries sharing no peak leave their sample out", {
    skip_if_not_installed("consensusRegions")

    sheet <- technicalSheet()
    sheet$peaks[sheet$sample == "A_1b"] <- writeTechnicalPeaks("technical_empty.narrowPeak", starts = 900000L)

    # consensusRegions reports the empty consensus of A_R1 first, then the sample is left out
    raisedWarnings <- testthat::capture_warnings(
        technicalRegions <- loadConsensusPeaks(sheet, groupBy = "condition", technicalBy = "bioRep", verbose = FALSE)
    )
    expect_true(any(grepl("share no peak", raisedWarnings)))
    # A_R2 is left alone in its group, and its ten shared peaks stand for it
    expect_length(consensusData(technicalRegions)$groups$A, 10L)
    expect_true(isTRUE(S4Vectors::mcols(regionRanges(technicalRegions)$consensus)$peak.A[1]))
})
