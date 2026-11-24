# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

## The CRAM example holds the alignments of the Rsamtools BAM example, so every number read from it must
## be the one read from the BAM file. No reference genome is reachable during these tests: the reference
## paths point to an empty folder, and the header of the file carries no UR tag.

noReference <- function() {
    previous <- Sys.getenv(c("REF_PATH", "REF_CACHE"), unset = NA)
    emptyFolder <- file.path(tempdir(), "noReference")
    dir.create(emptyFolder, showWarnings = FALSE)
    Sys.setenv(REF_PATH = file.path(emptyFolder, "%2s", "%s"), REF_CACHE = file.path(emptyFolder, "cache"))
    return(previous)
}


restoreReference <- function(previous) {
    for (variable in names(previous)) {
        if (is.na(previous[[variable]])) {
            Sys.unsetenv(variable)
        } else {
            do.call(Sys.setenv, stats::setNames(list(previous[[variable]]), variable))
        }
    }
}


test_that("the header, the index and the layout of a CRAM file are read", {
    previous <- noReference()

    expect_true(file.exists(toyCramFile()))
    expect_identical(RegionSetDE:::.bamIndexPath(toyCramFile()), paste0(toyCramFile(), ".crai"))
    expect_identical(
        RegionSetDE:::.alignmentTargets(toyCramFile()),
        Rsamtools::scanBamHeader(toyBamFile())[[1]]$targets
    )
    expect_true(RegionSetDE:::.bamIsPairedEnd(toyCramFile()))

    restoreReference(previous)
})


test_that("the reads of a CRAM file are those of its BAM file", {
    previous <- noReference()
    which <- GenomicRanges::GRanges(c("seq1", "seq2"), IRanges::IRanges(c(1L, 200L), c(1575L, 1200L)))

    # Single-end view of both files against Rsamtools on the BAM file: positions, strands and widths on the genome
    cigarWidth <- function(cigar) {
        vapply(regmatches(cigar, gregexpr("[0-9]+[MDN=X]", cigar)), function(operations) {
            sum(as.integer(sub("[MDN=X]$", "", operations)))
        }, integer(1))
    }
    rsamtoolsReads <- Rsamtools::scanBam(toyBamFile(), param = Rsamtools::ScanBamParam(
        what = c("pos", "strand", "cigar"),
        which = which,
        mapqFilter = 20L,
        flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isDuplicate = FALSE)
    ))
    for (alignmentFile in c(toyBamFile(), toyCramFile())) {
        reads <- RegionSetDE:::.scanAlignments(alignmentFile, which = which, isPairedEnd = FALSE, minMapq = 20, removeDuplicates = TRUE, maxFragmentLength = 1000)
        for (stretch in names(rsamtoolsReads)) {
            expect_identical(reads[[stretch]]$pos, rsamtoolsReads[[stretch]]$pos)
            expect_identical(as.character(reads[[stretch]]$strand), as.character(rsamtoolsReads[[stretch]]$strand))
            expect_identical(reads[[stretch]]$width, cigarWidth(rsamtoolsReads[[stretch]]$cigar))
        }
    }

    # Paired-end: the same fragments from both files
    cramPairs <- RegionSetDE:::.scanAlignments(toyCramFile(), which = which, isPairedEnd = TRUE, minMapq = 0, removeDuplicates = TRUE, maxFragmentLength = 1000)
    bamPairs <- RegionSetDE:::.scanAlignments(toyBamFile(), which = which, isPairedEnd = TRUE, minMapq = 0, removeDuplicates = TRUE, maxFragmentLength = 1000)
    for (stretch in names(bamPairs)) {
        expect_gt(length(bamPairs[[stretch]]$pos), 0)
        expect_identical(
            sort(paste(cramPairs[[stretch]]$pos, cramPairs[[stretch]]$mpos, cramPairs[[stretch]]$fragment.start, cramPairs[[stretch]]$isize, cramPairs[[stretch]]$tag$MQ)),
            sort(paste(bamPairs[[stretch]]$pos, bamPairs[[stretch]]$mpos, bamPairs[[stretch]]$fragment.start, bamPairs[[stretch]]$isize, bamPairs[[stretch]]$tag$MQ))
        )
    }

    restoreReference(previous)
})


test_that("the mapping quality of the mate comes from the mate itself", {
    previous <- noReference()
    which <- GenomicRanges::GRanges("seq1", IRanges::IRanges(1L, 1575L))

    cramPairs <- RegionSetDE:::.scanAlignments(toyCramFile(), which = which, isPairedEnd = TRUE, minMapq = 0, removeDuplicates = TRUE, maxFragmentLength = 1000)[[1]]

    # The mapping quality of the second mate of every pair, looked up by read name in the BAM file
    firstMates <- Rsamtools::scanBam(toyBamFile(), param = Rsamtools::ScanBamParam(
        what = c("qname", "pos"),
        which = which,
        flag = Rsamtools::scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isFirstMateRead = TRUE, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isUnmappedQuery = FALSE)
    ))[[1]]
    secondMates <- Rsamtools::scanBam(toyBamFile(), param = Rsamtools::ScanBamParam(
        what = c("qname", "mapq"),
        flag = Rsamtools::scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isSecondMateRead = TRUE, isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isUnmappedQuery = FALSE)
    ))[[1]]
    expectedMapq <- secondMates$mapq[match(firstMates$qname, secondMates$qname)]

    expect_gt(length(cramPairs$pos), 0)
    expect_identical(sort(cramPairs$tag$MQ), sort(expectedMapq[!is.na(expectedMapq)]))

    restoreReference(previous)
})


test_that("the fragment runs from the 5' end of the forward mate to the 5' end of the reverse mate", {
    previous <- noReference()

    # Six designed pairs: a normal one, a read running one base into the adapter, two mates run past
    # each other, a complete overrun, two mates on the same strand and a first mate on the reverse strand
    samFile <- system.file("extdata", "cram", "overrunPairs.sam", package = "RegionSetDE")
    bamFile <- Rsamtools::asBam(samFile, destination = file.path(tempdir(), "overrunPairs"), overwrite = TRUE, indexDestination = TRUE)
    cramFile <- system.file("extdata", "cram", "overrunPairs.cram", package = "RegionSetDE")
    expectedFragments <- IRanges::IRanges(start = c(101L, 501L, 1063L, 1701L), end = c(300L, 543L, 1071L, 1850L))

    for (alignmentFile in c(bamFile, cramFile)) {
        fragments <- RegionSetDE:::.windowFragments(
            bamFile = alignmentFile,
            windows = GenomicRanges::GRanges("chrT", IRanges::IRanges(1L, 2000L)),
            isPairedEnd = TRUE,
            fragmentLength = 150,
            maxFragmentLength = 1000,
            minMapq = 0,
            removeDuplicates = TRUE,
            chromosomeLengths = c(chrT = 2000L)
        )
        expect_identical(sort(IRanges::ranges(fragments)), expectedFragments)
    }

    restoreReference(previous)
})


test_that("countReads gives the counts of the BAM file from its CRAM file", {
    previous <- noReference()

    bamCounts <- toyCounts()
    cramCounts <- RegionSetDE::countReads(toyRegionSet(),
        bamFiles = toyCramFile(),
        sampleNames = "example",
        verbose = FALSE
    )

    expect_identical(SummarizedExperiment::assay(cramCounts, "counts"), SummarizedExperiment::assay(bamCounts, "counts"))
    expect_identical(cramCounts$library.size, bamCounts$library.size)
    expect_true(cramCounts$paired.end)

    # Single-end counting of the same file, with the mapping quality threshold
    bamSingle <- toyCounts(pairedEnd = FALSE, fragmentLength = 100)
    cramSingle <- RegionSetDE::countReads(toyRegionSet(),
        bamFiles = toyCramFile(),
        sampleNames = "example",
        pairedEnd = FALSE,
        fragmentLength = 100,
        verbose = FALSE
    )
    expect_identical(SummarizedExperiment::assay(cramSingle, "counts"), SummarizedExperiment::assay(bamSingle, "counts"))

    restoreReference(previous)
})


test_that("the mate is checked alike in BAM and CRAM files, without the MQ tag", {
    previous <- noReference()

    # The BAM example carries no MQ tag: the mapping quality of the mate is taken from its own record in both formats
    bamCounts <- toyCounts(minMapq = 30)
    cramCounts <- RegionSetDE::countReads(toyRegionSet(),
        bamFiles = toyCramFile(),
        sampleNames = "example",
        minMapq = 30,
        verbose = FALSE
    )

    expect_identical(SummarizedExperiment::assay(cramCounts, "counts"), SummarizedExperiment::assay(bamCounts, "counts"))
    expect_lt(bamCounts$library.size, toyCounts(minMapq = 0)$library.size)

    restoreReference(previous)
})


test_that("BAM and CRAM files are counted together", {
    previous <- noReference()

    counts <- RegionSetDE::countReads(toyRegionSet(),
        bamFiles = c(toyBamFile(), toyCramFile()),
        sampleNames = c("bam", "cram"),
        verbose = FALSE
    )
    countMatrix <- SummarizedExperiment::assay(counts, "counts")

    expect_identical(unname(countMatrix[, "bam"]), unname(countMatrix[, "cram"]))
    expect_identical(counts$library.size[1], counts$library.size[2])

    restoreReference(previous)
})


test_that("libInfo counts the records of a CRAM file by reading it", {
    previous <- noReference()

    bamInfo <- RegionSetDE::libInfo(toyCounts())
    cramInfo <- RegionSetDE::libInfo(RegionSetDE::countReads(toyRegionSet(),
        bamFiles = toyCramFile(),
        sampleNames = "example",
        verbose = FALSE
    ))

    expect_identical(cramInfo$bam.reads, bamInfo$bam.reads)
    expect_identical(cramInfo$bam.mapped, bamInfo$bam.mapped)
    expect_identical(cramInfo$FRiP, bamInfo$FRiP)

    restoreReference(previous)
})


test_that("a CRAM file without index is refused", {
    unindexed <- file.path(tempdir(), "unindexed.cram")
    file.copy(toyCramFile(), unindexed, overwrite = TRUE)

    expect_error(
        RegionSetDE::countReads(toyRegionSet(), bamFiles = unindexed, sampleNames = "example", verbose = FALSE),
        "not indexed"
    )
    expect_error(.Call(RegionSetDE:::C_cramTargets, toyBamFile()), "not a CRAM file")
})
