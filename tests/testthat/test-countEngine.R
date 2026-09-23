## The counting engine reads the first mate of each proper pair and rebuilds the
## fragment from its template length. These tests write small BAM files where the
## expected counts are known by construction: every pair below is placed by hand.

# One proper pair in the forward-reverse orientation, the first mate on the left
pairLines <- function(name, chromosome, leftStart, fragmentLength, readLength = 50L,
                      mapq = c(42L, 42L), duplicate = FALSE, proper = TRUE, mateTag = TRUE) {

  rightStart <- leftStart + fragmentLength - readLength
  extraFlag <- (if (proper) {2L} else {0L}) + (if (duplicate) {1024L} else {0L})

  firstTag <- if (mateTag) {paste0("\tMQ:i:", mapq[2])} else {""}
  secondTag <- if (mateTag) {paste0("\tMQ:i:", mapq[1])} else {""}

  c(paste0(name, "\t", 1L + 32L + 64L + extraFlag, "\t", chromosome, "\t", leftStart, "\t", mapq[1], "\t", readLength, "M\t=\t",
           rightStart, "\t", fragmentLength, "\t*\t*", firstTag),
    paste0(name, "\t", 1L + 16L + 128L + extraFlag, "\t", chromosome, "\t", rightStart, "\t", mapq[2], "\t", readLength, "M\t=\t",
           leftStart, "\t", -fragmentLength, "\t*\t*", secondTag))
}


# One single-end read, on either strand
singleLine <- function(name, chromosome, start, reverse = FALSE, cigar = "50M", mapq = 42L) {
  paste0(name, "\t", if (reverse) {16L} else {0L}, "\t", chromosome, "\t", start, "\t", mapq, "\t", cigar, "\t*\t0\t0\t*\t*")
}


# Sorted and indexed by Rsamtools, which is all the counting needs
writeToyBam <- function(recordLines, chromosomeLengths = c(chr1 = 20000L, chr2 = 5000L)) {
  samFile <- tempfile(fileext = ".sam")
  writeLines(c("@HD\tVN:1.6\tSO:unsorted",
               paste0("@SQ\tSN:", names(chromosomeLengths), "\tLN:", chromosomeLengths),
               recordLines),
             samFile)

  return(Rsamtools::asBam(samFile, destination = tempfile(), overwrite = TRUE, indexDestination = TRUE))
}


toyPairedBam <- function(mateTag = TRUE) {
  writeToyBam(c(pairLines("inside", "chr1", 1001L, 100L, mateTag = mateTag),
                pairLines("spanning", "chr1", 5001L, 400L, mateTag = mateTag),
                pairLines("duplicate", "chr1", 1001L, 100L, duplicate = TRUE, mateTag = mateTag),
                pairLines("lowMate", "chr1", 1011L, 100L, mapq = c(42L, 5L), mateTag = mateTag),
                pairLines("tooLong", "chr1", 9001L, 1500L, mateTag = mateTag),
                pairLines("improper", "chr1", 12001L, 300L, proper = FALSE, mateTag = mateTag),
                pairLines("otherChromosome", "chr2", 1001L, 200L, mateTag = mateTag)))
}


toyEngineRegions <- function() {
  regionRanges <- GenomicRanges::GRanges(seqnames = "chr1",
                                         ranges = IRanges::IRanges(start = c(1001L, 5201L, 15001L), width = c(100L, 50L, 100L)))
  names(regionRanges) <- c("inside", "spanned", "empty")
  regionRanges$setName <- "toy"

  RegionSetDE::splitLoadRegions(regionRanges, splitBy = "setName", seqlevelsStyle = NULL, verbose = FALSE)
}


countOf <- function(counts, regionName, column = 1) {
  as.numeric(SummarizedExperiment::assay(counts, "counts")[paste0("toy|", regionName), column])
}


test_that("a paired-end fragment is counted whole, gap between the mates included", {

  counts <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = toyPairedBam(), sampleNames = "toy", verbose = FALSE)

  expect_equal(countOf(counts, "inside"), 1)

  # Both mates lie outside this region, only the fragment between them crosses it
  expect_equal(countOf(counts, "spanned"), 1)
  expect_equal(countOf(counts, "empty"), 0)
})


test_that("the read filters apply to the counts and to the library sizes alike", {

  bamFile <- toyPairedBam()

  # Kept: inside, spanning, otherChromosome. Out: duplicate, lowMate, tooLong, improper.
  filtered <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy", verbose = FALSE)
  expect_equal(filtered$library.size, 3)

  withDuplicates <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy",
                                            removeDuplicates = FALSE, verbose = FALSE)
  expect_equal(withDuplicates$library.size, 4)
  expect_equal(countOf(withDuplicates, "inside"), 2)

  withLowQuality <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy",
                                            minMapq = 0, verbose = FALSE)
  expect_equal(withLowQuality$library.size, 4)

  withLongFragments <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy",
                                               maxFragmentLength = 2000, verbose = FALSE)
  expect_equal(withLongFragments$library.size, 4)
})


test_that("without the MQ tag the second mate passes whatever its quality", {

  counts <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = toyPairedBam(mateTag = FALSE), sampleNames = "toy", verbose = FALSE)

  expect_equal(counts$library.size, 4)

  countingMessages <- testthat::capture_messages(RegionSetDE::countReads(toyEngineRegions(), bamFiles = toyPairedBam(mateTag = FALSE),
                                                                         sampleNames = "toy"))
  expect_true(any(grepl("No MQ tag", countingMessages)))
})


test_that("fullLibrarySize = FALSE only reads the chromosomes carrying regions", {

  bamFile <- toyPairedBam()

  fullCounts <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy", verbose = FALSE)
  partialCounts <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy",
                                           fullLibrarySize = FALSE, verbose = FALSE)

  expect_equal(SummarizedExperiment::assay(partialCounts, "counts"), SummarizedExperiment::assay(fullCounts, "counts"))
  expect_equal(fullCounts$library.size, 3)
  expect_equal(partialCounts$library.size, 2)
})


test_that("excludeChromosomes changes the library sizes, never the counts", {

  bamFile <- toyPairedBam()
  allChromosomes <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy", verbose = FALSE)

  # chr2 carries no region, only its fragment leaves the library size
  withoutSecond <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy",
                                           excludeChromosomes = "chr2", verbose = FALSE)

  expect_equal(SummarizedExperiment::assay(withoutSecond, "counts"), SummarizedExperiment::assay(allChromosomes, "counts"))
  expect_equal(withoutSecond$library.size, 2)

  # chr1 carries the regions: they are still counted, but its fragments leave the library size
  withoutFirst <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy",
                                          excludeChromosomes = "chr1", verbose = FALSE)

  expect_equal(SummarizedExperiment::assay(withoutFirst, "counts"), SummarizedExperiment::assay(allChromosomes, "counts"))
  expect_equal(withoutFirst$library.size, 1)

  countingMessages <- testthat::capture_messages(suppressWarnings(
    RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy",
                            excludeChromosomes = c("chr1", "chrZ"))))
  expect_true(any(grepl("they are counted", countingMessages)))

  # A name that reaches no chromosome leaves it inside the library size, which the normalisation
  # is then built on, so it is a warning rather than a message a quiet call would never show
  expect_warning(RegionSetDE::countReads(toyEngineRegions(), bamFiles = bamFile, sampleNames = "toy",
                                         excludeChromosomes = c("chr1", "chrZ"), verbose = FALSE),
                 "exclude nothing")
})


test_that("a library size of zero raises a warning", {

  # Only chr1 is read, and chr1 is excluded: nothing is left for the library size
  expect_warning(RegionSetDE::countReads(toyEngineRegions(), bamFiles = toyPairedBam(), sampleNames = "toy",
                                         excludeChromosomes = "chr1", fullLibrarySize = FALSE, verbose = FALSE),
                 "No fragment entered the library size")
})


test_that("discardRegions drops the fragments with a read starting inside them", {

  blacklist <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 5301L, end = 5400L))

  counts <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = toyPairedBam(), sampleNames = "toy",
                                    discardRegions = blacklist, verbose = FALSE)

  # The second mate of the spanning pair starts at 5351
  expect_equal(countOf(counts, "spanned"), 0)
  expect_equal(countOf(counts, "inside"), 1)
  expect_equal(counts$library.size, 2)
})


test_that("single-end reads are extended from their 5' end", {

  bamFile <- writeToyBam(c(singleLine("forward", "chr1", 1001L),
                           singleLine("reverse", "chr1", 2951L, reverse = TRUE),
                           singleLine("clipped", "chr1", 7001L, reverse = TRUE, cigar = "10S40M"),
                           singleLine("lowQuality", "chr1", 4001L, mapq = 5L)))

  regionRanges <- GenomicRanges::GRanges(seqnames = "chr1",
                                         ranges = IRanges::IRanges(start = c(1150L, 2810L, 1300L, 6850L, 4050L), width = 10L))
  names(regionRanges) <- c("forward", "reverse", "beyond", "clipped", "lowQuality")
  regionRanges$setName <- "toy"
  regionSet <- RegionSetDE::splitLoadRegions(regionRanges, splitBy = "setName", seqlevelsStyle = NULL, verbose = FALSE)

  counts <- RegionSetDE::countReads(regionSet, bamFiles = bamFile, sampleNames = "toy",
                                    pairedEnd = FALSE, fragmentLength = 200, verbose = FALSE)

  # Forward read at 1001 covers 1001-1200, reverse read ending at 3000 covers 2801-3000
  expect_equal(countOf(counts, "forward"), 1)
  expect_equal(countOf(counts, "reverse"), 1)
  expect_equal(countOf(counts, "beyond"), 0)

  # The clipped read ends at 7040, its fragment covers 6841-7040
  expect_equal(countOf(counts, "clipped"), 1)
  expect_equal(countOf(counts, "lowQuality"), 0)
  expect_equal(counts$library.size, 3)
})


test_that("the layout is read from the files", {

  pairedCounts <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = toyPairedBam(), sampleNames = "toy", verbose = FALSE)
  expect_true(pairedCounts$paired.end)

  singleBam <- writeToyBam(singleLine("forward", "chr1", 1001L))
  singleCounts <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = singleBam, sampleNames = "toy", verbose = FALSE)
  expect_false(singleCounts$paired.end)
})


test_that("countBackground counts every fragment once, even across bins", {

  counts <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = toyPairedBam(), sampleNames = "toy", verbose = FALSE)

  # Bins much narrower than the fragments, a fragment counted by overlap would appear in several of them
  counts <- RegionSetDE::countBackground(counts, binSize = 100L, excludeRegions = FALSE, minCount = 0, verbose = FALSE)
  backgroundBins <- S4Vectors::metadata(counts)$background

  expect_equal(sum(SummarizedExperiment::assay(backgroundBins, "counts")), 3)
  expect_equal(backgroundBins$totals, 3)

  # The last bin of each chromosome stops at the chromosome end
  expect_equal(max(BiocGenerics::end(backgroundBins)), 20000)
})


test_that("the counts do not depend on the number of threads", {

  skip_on_os("windows")

  bamFile <- toyPairedBam()
  oneThread <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = c(bamFile, bamFile), sampleNames = c("a", "b"),
                                       nThreads = 1, verbose = FALSE)
  twoThreads <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = c(bamFile, bamFile), sampleNames = c("a", "b"),
                                        nThreads = 2, verbose = FALSE)

  expect_equal(SummarizedExperiment::assay(twoThreads, "counts"), SummarizedExperiment::assay(oneThread, "counts"))
  expect_equal(twoThreads$library.size, oneThread$library.size)
})


test_that("normalizeCounts warns when the library sizes are partial", {

  partialCounts <- RegionSetDE::countReads(toyEngineRegions(), bamFiles = toyPairedBam(), sampleNames = "toy",
                                           fullLibrarySize = FALSE, verbose = FALSE)

  expect_warning(RegionSetDE::normalizeCounts(partialCounts, method = "librarySize", verbose = FALSE),
                 "fullLibrarySize")
})


test_that("the CIGAR widths add up the reference operations only", {

  expect_equal(RegionSetDE:::.cigarReferenceWidth(c("50M", "10S40M", "20M2D30M", "25M3I22M", "25M100N25M")),
               c(50L, 40L, 52L, 47L, 150L))
})
