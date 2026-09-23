## Fixtures shared by the test files. testthat sources every helper-*.R before
## running the tests, so these are available everywhere.


# A handful of regions on the contigs of the Rsamtools example alignment, which
# is the only BAM file available at check time.
toyRegions <- function(width = 300L) {

  regionRanges <- GenomicRanges::GRanges(
    seqnames = rep(c("seq1", "seq2"), each = 3),
    ranges = IRanges::IRanges(start = rep(c(1L, 500L, 1000L), 2), width = width)
  )

  regionRanges$setName <- rep(c("firstSet", "secondSet"), each = 3)

  return(regionRanges)
}


toyRegionSet <- function(...) {
  RegionSetDE::splitLoadRegions(toyRegions(...),
                                splitBy = "setName",
                                seqlevelsStyle = NULL,
                                verbose = FALSE)
}


toyBamFile <- function() {
  system.file("extdata", "ex1.bam", package = "Rsamtools")
}


toyCounts <- function(...) {
  RegionSetDE::countReads(toyRegionSet(),
                          bamFiles = toyBamFile(),
                          sampleNames = "example",
                          verbose = FALSE,
                          ...)
}


# The packaged objects, wrapped so that the verbose argument is not repeated in
# every test.
exampleCounts <- function() {
  RegionSetDE::loadExampleData("counts", verbose = FALSE)
}


exampleFit <- function() {
  RegionSetDE::loadExampleData("fit", verbose = FALSE)
}


exampleContrast <- function() {
  c("condition", "SHR", "BN")
}


exampleResults <- function() {
  RegionSetDE::testRegions(exampleFit(),
                           contrast = exampleContrast(),
                           verbose = FALSE)
}


exampleSetResults <- function() {
  RegionSetDE::testRegionSets(exampleFit(),
                              contrast = exampleContrast(),
                              verbose = FALSE)
}



# A single-end BAM file on one contig: reads spread at random, plus an optional pile of reads over 1 kb.
# An empty second contig can be added, standing in for a chromosome without reads such as chrY in a female sample.
syntheticInput <- function(fileName,
                           pileStart = NULL,
                           emptyContig = FALSE,
                           seed = 1L,
                           chromosome = "chrT",
                           contigLength = 200000L,
                           readCount = 20000L) {
  set.seed(seed)

  positions <- sample.int(contigLength - 100L, readCount, replace = TRUE)
  if (!is.null(pileStart)) {positions <- c(positions, pileStart + sample.int(900L, 3000L, replace = TRUE))}
  positions <- sort(positions)

  samHeader <- c("@HD\tVN:1.6\tSO:coordinate", paste0("@SQ\tSN:", chromosome, "\tLN:", contigLength))
  if (isTRUE(emptyContig)) {samHeader <- c(samHeader, "@SQ\tSN:chrU\tLN:100000")}

  samRecords <- paste(paste0("read", seq_along(positions)), sample(c(0L, 16L), length(positions), replace = TRUE),
                      chromosome, positions, 60, "50M", "*", 0, 0, "*", "*", sep = "\t")

  samFile <- file.path(tempdir(), paste0(fileName, ".sam"))
  writeLines(c(samHeader, samRecords), samFile)

  return(Rsamtools::asBam(samFile, destination = file.path(tempdir(), fileName), overwrite = TRUE, indexDestination = TRUE))
}
