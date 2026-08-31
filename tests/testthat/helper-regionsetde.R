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
