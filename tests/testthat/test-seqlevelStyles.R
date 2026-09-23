# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

## Every place where two sets of coordinates meet is a place where one of them can be written chr19
## and the other 19. The overlap is then empty, which is not an error and produces no message: the
## blacklist removes nothing, the occupancy counts nobody, the excluded chromosome stays in the
## library size. These tests give each junction the same input twice, once in each style, and ask
## for the same answer.

testthat::skip_if_not_installed("consensusRegions")

sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)

consensusEnsembl <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)
consensusUCSC <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "UCSC", verbose = FALSE)

# Half of the counted window, so that the answer cannot be "nothing overlapped anyway"
halfWindowEnsembl <- GenomicRanges::GRanges("19", IRanges::IRanges(46e6, 52e6))
halfWindowUCSC <- GenomicRanges::GRanges("chr19", IRanges::IRanges(46e6, 52e6))


test_that("the regions are counted whichever style they carry", {

  ensemblCounts <- countTable(countReads(consensusEnsembl, sampleSheet = sampleSheet, pairedEnd = TRUE, verbose = FALSE),
                              format = "matrix")
  ucscCounts <- countTable(countReads(consensusUCSC, sampleSheet = sampleSheet, pairedEnd = TRUE, verbose = FALSE),
                           format = "matrix")

  expect_identical(unname(ensemblCounts), unname(ucscCounts))

  # The counting renames its own copy, the object keeps what it was built with
  expect_identical(GenomeInfoDb::seqlevels(consensusUCSC@regions$consensus), "chr19")
})


test_that("the regions removed by a filter do not depend on the style of the filter", {

  ensemblFiltered <- applyBlacklist(consensusEnsembl, blacklist = halfWindowEnsembl, verbose = FALSE)
  ucscFiltered <- applyBlacklist(consensusEnsembl, blacklist = halfWindowUCSC, verbose = FALSE)

  expect_identical(length(ensemblFiltered@regions$consensus), length(ucscFiltered@regions$consensus))
  expect_lt(length(ensemblFiltered@regions$consensus), length(consensusEnsembl@regions$consensus))
})


test_that("a filter sharing no chromosome with the regions is refused", {

  unknownAssembly <- GenomicRanges::GRanges("scaffold_1", IRanges::IRanges(1, 1e6))

  expect_error(applyBlacklist(consensusEnsembl, blacklist = unknownAssembly, verbose = FALSE),
               "shares no chromosome name")
})


test_that("the peaks excluded before the consensus do not depend on the style of the exclusion list", {

  ensemblConsensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl",
                                         excludeRegions = halfWindowEnsembl, verbose = FALSE)
  ucscConsensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl",
                                      excludeRegions = halfWindowUCSC, verbose = FALSE)

  expect_identical(length(consensusData(ensemblConsensus)$total), length(consensusData(ucscConsensus)$total))
  expect_lt(length(consensusData(ensemblConsensus)$total), length(consensusData(consensusEnsembl)$total))
})


test_that("the reads discarded by region do not depend on the style of the regions", {

  countInWindow <- function(discarded) {
    sum(countTable(countReads(consensusEnsembl, sampleSheet = sampleSheet, pairedEnd = TRUE,
                              discardRegions = discarded, verbose = FALSE),
                   format = "matrix"))
  }

  expect_identical(countInWindow(halfWindowEnsembl), countInWindow(halfWindowUCSC))
  expect_lt(countInWindow(halfWindowEnsembl), countInWindow(NULL))
})


test_that("the greenlist is counted whichever style it carries", {

  counts <- countReads(consensusEnsembl, sampleSheet = sampleSheet, pairedEnd = TRUE, verbose = FALSE)

  binsEnsembl <- GenomicRanges::GRanges("19", IRanges::IRanges(seq(46e6, 57.9e6, by = 1e5), width = 1e5))
  binsUCSC <- GenomicRanges::GRanges("chr19", IRanges::IRanges(seq(46e6, 57.9e6, by = 1e5), width = 1e5))

  ensemblGreenlist <- countGreenlist(counts, greenlist = binsEnsembl, bamFiles = sampleSheet$bam,
                                     pairedEnd = TRUE, verbose = FALSE)
  ucscGreenlist <- countGreenlist(counts, greenlist = binsUCSC, bamFiles = sampleSheet$bam,
                                  pairedEnd = TRUE, verbose = FALSE)

  expect_identical(sum(RegionSetDE:::.greenlistMatrix(ensemblGreenlist)),
                   sum(RegionSetDE:::.greenlistMatrix(ucscGreenlist)))
})


test_that("a chromosome to exclude is recognised in either style", {

  librarySizeWithout <- function(excluded) {
    counts <- countReads(consensusEnsembl, sampleSheet = sampleSheet, pairedEnd = TRUE,
                         excludeChromosomes = excluded, verbose = FALSE)
    SummarizedExperiment::colData(counts)$library.size[1]
  }

  # The files name it 19, and chr19 has to reach the same place
  expect_identical(suppressWarnings(librarySizeWithout("19")), suppressWarnings(librarySizeWithout("chr19")))
  expect_gt(librarySizeWithout(NULL), 0)

  # A name that matches nothing after the conversion is worth a warning, since it silently changes
  # the library sizes that the normalisation is built on
  expect_warning(countReads(consensusEnsembl, sampleSheet = sampleSheet, pairedEnd = TRUE,
                            excludeChromosomes = "chrY", verbose = FALSE),
                 "exclude nothing")
})


test_that("the mitochondrion crosses between the two naming styles", {

  expect_identical(RegionSetDE:::.matchChromosomeNames(c("chr19", "chrM"), c("19", "MT")), c("19", "MT"))
  expect_identical(RegionSetDE:::.matchChromosomeNames(c("19", "MT"), c("chr19", "chrM")), c("chr19", "chrM"))

  mitochondrial <- GenomicRanges::GRanges("chrM", IRanges::IRanges(1, 1000))
  expect_identical(GenomeInfoDb::seqlevels(RegionSetDE:::.matchSeqlevels(mitochondrial, c("19", "MT"), verbose = FALSE)),
                   "MT")
})
