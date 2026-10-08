# Two libraries piling fragments on nine sites of one contig, and an input piling reads on the third of them
listSamples <- function() {
  c(syntheticFragments("lists_sample_1"), syntheticFragments("lists_sample_2", seed = 2L))
}


listInput <- function() {
  syntheticInput("lists_input", pileStart = 60000L)
}


listRegions <- function() {
  syntheticSiteRegions(seq(20000L, 180000L, by = 20000L))
}


listCounts <- function(regions = listRegions(), ...) {
  countReads(regions, bamFiles = listSamples(), sampleNames = c("first", "second"), pairedEnd = FALSE, verbose = FALSE, ...)
}


test_that("a second blacklist adds to the stored one and a second whitelist restricts it", {

  firstList <- GenomicRanges::GRanges("seq1", IRanges::IRanges(1, 100))
  secondList <- GenomicRanges::GRanges("seq2", IRanges::IRanges(480, 520))

  once <- applyBlacklist(toyRegionSet(), blacklist = firstList, verbose = FALSE)
  twice <- applyBlacklist(once, blacklist = secondList, verbose = FALSE)

  # Both lists are in the slot, and each call left its own rows in the log
  expect_length(once@blacklist, 1)
  expect_length(twice@blacklist, 2)
  expect_true(all(IRanges::overlapsAny(firstList, twice@blacklist)))
  expect_true(all(IRanges::overlapsAny(secondList, twice@blacklist)))
  expect_identical(nrow(twice@filtering.log), 4L)
  expect_identical(sum(twice@filtering.log$n.removed), 2)

  # Overlapping lists are merged rather than stacked
  widened <- applyBlacklist(once, blacklist = GenomicRanges::GRanges("seq1", IRanges::IRanges(50, 150)), verbose = FALSE)
  expect_length(widened@blacklist, 1)
  expect_identical(BiocGenerics::width(widened@blacklist), 150L)

  # The greylist still leaves the stored blacklist alone
  expect_length(applyGreylist(once, greylist = secondList, verbose = FALSE)@blacklist, 1)

  # Two whitelists leave the stretches they share
  wideList <- GenomicRanges::GRanges(c("seq1", "seq2"), IRanges::IRanges(c(1, 1), c(1200, 2000)))
  narrowList <- GenomicRanges::GRanges("seq1", IRanges::IRanges(400, 2000))

  restricted <- applyWhitelist(applyWhitelist(toyRegionSet(), whitelist = wideList, verbose = FALSE),
                               whitelist = narrowList, emptySets = "keep", verbose = FALSE)

  expect_length(restricted@whitelist, 1)
  expect_identical(BiocGenerics::start(restricted@whitelist), 400L)
  expect_identical(BiocGenerics::end(restricted@whitelist), 1200L)
  expect_identical(as.integer(lengths(restricted@regions)), c(2L, 0L))
})


test_that("countReads removes the blacklisted regions and ignores their reads", {

  # The fifth site lies on the blacklist
  blacklist <- GenomicRanges::GRanges("chrT", IRanges::IRanges(99000, 101000))

  plain <- listCounts()
  listed <- listCounts(blacklist = blacklist)

  expect_identical(nrow(plain), 9L)
  expect_identical(nrow(listed), 8L)
  expect_false(any(IRanges::overlapsAny(SummarizedExperiment::rowRanges(listed), blacklist)))

  # The fragments taken out are counted, and they are what the library sizes lost
  expect_true(all(listed$discarded.reads > 1000))
  expect_equal(listed$library.size + listed$discarded.reads, plain$library.size)
  expect_false("discarded.reads" %in% colnames(SummarizedExperiment::colData(plain)))

  # The list, the step and the discarded regions follow the object
  expect_length(listed@blacklist, 1)
  expect_identical(listed@filtering.log$step, "blacklist")
  expect_identical(listed@filtering.log$n.removed, 1)
  expect_true(all(IRanges::overlapsAny(blacklist, S4Vectors::metadata(listed)$discard.regions)))
  expect_true(listed@parameters$countReads$discardListedReads)
  expect_identical(listed@parameters$countReads$n.discard.regions, 1L)

  # libInfo reports them beside the library sizes
  infoTable <- libInfo(listed)
  expect_identical(infoTable$discarded.reads, as.numeric(listed$discarded.reads))
  expect_false("discarded.reads" %in% colnames(libInfo(plain)))

  # Without the read filter the regions go and the library sizes stay
  kept <- listCounts(blacklist = blacklist, discardListedReads = FALSE)
  expect_identical(nrow(kept), 8L)
  expect_equal(kept$library.size, plain$library.size)
  expect_null(S4Vectors::metadata(kept)$discard.regions)
  expect_length(kept@blacklist, 1)

  # The remaining regions are far from the list, their counts do not move
  expect_identical(countTable(kept, format = "matrix"), countTable(listed, format = "matrix"))

  expect_error(listCounts(blacklist = blacklist, discardListedReads = NA), "'discardListedReads' parameter must be TRUE or FALSE")
})


test_that("a blacklist given to countReads merges with the one stored in the regions", {

  firstList <- GenomicRanges::GRanges("chrT", IRanges::IRanges(99000, 101000))
  secondList <- GenomicRanges::GRanges("chrT", IRanges::IRanges(139000, 141000))

  blacklisted <- applyBlacklist(listRegions(), blacklist = firstList, verbose = FALSE)

  # The stored list alone already takes its reads out
  stored <- listCounts(regions = blacklisted)
  expect_identical(nrow(stored), 8L)
  expect_true(all(stored$discarded.reads > 1000))
  expect_length(S4Vectors::metadata(stored)$discard.regions, 1)

  # A second list in the call adds its regions and its reads
  merged <- listCounts(regions = blacklisted, blacklist = secondList)
  expect_identical(nrow(merged), 7L)
  expect_length(merged@blacklist, 2)
  expect_length(S4Vectors::metadata(merged)$discard.regions, 2)
  expect_identical(merged@filtering.log$n.removed, c(1, 1))
  expect_true(all(merged$discarded.reads > stored$discarded.reads))

  # The same object as with both lists applied beforehand, or both given in the call
  beforehand <- listCounts(regions = applyBlacklist(blacklisted, blacklist = secondList, verbose = FALSE))
  together <- listCounts(blacklist = list(firstList, secondList))

  expect_identical(countTable(merged, format = "matrix"), countTable(beforehand, format = "matrix"))
  expect_identical(countTable(merged, format = "matrix"), countTable(together, format = "matrix"))
  expect_equal(merged$library.size, beforehand$library.size)
  expect_equal(merged$library.size, together$library.size)

  # 'discardRegions' is pooled with the lists, and removes no region
  pooled <- listCounts(regions = blacklisted, discardRegions = secondList)
  expect_identical(nrow(pooled), 8L)
  expect_length(S4Vectors::metadata(pooled)$discard.regions, 2)
  expect_equal(pooled$library.size, merged$library.size)

  # With the read filter off, 'discardRegions' is the only list left
  onlyGiven <- listCounts(regions = blacklisted, discardRegions = secondList, discardListedReads = FALSE)
  expect_length(S4Vectors::metadata(onlyGiven)$discard.regions, 1)
  expect_true(all(IRanges::overlapsAny(secondList, S4Vectors::metadata(onlyGiven)$discard.regions)))

  # Lists lying on different chromosomes are pooled without a word
  expect_warning(apart <- toyCounts(blacklist = GenomicRanges::GRanges("seq1", IRanges::IRanges(1, 100)),
                                    discardRegions = GenomicRanges::GRanges("seq2", IRanges::IRanges(480, 520))),
                 NA)
  expect_identical(sort(as.character(GenomeInfoDb::seqnames(S4Vectors::metadata(apart)$discard.regions))), c("seq1", "seq2"))
})


test_that("greylist = TRUE builds the greylist from the inputs and applies it", {

  inputFile <- listInput()
  plain <- listCounts(inputFiles = inputFile)

  expect_message(greylisted <- countReads(listRegions(), bamFiles = listSamples(), sampleNames = c("first", "second"),
                                          inputFiles = inputFile, pairedEnd = FALSE, greylist = TRUE),
                 "Greylist: 1 regions")

  # The site under the pile of the input is gone, with the reads the samples have there
  builtGreylist <- S4Vectors::metadata(greylisted)$greylist
  expect_s4_class(builtGreylist, "GRanges")
  expect_true(all(IRanges::overlapsAny(GenomicRanges::GRanges("chrT", IRanges::IRanges(60000, 60900)), builtGreylist)))
  expect_identical(S4Vectors::metadata(builtGreylist)$thresholds$input, "lists_input")

  expect_identical(nrow(greylisted), 8L)
  expect_false(any(IRanges::overlapsAny(SummarizedExperiment::rowRanges(greylisted), builtGreylist)))
  expect_identical(greylisted@filtering.log$step, "greylist")
  expect_identical(greylisted@parameters$greylist$n.regions, 1L)
  expect_identical(greylisted@parameters$countReads$greylistSource, "inputs")
  expect_null(greylisted@blacklist)

  expect_equal(greylisted$library.size + greylisted$discarded.reads, plain$library.size)
  expect_true(all(greylisted$discarded.reads > 1000))

  # The input loses the pile that made the greylist
  expect_true(all(plain$input.library.size - greylisted$input.library.size >= 3000))

  # A greylist already built gives the same object without reading the inputs again
  given <- listCounts(inputFiles = inputFile, greylist = builtGreylist)
  expect_identical(countTable(given, format = "matrix"), countTable(greylisted, format = "matrix"))
  expect_equal(given$library.size, greylisted$library.size)
  expect_identical(given@parameters$countReads$greylistSource, "given")
  expect_null(S4Vectors::metadata(given)$greylist)

  # The two lists together, each with its own step
  both <- listCounts(inputFiles = inputFile, greylist = builtGreylist,
                     blacklist = GenomicRanges::GRanges("chrT", IRanges::IRanges(99000, 101000)))
  expect_identical(nrow(both), 7L)
  expect_identical(both@filtering.log$step, c("blacklist", "greylist"))
  expect_length(S4Vectors::metadata(both)$discard.regions, 2)

  expect_error(listCounts(greylist = TRUE), "no sample has one")
  expect_error(listCounts(greylist = NA), "'greylist' parameter must be TRUE, FALSE")
})


test_that("the background and the greenlist leave the discarded reads out", {

  blacklist <- GenomicRanges::GRanges("chrT", IRanges::IRanges(99000, 101000))

  plain <- listCounts()
  listed <- listCounts(blacklist = blacklist)

  # The bins cover the whole contig, so their totals are the library sizes
  plainBins <- S4Vectors::metadata(countBackground(plain, binSize = 10000, verbose = FALSE))$background
  listedBins <- S4Vectors::metadata(countBackground(listed, binSize = 10000, verbose = FALSE))$background

  expect_equal(plainBins$totals, plain$library.size)
  expect_equal(listedBins$totals, listed$library.size)

  # A greenlist region lying on the blacklist holds no read once they are discarded
  greenlist <- GenomicRanges::GRanges("chrT", IRanges::IRanges(c(99500, 150000), width = 1000))

  plainGreenlist <- suppressWarnings(S4Vectors::metadata(countGreenlist(plain, greenlist = greenlist, excludeCounted = FALSE, minCount = 0, verbose = FALSE))$greenlist)
  listedGreenlist <- suppressWarnings(S4Vectors::metadata(countGreenlist(listed, greenlist = greenlist, excludeCounted = FALSE, minCount = 0, verbose = FALSE))$greenlist)

  expect_true(all(SummarizedExperiment::assay(plainGreenlist, "counts")[1, ] > 1000))
  expect_true(all(SummarizedExperiment::assay(listedGreenlist, "counts")[1, ] == 0))
  expect_identical(SummarizedExperiment::assay(plainGreenlist, "counts")[2, ], SummarizedExperiment::assay(listedGreenlist, "counts")[2, ])
})


test_that("the lists are recorded for plain ranges, and the emptied sets leave with a warning", {

  blacklist <- GenomicRanges::GRanges("chrT", IRanges::IRanges(99000, 101000))
  siteRanges <- listRegions()@regions$sites

  # A GRangesList carries no history, the counting writes the step itself
  fromList <- listCounts(regions = GenomicRanges::GRangesList(sites = siteRanges), blacklist = blacklist)
  expect_identical(nrow(fromList), 8L)
  expect_length(fromList@blacklist, 1)
  expect_identical(fromList@filtering.log$step, "blacklist")
  expect_identical(fromList@filtering.log$n.before, 9)
  expect_identical(fromList@filtering.log$n.after, 8)

  # A single GRanges goes through loadRegions first, and the list given under the other naming style still finds it
  fromRanges <- listCounts(regions = siteRanges, blacklist = GenomicRanges::GRanges("T", IRanges::IRanges(99000, 101000)))
  expect_identical(nrow(fromRanges), 8L)
  expect_true(all(fromRanges$discarded.reads > 1000))

  expect_error(listCounts(regions = GenomicRanges::GRangesList(siteRanges), blacklist = blacklist), "must be named")

  # One set emptied: it is dropped and said. Every set emptied: nothing is left to count
  twoSets <- loadRegions(list(left = siteRanges[1:3], right = siteRanges[7:9]), seqlevelsStyle = NULL, verbose = FALSE)

  expect_warning(halved <- listCounts(regions = twoSets, blacklist = GenomicRanges::GRanges("chrT", IRanges::IRanges(1, 90000))),
                 "not counted: left")
  expect_identical(unique(SummarizedExperiment::rowData(halved)$region.set), "right")

  expect_error(listCounts(blacklist = GenomicRanges::GRanges("chrT", IRanges::IRanges(1, 200000))), "left no region to count")

  # An empty list changes nothing
  expect_identical(nrow(listCounts(blacklist = GenomicRanges::GRanges())), 9L)
})


test_that("a list of another assembly is refused, and an alias of the same assembly is not", {

  siteRanges <- listRegions()@regions$sites
  declaredRegions <- loadRegions(list(sites = siteRanges), seqlevelsStyle = NULL, genomeAssembly = "GRCh38", verbose = FALSE)

  blacklist <- GenomicRanges::GRanges("chrT", IRanges::IRanges(99000, 101000))

  GenomeInfoDb::genome(blacklist) <- "hg38"
  expect_identical(nrow(listCounts(regions = declaredRegions, blacklist = blacklist)), 8L)
  expect_length(applyBlacklist(declaredRegions, blacklist = blacklist, verbose = FALSE)@blacklist, 1)

  GenomeInfoDb::genome(blacklist) <- "mm10"
  expect_error(listCounts(regions = declaredRegions, blacklist = blacklist), "blacklist was built for mm10")
  expect_error(listCounts(regions = declaredRegions, greylist = blacklist), "greylist was built for mm10")
  expect_error(applyBlacklist(declaredRegions, blacklist = blacklist, verbose = FALSE), "blacklist was built for mm10")
})


test_that("countBigwig removes the blacklisted regions", {

  bamFile <- listSamples()[1]
  chromosomeLengths <- Rsamtools::scanBamHeader(bamFile)[[1]]$targets
  sampleCoverage <- RegionSetDE:::.windowCoverage(bamFile = bamFile,
                                                  windows = GenomicRanges::GRanges("chrT", IRanges::IRanges(1L, chromosomeLengths[["chrT"]])),
                                                  chromosomeLengths = chromosomeLengths, isPairedEnd = FALSE, fragmentLength = 150,
                                                  maxFragmentLength = 1000, minMapq = 20, removeDuplicates = TRUE)$coverage

  bigwigFile <- file.path(tempdir(), "lists_sample_1.bw")
  rtracklayer::export.bw(sampleCoverage, bigwigFile)

  # rtracklayer writes bigWig files on Windows but cannot read them back from an absolute path
  bigwigReadable <- suppressWarnings(try(GenomeInfoDb::seqlengths(rtracklayer::BigWigFile(bigwigFile)), silent = TRUE))
  skip_if(inherits(bigwigReadable, "try-error"), "bigWig reading is not supported here")

  blacklist <- GenomicRanges::GRanges("chrT", IRanges::IRanges(99000, 101000))

  plain <- countBigwig(listRegions(), bigwigFiles = bigwigFile, sampleNames = "first", verbose = FALSE)
  listed <- countBigwig(listRegions(), bigwigFiles = bigwigFile, sampleNames = "first", blacklist = blacklist, verbose = FALSE)

  expect_identical(nrow(plain), 9L)
  expect_identical(nrow(listed), 8L)
  expect_length(listed@blacklist, 1)
  expect_identical(listed@filtering.log$step, "blacklist")

  # The regions left carry the signal they had
  expect_equal(SummarizedExperiment::assay(listed, "counts")[, 1],
               SummarizedExperiment::assay(plain, "counts")[rownames(listed), 1])
})
