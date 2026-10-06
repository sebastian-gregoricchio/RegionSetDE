# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

## Two things are tested here. A single GRanges handed to the counting functions has to behave as a
## region set of one. And files that name the chromosomes of one assembly in different ways, chrT in
## one header and T in the next, have to be counted together: before, the counting stopped at the
## first header that differed. The twin files below hold the same reads under the two names, so the
## answer is known in advance: their columns must be equal.

twinUCSC <- syntheticInput("twin_ucsc", pileStart = 60000L, chromosome = "chrT", seed = 1L)
twinEnsembl <- syntheticInput("twin_ensembl", pileStart = 60000L, chromosome = "T", seed = 1L)
otherUCSC <- syntheticInput("other_ucsc", pileStart = 120000L, chromosome = "chrT", seed = 2L)

siteRanges <- GenomicRanges::GRanges("chrT", IRanges::IRanges(start = seq(10000L, 190000L, by = 10000L), width = 2000L))
siteRegions <- loadRegions(list(sites = siteRanges), seqlevelsStyle = NULL, verbose = FALSE)


test_that("the same chromosome is found under the names the styles give it", {

  translate <- RegionSetDE:::.translateChromosomeNames

  expect_identical(translate(c("chr1", "chrX", "chrM"), c("1", "X", "MT")), c("1", "X", "MT"))
  expect_identical(translate(c("1", "X", "MT"), c("chr1", "chrX", "chrM")), c("chr1", "chrX", "chrM"))

  # Each name is settled on its own, so a vector mixing the styles is no harder than a uniform one
  expect_identical(translate(c("chr1", "2", "scaffold_7"), c("chr1", "chr2")), c("chr1", "chr2", NA))

  # The mitochondrion goes under four names
  expect_identical(translate(c("chrM", "M", "MT"), c("chr1", "chrMT")), rep("chrMT", 3))

  expect_identical(translate(character(0), "chr1"), character(0))
})


test_that("ranges are renamed chromosome by chromosome", {

  mixedRanges <- GenomicRanges::GRanges(c("chr1", "1", "chr2", "scaffold_7"), IRanges::IRanges(c(1, 50, 5, 7), width = 10), score = 1:4)
  names(mixedRanges) <- letters[1:4]

  # chr1 and 1 are one chromosome for a file that calls it 1, and the scaffold the file lacks keeps its name
  renamedRanges <- RegionSetDE:::.matchSeqlevels(mixedRanges, c("1", "2", "MT"), verbose = FALSE)

  expect_identical(as.character(GenomeInfoDb::seqnames(renamedRanges)), c("1", "1", "2", "scaffold_7"))
  expect_identical(names(renamedRanges), letters[1:4])
  expect_identical(renamedRanges$score, 1:4)

  # Nothing in common under any name is still an error
  expect_error(RegionSetDE:::.matchSeqlevels(GenomicRanges::GRanges("seq1", IRanges::IRanges(1, 10)), c("1", "2"), verbose = FALSE),
               "cannot be reconciled")
})


test_that("a single GRanges is counted as a region set of one", {

  expect_message(rangeCounts <- countReads(siteRanges, bamFiles = twinUCSC, sampleNames = "twin"),
                 "loaded with loadRegions\\(\\) as a single set named 'regions'")

  setCounts <- countReads(siteRegions, bamFiles = twinUCSC, sampleNames = "twin", verbose = FALSE)

  expect_s4_class(rangeCounts, "RegionSetDE.counts")
  expect_identical(unname(countTable(rangeCounts, format = "matrix")), unname(countTable(setCounts, format = "matrix")))
  expect_identical(unique(SummarizedExperiment::rowData(rangeCounts)$region.set), "regions")

  # The rows follow the ranges as they were given, which is what lets them be matched back
  reversedCounts <- countReads(rev(siteRanges), bamFiles = twinUCSC, sampleNames = "twin", verbose = FALSE)
  expect_identical(BiocGenerics::start(SummarizedExperiment::rowRanges(reversedCounts)), rev(BiocGenerics::start(siteRanges)))
  expect_identical(unname(countTable(reversedCounts, format = "matrix")), unname(countTable(setCounts, format = "matrix"))[rev(seq_along(siteRanges)), , drop = FALSE])

  # Ranges given twice would be two rows with one name
  expect_identical(nrow(countReads(c(siteRanges[1:3], siteRanges[1:5]), bamFiles = twinUCSC, verbose = FALSE)), 5L)

  expect_error(countReads(siteRanges[0], bamFiles = twinUCSC, verbose = FALSE), "holds no region")

  # A table of counts is attached to a GRanges the same way
  tableCounts <- loadCounts(siteRanges,
                            counts = data.frame(chr = "chrT", start = BiocGenerics::start(siteRanges), end = BiocGenerics::end(siteRanges),
                                                twin = countTable(setCounts, format = "matrix")[, 1]),
                            verbose = FALSE)
  expect_identical(unname(countTable(tableCounts, format = "matrix")), unname(countTable(setCounts, format = "matrix")))
})


test_that("BAM files naming the chromosomes differently are counted together", {

  expect_message(mixedCounts <- countReads(siteRegions, bamFiles = c(twinUCSC, twinEnsembl, otherUCSC), sampleNames = c("ucsc", "ensembl", "other")),
                 "do not name their chromosomes alike")

  mixedMatrix <- countTable(mixedCounts, format = "matrix")

  # Same reads, two names: the same counts and the same library
  expect_identical(mixedMatrix[, "ucsc"], mixedMatrix[, "ensembl"])
  expect_identical(mixedCounts$library.size[1], mixedCounts$library.size[2])
  expect_gt(sum(mixedMatrix[, "ensembl"]), 0)

  # Each file alone gives what it gives in the mix
  aloneCounts <- countReads(siteRegions, bamFiles = c(twinUCSC, otherUCSC), sampleNames = c("ucsc", "other"), verbose = FALSE)
  expect_identical(mixedMatrix[, c("ucsc", "other")], countTable(aloneCounts, format = "matrix"))
  expect_identical(mixedCounts$library.size[c(1, 3)], aloneCounts$library.size)

  # Which file comes first decides the names the regions are converted to, not the counts
  reorderedCounts <- countReads(siteRegions, bamFiles = c(twinEnsembl, otherUCSC, twinUCSC), sampleNames = c("ensembl", "other", "ucsc"), verbose = FALSE)
  expect_identical(countTable(reorderedCounts, format = "matrix")[, colnames(mixedMatrix)], mixedMatrix)

  # The object keeps the names of the regions it was given
  expect_identical(GenomeInfoDb::seqlevels(mixedCounts), "chrT")

  # The pieces of the files are shared among the threads under the same names
  threadedCounts <- countReads(siteRegions, bamFiles = c(twinUCSC, twinEnsembl, otherUCSC), sampleNames = c("ucsc", "ensembl", "other"),
                               nThreads = 2, verbose = FALSE)
  expect_identical(countTable(threadedCounts, format = "matrix"), mixedMatrix)
})


test_that("what is given beside the regions follows the files as well", {

  discardedUCSC <- GenomicRanges::GRanges("chrT", IRanges::IRanges(55000L, 70000L))
  discardedEnsembl <- GenomicRanges::GRanges("T", IRanges::IRanges(55000L, 70000L))

  countWithout <- function(discarded, files) {
    countTable(countReads(siteRegions, bamFiles = files, sampleNames = c("first", "second"), discardRegions = discarded, verbose = FALSE),
               format = "matrix")
  }

  # Reads discarded in one style are discarded in the files of the other
  discardedMatrix <- countWithout(discardedUCSC, c(twinUCSC, twinEnsembl))
  expect_identical(discardedMatrix[, 1], discardedMatrix[, 2])
  expect_identical(unname(discardedMatrix), unname(countWithout(discardedEnsembl, c(twinEnsembl, twinUCSC))))
  expect_lt(sum(discardedMatrix), sum(countTable(countReads(siteRegions, bamFiles = c(twinUCSC, twinEnsembl), verbose = FALSE), format = "matrix")))

  # A chromosome excluded under one name leaves the library sizes of every file, here the only one they have
  expect_warning(excludedCounts <- countReads(siteRegions, bamFiles = c(twinUCSC, twinEnsembl), excludeChromosomes = "T", verbose = FALSE),
                 "No fragment entered the library size")
  expect_identical(excludedCounts$library.size, c(0, 0))

  # An input named the other way is counted over the same rows
  inputCounts <- countReads(siteRegions, bamFiles = c(twinUCSC, otherUCSC), sampleNames = c("ucsc", "other"), inputFiles = twinEnsembl, verbose = FALSE)
  expect_identical(SummarizedExperiment::assay(inputCounts, "input")[, "ucsc"], SummarizedExperiment::assay(inputCounts, "counts")[, "ucsc"])
  expect_identical(inputCounts$input.library.size[1], inputCounts$library.size[1])
})


test_that("the functions reading the files after the counting follow them too", {

  mixedCounts <- countReads(siteRegions, bamFiles = c(twinUCSC, twinEnsembl), sampleNames = c("ucsc", "ensembl"), verbose = FALSE)

  # Background bins, named after the regions of the object
  backgroundBins <- S4Vectors::metadata(countBackground(mixedCounts, binSize = 5000, excludeRegions = FALSE, verbose = FALSE))$background
  backgroundMatrix <- SummarizedExperiment::assay(backgroundBins, "counts")

  expect_identical(backgroundMatrix[, 1], backgroundMatrix[, 2])
  expect_gt(sum(backgroundMatrix), 0)
  expect_identical(GenomeInfoDb::seqlevels(backgroundBins), "chrT")

  # Summits, placed on the pileup of both files
  summitCounts <- countReads(siteRegions, bamFiles = c(twinUCSC, twinEnsembl), summits = 0, verbose = FALSE)
  aloneSummits <- countReads(siteRegions, bamFiles = twinUCSC, summits = 0, verbose = FALSE)
  expect_identical(SummarizedExperiment::rowData(summitCounts)$summit, SummarizedExperiment::rowData(aloneSummits)$summit)

  # Fragment lengths, estimated over regions written in the style of the first file only
  lengthTable <- suppressWarnings(estimateFragmentLength(c(twinUCSC, twinEnsembl), regions = siteRanges, verbose = FALSE))$table
  expect_identical(lengthTable$n.reads[1], lengthTable$n.reads[2])
  expect_gt(lengthTable$n.reads[2], 0)

  # The greylist does not change when one of the inputs is named the other way
  mixedGreylist <- suppressWarnings(makeGreylist(c(otherUCSC, twinEnsembl), verbose = FALSE))
  uniformGreylist <- suppressWarnings(makeGreylist(c(otherUCSC, twinUCSC), verbose = FALSE))
  expect_identical(BiocGenerics::start(mixedGreylist), BiocGenerics::start(uniformGreylist))
  expect_gt(length(mixedGreylist), 0)
})


test_that("contigs some files lack are no reason to stop, another assembly is", {

  # chrU is in one header only, and holds no read
  withContig <- syntheticInput("twin_contig", pileStart = 60000L, chromosome = "T", emptyContig = TRUE, seed = 1L)

  expect_message(contigCounts <- countReads(siteRegions, bamFiles = c(twinUCSC, withContig), sampleNames = c("plain", "contig")),
                 "missing from at least one BAM file")

  contigMatrix <- countTable(contigCounts, format = "matrix")
  expect_identical(contigMatrix[, "plain"], contigMatrix[, "contig"])
  expect_identical(contigCounts$library.size[1], contigCounts$library.size[2])

  # The same chromosome with another length is another assembly
  shorterAssembly <- syntheticInput("twin_short", chromosome = "T", contigLength = 150000L, seed = 1L)
  expect_error(countReads(siteRegions, bamFiles = c(twinUCSC, shorterAssembly), verbose = FALSE), "not aligned to the same assembly")
  expect_error(countReads(siteRegions, bamFiles = twinUCSC, inputFiles = shorterAssembly, verbose = FALSE), "not aligned to the same assembly")

  # And a file sharing no chromosome under any name has nothing to do with the others
  expect_error(countReads(siteRegions, bamFiles = c(twinUCSC, toyBamFile()), verbose = FALSE), "cannot be reconciled")
})


test_that("a filter follows the names of sets that declare no style", {

  ensemblRanges <- GenomicRanges::GRanges("T", IRanges::ranges(siteRanges))
  ensemblRegions <- loadRegions(list(sites = ensemblRanges), seqlevelsStyle = NULL, verbose = FALSE)

  firstHalfUCSC <- GenomicRanges::GRanges("chrT", IRanges::IRanges(1L, 100000L))
  firstHalfEnsembl <- GenomicRanges::GRanges("T", IRanges::IRanges(1L, 100000L))

  # The list written as the regions used to be renamed to UCSC and then refused
  ensemblFiltered <- applyBlacklist(ensemblRegions, blacklist = firstHalfEnsembl, verbose = FALSE)
  ucscFiltered <- applyBlacklist(ensemblRegions, blacklist = firstHalfUCSC, verbose = FALSE)

  expect_identical(length(ensemblFiltered@regions$sites), length(ucscFiltered@regions$sites))
  expect_lt(length(ucscFiltered@regions$sites), length(ensemblRanges))
  expect_identical(GenomeInfoDb::seqlevels(ucscFiltered@regions$sites), "T")

  # Plain ranges declare no style either
  expect_identical(length(applyBlacklist(ensemblRanges, blacklist = firstHalfUCSC, verbose = FALSE)), length(ucscFiltered@regions$sites))
  expect_identical(length(applyWhitelist(GenomicRanges::GRangesList(sites = ensemblRanges), whitelist = firstHalfUCSC, verbose = FALSE)$sites),
                   length(ensemblRanges) - length(ucscFiltered@regions$sites))
})


test_that("a count table written with other chromosome names finds its regions", {

  setCounts <- countReads(siteRegions, bamFiles = c(twinUCSC, otherUCSC), sampleNames = c("twin", "other"), verbose = FALSE)

  ensemblTable <- data.frame(chr = "T", start = BiocGenerics::start(siteRanges), end = BiocGenerics::end(siteRanges),
                             countTable(setCounts, format = "matrix"))

  tableCounts <- loadCounts(siteRegions, counts = ensemblTable, verbose = FALSE)

  expect_identical(unname(countTable(tableCounts, format = "matrix")), unname(countTable(setCounts, format = "matrix")))
})


test_that("bigWig files naming the chromosomes differently are read together", {

  # The coverage of one file, written once under each name
  chromosomeLengths <- Rsamtools::scanBamHeader(twinUCSC)[[1]]$targets
  twinCoverage <- RegionSetDE:::.windowCoverage(bamFile = twinUCSC,
                                                windows = GenomicRanges::GRanges("chrT", IRanges::IRanges(1L, chromosomeLengths[["chrT"]])),
                                                chromosomeLengths = chromosomeLengths, isPairedEnd = FALSE, fragmentLength = 150,
                                                maxFragmentLength = 1000, minMapq = 20, removeDuplicates = TRUE)$coverage

  ucscBigwig <- file.path(tempdir(), "twin_ucsc.bw")
  ensemblBigwig <- file.path(tempdir(), "twin_ensembl.bw")

  rtracklayer::export.bw(twinCoverage, ucscBigwig)
  names(twinCoverage) <- "T"
  rtracklayer::export.bw(twinCoverage, ensemblBigwig)

  # rtracklayer writes bigWig files on Windows but cannot read them back from an absolute path
  bigwigReadable <- suppressWarnings(try(GenomeInfoDb::seqlengths(rtracklayer::BigWigFile(ucscBigwig)), silent = TRUE))
  skip_if(inherits(bigwigReadable, "try-error"), "bigWig reading is not supported here")

  expect_message(mixedSignal <- countBigwig(siteRanges, bigwigFiles = c(ensemblBigwig, ucscBigwig), sampleNames = c("ensembl", "ucsc")),
                 "do not name their chromosomes alike")

  signalMatrix <- SummarizedExperiment::assay(mixedSignal, "counts")

  expect_identical(signalMatrix[, "ensembl"], signalMatrix[, "ucsc"])
  expect_gt(sum(signalMatrix), 0)
  expect_identical(GenomeInfoDb::seqlevels(mixedSignal), "chrT")

  # The profiles read the files the same way
  setCounts <- countReads(siteRegions, bamFiles = c(twinEnsembl, twinUCSC), sampleNames = c("ensembl", "ucsc"), verbose = FALSE)
  profileWindows <- GenomicRanges::resize(SummarizedExperiment::rowRanges(setCounts), width = 2000, fix = "center")

  bigwigProfiles <- RegionSetDE:::.binnedSignal(files = c(ensemblBigwig, ucscBigwig), type = "bigwig", windows = profileWindows, binWidth = 100, counts = setCounts)
  bamProfiles <- RegionSetDE:::.binnedSignal(files = c(twinEnsembl, twinUCSC), type = "bam", windows = profileWindows, binWidth = 100, counts = setCounts)

  expect_identical(bigwigProfiles$matrices[[1]], bigwigProfiles$matrices[[2]])
  expect_identical(bamProfiles$matrices[[1]], bamProfiles$matrices[[2]])
  expect_gt(sum(bamProfiles$matrices[[1]]), 0)
})
