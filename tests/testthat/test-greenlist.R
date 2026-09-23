# Two libraries over the contigs of the example alignment, and a greenlist beside the regions they carry
greenlistCountsObject <- function() {
  RegionSetDE::countReads(toyRegionSet(),
                          bamFiles = rep(toyBamFile(), 2),
                          sampleNames = c("first", "second"),
                          sampleMetadata = data.frame(sample = c("first", "second"),
                                                      condition = c("A", "B")),
                          pairedEnd = FALSE,
                          verbose = FALSE)
}


toyGreenlist <- function() {
  GenomicRanges::GRanges(seqnames = rep(c("seq1", "seq2"), each = 5),
                         ranges = IRanges::IRanges(start = rep(seq(100L, 1300L, by = 300L), 2), width = 200))
}


test_that("the lists shipped with the package are indexed and readable", {

  listIndex <- availableRegionLists()

  expect_s3_class(listIndex, "data.frame")
  expect_true(all(c("type", "genome", "assay", "source", "version", "n.regions", "covered.bp") %in% colnames(listIndex)))
  expect_setequal(unique(listIndex$type), c("blacklist", "greenlist"))
  expect_false("file" %in% colnames(listIndex))

  expect_identical(unique(availableRegionLists(type = "greenlist")$type), "greenlist")
  expect_identical(unique(availableRegionLists(genome = "GRCh38")$genome), "hg38")

  expect_error(availableRegionLists(type = "greylist"), "'blacklist' or 'greenlist'")
})


test_that("the blacklists come back as GRanges carrying their provenance", {

  blacklist <- loadBlacklist("hg38", verbose = FALSE)

  expect_s4_class(blacklist, "GRanges")
  expect_length(blacklist, 636L)
  expect_true(all(c("High Signal Region", "Low Mappability") %in% blacklist$name))
  expect_identical(unname(GenomeInfoDb::genome(blacklist))[1], "hg38")

  listMetadata <- S4Vectors::metadata(blacklist)
  expect_identical(listMetadata$source, "ENCODE")
  expect_identical(listMetadata$assay, "any")
  expect_match(listMetadata$reference, "Amemiya")

  # The assembly aliases reach the same file, and the style of the names can be changed on the way out
  expect_identical(length(loadBlacklist("GRCh38", verbose = FALSE)), length(blacklist))
  expect_false(any(grepl("^chr", GenomeInfoDb::seqlevels(loadBlacklist("hg38", seqlevelsStyle = "Ensembl", verbose = FALSE)))))

  # Naming an assay gives the high signal regions of that protocol instead of the ENCODE list
  cutrunBlacklist <- loadBlacklist("hg38", assay = "cutrun", verbose = FALSE)
  expect_identical(S4Vectors::metadata(cutrunBlacklist)$source, "deMello")
  expect_false(length(cutrunBlacklist) == length(blacklist))

  # A genome carrying several lists is told apart by the source
  expect_identical(length(loadBlacklist("hg38", source = "ENCODE", verbose = FALSE)), length(blacklist))
  expect_error(loadBlacklist("hg38", source = "kundaje"), "Available: ENCODE, deMello")

  # T2T answers to every name it goes by, and ENCODE never covered it, so its list comes from elsewhere
  t2tBlacklist <- loadBlacklist("T2T-CHM13v2.0", verbose = FALSE)

  expect_length(t2tBlacklist, 3565L)
  expect_identical(S4Vectors::metadata(t2tBlacklist)$source, "excluderanges")
  expect_identical(unname(GenomeInfoDb::genome(t2tBlacklist))[1], "hs1")
  expect_identical(length(loadBlacklist("chm13", verbose = FALSE)), length(t2tBlacklist))
  expect_identical(length(loadBlacklist("hs1", verbose = FALSE)), length(t2tBlacklist))
})


test_that("the greenlists are named per assay and refuse to be guessed", {

  greenlist <- loadGreenlist("hg38", assay = "cutrun", verbose = FALSE)

  expect_s4_class(greenlist, "GRanges")
  expect_length(greenlist, 869L)
  expect_identical(S4Vectors::metadata(greenlist)$assay, "cutrun")

  # 'CUT&RUN', 'CUTnRUN' and 'cut_run' all name the same list
  expect_identical(length(loadGreenlist("hg38", assay = "CUT&RUN", verbose = FALSE)), length(greenlist))
  expect_false(length(loadGreenlist("hg38", assay = "CUT&Tag", verbose = FALSE)) == length(greenlist))

  expect_error(loadGreenlist("hg38"), "Name the assay")
  expect_error(loadGreenlist("mm10", assay = "cutrun"), "No greenlist is shipped for the genome 'mm10'")
  expect_error(loadGreenlist("hg38", assay = "atac"), "the assay 'atac'")
  expect_message(loadGreenlist("hg38", assay = "cutrun"), "869 regions")
})


test_that("a list built for another assembly is refused", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
                              splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

  # Rat and human both have a chr1, so the overlap would run and remove regions for no reason
  expect_error(applyBlacklist(regions, blacklist = loadBlacklist("hg38", verbose = FALSE), verbose = FALSE),
               "was built for hg38")

  forcedBlacklist <- loadBlacklist("hg38", verbose = FALSE)
  GenomeInfoDb::genome(forcedBlacklist) <- NA

  expect_s4_class(applyBlacklist(regions, blacklist = forcedBlacklist, verbose = FALSE), "RegionSetDE")
})


test_that("countGreenlist stores the counts beside the regions", {

  counts <- greenlistCountsObject()
  counts <- countGreenlist(counts, greenlist = toyGreenlist(), pairedEnd = FALSE, verbose = FALSE)

  greenlistCounts <- S4Vectors::metadata(counts)$greenlist

  expect_s4_class(greenlistCounts, "RangedSummarizedExperiment")
  expect_identical(colnames(greenlistCounts), colnames(counts))
  expect_identical(nrow(greenlistCounts), 10L)
  expect_true(all(SummarizedExperiment::assay(greenlistCounts, "counts") > 0))

  # The two libraries are the same file, so every region carries the same count twice
  countMatrix <- SummarizedExperiment::assay(greenlistCounts, "counts")
  expect_equal(countMatrix[, 1], countMatrix[, 2], ignore_attr = TRUE)

  countingParameters <- counts@parameters$countGreenlist
  expect_identical(countingParameters$n.regions, 10L)
  expect_identical(countingParameters$covered.bp, 2000)

  expect_error(countGreenlist(counts, greenlist = GenomicRanges::GRanges("chrZZ", IRanges::IRanges(1, 100)), verbose = FALSE),
               "cannot be reconciled")
  expect_error(countGreenlist(counts, greenlist = 42, verbose = FALSE), "GRanges or the path to a BED file")
})


test_that("the greenlist factors follow the median of ratios", {

  counts <- greenlistCountsObject()
  countedCounts <- countGreenlist(counts, greenlist = toyGreenlist(), pairedEnd = FALSE, verbose = FALSE)

  # Two copies of one library have nothing to scale between them
  normalized <- normalizeCounts(countedCounts, method = "greenlist", verbose = FALSE)
  expect_equal(normalized$scaling.factor, c(1, 1), ignore_attr = TRUE)
  expect_identical(S4Vectors::metadata(normalized)$normalization$method, "greenlist")

  # A sample twice as deep over the whole list gets twice the factor
  greenlistMatrix <- matrix(c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
                              20, 40, 60, 80, 100, 120, 140, 160, 180, 200),
                            ncol = 2, dimnames = list(NULL, c("first", "second")))

  suppliedCounts <- normalizeCounts(counts, method = "greenlist", greenlistCounts = greenlistMatrix, verbose = FALSE)
  expect_equal(suppliedCounts$scaling.factor, c(2/3, 4/3), ignore_attr = TRUE)

  # One total per sample leaves the sum as the only estimator, whatever was asked for
  suppliedTotals <- normalizeCounts(counts, method = "greenlist",
                                    greenlistCounts = c(first = 1000, second = 3000), verbose = FALSE)
  expect_equal(suppliedTotals$scaling.factor, c(0.5, 1.5), ignore_attr = TRUE)

  for (estimator in c("medianRatio", "TMM", "sum")) {
    estimated <- normalizeCounts(countedCounts, method = "greenlist", greenlistEstimator = estimator, verbose = FALSE)
    expect_equal(estimated$scaling.factor, c(1, 1), ignore_attr = TRUE, info = estimator)
  }

  expect_error(normalizeCounts(counts, method = "greenlist", verbose = FALSE), "run 'countGreenlist'")
  expect_error(normalizeCounts(countedCounts, method = "greenlist", greenlistEstimator = "median", verbose = FALSE),
               "'medianRatio', 'TMM' or 'sum'")
})


test_that("the greenlist joins the other methods in the normalisation comparison", {

  counts <- countGreenlist(greenlistCountsObject(), greenlist = toyGreenlist(), pairedEnd = FALSE, verbose = FALSE)

  comparisonTable <- plotNormComparison(counts, methods = c("librarySize", "greenlist"), returnData = TRUE)

  expect_s3_class(comparisonTable, "data.frame")
  expect_setequal(unique(comparisonTable$method), c("librarySize", "greenlist"))
})
