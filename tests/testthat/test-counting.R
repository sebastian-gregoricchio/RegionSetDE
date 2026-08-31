test_that("countReads returns one row per region and one column per sample", {

  counts <- toyCounts()

  expect_s4_class(counts, "RegionSetDE.counts")
  expect_equal(nrow(counts), 6)
  expect_equal(ncol(counts), 1)
  expect_identical(colnames(counts), "example")

  countMatrix <- SummarizedExperiment::assay(counts, "counts")
  expect_true(all(countMatrix >= 0))
  expect_true(is.numeric(countMatrix))
})


test_that("countReads keeps the region set membership", {

  counts <- toyCounts()
  membership <- SummarizedExperiment::rowData(counts)$region.set

  expect_setequal(unique(as.character(membership)), c("firstSet", "secondSet"))
})


test_that("countReads tiles the regions when asked", {

  plainCounts <- toyCounts()
  tiledCounts <- toyCounts(tileWidth = 100L)

  expect_gt(nrow(tiledCounts), nrow(plainCounts))
})


test_that("countReads carries the sample metadata into colData", {

  metadataTable <- data.frame(sample = "example",
                              condition = "ctrl",
                              stringsAsFactors = FALSE)

  counts <- RegionSetDE::countReads(toyRegionSet(),
                                    bamFiles = toyBamFile(),
                                    sampleNames = "example",
                                    sampleMetadata = metadataTable,
                                    verbose = FALSE)

  expect_true("condition" %in% colnames(SummarizedExperiment::colData(counts)))
  expect_identical(as.character(SummarizedExperiment::colData(counts)$condition), "ctrl")
})


test_that("countReads errors when the file is missing", {

  expect_error(
    RegionSetDE::countReads(toyRegionSet(),
                            bamFiles = "there/is/no/such.bam",
                            verbose = FALSE)
  )
})


test_that("countReads errors when sampleNames does not match the files", {

  expect_error(
    RegionSetDE::countReads(toyRegionSet(),
                            bamFiles = toyBamFile(),
                            sampleNames = c("one", "two"),
                            verbose = FALSE)
  )
})


test_that("countBigwig reads a bigWig into the same shape", {

  skip_if_not_installed("rtracklayer")

  bigwigFile <- file.path(system.file("tests", package = "rtracklayer"), "test.bw")
  skip_if(!file.exists(bigwigFile), "the rtracklayer test bigWig is not available")

  bigwigRanges <- try(rtracklayer::import(bigwigFile), silent = TRUE)
  skip_if(inherits(bigwigRanges, "try-error"), "bigWig reading is not supported here")

  exampleRegions <- GenomicRanges::reduce(bigwigRanges)
  exampleRegions$setName <- "covered"

  exampleSets <- splitLoadRegions(exampleRegions, splitBy = "setName",
                                  seqlevelsStyle = NULL, verbose = FALSE)

  signal <- countBigwig(exampleSets, bigwigFiles = bigwigFile,
                        sampleNames = "example", verbose = FALSE)

  expect_s4_class(signal, "RegionSetDE.counts")
  expect_equal(ncol(signal), 1)
})


test_that("the packaged counts already carry their background bins", {

  counts <- exampleCounts()
  backgroundBins <- S4Vectors::metadata(counts)$background

  expect_false(is.null(backgroundBins))
  expect_gt(nrow(backgroundBins), 0)
  expect_equal(ncol(backgroundBins), ncol(counts))
})


test_that("loadCounts rebuilds an object from a plain table", {

  counts <- exampleCounts()
  regionRanges <- SummarizedExperiment::rowRanges(counts)

  countTable <- data.frame(
    seqnames = as.character(GenomicRanges::seqnames(regionRanges)),
    start = GenomicRanges::start(regionRanges),
    end = GenomicRanges::end(regionRanges),
    SummarizedExperiment::assay(counts, "counts"),
    check.names = FALSE)

  regions <- splitLoadRegions(regionRanges, splitBy = "region.set",
                              genomeAssembly = "rn4", verbose = FALSE)

  reloaded <- loadCounts(regions, counts = countTable, verbose = FALSE)

  expect_s4_class(reloaded, "RegionSetDE.counts")
  expect_equal(nrow(reloaded), nrow(counts))
  expect_equal(ncol(reloaded), ncol(counts))
  expect_equal(unname(colSums(SummarizedExperiment::assay(reloaded, "counts"))),
               unname(colSums(SummarizedExperiment::assay(counts, "counts"))))
})


test_that("selectSamples keeps the requested samples only", {

  counts <- exampleCounts()
  brownNorway <- selectSamples(counts, condition == "BN", verbose = FALSE)

  expect_equal(ncol(brownNorway), 2)
  expect_setequal(
    as.character(SummarizedExperiment::colData(brownNorway)$condition), "BN")
  expect_equal(nrow(brownNorway), nrow(counts))
})


test_that("splitSamples returns one object per level", {

  byStrain <- splitSamples(exampleCounts(), by = "condition", verbose = FALSE)

  expect_type(byStrain, "list")
  expect_setequal(names(byStrain), c("BN", "SHR"))
  expect_s4_class(byStrain[[1]], "RegionSetDE.counts")
})


test_that("asDGEList carries the library sizes across", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)
  dgeList <- asDGEList(counts)

  expect_s4_class(dgeList, "DGEList")
  expect_equal(nrow(dgeList$samples), ncol(counts))
})
