## Loading, export and coercion options, plus the show methods, which are printed
## by every example but exercised by no assertion.


test_that("loadRegions returns each of its output formats", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  regionList <- split(regionTable[, c("seqnames", "start", "end")], regionTable$setName)

  asObject <- loadRegions(regionList, genomeAssembly = "rn4", verbose = FALSE)
  asGRangesList <- loadRegions(regionList, genomeAssembly = "rn4",
                               outputFormat = "GRangesList", verbose = FALSE)
  asList <- loadRegions(regionList, genomeAssembly = "rn4",
                        outputFormat = "list", verbose = FALSE)

  expect_s4_class(asObject, "RegionSetDE")
  expect_s4_class(asGRangesList, "GRangesList")
  expect_type(asList, "list")

  expect_equal(length(asGRangesList), length(regionList))
})


test_that("loadRegions handles identical sets three ways", {

  regionRanges <- toyRegions()

  expect_error(
    suppressWarnings(loadRegions(list(first = regionRanges, second = regionRanges),
                                 seqlevelsStyle = NULL, verbose = FALSE))
  )

  removed <- suppressWarnings(
    loadRegions(list(first = regionRanges, second = regionRanges),
                duplicatedSets = "remove", seqlevelsStyle = NULL, verbose = FALSE))

  kept <- suppressWarnings(
    loadRegions(list(first = regionRanges, second = regionRanges),
                duplicatedSets = "keep", seqlevelsStyle = NULL, verbose = FALSE))

  expect_length(removed@regions, 1)
  expect_length(kept@regions, 2)
})


test_that("loadRegions can reduce and can leave the metadata behind", {

  regionRanges <- toyRegions()

  reduced <- loadRegions(list(firstSet = regionRanges), reduceRegions = TRUE,
                         seqlevelsStyle = NULL, verbose = FALSE)
  withoutMetadata <- loadRegions(list(firstSet = regionRanges), keepMetadata = FALSE,
                                 seqlevelsStyle = NULL, verbose = FALSE)
  unsorted <- loadRegions(list(firstSet = regionRanges), sortRegions = FALSE,
                          seqlevelsStyle = NULL, verbose = FALSE)

  expect_s4_class(reduced, "RegionSetDE")
  expect_null(withoutMetadata@regions$firstSet$setName)
  expect_s4_class(unsorted, "RegionSetDE")
})


test_that("splitLoadRegions can pick and cap the sets", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  regionRanges <- GenomicRanges::makeGRangesFromDataFrame(regionTable,
                                                          keep.extra.columns = TRUE)

  selected <- splitLoadRegions(regionRanges, splitBy = "setName",
                               selectedSets = c("promoterCpG", "intergenic"),
                               genomeAssembly = "rn4", verbose = FALSE)

  expect_setequal(regionSetNames(selected), c("promoterCpG", "intergenic"))

  expect_error(
    splitLoadRegions(regionRanges, splitBy = "setName", maxSets = 2,
                     genomeAssembly = "rn4", verbose = FALSE)
  )

  withColumn <- splitLoadRegions(regionRanges, splitBy = "setName",
                                 keepSplitColumn = TRUE,
                                 genomeAssembly = "rn4", verbose = FALSE)

  expect_false(is.null(withColumn@regions[[1]]$setName))
})


test_that("applyBlacklist can trim rather than drop", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)

  regions <- splitLoadRegions(
    GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
    splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

  trimmed <- applyBlacklist(regions, blacklist = exclusionRegions,
                            trimRegions = TRUE, verbose = FALSE)
  dropped <- applyBlacklist(regions, blacklist = exclusionRegions,
                            trimRegions = FALSE, verbose = FALSE)

  expect_gte(sum(lengths(trimmed@regions)), sum(lengths(dropped@regions)))
})


test_that("applyBlacklist honours the overlap thresholds", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)

  regions <- splitLoadRegions(
    GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
    splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

  anyOverlap <- applyBlacklist(regions, blacklist = exclusionRegions,
                               minOverlapBp = 1, verbose = FALSE)
  substantialOverlap <- applyBlacklist(regions, blacklist = exclusionRegions,
                                       minOverlapFraction = 0.9, verbose = FALSE)

  expect_lte(sum(lengths(anyOverlap@regions)),
             sum(lengths(substantialOverlap@regions)))
})


test_that("loadCounts matches on coordinates and on identifiers", {

  counts <- exampleCounts()
  regionRanges <- SummarizedExperiment::rowRanges(counts)

  countTable <- data.frame(
    seqnames = as.character(GenomicRanges::seqnames(regionRanges)),
    start = GenomicRanges::start(regionRanges),
    end = GenomicRanges::end(regionRanges),
    region.id = names(regionRanges),
    SummarizedExperiment::assay(counts, "counts"),
    check.names = FALSE)

  regions <- splitLoadRegions(regionRanges, splitBy = "region.set",
                              genomeAssembly = "rn4", verbose = FALSE)

  byCoordinates <- loadCounts(regions, counts = countTable,
                              matchBy = "coordinates", verbose = FALSE)

  expect_s4_class(byCoordinates, "RegionSetDE.counts")
  expect_equal(nrow(byCoordinates), nrow(counts))
})


test_that("loadCounts refuses a table missing rows it needs", {

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

  expect_error(
    loadCounts(regions, counts = utils::head(countTable, 50),
               missingRegions = "stop", verbose = FALSE)
  )
})


test_that("selectSamples takes names as well as an expression", {

  counts <- exampleCounts()

  byExpression <- selectSamples(counts, condition == "BN", verbose = FALSE)
  byName <- selectSamples(counts, samples = colnames(counts)[1:2], verbose = FALSE)

  expect_equal(ncol(byExpression), 2)
  expect_equal(ncol(byName), 2)
})


test_that("selectSamples refuses a sample that is not there", {

  expect_error(selectSamples(exampleCounts(), samples = "notASample", verbose = FALSE))
})


test_that("splitSamples honours minSamples", {

  counts <- exampleCounts()

  bySex <- splitSamples(counts, by = "sex", minSamples = 1, verbose = FALSE)
  expect_gte(length(bySex), 1)

  # The female group holds a single library, so a floor of two drops it
  strictSplit <- splitSamples(counts, by = "sex", minSamples = 2, verbose = FALSE)
  expect_lte(length(strictSplit), length(bySex))
})


test_that("splitSamples refuses a column that is not there", {

  expect_error(splitSamples(exampleCounts(), by = "notAColumn", verbose = FALSE))
})


test_that("asDGEList works with and without the offsets", {

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  withOffsets <- asDGEList(counts, useOffsets = TRUE, verbose = FALSE)
  withoutOffsets <- asDGEList(counts, useOffsets = FALSE, verbose = FALSE)

  expect_equal(nrow(withOffsets$samples), ncol(counts))
  expect_equal(nrow(withoutOffsets$samples), ncol(counts))

  # The registered coercion reaches the same place
  coerced <- methods::as(counts, "DGEList")
  expect_equal(nrow(coerced$samples), ncol(counts))
})


test_that("asDESeqDataSet builds a dataset", {

  skip_if_not_installed("DESeq2")

  counts <- normalizeCounts(exampleCounts(), method = "background", verbose = FALSE)

  deseqData <- asDESeqDataSet(counts, design = ~ condition, verbose = FALSE)

  expect_s4_class(deseqData, "DESeqDataSet")
  expect_equal(ncol(deseqData), ncol(counts))
})


test_that("exportResults writes each of its score mappings", {

  results <- exampleResults()

  for (scoreColumn in c("FDR", "log2FC", "none")) {

    outputDirectory <- file.path(tempdir(), paste0("regionSetDE-score-", scoreColumn))
    dir.create(outputDirectory, showWarnings = FALSE)
    on.exit(unlink(outputDirectory, recursive = TRUE), add = TRUE)

    exportResults(results, path = outputDirectory, prefix = "score",
                  bedScore = scoreColumn, verbose = FALSE)

    expect_gt(length(list.files(outputDirectory)), 0)
  }
})


test_that("exportResults honours its layout arguments", {

  results <- exampleResults()

  outputDirectory <- file.path(tempdir(), "regionSetDE-layout")
  dir.create(outputDirectory, showWarnings = FALSE)
  on.exit(unlink(outputDirectory, recursive = TRUE), add = TRUE)

  exportResults(results, path = outputDirectory, prefix = "plain",
                colourByStatus = FALSE, provenance = FALSE, compress = FALSE,
                verbose = FALSE)

  exportResults(results, path = outputDirectory, prefix = "split",
                splitByDirection = TRUE, verbose = FALSE)

  exportResults(results, path = outputDirectory, prefix = "oneSet",
                set = "promoterCpG", verbose = FALSE)

  writtenFiles <- list.files(outputDirectory)

  expect_true(any(grepl("^plain", writtenFiles)))
  expect_true(any(grepl("^split", writtenFiles)))
  expect_true(any(grepl("^oneSet", writtenFiles)))

  # provenance = FALSE means no parameters file for that prefix
  expect_false(any(grepl("^plain.*parameters", writtenFiles)))
})


test_that("exportResults can write a results list", {

  resultsList <- testRegions(exampleFit(),
                             contrast = list(strainEffect = c("condition", "SHR", "BN"),
                                             reversed = c("condition", "BN", "SHR")),
                             verbose = FALSE)

  outputDirectory <- file.path(tempdir(), "regionSetDE-list-export")
  dir.create(outputDirectory, showWarnings = FALSE)
  on.exit(unlink(outputDirectory, recursive = TRUE), add = TRUE)

  exportResults(resultsList, path = outputDirectory, prefix = "both", verbose = FALSE)

  expect_gt(length(list.files(outputDirectory)), 1)
})


test_that("every class prints without erroring", {

  regionTable <- loadExampleData("regions", verbose = FALSE)

  regions <- splitLoadRegions(
    GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
    splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)

  resultsList <- testRegions(exampleFit(),
                             contrast = list(one = c("condition", "SHR", "BN")),
                             verbose = FALSE)

  printedObjects <- list(regions,
                         exampleCounts(),
                         exampleFit(),
                         exampleResults(),
                         exampleSetResults(),
                         resultsList,
                         makeSetUniverse(exampleFit(), verbose = FALSE))

  for (printedObject in printedObjects) {
    expect_output(print(printedObject))
  }
})
