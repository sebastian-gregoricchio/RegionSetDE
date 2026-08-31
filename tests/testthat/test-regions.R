test_that("loadRegions builds an object from a named list of data.frames", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  regionList <- split(regionTable[, c("seqnames", "start", "end")], regionTable$setName)

  regions <- loadRegions(regionList, genomeAssembly = "rn4", verbose = FALSE)

  expect_s4_class(regions, "RegionSetDE")
  expect_setequal(regionSetNames(regions), names(regionList))
  expect_equal(sum(lengths(regions@regions)), nrow(regionTable))
})


test_that("loadRegions accepts GRanges and keeps the assembly", {

  regions <- loadRegions(list(firstSet = toyRegions()),
                         seqlevelsStyle = NULL,
                         genomeAssembly = "toy",
                         verbose = FALSE)

  expect_s4_class(regions, "RegionSetDE")
  expect_identical(regions@genome.assembly, "toy")
  expect_length(regions@regions, 1)
})


test_that("loadRegions rejects duplicated set names", {

  regionRanges <- toyRegions()

  expect_error(
    suppressWarnings(
      loadRegions(list(sameName = regionRanges, sameName = regionRanges),
                  seqlevelsStyle = NULL, verbose = FALSE))
  )
})


test_that("splitLoadRegions splits on a metadata column", {

  regions <- toyRegionSet()

  expect_s4_class(regions, "RegionSetDE")
  expect_setequal(regionSetNames(regions), c("firstSet", "secondSet"))
  expect_equal(unname(lengths(regions@regions)), c(3L, 3L))
})


test_that("splitLoadRegions drops sets below minRegionsPerSet", {

  regionRanges <- toyRegions()
  regionRanges$setName <- c(rep("bigSet", 5), "tinySet")

  expect_warning(
    regions <- splitLoadRegions(regionRanges,
                                splitBy = "setName",
                                minRegionsPerSet = 2,
                                seqlevelsStyle = NULL,
                                verbose = FALSE),
    "discarded")

  expect_setequal(regionSetNames(regions), "bigSet")
})


test_that("splitLoadRegions errors on a column that is not there", {

  expect_error(
    splitLoadRegions(toyRegions(), splitBy = "notAColumn",
                     seqlevelsStyle = NULL, verbose = FALSE)
  )
})


# The filter tests run on the packaged rn4 regions rather than on the toy
# contigs, so that the chromosome naming style is never in question.
exampleRegionSet <- function() {
  regionTable <- loadExampleData("regions", verbose = FALSE)

  splitLoadRegions(
    GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
    splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)
}


test_that("applyBlacklist removes the overlapping regions and logs it", {

  exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)
  regions <- exampleRegionSet()

  before <- sum(lengths(regions@regions))
  filtered <- applyBlacklist(regions, blacklist = exclusionRegions, verbose = FALSE)
  after <- sum(lengths(filtered@regions))

  expect_lt(after, before)
  expect_s4_class(filtered, "RegionSetDE")
  expect_gt(nrow(filtered@filtering.log), 0)
})


test_that("applyBlacklist leaves the regions alone when nothing overlaps", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  regionRanges <- GenomicRanges::makeGRangesFromDataFrame(regionTable)

  # The widest gap between two regions is guaranteed to overlap nothing
  gapRanges <- GenomicRanges::setdiff(range(regionRanges), regionRanges)
  emptyBlacklist <- gapRanges[which.max(GenomicRanges::width(gapRanges))]

  regions <- exampleRegionSet()
  filtered <- applyBlacklist(regions, blacklist = emptyBlacklist, verbose = FALSE)

  expect_equal(lengths(filtered@regions), lengths(regions@regions))
})


test_that("applyWhitelist keeps only what falls inside", {

  regions <- exampleRegionSet()

  whitelist <- GenomicRanges::GRanges(
    seqnames = "chr12",
    ranges = IRanges::IRanges(start = 1, end = 25e6))

  filtered <- applyWhitelist(regions, whitelist = whitelist, verbose = FALSE)

  keptRanges <- unlist(GenomicRanges::GRangesList(filtered@regions))

  expect_true(all(GenomicRanges::end(keptRanges) <= 25e6))
  expect_lt(length(keptRanges), sum(lengths(regions@regions)))

  # Every set has to survive, otherwise the default emptySets would have stopped
  expect_setequal(regionSetNames(filtered), regionSetNames(regions))
})


test_that("emptySets stops when a set is wiped out", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  regions <- exampleRegionSet()

  # The sets are disjoint, so whitelisting one of them empties the other three
  whitelist <- GenomicRanges::makeGRangesFromDataFrame(
    regionTable[regionTable$setName == "promoterCpG", ])

  expect_error(
    applyWhitelist(regions, whitelist = whitelist, emptySets = "stop", verbose = FALSE)
  )
})


test_that("emptySets removes them instead when asked", {

  regionTable <- loadExampleData("regions", verbose = FALSE)
  regions <- exampleRegionSet()

  whitelist <- GenomicRanges::makeGRangesFromDataFrame(
    regionTable[regionTable$setName == "promoterCpG", ])

  expect_warning(
    filtered <- applyWhitelist(regions, whitelist = whitelist,
                               emptySets = "remove", verbose = FALSE)
  )

  expect_setequal(regionSetNames(filtered), "promoterCpG")
})


test_that("renameBedColumns gives the first three columns usable names", {

  bedTable <- data.frame(V1 = "chr12",
                         V2 = c(1000, 5000),
                         V3 = c(2000, 6000),
                         V4 = c("peak_1", "peak_2"),
                         V5 = c(120, 340),
                         V6 = c("+", "-"))

  renamed <- renameBedColumns(bedTable, bedFormat = 6)

  expect_equal(ncol(renamed), 6)
  expect_false(any(colnames(renamed)[1:3] == colnames(bedTable)[1:3]))

  # The point of the renaming is that the table becomes a GRanges
  expect_s4_class(GenomicRanges::makeGRangesFromDataFrame(renamed), "GRanges")
})


test_that("renameBedColumns leaves the columns beyond the format alone", {

  bedTable <- data.frame(V1 = "chr12",
                         V2 = c(1000, 5000),
                         V3 = c(2000, 6000),
                         V4 = c("peak_1", "peak_2"),
                         V5 = c(120, 340),
                         V6 = c("+", "-"))

  partial <- renameBedColumns(bedTable, bedFormat = 3)

  expect_equal(colnames(partial)[4:6], colnames(bedTable)[4:6])
})


test_that("renameBedColumns rejects a format that is not BED", {
  bedTable <- data.frame(V1 = "chr12", V2 = 1000, V3 = 2000)
  expect_error(renameBedColumns(bedTable, bedFormat = 4))
})
