## The region to gene links, the gene sets built from them, and the sets a set level test runs on.


# A counts object with the same coordinates held by two sets, which is the case the keys exist for
sharedSetCounts <- function() {

  counts <- exampleCounts()
  regions <- SummarizedExperiment::rowRanges(counts)

  sharedRegion <- regions[1]
  S4Vectors::mcols(sharedRegion)$region.set <- "sharedSet"

  sharedCounts <- SummarizedExperiment::rbind(counts, counts[1, ])
  SummarizedExperiment::rowRanges(sharedCounts)[nrow(sharedCounts)] <- sharedRegion

  return(sharedCounts)
}


# The coordinates of the first regions, the shape a contact table carries them in
linkCoordinates <- function(counts, rows) {
  as.data.frame(SummarizedExperiment::rowRanges(counts))[rows, c("seqnames", "start", "end")]
}



test_that("the links are stored whole, and the columns next to the regions are a view of them", {

  counts <- exampleCounts()
  coordinates <- linkCoordinates(counts, 1:4)

  # One region linked to two genes, which is what a bidirectional promoter looks like
  linkTable <- data.frame(rbind(coordinates, coordinates[1, ]),
                          gene.id = c("near", "second", "third", "fourth", "far"),
                          distance.to.tss = c(300L, 1200L, 4000L, 8000L, 95000L))

  annotated <- annotateRegions(counts, linkTable = linkTable, verbose = FALSE)

  expect_equal(nrow(regionLinks(annotated)), 5)
  expect_true(all(c("region", "gene.id", "distance.to.tss") %in% colnames(regionLinks(annotated))))

  # The coordinate columns are consumed into the key rather than kept alongside it
  expect_false(any(c("seqnames", "start", "end") %in% colnames(regionLinks(annotated))))

  # The representative link is the nearest one, and the count says how many were collapsed
  regionTable <- SummarizedExperiment::rowData(annotated)
  expect_identical(regionTable$gene.id[1], "near")
  expect_identical(regionTable$n.links[1], 2L)
  expect_identical(regionTable$n.links[5], 0L)

  expect_identical(annotated@parameters$annotateRegions$source, "table")
  expect_identical(annotated@parameters$annotateRegions$n.genes, 5L)
})


test_that("the regions are named by coordinates, by key or by identifier, and by a GRanges", {

  counts <- exampleCounts()
  coordinates <- linkCoordinates(counts, 1:4)
  geneIds <- c("one", "two", "three", "four")

  byCoordinates <- annotateRegions(counts, linkTable = data.frame(coordinates, gene.id = geneIds), verbose = FALSE)

  byKey <- annotateRegions(counts,
                           linkTable = data.frame(region = regionLinks(byCoordinates)$region, gene.id = geneIds),
                           verbose = FALSE)

  # The identifiers the lookup uses, taken the same way the function takes them, so the test
  # cannot disagree with it about which column holds them
  regionIdentifiers <- RegionSetDE:::.regionIdentifiers(SummarizedExperiment::rowRanges(counts))

  byIdentifier <- annotateRegions(counts,
                                  linkTable = data.frame(region.id = regionIdentifiers[1:4],
                                                         gene.id = geneIds),
                                  verbose = FALSE)

  # ChIPseeker names its columns differently, and returns a GRanges
  asRanges <- GenomicRanges::makeGRangesFromDataFrame(data.frame(coordinates, geneId = geneIds), keep.extra.columns = TRUE)
  byRanges <- annotateRegions(counts, linkTable = asRanges, verbose = FALSE)

  for (annotated in list(byKey, byIdentifier, byRanges)) {
    expect_identical(regionLinks(annotated)$gene.id, regionLinks(byCoordinates)$gene.id)
  }

  # Identifiers from somewhere else are a different failure from coordinates that do not match,
  # and the message has to name the identifiers rather than keys it does not have
  expect_error(annotateRegions(counts,
                               linkTable = data.frame(region.id = paste0("elsewhere", 1:4),
                                                      gene.id = geneIds),
                               verbose = FALSE),
               "Not one 'region.id'")
})


test_that("the links can be built by proximity to gene ranges supplied directly", {

  counts <- exampleCounts()
  regions <- SummarizedExperiment::rowRanges(counts)

  # Start sites a little upstream and downstream of the first regions, on both strands
  geneRanges <- GenomicRanges::GRanges(
    seqnames = GenomeInfoDb::seqnames(regions)[1:3],
    ranges = IRanges::IRanges(start = BiocGenerics::start(regions)[1:3] + 100L, width = 200L),
    strand = c("+", "-", "+"))

  S4Vectors::mcols(geneRanges)$gene.id <- c("first", "second", "third")
  S4Vectors::mcols(geneRanges)$gene.symbol <- c("FIRST", "SECOND", "THIRD")

  annotated <- annotateRegions(counts, genes = geneRanges, verbose = FALSE)
  links <- regionLinks(annotated)

  expect_equal(nrow(links), nrow(counts))
  expect_identical(annotated@parameters$annotateRegions$source, "nearestGene")
  expect_true(all(c("distance.to.tss", "gene.symbol") %in% colnames(links)))
  expect_true(all(links$gene.id %in% c("first", "second", "third")))

  # The start site of a minus strand gene is its 5' end, the higher coordinate
  minusSite <- BiocGenerics::end(geneRanges)[2]
  plusSite <- BiocGenerics::start(geneRanges)[1]
  expect_gt(minusSite, BiocGenerics::start(geneRanges)[2])
  expect_equal(plusSite, BiocGenerics::start(geneRanges)[1])

  # The identifier may be the names of the ranges instead of a column
  namedRanges <- geneRanges
  S4Vectors::mcols(namedRanges) <- NULL
  names(namedRanges) <- c("first", "second", "third")
  expect_equal(nrow(regionLinks(annotateRegions(counts, genes = namedRanges, verbose = FALSE))),
               nrow(counts))

  # Beyond maxDistance the regions stay in the object with no gene
  limited <- annotateRegions(counts, genes = geneRanges, maxDistance = 500, verbose = FALSE)
  expect_equal(nrow(regionLinks(limited)), nrow(counts))
  expect_gt(sum(is.na(regionLinks(limited)$gene.id)), 0)
  expect_equal(nrow(limited), nrow(counts))

  # One source at a time, and the genes have to be addressable
  bareRanges <- geneRanges
  S4Vectors::mcols(bareRanges) <- NULL
  expect_error(annotateRegions(counts, genes = bareRanges, verbose = FALSE), "carry no identifier")
  expect_error(annotateRegions(counts, genes = geneRanges,
                               linkTable = data.frame(linkCoordinates(counts, 1), gene.id = "one"),
                               verbose = FALSE),
               "cannot be used together")
  expect_error(annotateRegions(counts, genes = as.data.frame(geneRanges), verbose = FALSE),
               "must be a GRanges")
})


test_that("a flank links every gene within reach, and tssRegion labels what is a promoter", {

  counts <- exampleCounts()
  regions <- SummarizedExperiment::rowRanges(counts)
  firstCentre <- (BiocGenerics::start(regions)[1] + BiocGenerics::end(regions)[1]) %/% 2L

  # Two start sites straddling the first region, which is the bidirectional promoter case
  # leftGene sits 50 bp from the centre, rightGene is on the minus strand with its start site at
  # its 5' end, 1999 bp away, so a window of 500 separates the two
  geneRanges <- GenomicRanges::GRanges(
    seqnames = rep(GenomeInfoDb::seqnames(regions)[1], 2),
    ranges = IRanges::IRanges(start = c(firstCentre - 50L, firstCentre + 1900L), width = 100L),
    strand = c("+", "-"))
  S4Vectors::mcols(geneRanges)$gene.id <- c("leftGene", "rightGene")

  # The nearest route takes one of them, the flank takes both
  nearestLinks <- regionLinks(annotateRegions(counts, genes = geneRanges, verbose = FALSE))
  firstKey <- unique(nearestLinks$region)[1]
  expect_equal(sum(nearestLinks$region == firstKey), 1)

  flanked <- annotateRegions(counts, genes = geneRanges, flankDistance = 5000, verbose = FALSE)
  flankLinks <- regionLinks(flanked)

  expect_setequal(stats::na.omit(flankLinks$gene.id[flankLinks$region == firstKey]),
                  c("leftGene", "rightGene"))
  expect_equal(flanked@parameters$annotateRegions$flank.distance, 5000)

  # Regions the flank reached nothing on stay in the table, unassigned, so the universe is intact
  expect_equal(length(unique(flankLinks$region)), nrow(counts))
  expect_gt(sum(is.na(flankLinks$gene.id)), 0)

  # The window decides the label and nothing else
  splitLinks <- regionLinks(annotateRegions(counts, genes = geneRanges, flankDistance = 5000,
                                            tssRegion = c(-500, 500), verbose = FALSE))
  firstLinks <- splitLinks[splitLinks$region == firstKey, , drop = FALSE]

  expect_identical(firstLinks$link.type[firstLinks$gene.id == "leftGene"], "promoter")
  expect_identical(firstLinks$link.type[firstLinks$gene.id == "rightGene"], "distal")

  wideLinks <- regionLinks(annotateRegions(counts, genes = geneRanges, flankDistance = 5000,
                                           tssRegion = c(-5000, 5000), verbose = FALSE))
  expect_true(all(wideLinks$link.type[wideLinks$region == firstKey] == "promoter"))

  tightLinks <- regionLinks(annotateRegions(counts, genes = geneRanges, flankDistance = 5000,
                                            tssRegion = c(-10, 10), verbose = FALSE))
  expect_true(all(tightLinks$link.type[tightLinks$region == firstKey] == "distal"))

  # A gene with two start sites, alternative promoters, reached through both by a wide flank:
  # one link per pair, through the nearer site, rather than the gene counting twice
  twoSited <- geneRanges
  S4Vectors::mcols(twoSited)$gene.id <- c("oneGene", "oneGene")

  pairedLinks <- regionLinks(annotateRegions(counts, genes = twoSited, flankDistance = 5000,
                                             verbose = FALSE))
  firstPaired <- pairedLinks[pairedLinks$region == firstKey, , drop = FALSE]

  expect_equal(nrow(firstPaired), 1)
  expect_equal(abs(firstPaired$distance.to.tss), 50)
  expect_equal(anyDuplicated(pairedLinks[, c("region", "gene.id")]), 0)

  # The links stay in the order of the regions rather than being sorted by their key
  expect_identical(unique(pairedLinks$region)[1], firstKey)

  expect_error(annotateRegions(counts, genes = geneRanges, maxDistance = 1000,
                               flankDistance = 1000, verbose = FALSE),
               "cannot be used together")
  expect_error(annotateRegions(counts, genes = geneRanges, tssRegion = 3000, verbose = FALSE),
               "must be two numbers")
})


test_that("the link type splits the gene sets and the gene table", {

  fit <- exampleFit()
  regions <- SummarizedExperiment::rowRanges(fit@counts)

  # One start site on every region, half of them inside a tight promoter window
  geneCentres <- (BiocGenerics::start(regions) + BiocGenerics::end(regions)) %/% 2L
  siteOffset <- rep(c(0L, 20000L), length.out = length(regions))

  geneRanges <- GenomicRanges::GRanges(
    seqnames = GenomeInfoDb::seqnames(regions),
    ranges = IRanges::IRanges(start = geneCentres + siteOffset, width = 100L),
    strand = "+")
  geneIdentifiers <- rep(paste0("gene", seq_len(ceiling(length(regions) / 2))), each = 2)
  S4Vectors::mcols(geneRanges)$gene.id <- geneIdentifiers[seq_along(regions)]

  annotatedFit <- annotateRegions(fit, genes = geneRanges, flankDistance = 50000,
                                  tssRegion = c(-1000, 1000), verbose = FALSE)

  storedTypes <- stats::na.omit(regionLinks(annotatedFit)$link.type)
  expect_true(all(c("promoter", "distal") %in% storedTypes))

  geneSets <- list(everything = unique(S4Vectors::mcols(geneRanges)$gene.id))

  allSets <- geneSetsToRegionSets(annotatedFit, geneSets = geneSets, verbose = FALSE)
  promoterSets <- geneSetsToRegionSets(annotatedFit, geneSets = geneSets,
                                       linkType = "promoter", verbose = FALSE)

  expect_lte(length(promoterSets$everything), length(allSets$everything))
  expect_gt(length(promoterSets$everything), 0)

  expect_error(geneSetsToRegionSets(annotatedFit, geneSets = geneSets,
                                    linkType = "enhancer", verbose = FALSE),
               "carry no enhancer link")

  # And the same split on the gene table
  results <- testRegions(annotatedFit, contrast = exampleContrast(), verbose = FALSE)

  allGenes <- resultsTable(aggregateToGenes(results, verbose = FALSE))
  promoterGenes <- resultsTable(aggregateToGenes(results, linkType = "promoter", verbose = FALSE))

  expect_lte(sum(promoterGenes$n.regions), sum(allGenes$n.regions))

  # Links with no type at all cannot be filtered on
  byTable <- annotateRegions(fit, linkTable = data.frame(linkCoordinates(fit@counts, 1:4),
                                                         gene.id = letters[1:4]),
                             verbose = FALSE)
  tableResults <- testRegions(byTable, contrast = exampleContrast(), verbose = FALSE)
  expect_error(aggregateToGenes(tableResults, linkType = "promoter", verbose = FALSE),
               "carry no 'link.type' column")
})


test_that("the chromosome names of the links follow the ones of the regions", {

  counts <- exampleCounts()
  coordinates <- linkCoordinates(counts, 1:4)

  # The same links written the Ensembl way, as a contact table computed on another annotation would be
  ensemblCoordinates <- dplyr::mutate(coordinates, seqnames = sub("^chr", "", .data$seqnames))

  ucscLinks <- annotateRegions(counts, linkTable = data.frame(coordinates, gene.id = letters[1:4]), verbose = FALSE)
  ensemblLinks <- annotateRegions(counts, linkTable = data.frame(ensemblCoordinates, gene.id = letters[1:4]), verbose = FALSE)

  expect_identical(regionLinks(ensemblLinks)$region, regionLinks(ucscLinks)$region)
})


test_that("links outside the object are dropped while regions without a link are kept", {

  counts <- exampleCounts()
  coordinates <- linkCoordinates(counts, 1:2)

  # A contact map is routinely computed over more peaks than the ones being tested
  linkTable <- rbind(data.frame(coordinates, gene.id = c("one", "two")),
                     data.frame(seqnames = "chr1", start = 1L, end = 500L, gene.id = "elsewhere"))

  expect_message(annotateRegions(counts, linkTable = linkTable), "absent from the object")

  annotated <- annotateRegions(counts, linkTable = linkTable, verbose = FALSE)

  expect_equal(nrow(regionLinks(annotated)), 2)
  expect_equal(annotated@parameters$annotateRegions$n.regions.linked, 2)

  # Dropping the unlinked regions would shrink the universe of every enrichment run afterwards
  expect_equal(annotated@parameters$annotateRegions$n.regions.unlinked, nrow(counts) - 2)
  expect_equal(nrow(annotated), nrow(counts))
})


test_that("an ambiguous link table is refused rather than collapsed", {

  counts <- exampleCounts()
  coordinates <- linkCoordinates(counts, 1:3)
  linkTable <- data.frame(coordinates, gene.id = c("one", "two", "three"))

  # A duplicated pair would count its gene twice in every aggregation built on the links
  expect_error(annotateRegions(counts, linkTable = rbind(linkTable, linkTable[1, ]), verbose = FALSE),
               "duplicated region-gene pairs")

  # A missing weight is ambiguous between "not measured" and "zero"
  expect_error(annotateRegions(counts, linkTable = data.frame(linkTable, link.weight = c(0.5, NA, 0.5)), verbose = FALSE),
               "negative or missing")
  expect_error(annotateRegions(counts, linkTable = data.frame(linkTable, link.weight = c(0.5, -1, 0.5)), verbose = FALSE),
               "negative or missing")

  expect_error(annotateRegions(counts, linkTable = coordinates, verbose = FALSE), "requires a 'gene.id' column")
  expect_error(annotateRegions(counts, linkTable = data.frame(gene.id = "one"), verbose = FALSE), "must name the regions")
  expect_error(annotateRegions(counts, verbose = FALSE), "must be given")
  expect_error(annotateRegions(counts, txdb = "any", linkTable = linkTable, verbose = FALSE), "cannot be used together")

  expect_error(annotateRegions(counts,
                               linkTable = data.frame(seqnames = "chr1", start = 1L, end = 2L, gene.id = "elsewhere"),
                               verbose = FALSE),
               "Not one link")
})


test_that("links carry a weight or carry none, and the source of it is recorded", {

  counts <- exampleCounts()
  coordinates <- linkCoordinates(counts, 1:4)
  linkTable <- data.frame(coordinates, gene.id = c("one", "two", "three", "four"))

  # Supplying no weight is the normal case, not a degraded one
  unweighted <- annotateRegions(counts, linkTable = linkTable, verbose = FALSE)
  expect_false("link.weight" %in% colnames(regionLinks(unweighted)))
  expect_false(unweighted@parameters$annotateRegions$weighted)

  # Whatever the weight measures, it is stored as it is and named rather than rescaled
  weighted <- annotateRegions(counts, linkTable = data.frame(linkTable, link.weight = c(0.9, 0.5, 0.2, 0)),
                              weightSource = "Hi-C loop score", verbose = FALSE)

  expect_true(weighted@parameters$annotateRegions$weighted)
  expect_identical(weighted@parameters$annotateRegions$weight.source, "Hi-C loop score")
  expect_identical(regionLinks(weighted)$link.weight, c(0.9, 0.5, 0.2, 0))
})


test_that("the links can be attached to an object that never saw an annotation package", {

  counts <- exampleCounts()
  expect_null(regionLinks(counts))

  regionLinks(counts) <- data.frame(linkCoordinates(counts, 1:3), gene.id = c("one", "two", "three"))
  expect_equal(nrow(regionLinks(counts)), 3)

  # The replacement form runs the same checks the function does
  expect_error(regionLinks(counts) <- data.frame(gene.id = "one"), "must name the regions")
})


test_that("the links reach a fit and the results, and are read back from either", {

  fit <- exampleFit()
  coordinates <- linkCoordinates(fit@counts, 1:4)

  annotatedFit <- annotateRegions(fit, linkTable = data.frame(coordinates, gene.id = c("one", "two", "three", "four")),
                                  verbose = FALSE)

  expect_s4_class(annotatedFit, "RegionSetDE.fit")
  expect_equal(nrow(regionLinks(annotatedFit)), 4)

  results <- testRegions(annotatedFit, contrast = exampleContrast(), verbose = FALSE)
  expect_equal(nrow(regionLinks(results)), 4)
  expect_true("gene.id" %in% colnames(S4Vectors::mcols(resultRanges(results))))

  # The links are built from the gene sets of a fit as readily as from counts
  annotatedFit <- annotateRegions(annotatedFit,
                                  linkTable = data.frame(coordinates, gene.id = rep("oneGene", 4)),
                                  verbose = FALSE)
  expect_length(geneSetsToRegionSets(annotatedFit, geneSets = list(onePathway = "oneGene"),
                                     verbose = FALSE)$onePathway, 4)

  # A set score object holds no counts slot at all, so it is refused rather than reached into
  expect_error(annotateRegions(methods::new("RegionSetDE.setScores"),
                               linkTable = data.frame(coordinates, gene.id = letters[1:4]), verbose = FALSE),
               "must be a RegionSetDE.counts")
})


test_that("gene sets become region sets, and only the ones a test can run on come back", {

  counts <- exampleCounts()
  coordinates <- linkCoordinates(counts, 1:6)

  # Two regions per gene, which is the situation the whole exercise exists for
  annotated <- annotateRegions(counts,
                               linkTable = data.frame(coordinates, gene.id = rep(c("geneA", "geneB", "geneC"), each = 2)),
                               verbose = FALSE)

  geneSets <- list(pathwayOne = c("geneA", "geneB"),
                   pathwayTwo = c("geneC", "notLinked"),
                   pathwayThree = "notLinked")

  regionSets <- geneSetsToRegionSets(annotated, geneSets = geneSets, verbose = FALSE)

  expect_identical(sort(names(regionSets)), c("pathwayOne", "pathwayTwo"))
  expect_length(regionSets$pathwayOne, 4)
  expect_true(all(grepl("|", regionSets$pathwayOne, fixed = TRUE)))

  summaryTable <- attr(regionSets, "geneSetSummary")
  expect_identical(summaryTable$n.genes.linked[summaryTable$set.name == "pathwayTwo"], 1L)
  expect_identical(summaryTable$n.genes[summaryTable$set.name == "pathwayTwo"], 2L)

  # The two shapes gene sets arrive in
  termToGene <- data.frame(term = rep("pathwayOne", 2), gene = c("geneA", "geneB"))
  expect_identical(lengths(geneSetsToRegionSets(annotated, geneSets = termToGene, verbose = FALSE)),
                   lengths(regionSets["pathwayOne"]))

  expect_error(geneSetsToRegionSets(counts, geneSets = geneSets, verbose = FALSE), "carries no region to gene links")
  expect_error(geneSetsToRegionSets(annotated, geneSets = list(c("geneA")), verbose = FALSE), "must be named")
  expect_error(geneSetsToRegionSets(annotated, geneSets = list(empty = "notLinked"), verbose = FALSE),
               "Not one gene set reaches")
})


test_that("a weight filters the links it is given, and asking for one that is not there stops", {

  counts <- exampleCounts()
  coordinates <- linkCoordinates(counts, 1:4)
  linkTable <- data.frame(coordinates, gene.id = rep(c("geneA", "geneB"), each = 2))

  unweighted <- annotateRegions(counts, linkTable = linkTable, verbose = FALSE)
  weighted <- annotateRegions(counts, linkTable = data.frame(linkTable, link.weight = c(0.9, 0.1, 0.8, 0.2)), verbose = FALSE)

  geneSets <- list(pathwayOne = c("geneA", "geneB"))

  expect_length(geneSetsToRegionSets(weighted, geneSets = geneSets, verbose = FALSE)$pathwayOne, 4)
  expect_length(geneSetsToRegionSets(weighted, geneSets = geneSets, minWeight = 0.5, verbose = FALSE)$pathwayOne, 2)

  # A filter that quietly does nothing is worse than one that is not there
  expect_error(geneSetsToRegionSets(unweighted, geneSets = geneSets, minWeight = 0.5, verbose = FALSE),
               "carry no 'link.weight'")
})


test_that("one link reaches every set holding the same coordinates", {

  counts <- sharedSetCounts()
  annotated <- annotateRegions(counts, linkTable = data.frame(linkCoordinates(counts, 1), gene.id = "geneA"),
                               verbose = FALSE)

  # The coordinates sit in two sets, so the single link has to reach both rows
  expect_equal(nrow(regionLinks(annotated)), 1)
  expect_length(geneSetsToRegionSets(annotated, geneSets = list(pathwayOne = "geneA"), verbose = FALSE)$pathwayOne, 2)

  # Restricting to one set is how an analysis is kept to one class of element
  expect_error(geneSetsToRegionSets(annotated, geneSets = list(pathwayOne = "geneA"),
                                    regionSets = "sharedSet", verbose = FALSE),
               "Not one gene set reaches")
})


test_that("the sets a test runs on come from the column or from a supplied list", {

  regionStats <- data.frame(region.set = rep(c("firstSet", "secondSet"), each = 3),
                            region.id = paste0("region", 1:6),
                            stringsAsFactors = FALSE)

  # The sets of the object are a partition
  columnSets <- RegionSetDE:::.setIndexList(regionStats = regionStats)
  expect_identical(columnSets$firstSet, 1:3)
  expect_identical(names(RegionSetDE:::.setIndexList(regionStats, regionSets = "secondSet")), "secondSet")
  expect_error(RegionSetDE:::.setIndexList(regionStats, regionSets = "absentSet"), "absent from the object")

  # Supplied sets overlap each other, which is what gene sets do
  suppliedSets <- list(pathwayOne = c("firstSet|region1", "secondSet|region4"),
                       pathwayTwo = c("firstSet|region1", "firstSet|region2"))

  resolvedSets <- RegionSetDE:::.setIndexList(regionStats, regionSets = suppliedSets)
  expect_identical(resolvedSets$pathwayOne, c(1L, 4L))
  expect_identical(resolvedSets$pathwayTwo, c(1L, 2L))

  expect_error(RegionSetDE:::.setIndexList(regionStats, regionSets = list(pathwayOne = "nothing|here")),
               "Not one region of the supplied sets")
  expect_warning(RegionSetDE:::.setIndexList(regionStats, regionSets = c(suppliedSets, list(empty = "nothing|here"))),
                 "match no row of the object")
})


test_that("a set level test runs on supplied sets, against every other row", {

  fit <- exampleFit()
  coordinates <- linkCoordinates(fit@counts, 1:8)

  annotatedFit <- annotateRegions(fit,
                                  linkTable = data.frame(coordinates, gene.id = rep(c("geneA", "geneB"), each = 4)),
                                  verbose = FALSE)

  geneSets <- geneSetsToRegionSets(annotatedFit,
                                   geneSets = list(pathwayOne = "geneA", pathwayTwo = "geneB"),
                                   verbose = FALSE)

  setResults <- testRegionSets(annotatedFit, contrast = exampleContrast(), regionSets = geneSets, verbose = FALSE)
  setTable <- resultsTable(setResults)

  expect_identical(sort(setTable$region.set), c("pathwayOne", "pathwayTwo"))
  expect_true(all(setTable$n.regions > 1))
  expect_true(all(setTable$n.regions <= 4))

  # Each set is compared against everything else that was tested rather than against the other
  # sets, give or take the comparison rows overlapPolicy drops for overlapping the set
  expect_true(all(setTable$n.comparison <= nrow(annotatedFit@counts) - setTable$n.regions))
  expect_true(all(setTable$n.comparison > 0))
  expect_identical(setResults@universe@type, "allRows")

  # The two questions fry answers are both reported
  expect_true(all(c("fry.p", "fry.mixed.p", "n.up", "n.down") %in% colnames(setTable)))
  expect_true(all(setTable$fry.mixed.p >= 0 & setTable$fry.mixed.p <= 1))
})


test_that("the sets of the object are tested as before", {

  fit <- exampleFit()

  byName <- testRegionSets(fit, contrast = exampleContrast(), verbose = FALSE)
  byDefault <- testRegionSets(fit, contrast = exampleContrast(),
                              regionSets = unique(as.character(SummarizedExperiment::rowData(fit@counts)$region.set)),
                              verbose = FALSE)

  expect_equal(resultsTable(byName)$camera.p, resultsTable(byDefault)$camera.p)
  expect_equal(resultsTable(byName)$fry.p, resultsTable(byDefault)$fry.p)

  # Both questions are answered, and the mixed p-value is a p-value
  setTable <- resultsTable(byName)
  expect_true(all(c("fry.p", "fry.mixed.p", "fry.mixed.FDR") %in% colnames(setTable)))
  expect_true(all(setTable$fry.mixed.p >= 0 & setTable$fry.mixed.p <= 1))
  expect_true(all(setTable$n.up + setTable$n.down <= setTable$n.regions))
})
