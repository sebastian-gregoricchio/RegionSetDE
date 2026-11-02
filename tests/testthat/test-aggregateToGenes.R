## The per-gene table built from a per-region result: the directions, the two aggregations and the weighting.


# A result whose regions have known directions, with a known number of regions per gene
annotatedResults <- function(weights = NULL, distances = NULL, genesPerRegion = NULL) {

  fit <- exampleFit()
  regionTable <- as.data.frame(SummarizedExperiment::rowRanges(fit@counts))[, c("seqnames", "start", "end")]

  # Four regions per gene unless the caller says otherwise
  geneIds <- if (is.null(genesPerRegion)) {
    rep(paste0("gene", seq_len(ceiling(nrow(regionTable) / 4))), each = 4)[seq_len(nrow(regionTable))]
  } else {
    genesPerRegion
  }

  linkTable <- data.frame(regionTable, gene.id = geneIds, stringsAsFactors = FALSE)
  if (!is.null(weights)) {linkTable$link.weight <- weights}
  if (!is.null(distances)) {linkTable$distance.to.tss <- distances}

  annotatedFit <- annotateRegions(fit, linkTable = linkTable, verbose = FALSE)

  return(testRegions(annotatedFit, contrast = exampleContrast(), verbose = FALSE))
}



test_that("the regions of a gene are counted and their directions reported rather than resolved", {

  results <- annotatedResults()
  geneResults <- aggregateToGenes(results, verbose = FALSE)

  expect_s4_class(geneResults, "RegionSetDE.geneResults")
  expect_true(methods::validObject(geneResults))

  geneTable <- resultsTable(geneResults)
  expect_equal(anyDuplicated(geneTable$gene.id), 0)
  expect_true(all(geneTable$n.regions > 0))
  expect_true(all(geneTable$direction %in% c("up", "down", "mixed", "unchanged")))

  # The four classes are mutually exclusive and follow the two counts
  expect_identical(geneTable$direction == "mixed", geneTable$n.up > 0 & geneTable$n.down > 0)
  expect_identical(geneTable$direction == "unchanged", geneTable$n.up == 0 & geneTable$n.down == 0)
  expect_true(all(geneTable$n.up + geneTable$n.down <= geneTable$n.regions))
})


test_that("the directions are read from diff.status, so updateThresholds travels into the gene table", {

  results <- annotatedResults()
  strictTable <- resultsTable(aggregateToGenes(results, verbose = FALSE))

  # Loosening the labels of the regions has to loosen the labels of the genes
  loosened <- updateThresholds(results, FDR = 1, verbose = FALSE)
  loosenedTable <- resultsTable(aggregateToGenes(loosened, verbose = FALSE))

  expect_equal(nrow(strictTable), nrow(loosenedTable))
  expect_gte(sum(loosenedTable$n.up + loosenedTable$n.down),
             sum(strictTable$n.up + strictTable$n.down))
  expect_lte(sum(loosenedTable$direction == "unchanged"), sum(strictTable$direction == "unchanged"))
})


test_that("the Simes combination matches its definition and carries the right fold change", {

  results <- annotatedResults()
  geneResults <- aggregateToGenes(results, verbose = FALSE)

  geneTable <- resultsTable(geneResults)
  regionTable <- resultsTable(results)
  links <- regionLinks(geneResults)

  regionKeys <- RegionSetDE:::.regionKeys(regionTable)

  for (geneId in utils::head(geneTable$gene.id, 5)) {
    regionIndex <- match(links$region[links$gene.id == geneId], regionKeys)
    regionP <- regionTable$p.value[regionIndex]

    simesValue <- min(1, min(length(regionP) / seq_along(regionP) * sort(regionP)))
    expect_equal(geneTable$p.value[geneTable$gene.id == geneId], simesValue, tolerance = 1e-10)

    # The fold change beside the p-value belongs to the region that produced it, not to an average
    expect_true(geneTable$log2FC[geneTable$gene.id == geneId] %in% regionTable$log2FC[regionIndex])
  }

  expect_true(all(geneTable$FDR >= geneTable$p.value - 1e-12))
  expect_false(is.unsorted(geneTable$p.value))
})


test_that("the potential score cancels on a mixed gene and its magnitude does not", {

  results <- annotatedResults()
  geneResults <- aggregateToGenes(results, method = "potential", verbose = FALSE)
  geneTable <- resultsTable(geneResults)

  expect_true(all(c("score", "score.magnitude") %in% colnames(geneTable)))
  expect_false("p.value" %in% colnames(geneTable))

  # The unsigned sum is never smaller than the signed one, and they part ways exactly on the mixed genes
  expect_true(all(geneTable$score.magnitude >= abs(geneTable$score) - 1e-12))

  mixedGenes <- geneTable$direction == "mixed"
  if (any(mixedGenes)) {
    expect_true(all(geneTable$score.magnitude[mixedGenes] > abs(geneTable$score[mixedGenes])))
  }
})


test_that("the weighting falls back from the link weight to the distance to nothing", {

  regionNumber <- nrow(exampleFit()@counts)

  # Nothing to weight with
  plainResults <- aggregateToGenes(annotatedResults(), method = "potential", verbose = FALSE)
  expect_identical(plainResults@weighting, "none")

  # A distance, turned into exp(-d / decay)
  distances <- rep(c(0L, 10000L, 20000L, 30000L), length.out = regionNumber)
  distanceResults <- aggregateToGenes(annotatedResults(distances = distances),
                                      method = "potential", decay = 5000, verbose = FALSE)
  expect_identical(distanceResults@weighting, "distance")

  # A weight of its own, which beats the distance whatever that says
  bothResults <- aggregateToGenes(annotatedResults(weights = rep(0.5, regionNumber), distances = distances),
                                  method = "potential", verbose = FALSE)
  expect_identical(bothResults@weighting, "link.weight")

  # A region weighted zero cannot reach the score of its gene, so every score is the fold change
  # of the one region of that gene that kept a weight
  zeroWeights <- rep(c(1, 0, 0, 0), length.out = regionNumber)
  zeroSource <- annotatedResults(weights = zeroWeights)
  zeroResults <- aggregateToGenes(zeroSource, method = "potential", verbose = FALSE)

  regionTable <- resultsTable(zeroSource)
  regionKeys <- RegionSetDE:::.regionKeys(regionTable)
  links <- regionLinks(zeroResults)

  for (geneId in utils::head(resultsTable(zeroResults)$gene.id, 5)) {
    weightedRegion <- links$region[links$gene.id == geneId & links$link.weight > 0]
    expect_length(weightedRegion, 1)

    expect_equal(resultsTable(zeroResults)$score[resultsTable(zeroResults)$gene.id == geneId],
                 regionTable$log2FC[match(weightedRegion, regionKeys)],
                 tolerance = 1e-10)
  }
})


test_that("both aggregations can run together and the object records what travelled", {

  results <- annotatedResults()
  geneResults <- aggregateToGenes(results, method = c("simes", "potential"), verbose = FALSE)

  expect_identical(geneResults@method, c("simes", "potential"))
  expect_true(all(c("p.value", "FDR", "score", "score.magnitude") %in% colnames(resultsTable(geneResults))))

  expect_identical(geneResults@contrast, results@contrast)
  expect_identical(geneResults@engine, results@engine)
  expect_identical(geneResults@genome.assembly, results@genome.assembly)
  expect_equal(geneResults@thresholds$FDR, results@thresholds$FDR)
  expect_equal(nrow(regionLinks(geneResults)), nrow(regionLinks(results)))
  expect_equal(geneResults@parameters$aggregateToGenes$n.genes, nrow(resultsTable(geneResults)))

  expect_output(show(geneResults), "RegionSetDE.geneResults")
  expect_output(show(geneResults), "aggregation")
})


test_that("the aggregation refuses what it cannot aggregate", {

  results <- annotatedResults()

  expect_error(aggregateToGenes(results, method = "best", verbose = FALSE), "must contain only")
  expect_error(aggregateToGenes(results, decay = 0, verbose = FALSE), "must be a positive number")
  expect_error(aggregateToGenes(results, regionSets = "absentSet", verbose = FALSE), "absent from the object")

  # Without links there is nothing to aggregate through
  bareResults <- testRegions(exampleFit(), contrast = exampleContrast(), verbose = FALSE)
  expect_error(aggregateToGenes(bareResults, verbose = FALSE), "carries no region to gene links")

  # Rows that are tiles would weight a gene by how finely its regions were cut
  tiledResults <- results
  tiledResults@counting.level <- "tile"
  tiledResults@combination <- list(applied = FALSE, method = NA_character_)
  expect_error(aggregateToGenes(tiledResults, verbose = FALSE), "rows of this result are tiles")

  expect_error(aggregateToGenes(exampleFit(), verbose = FALSE), "must be a RegionSetDE.results object")
})


test_that("restricting to one class of element drops the genes that only had the others", {

  results <- annotatedResults()
  setNames <- unique(as.character(resultsTable(results)$region.set))

  allGenes <- resultsTable(aggregateToGenes(results, verbose = FALSE))
  someGenes <- resultsTable(aggregateToGenes(results, regionSets = setNames[1], verbose = FALSE))

  expect_lte(nrow(someGenes), nrow(allGenes))
  expect_true(all(someGenes$gene.id %in% allGenes$gene.id))

  # A gene kept in both can only have lost regions
  sharedGenes <- intersect(someGenes$gene.id, allGenes$gene.id)
  expect_true(all(someGenes$n.regions[match(sharedGenes, someGenes$gene.id)] <=
                    allGenes$n.regions[match(sharedGenes, allGenes$gene.id)]))
})
