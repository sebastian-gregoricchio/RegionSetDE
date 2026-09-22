## The brackets of plotRegion: the tests on the plotted values, the pairing by a
## metadata column, the statistics of the fit, and the checks that stop a
## comparison the data cannot hold.


# The example has no matched replicates, so the pairs are made up: each BN sample gets one SHR partner
pairedCounts <- function() {
  counts <- resultCounts(exampleResults())
  counts$pair <- c("p1", "p2", "p1", "p2")
  return(counts)
}


# Eight libraries in three groups, refitted: A and C hold the BN samples and a copy of them, B the SHR samples twice
threeGroupCounts <- function() {
  counts <- exampleCounts()

  # Repeating the columns gives eight libraries without binding two objects, which depends on the Bioconductor release
  mergedCounts <- counts[, rep(seq_len(ncol(counts)), times = 2)]
  colnames(mergedCounts) <- c(colnames(counts), paste0(colnames(counts), "-copy"))
  mergedCounts$group3 <- c("A", "A", "B", "B", "C", "C", "B", "B")

  return(normalizeCounts(mergedCounts, method = "TMM", verbose = FALSE))
}


# One contrast per form: the three-element one, and an expression that reduces to two levels
threeGroupResults <- function() {
  fit <- fitRegions(threeGroupCounts(), design = ~ group3, verbose = FALSE)

  return(testRegions(fit,
                     contrast = list(BvsA = c("group3", "B", "A"),
                                     CvsB = "group3C - group3B"),
                     verbose = FALSE))
}


topRegionOf <- function(results) {
  return(topRegions(results, n = 1, FDR = 1)$region.id)
}


# The values of one region as plotRegion draws them
regionValues <- function(counts, regionId, assay = "norm.counts") {
  rowIndex <- which(SummarizedExperiment::rowData(counts)$region.id == regionId)
  return(log2(as.numeric(SummarizedExperiment::assay(counts, assay)[rowIndex, ]) + 1))
}


# The table behind the labels of the brackets, NULL when the plot has none
bracketData <- function(regionPlot) {
  isLabel <- vapply(regionPlot$layers, function(layer) {inherits(layer$geom, "GeomRichText")}, logical(1))

  if (!any(isLabel)) {
    return(NULL)
  }

  return(regionPlot$layers[[which(isLabel)[1]]]$data)
}


# Building the grobs is where a label that the markdown parser cannot read fails
renderPlot <- function(regionPlot) {
  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())
  ggplot2::ggplotGrob(regionPlot)
  return(invisible(TRUE))
}




test_that("the t-test bracket carries the p-value of stats::t.test on the plotted values", {

  results <- exampleResults()
  topRegion <- topRegionOf(results)

  regionPlot <- plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "t.test")
  labelTable <- bracketData(regionPlot)

  values <- regionValues(resultCounts(results), topRegion)
  groups <- as.character(resultCounts(results)$condition)

  expect_s3_class(regionPlot, "ggplot")
  expect_identical(nrow(labelTable), 1L)
  expect_equal(labelTable$p.value, stats::t.test(values[groups == "BN"], values[groups == "SHR"])$p.value)
  expect_match(labelTable$label, "^p = ")
  expect_match(regionPlot$labels$caption, "^Welch t-test on log<sub>2</sub> norm.counts")
  expect_no_error(renderPlot(regionPlot))

  # Without brackets the plot keeps no caption
  expect_null(plotRegion(results, region = topRegion, groupBy = "condition")$labels$caption)
})


test_that("pairBy matches the samples on the column, whatever their order", {

  counts <- pairedCounts()
  topRegion <- topRegionOf(exampleResults())

  pairedPlot <- plotRegion(counts, region = topRegion, groupBy = "condition",
                           pairwiseTest = "t.test", pairBy = "pair")
  reversedPlot <- plotRegion(counts[, rev(seq_len(ncol(counts)))], region = topRegion, groupBy = "condition",
                             pairwiseTest = "t.test", pairBy = "pair")

  # BN holds the first two samples and SHR the last two, both in the order p1, p2
  values <- regionValues(counts, topRegion)
  pairedP <- stats::t.test(values[1:2], values[3:4], paired = TRUE)$p.value

  expect_equal(bracketData(pairedPlot)$p.value, pairedP)
  expect_equal(bracketData(reversedPlot)$p.value, pairedP)
  expect_match(pairedPlot$labels$caption, "^Paired t-test .*, paired by pair")
})


test_that("samples without a partner are left out, and fewer than two pairs stop the test", {

  results <- exampleResults()
  topRegion <- topRegionOf(results)

  # bio2 is the only replicate present in both conditions
  expect_error(plotRegion(results, region = topRegion, groupBy = "condition",
                          pairwiseTest = "t.test", pairBy = "biologicalReplicate"),
               "Fewer than two pairs")

  counts <- threeGroupCounts()
  counts$pair <- c("p1", "p2", "p1", "p2", "p1", "p2", "p3", "p4")

  expect_warning(regionPlot <- plotRegion(counts, region = topRegion, groupBy = "group3",
                                          pairwiseTest = "t.test", pairBy = "pair",
                                          comparisons = list(c("A", "B"))),
                 "No partner in 'pair'")
  expect_identical(bracketData(regionPlot)$n.first, 2L)

  # A value repeated inside one group gives a sample two possible partners
  counts$pair <- c("p1", "p1", "p1", "p2", "p1", "p2", "p3", "p4")
  expect_error(plotRegion(counts, region = topRegion, groupBy = "group3",
                          pairwiseTest = "t.test", pairBy = "pair",
                          comparisons = list(c("A", "B"))),
               "appears more than once")
})


test_that("the Wilcoxon bracket says when the groups are too small for the exact test", {

  results <- exampleResults()
  topRegion <- topRegionOf(results)

  expect_message(regionPlot <- plotRegion(results, region = topRegion, groupBy = "condition",
                                          pairwiseTest = "wilcox.test", pLabel = "stars"),
                 "cannot go below p = 0.33")

  expect_equal(bracketData(regionPlot)$p.value, 1 / 3)
  expect_identical(bracketData(regionPlot)$label, "ns")
  expect_match(regionPlot$labels$caption, "^Wilcoxon rank-sum test")
})


test_that("a test on raw counts warns about the sequencing depth", {

  topRegion <- topRegionOf(exampleResults())

  expect_warning(plotRegion(exampleCounts(), region = topRegion, groupBy = "condition", pairwiseTest = "t.test"),
                 "not normalised")
  expect_no_warning(plotRegion(exampleResults(), region = topRegion, groupBy = "condition", pairwiseTest = "t.test"))
})


test_that("three groups get one bracket per pair, adjusted across the plot when asked", {

  counts <- threeGroupCounts()
  regionTable <- resultsTable(threeGroupResults()@results$BvsA)
  topRegion <- regionTable$region.id[which.min(regionTable$p.value)]

  regionPlot <- plotRegion(counts, region = topRegion, groupBy = "group3",
                           pairwiseTest = "t.test", pAdjustMethod = "holm")
  labelTable <- bracketData(regionPlot)

  expect_identical(nrow(labelTable), 3L)
  expect_equal(labelTable$p.adjusted, stats::p.adjust(labelTable$p.value, method = "holm"))
  expect_true(all(grepl("^p<sub>adj</sub>", labelTable$label)))
  expect_match(regionPlot$labels$caption, "holm adjustment over 3 comparisons")

  # Brackets sharing a group cannot share a level, the widest one goes on top
  expect_identical(length(unique(labelTable$y.bracket)), 3L)
  expect_equal(labelTable$y.bracket[labelTable$group1 == "A" & labelTable$group2 == "C"], max(labelTable$y.bracket))

  onePair <- plotRegion(counts, region = topRegion, groupBy = "group3",
                        pairwiseTest = "t.test", comparisons = list(c("C", "A")))
  expect_identical(nrow(bracketData(onePair)), 1L)
  expect_no_error(renderPlot(regionPlot))
})


test_that("the model brackets carry the statistics of every contrast on the grouping column", {

  results <- threeGroupResults()
  regionTable <- resultsTable(results@results$BvsA)
  topRegion <- regionTable$region.id[which.min(regionTable$p.value)]

  modelPlot <- plotRegion(results, region = topRegion, groupBy = "group3", pairwiseTest = "model")
  labelTable <- bracketData(modelPlot)

  expect_identical(sort(labelTable$contrast), c("BvsA", "CvsB"))
  expect_equal(labelTable$FDR[labelTable$contrast == "BvsA"],
               regionTable$FDR[regionTable$region.id == topRegion])
  expect_equal(labelTable$log2FC[labelTable$contrast == "BvsA"],
               regionTable$log2FC[regionTable$region.id == topRegion])
  expect_true(all(grepl("^log<sub>2</sub>FC = .*<br>FDR ", labelTable$label)))
  expect_match(modelPlot$labels$caption, "from the edgeR fit")

  # Several contrasts on the brackets leave the subtitle empty
  expect_null(modelPlot$labels$subtitle)

  # A named contrast keeps its own bracket, and its statistics in the subtitle
  singlePlot <- plotRegion(results, region = topRegion, groupBy = "group3", pairwiseTest = "model", contrast = "CvsB")
  expect_identical(bracketData(singlePlot)$contrast, "CvsB")
  expect_match(singlePlot$labels$subtitle, "group3C - group3B")

  narrowPlot <- plotRegion(results, region = topRegion, groupBy = "group3", pairwiseTest = "model",
                           comparisons = list(c("A", "B")))
  expect_identical(bracketData(narrowPlot)$contrast, "BvsA")
  expect_no_error(renderPlot(modelPlot))
})


test_that("the model brackets need a results object and a contrast on the grouping column", {

  results <- exampleResults()
  topRegion <- topRegionOf(results)

  expect_error(plotRegion(resultCounts(results), region = topRegion, groupBy = "condition", pairwiseTest = "model"),
               "testRegions")
  expect_error(plotRegion(results, region = topRegion, groupBy = "sex", pairwiseTest = "model"),
               "None of the contrasts compares two levels of 'sex'")
  expect_warning(plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "model",
                            pairBy = "biologicalReplicate"),
                 "design of the fit")
  expect_error(plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "model",
                          comparisons = list(c("SHR", "SHR"))),
               "two different groups")
})


test_that("the brackets refuse what they cannot draw", {

  results <- exampleResults()
  topRegion <- topRegionOf(results)

  expect_error(plotRegion(results, region = topRegion, pairwiseTest = "t.test"), "needs a 'groupBy' column")
  expect_error(plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "anova"), "must be one of")
  expect_error(plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "t.test",
                          pLabel = "numbers"), "'pLabel'")
  expect_error(plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "t.test",
                          pAdjustMethod = "sidak"), "'pAdjustMethod'")
  expect_error(plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "t.test",
                          pDecimals = -1), "'pDecimals'")
  expect_error(plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "t.test",
                          comparisons = list(c("BN", "WKY"))), "absent from the object: WKY")
  expect_error(plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "t.test",
                          comparisons = list("BN")), "length two")

  # One female sample leaves the t-test without a variance on that side
  expect_error(plotRegion(results, region = topRegion, groupBy = "sex", pairwiseTest = "t.test"),
               "at least two samples per group")
  expect_error(plotRegion(results, region = topRegion, groupBy = "paired.end", pairwiseTest = "t.test"),
               "single group")

  expect_warning(plotRegion(results, region = topRegion, groupBy = "condition", pairBy = "biologicalReplicate"),
                 "have been ignored")

  tiledCounts <- toyCounts(tileWidth = 100L)
  expect_error(plotRegion(tiledCounts, region = "seq1:1-300", groupBy = "sample", pairwiseTest = "t.test"),
               "single row")
})


test_that("the labels are written as numbers or as symbols", {

  expect_identical(RegionSetDE:::.formatPvalue(c(0.5, 0.032, 1e-10, NA)),
                   c("0.50", "3.20&times;10<sup>-2</sup>", "1.00&times;10<sup>-10</sup>", "NA"))

  # A mantissa rounded up to 10 moves to the next power, and back to plain decimals from 0.1
  expect_identical(RegionSetDE:::.formatPvalue(c(0.009999, 0.09999)), c("1.00&times;10<sup>-2</sup>", "0.10"))
  expect_identical(RegionSetDE:::.formatPvalue(0, decimals = 1), "&lt; 2.2&times;10<sup>-308</sup>")

  expect_identical(RegionSetDE:::.bracketLabels(c(0.2, 0.03, 0.004, 5e-4, 5e-5, NA), pLabel = "stars"),
                   c("ns", "\\*", "\\*\\*", "\\*\\*\\*", "\\*\\*\\*\\*", "NA"))
  expect_identical(RegionSetDE:::.bracketLabels(0, prefix = "FDR"), "FDR &lt; 2.23&times;10<sup>-308</sup>")
  expect_identical(RegionSetDE:::.bracketLabels(0.5, prefix = "FDR", log2FC = 1.234),
                   "log<sub>2</sub>FC = 1.23<br>FDR = 0.50")
})


test_that("testRegions records the two levels a contrast compares", {

  fit <- exampleFit()

  byLevels <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
  byCoefficient <- testRegions(fit, contrast = "conditionSHR", verbose = FALSE)

  expect_identical(byLevels@contrast.groups, list(column = "condition", groups = c("SHR", "BN")))

  # The sample names separate any two samples too, the variable of the design is the one that counts
  expect_identical(byCoefficient@contrast.groups, byLevels@contrast.groups)

  setResults <- testRegionSets(fit, contrast = "conditionSHR", verbose = FALSE)
  expect_identical(setResults@contrast.groups$column, "condition")

  # With one sample per level the design variable looks like a column of sample labels, and still wins
  sampleTable <- data.frame(sample = c("s1", "s2"), condition = c("A", "B"))
  design <- stats::model.matrix(~ condition, data = sampleTable)

  expect_identical(RegionSetDE:::.contrastGroups(contrastVector = c(0, 1), design = design, colData = sampleTable),
                   list(column = "condition", groups = c("B", "A")))
})
