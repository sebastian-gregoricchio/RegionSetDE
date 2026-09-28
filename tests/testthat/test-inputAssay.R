# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)


test_that("countReads counts the inputs of the sample sheet into an input assay", {

  counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)

  expect_identical(SummarizedExperiment::assayNames(counts), c("counts", "input"))
  inputMatrix <- SummarizedExperiment::assay(counts, "input")
  expect_identical(dim(inputMatrix), dim(counts))
  expect_identical(dimnames(inputMatrix), dimnames(SummarizedExperiment::assay(counts, "counts")))

  # One input serves every sample, so the columns are the same, and equal to the input counted as a sample
  expect_true(all(inputMatrix == inputMatrix[, 1]))
  inputAsSample <- countReads(peakRegions, bamFiles = unique(sampleSheet$input), sampleNames = "input", verbose = FALSE)
  expect_equal(as.numeric(inputMatrix[, 1]), as.numeric(SummarizedExperiment::assay(inputAsSample, "counts")[, 1]))

  sampleTable <- SummarizedExperiment::colData(counts)
  expect_true(all(c("input", "input.id", "input.library.size") %in% colnames(sampleTable)))
  expect_true(all(sampleTable$input.library.size == SummarizedExperiment::colData(inputAsSample)$library.size))
  expect_identical(counts@parameters$countReads$inputFiles, sampleSheet$input)
})


test_that("the inputs can be skipped, given directly, or missing for some samples", {

  withoutInput <- countReads(peakRegions, sampleSheet = sampleSheet, countInput = FALSE, verbose = FALSE)
  expect_identical(SummarizedExperiment::assayNames(withoutInput), "counts")

  bamFiles <- sampleSheet$bam[1:3]
  partialInput <- countReads(peakRegions, bamFiles = bamFiles, sampleNames = sampleSheet$sample[1:3],
                             inputFiles = c(sampleSheet$input[1], NA, sampleSheet$input[1]), verbose = FALSE)

  inputMatrix <- SummarizedExperiment::assay(partialInput, "input")
  expect_true(all(is.na(inputMatrix[, 2])))
  expect_false(anyNA(inputMatrix[, c(1, 3)]))
  expect_identical(SummarizedExperiment::colData(partialInput)$input.id, c("input", NA, "input"))

  noInputAtAll <- countReads(peakRegions, bamFiles = bamFiles, verbose = FALSE)
  expect_identical(SummarizedExperiment::assayNames(noInputAtAll), "counts")

  expect_error(countReads(peakRegions, bamFiles = bamFiles, inputFiles = c("no/such/input.bam", NA, NA), verbose = FALSE),
               "input files do not exist")
  expect_error(countReads(peakRegions, bamFiles = bamFiles, inputFiles = sampleSheet$input[1:2], verbose = FALSE),
               "one path per sample")
})


test_that("countTable returns the inputs, and they follow the samples through the analysis", {

  counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)

  inputTable <- countTable(counts, input = TRUE, format = "matrix")
  expect_equal(inputTable, as.matrix(SummarizedExperiment::assay(counts, "input")), ignore_attr = TRUE)

  longTable <- countTable(counts, input = TRUE, format = "long", extraColumns = FALSE)
  expect_true("input" %in% colnames(longTable))

  expect_error(countTable(counts, input = TRUE, normalized = TRUE), "not normalised")
  expect_error(countTable(countReads(peakRegions, sampleSheet = sampleSheet, countInput = FALSE, verbose = FALSE), input = TRUE),
               "carries no input counts")

  # Subsetting and fitting leave the assay where it belongs
  selectedCounts <- selectSamples(counts, samples = 1:4, verbose = FALSE)
  expect_identical(dim(SummarizedExperiment::assay(selectedCounts, "input")), c(nrow(counts), 4L))

  normalisedCounts <- normalizeCounts(counts, method = "TMM", verbose = FALSE)
  expect_s4_class(fitRegions(normalisedCounts, design = ~ condition, verbose = FALSE), "RegionSetDE.fit")
})
