# A project folder with indexed BAM files copied from the Rsamtools example, a peak file and an unindexed BAM
sheetProject <- function() {
  projectPath <- tempfile(pattern = "sheetProject_")
  dir.create(file.path(projectPath, "reads"), recursive = TRUE)
  dir.create(file.path(projectPath, "peaks"))

  for (fileName in c("s1", "s2", "s3", "inputA", "inputB")) {
    file.copy(toyBamFile(), file.path(projectPath, "reads", paste0(fileName, ".bam")))
    file.copy(paste0(toyBamFile(), ".bai"), file.path(projectPath, "reads", paste0(fileName, ".bam.bai")))
  }

  file.copy(toyBamFile(), file.path(projectPath, "reads", "unindexed.bam"))
  writeLines("seq1\t1\t100", file.path(projectPath, "peaks", "s1.bed"))

  return(projectPath)
}


test_that("a table with relative paths is resolved against the folder of the sheet", {

  projectPath <- sheetProject()
  sheetTable <- data.frame(sample = c("s1", "s2", "s3"),
                           bam = c("reads/s1.bam", "reads/s2.bam", "reads/s3.bam"),
                           peaks = c("peaks/s1.bed", NA, NA),
                           input = c("reads/inputA.bam", "reads/inputA.bam", ""),
                           condition = c("A", "A", "B"))

  sheetFile <- file.path(projectPath, "samples.csv")
  utils::write.csv(sheetTable, sheetFile, row.names = FALSE)

  sampleSheet <- loadSampleSheet(sheetFile, verbose = FALSE)

  expect_s3_class(sampleSheet, "data.frame")
  expect_identical(colnames(sampleSheet), c("sample", "bam", "peaks", "input", "input.id", "condition"))
  expect_true(all(file.exists(sampleSheet$bam)))
  expect_identical(sampleSheet$bam[1], file.path(normalizePath(projectPath), "reads", "s1.bam"))

  # A blank cell is a missing input, and a shared input keeps one name
  expect_identical(sampleSheet$input.id, c("inputA", "inputA", NA))
  expect_true(is.na(sampleSheet$input[3]))
  expect_true(is.na(sampleSheet$peaks[2]))
})


test_that("a tab separated sheet is read as well", {

  projectPath <- sheetProject()
  sheetTable <- data.frame(sample = c("s1", "s2"), bam = c("reads/s1.bam", "reads/s2.bam"))

  sheetFile <- file.path(projectPath, "samples.txt")
  utils::write.table(sheetTable, sheetFile, sep = "\t", row.names = FALSE, quote = FALSE)

  expect_identical(loadSampleSheet(sheetFile, verbose = FALSE)$sample, c("s1", "s2"))
})


test_that("a DiffBind sample sheet is read as it is", {

  projectPath <- sheetProject()
  diffbindTable <- data.frame(SampleID = c("s1", "s2"),
                              Condition = c("resistant", "responsive"),
                              Replicate = c(1, 1),
                              bamReads = c("reads/s1.bam", "reads/s2.bam"),
                              ControlID = c("inputs_1", "inputs_2"),
                              bamControl = c("reads/inputA.bam", "reads/inputB.bam"),
                              Peaks = c("peaks/s1.bed", "peaks/s1.bed"),
                              PeakCaller = "bed")

  expect_message(sampleSheet <- loadSampleSheet(diffbindTable, basePath = projectPath), "SampleID as sample")

  expect_identical(colnames(sampleSheet)[seq_len(5)], c("sample", "bam", "peaks", "input", "input.id"))
  expect_true(all(c("Condition", "Replicate", "PeakCaller") %in% colnames(sampleSheet)))
  expect_identical(sampleSheet$input.id, c("inputs_1", "inputs_2"))
})


test_that("the columns can be pointed to by hand", {

  sheetTable <- data.frame(library = c("a", "b"), file = c("a.bam", "b.bam"), group = c("x", "y"))
  sampleSheet <- loadSampleSheet(sheetTable, columns = c(sample = "library", bam = "file"),
                                 basePath = "/data", checkFiles = FALSE, verbose = FALSE)

  expect_identical(sampleSheet$sample, c("a", "b"))
  expect_true("group" %in% colnames(sampleSheet))

  expect_error(loadSampleSheet(sheetTable, columns = c(sampleName = "library"), checkFiles = FALSE, verbose = FALSE), "not standard fields")
  expect_error(loadSampleSheet(sheetTable, columns = c(sample = "absent"), checkFiles = FALSE, verbose = FALSE), "absent from the sample sheet")
})


test_that("inputs named after their files stay distinct when two folders share a file name", {

  sheetTable <- data.frame(sample = c("a", "b"), bam = c("a.bam", "b.bam"),
                           input = c("run1/input.bam", "run2/input.bam"))

  sampleSheet <- loadSampleSheet(sheetTable, basePath = "/data", checkFiles = FALSE, verbose = FALSE)

  expect_identical(sampleSheet$input.id, c("input", "input_1"))
})


test_that("the names of the inputs must match their files one to one", {

  sheetTable <- data.frame(sample = c("a", "b"), bam = c("a.bam", "b.bam"),
                           input = c("input1.bam", "input2.bam"), input.id = c("same", "same"))

  expect_error(loadSampleSheet(sheetTable, basePath = "/data", checkFiles = FALSE, verbose = FALSE), "one name")
})


test_that("the sample names and the signal files are required", {

  expect_error(loadSampleSheet(data.frame(bam = "a.bam"), checkFiles = FALSE, verbose = FALSE), "No sample column")
  expect_error(loadSampleSheet(data.frame(sample = "a", peaks = "a.bed"), checkFiles = FALSE, verbose = FALSE), "'bam' or a 'bigwig'")
  expect_error(loadSampleSheet(data.frame(sample = c("a", "a"), bam = c("a.bam", "b.bam")), checkFiles = FALSE, verbose = FALSE), "unique")
  expect_error(loadSampleSheet(data.frame(sample = c("a", "b"), bam = c("a.bam", NA)), checkFiles = FALSE, verbose = FALSE), "no signal file")
  expect_error(loadSampleSheet(file.path(tempdir(), "absent.csv")), "does not exist")
})


test_that("missing and unindexed files are reported together", {

  projectPath <- sheetProject()
  sheetTable <- data.frame(sample = c("s1", "s2", "s3"),
                           bam = c("reads/s1.bam", "reads/absent.bam", "reads/unindexed.bam"))

  expect_error(loadSampleSheet(sheetTable, basePath = projectPath, verbose = FALSE), "missing bam files.*unindexed bam files")
  expect_s3_class(loadSampleSheet(sheetTable, basePath = projectPath, checkFiles = FALSE, verbose = FALSE), "data.frame")
})


test_that("a library listed twice is flagged, and remote files are left alone", {

  sheetTable <- data.frame(sample = c("a", "b"), bam = c("same.bam", "same.bam"))
  expect_warning(loadSampleSheet(sheetTable, basePath = "/data", checkFiles = FALSE, verbose = FALSE), "more than one sample")

  remoteTable <- data.frame(sample = "a", bigwig = "https://example.org/a.bw")
  expect_identical(loadSampleSheet(remoteTable, verbose = FALSE)$bigwig, "https://example.org/a.bw")
})


test_that("the sample sheet feeds countReads", {

  projectPath <- sheetProject()
  sheetTable <- data.frame(sample = c("s1", "s2"), bam = c("reads/s1.bam", "reads/s2.bam"), condition = c("A", "B"))
  sampleSheet <- loadSampleSheet(sheetTable, basePath = projectPath, verbose = FALSE)

  counts <- countReads(toyRegionSet(),
                       bamFiles = sampleSheet$bam,
                       sampleNames = sampleSheet$sample,
                       sampleMetadata = sampleSheet,
                       verbose = FALSE)

  expect_identical(colnames(counts), c("s1", "s2"))
  expect_identical(as.character(counts$condition), c("A", "B"))
})
