## make-data.R
##
## Builds the example objects shipped in inst/extdata for RegionSetDE.
##
## Source data: liver ChIP-seq from the EURATRANS project, comparing the Brown
## Norway (BN) and spontaneously hypertensive (SHR) rat strains. The BAM files
## are distributed inside the chromstaRData package, aligned to rn4 and
## restricted to chromosome 12.
##
## Contrast: H3K4me3, SHR against BN, two biological replicates per strain, all
## tech1. This is the only mark in the package with two biological replicates on
## both sides. BN bio1 is female and the other three libraries are male, so sex
## is not confounded with strain but does add within-group variance. With two
## replicates per side there is no way to model it.
##
## The input libraries take no part in the counting. They are used only to build
## the exclusion list, which keeps that list independent of the H3K4me3
## libraries that get tested afterwards.
##
## Nothing here contacts ExperimentHub or AnnotationHub. The BAM files ship with
## the package and the annotation tables come from the UCSC download server.
##
## Run once from the package root. The script installs with the package under
## inst/scripts and should not be added to .Rbuildignore.
##
## Author: Sebastian Gregoricchio


# ---- Parameters -------------------------------------------------------------

randomSeed <- 20260101L
targetChromosomes <- "chr12"
genomeBuild <- "rn4"

regionWidth <- 1000L
promoterFlank <- 500L
geneBodyTileWidth <- 5000L
intergenicTileWidth <- 10000L
minimumDistanceFromTss <- 2000L
minimumDistanceFromGene <- 10000L
maxRegionsPerSet <- 1500L
minRegionsPerSet <- 100L

exclusionBinWidth <- 1000L
exclusionQuantile <- 0.9995
backgroundBinSize <- 10000L
maxCrossCorrelationDistance <- 500L
minimumMappingQuality <- 10L

normalizationMethod <- "background"
filterMethod <- "background"
filterFoldChange <- 2
fitEngine <- "edgeR"

outputDir <- file.path("inst", "extdata")
scriptDir <- file.path("inst", "scripts")
ucscBaseUrl <- paste0("https://hgdownload.soe.ucsc.edu/goldenPath/",
                      genomeBuild, "/database/")

set.seed(randomSeed)


# ---- Dependency check -------------------------------------------------------

# All of this runs offline, so none of these belong in DESCRIPTION.
requiredPackages <- c(
  "chromstaRData",
  "TxDb.Rnorvegicus.UCSC.rn4.ensGene",
  "GenomicFeatures",
  "GenomicRanges",
  "GenomicAlignments",
  "SummarizedExperiment",
  "GenomeInfoDb",
  "IRanges",
  "Rsamtools",
  "csaw",
  "dplyr",
  "RegionSetDE"
)

missingPackages <- requiredPackages[
  !vapply(requiredPackages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missingPackages) > 0) {
  stop(
    "Missing packages required to rebuild the example data: ",
    paste(missingPackages, collapse = ", ")
  )
}

if (!dir.exists("R") || !file.exists("DESCRIPTION")) {
  stop("Run this script from the package root.")
}

dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
dir.create(scriptDir, recursive = TRUE, showWarnings = FALSE)


# ---- UCSC helpers -----------------------------------------------------------

# Older UCSC assemblies store some tracks per chromosome, so a genome-wide table
# name may not exist. Returns NA_character_ when a file is absent, which the
# caller uses to pick the other layout.
downloadUcscFile <- function(tableName) {

  destinationFile <- tempfile(fileext = ".txt.gz")
  downloadUrl <- paste0(ucscBaseUrl, tableName, ".txt.gz")

  downloadStatus <- tryCatch(
    utils::download.file(downloadUrl, destinationFile,
                         mode = "wb", quiet = TRUE, method = "libcurl"),
    error = function(condition) 1L,
    warning = function(condition) 1L
  )

  if (downloadStatus != 0 || !file.exists(destinationFile) ||
      file.size(destinationFile) == 0) {
    return(NA_character_)
  }

  destinationFile
}

# UCSC database dumps carry no header and an optional leading bin column, so the
# chromosome column is located by content rather than by position. Coordinates
# are 0-based half-open and get shifted here.
readUcscTable <- function(tableName, chromosomes, chromosomeInfo) {

  # Genome-wide layout first, per-chromosome layout second.
  downloadedFiles <- downloadUcscFile(tableName)

  if (is.na(downloadedFiles)) {
    downloadedFiles <- vapply(
      paste0(chromosomes, "_", tableName), downloadUcscFile, character(1)
    )
  }

  if (all(is.na(downloadedFiles))) {
    stop("Could not download the UCSC ", tableName, " table for ", genomeBuild,
         " under either the genome-wide or the per-chromosome layout. ",
         "List the database directory to see which tables exist.")
  }

  rawTable <- do.call(
    rbind,
    lapply(
      downloadedFiles[!is.na(downloadedFiles)],
      function(filePath) {
        utils::read.delim(gzfile(filePath), header = FALSE,
                          stringsAsFactors = FALSE)
      }
    )
  )

  chromosomeColumn <- which(
    vapply(
      rawTable,
      function(column) {
        is.character(column) && all(grepl("^chr", utils::head(column, 100)))
      },
      logical(1)
    )
  )[1]

  if (is.na(chromosomeColumn)) {
    stop("No chromosome column found in the UCSC ", tableName, " table.")
  }

  regionTable <-
    data.frame(
      seqnames = rawTable[[chromosomeColumn]],
      start = as.integer(rawTable[[chromosomeColumn + 1L]]) + 1L,
      end = as.integer(rawTable[[chromosomeColumn + 2L]]),
      stringsAsFactors = FALSE
    ) |>
    dplyr::filter(seqnames %in% chromosomes)

  if (nrow(regionTable) == 0) {
    stop("The UCSC ", tableName, " table holds no rows on ",
         paste(chromosomes, collapse = ", "), ".")
  }

  GenomicRanges::makeGRangesFromDataFrame(regionTable, seqinfo = chromosomeInfo)
}


# ---- Sample sheet -----------------------------------------------------------

euratransDir <- system.file("extdata", "euratrans", package = "chromstaRData",
                            mustWork = TRUE)

bamFileNames <- dir(euratransDir, pattern = "\\.bam$")

# File names encode everything: tissue-mark-strain-sex-bioRep-techRep.
nameFields <- do.call(
  rbind,
  strsplit(sub("\\.bam$", "", bamFileNames), "-", fixed = TRUE)
)

if (ncol(nameFields) != 6) {
  stop("Expected six dash-separated fields in the chromstaRData file names, found ",
       ncol(nameFields), ".")
}

colnames(nameFields) <- c("tissue", "mark", "strain", "sex",
                          "biologicalReplicate", "technicalReplicate")

sampleSheet <-
  as.data.frame(nameFields, stringsAsFactors = FALSE) |>
  dplyr::mutate(
    sampleName = sub("\\.bam$", "", bamFileNames),
    bamPath = file.path(euratransDir, bamFileNames),
    isInput = mark == "input"
  ) |>
  dplyr::select(sampleName, tissue, mark, strain, sex, biologicalReplicate,
                technicalReplicate, isInput, bamPath) |>
  dplyr::arrange(mark, strain, biologicalReplicate, technicalReplicate)

# H3K4me3 tech1 only, which leaves two biological replicates per strain.
chipSamples <-
  sampleSheet |>
  dplyr::filter(mark == "H3K4me3", technicalReplicate == "tech1") |>
  dplyr::mutate(condition = strain)

inputSamples <- sampleSheet |> dplyr::filter(isInput)

replicateCount <- chipSamples |> dplyr::count(condition, name = "nReplicates")
print(replicateCount)

if (nrow(replicateCount) != 2 || any(replicateCount$nReplicates < 2)) {
  stop("The H3K4me3 contrast does not have two replicates in both strains.")
}

if (nrow(inputSamples) == 0) {
  stop("No input libraries found; the exclusion list needs them.")
}


# ---- Sequence information ---------------------------------------------------

# Chromosome lengths come from the BAM header, which avoids pulling in a
# BSgenome package worth close to a gigabyte for one field.
bamTargets <- Rsamtools::scanBamHeader(chipSamples$bamPath[1])[[1]]$targets

if (!all(targetChromosomes %in% names(bamTargets))) {
  stop("Chromosome ", paste(targetChromosomes, collapse = ", "),
       " is absent from the BAM header.")
}

chromosomeInfo <- GenomeInfoDb::Seqinfo(
  seqnames = targetChromosomes,
  seqlengths = as.integer(bamTargets[targetChromosomes]),
  isCircular = rep(FALSE, length(targetChromosomes)),
  genome = genomeBuild
)


# ---- Fragment extension length ----------------------------------------------

# Reads are single-end, so the extension length is estimated once here and
# recorded rather than recomputed when the vignette builds. Duplicates are not
# flagged in these files, so deduplication stays off throughout.
crossCorrelation <- csaw::correlateReads(
  chipSamples$bamPath[1],
  max.dist = maxCrossCorrelationDistance,
  param = csaw::readParam(
    dedup = FALSE,
    minq = minimumMappingQuality,
    restrict = targetChromosomes
  )
)

extensionLength <- csaw::maximizeCcf(crossCorrelation)
message("Estimated fragment extension length: ", extensionLength, " bp")


# ---- Exclusion regions ------------------------------------------------------

# rn4 has no curated blacklist. The exclusion list combines assembly gaps with
# bins carrying implausible input coverage, a stripped-down version of how the
# ENCODE lists were built.
gapRegions <- readUcscTable("gap", targetChromosomes, chromosomeInfo)

coverageBins <- GenomicRanges::tileGenome(
  chromosomeInfo,
  tilewidth = exclusionBinWidth,
  cut.last.tile.in.chrom = TRUE
)

inputBinCounts <- GenomicAlignments::summarizeOverlaps(
  features = coverageBins,
  reads = Rsamtools::BamFileList(inputSamples$bamPath),
  singleEnd = TRUE,
  ignore.strand = TRUE
)

pooledInputCounts <- rowSums(SummarizedExperiment::assay(inputBinCounts))

coverageThreshold <- stats::quantile(
  pooledInputCounts[pooledInputCounts > 0], exclusionQuantile
)

highCoverageRegions <-
  coverageBins[pooledInputCounts > coverageThreshold] |>
  GenomicRanges::reduce(min.gapwidth = exclusionBinWidth)

exclusionRegions <- GenomicRanges::reduce(c(gapRegions, highCoverageRegions))
GenomeInfoDb::genome(exclusionRegions) <- genomeBuild

if (length(exclusionRegions) == 0) {
  stop("The exclusion list came out empty. Check the gap download and the ",
       "coverage threshold.")
}

message("Exclusion regions: ", length(exclusionRegions), " intervals covering ",
        round(sum(GenomicRanges::width(exclusionRegions)) / 1e3), " kb")


# ---- Annotation -------------------------------------------------------------

txdb <- TxDb.Rnorvegicus.UCSC.rn4.ensGene::TxDb.Rnorvegicus.UCSC.rn4.ensGene

transcriptRanges <-
  GenomicFeatures::transcripts(txdb) |>
  GenomeInfoDb::keepSeqlevels(targetChromosomes, pruning.mode = "coarse")

geneRanges <-
  GenomicFeatures::genes(txdb) |>
  GenomeInfoDb::keepSeqlevels(targetChromosomes, pruning.mode = "coarse")

transcriptStarts <- GenomicRanges::resize(transcriptRanges, width = 1L, fix = "start")
GenomicRanges::strand(transcriptStarts) <- "*"
transcriptStarts <- unique(GenomicRanges::granges(transcriptStarts))

GenomicRanges::strand(geneRanges) <- "*"
geneRanges <- GenomicRanges::reduce(GenomicRanges::granges(geneRanges))

cpgIslands <- readUcscTable("cpgIslandExt", targetChromosomes, chromosomeInfo)


# ---- Candidate anchors per region set ---------------------------------------

# Each set gets one anchor point per candidate region. Windows are cut to a
# fixed width later so that counts stay comparable across sets.
promoterWindows <- GenomicRanges::resize(
  transcriptStarts, width = 2L * promoterFlank, fix = "center"
)

# Positions inside genes but clear of any transcription start site.
geneBodyAnchors <-
  unlist(GenomicRanges::tile(geneRanges, width = geneBodyTileWidth)) |>
  GenomicRanges::resize(width = 1L, fix = "center")

geneBodyAnchors <- geneBodyAnchors[
  !IRanges::overlapsAny(
    GenomicRanges::resize(geneBodyAnchors,
                          width = 2L * minimumDistanceFromTss, fix = "center"),
    transcriptStarts,
    ignore.strand = TRUE
  )
]

# Positions well away from anything annotated. This set is the low-signal control.
intergenicAnchors <-
  GenomicRanges::tileGenome(chromosomeInfo, tilewidth = intergenicTileWidth,
                            cut.last.tile.in.chrom = TRUE) |>
  GenomicRanges::resize(width = 1L, fix = "center")

intergenicAnchors <- intergenicAnchors[
  !IRanges::overlapsAny(
    GenomicRanges::resize(intergenicAnchors,
                          width = 2L * minimumDistanceFromGene, fix = "center"),
    geneRanges,
    ignore.strand = TRUE
  )
]

# Promoters split by CpG island overlap. H3K4me3 sits mostly on the CpG side, so
# the two halves give a signal contrast inside the same feature class.
promoterIsCpG <- IRanges::overlapsAny(promoterWindows, cpgIslands,
                                      ignore.strand = TRUE)

anchorSets <- list(
  promoterCpG = GenomicRanges::resize(promoterWindows[promoterIsCpG],
                                      width = 1L, fix = "center"),
  promoterNonCpG = GenomicRanges::resize(promoterWindows[!promoterIsCpG],
                                         width = 1L, fix = "center"),
  geneBody = geneBodyAnchors,
  intergenic = intergenicAnchors
)


# ---- Region universe --------------------------------------------------------

# Sets are resolved in order, so a window claimed by an earlier set is never
# reused by a later one. Promoters win any collision. The exclusion list is not
# applied here: applyBlacklist does that on the object, where it gets logged.
keptRegions <- GenomicRanges::GRanges(seqinfo = chromosomeInfo)

for (setName in names(anchorSets)) {

  candidateRegions <-
    anchorSets[[setName]] |>
    GenomicRanges::resize(width = regionWidth, fix = "center") |>
    GenomicRanges::trim() |>
    sort()

  # Windows clipped at a chromosome edge are dropped rather than kept ragged.
  candidateRegions <- candidateRegions[
    GenomicRanges::width(candidateRegions) == regionWidth
  ]

  # disjointBins puts overlapping ranges in different bins, so the first bin is
  # a non-overlapping subset.
  candidateRegions <- candidateRegions[
    GenomicRanges::disjointBins(candidateRegions) == 1L
  ]

  candidateRegions <- candidateRegions[
    !IRanges::overlapsAny(candidateRegions, keptRegions, ignore.strand = TRUE)
  ]

  if (length(candidateRegions) > maxRegionsPerSet) {
    candidateRegions <- sort(sample(candidateRegions, maxRegionsPerSet))
  }

  candidateRegions <- GenomicRanges::granges(candidateRegions)
  candidateRegions$setName <- setName

  keptRegions <- c(keptRegions, candidateRegions)
}

keptRegions <- sort(keptRegions)

# Strand stays unset throughout: two windows differing only by strand would draw
# their counts from the same reads.
GenomicRanges::strand(keptRegions) <- "*"


# ---- Region sets ------------------------------------------------------------

regionTable <-
  data.frame(
    seqnames = as.character(GenomicRanges::seqnames(keptRegions)),
    start = GenomicRanges::start(keptRegions),
    end = GenomicRanges::end(keptRegions),
    setName = keptRegions$setName,
    stringsAsFactors = FALSE
  ) |>
  dplyr::mutate(regionId = sprintf("region_%05d", dplyr::row_number())) |>
  dplyr::select(regionId, seqnames, start, end, setName)

regionsGRanges <- GenomicRanges::makeGRangesFromDataFrame(
  regionTable,
  keep.extra.columns = TRUE,
  seqinfo = chromosomeInfo
)
names(regionsGRanges) <- regionTable$regionId

regions <- RegionSetDE::splitLoadRegions(
  regions = regionsGRanges,
  splitBy = "setName",
  minRegionsPerSet = minRegionsPerSet,
  seqlevelsStyle = "UCSC",
  genomeAssembly = genomeBuild
)

# Filtering runs on the object so that it lands in filtering.log and the vignette
# can show what was removed.
regions <- RegionSetDE::applyBlacklist(
  regionSet = regions,
  blacklist = exclusionRegions,
  ignoreStrand = TRUE,
  emptySets = "stop"
)

setSizes <- data.frame(
  setName = names(regions@regions),
  nRegions = lengths(regions@regions),
  stringsAsFactors = FALSE
)

print(setSizes)

if (nrow(setSizes) != length(anchorSets) || any(setSizes$nRegions < minRegionsPerSet)) {
  stop("At least one region set holds fewer than ", minRegionsPerSet,
       " regions after blacklisting. Widen the promoter flank or loosen the ",
       "distance thresholds.")
}


# ---- Sample metadata --------------------------------------------------------

# BN is the reference level, so the contrast reads as SHR against BN.
sampleMetadata <-
  chipSamples |>
  dplyr::select(sample = sampleName, condition, sex, biologicalReplicate) |>
  dplyr::mutate(condition = factor(condition, levels = c("BN", "SHR")))


# ---- Counting ---------------------------------------------------------------

# Duplicates are not flagged in these BAM files, so removeDuplicates would be a
# no-op and is switched off to keep the record honest.
counts <- RegionSetDE::countReads(
  regionSet = regions,
  bamFiles = chipSamples$bamPath,
  sampleNames = chipSamples$sampleName,
  sampleMetadata = sampleMetadata,
  pairedEnd = FALSE,
  fragmentLength = extensionLength,
  minMapq = minimumMappingQuality,
  removeDuplicates = FALSE,
  restrictChromosomes = targetChromosomes,
  nThreads = 1,
  verbose = TRUE
)

# The bins go into metadata(counts)$background, so this returns the same object.
counts <- RegionSetDE::countBackground(
  counts = counts,
  binSize = backgroundBinSize,
  excludeRegions = TRUE,
  nThreads = 1,
  verbose = TRUE
)


# ---- Fit --------------------------------------------------------------------

# The shipped counts object stays raw so that the vignette can run the
# normalisation and the filtering itself. The fit is built on a separate copy.
fitCountsObject <- RegionSetDE::normalizeCounts(
  counts = counts,
  method = normalizationMethod,
  verbose = TRUE
)

fitCountsObject <- RegionSetDE::filterRegions(
  counts = fitCountsObject,
  method = filterMethod,
  foldChange = filterFoldChange,
  verbose = TRUE
)

# The formula is built against baseenv(): a formula created here would otherwise
# carry the whole script environment into the serialised object.
designFormula <- stats::as.formula("~ condition", env = baseenv())

fit <- RegionSetDE::fitRegions(
  counts = fitCountsObject,
  design = designFormula,
  engine = fitEngine,
  verbose = TRUE
)

# One test run so that a broken fit is caught here and not in the man pages.
testRun <- RegionSetDE::testRegions(
  fit = fit,
  contrast = c("condition", "SHR", "BN"),
  verbose = FALSE
)

print(utils::head(RegionSetDE::topRegions(testRun, n = 5)))


# ---- Save -------------------------------------------------------------------

# bamPath is dropped: those paths point inside the installed chromstaRData tree
# and mean nothing on another machine.
objectsToSave <- list(
  "euratrans_sampleSheet.rds" = sampleSheet |> dplyr::select(-bamPath),
  "euratrans_regions.rds" = regionTable,
  "euratrans_exclusionRegions.rds" = exclusionRegions,
  "euratrans_counts.rds" = counts,
  "euratrans_fit.rds" = fit
)

invisible(
  mapply(
    function(object, fileName) {
      saveRDS(object, file.path(outputDir, fileName), compress = "xz")
    },
    objectsToSave,
    names(objectsToSave)
  )
)


# ---- Provenance -------------------------------------------------------------

buildMetadata <- list(
  buildDate = Sys.Date(),
  randomSeed = randomSeed,
  genomeBuild = genomeBuild,
  targetChromosomes = targetChromosomes,
  mark = "H3K4me3",
  contrast = "SHR against BN, rat liver",
  regionWidth = regionWidth,
  promoterFlank = promoterFlank,
  minimumDistanceFromTss = minimumDistanceFromTss,
  minimumDistanceFromGene = minimumDistanceFromGene,
  fragmentLength = extensionLength,
  minMapq = minimumMappingQuality,
  backgroundBinSize = backgroundBinSize,
  exclusionSource = "UCSC rn4 gap track and high-coverage input bins",
  exclusionQuantile = exclusionQuantile,
  annotationSource = "TxDb.Rnorvegicus.UCSC.rn4.ensGene and UCSC cpgIslandExt",
  normalizationMethod = normalizationMethod,
  filterMethod = filterMethod,
  filterFoldChange = filterFoldChange,
  fitEngine = fitEngine,
  fitDesign = "~ condition",
  sourceProject = "EURATRANS",
  sourcePackage = paste0(
    "chromstaRData ",
    as.character(utils::packageVersion("chromstaRData"))
  )
)

saveRDS(buildMetadata, file.path(outputDir, "euratrans_buildMetadata.rds"),
        compress = "xz")

writeLines(
  utils::capture.output(utils::sessionInfo()),
  file.path(scriptDir, "make-data-sessionInfo.txt")
)


# ---- Size report ------------------------------------------------------------

sizeReport <-
  data.frame(
    file = list.files(outputDir, pattern = "^euratrans_"),
    stringsAsFactors = FALSE
  ) |>
  dplyr::mutate(sizeKb = round(file.size(file.path(outputDir, file)) / 1024, 1)) |>
  dplyr::arrange(dplyr::desc(sizeKb))

print(sizeReport)
message("Total: ", round(sum(sizeReport$sizeKb) / 1024, 2), " MB in ", outputDir)
