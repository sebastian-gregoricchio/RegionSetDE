## make-data.R
##
## Builds the example objects shipped in inst/extdata for RegionSetDE.
##
## Source data: NF-YA ChIP-seq in mouse embryonic stem cells and terminal
## neurons (GSE25532, Tiwari et al. 2012), aligned to mm10 and distributed as
## sorted, indexed, duplicate-marked BAM files by the chipseqDBData package.
## Two biological replicates per condition plus one input library.
##
## Analysis is restricted to chr17 and chr19 to keep the shipped objects small.
##
## Run once from the package root. The script is not part of the build and
## belongs in .Rbuildignore. Expect the first run to take a while: the BAM files
## are downloaded from ExperimentHub and the motif scan covers two whole
## chromosomes.
##
## Author: Sebastian Gregoricchio


# ---- Parameters -------------------------------------------------------------

randomSeed <- 20260101L
targetChromosomes <- c("chr17", "chr19")
regionWidth <- 300L
promoterUpstream <- 1000L
promoterDownstream <- 500L
motifPvalueCutoff <- 1e-05
nullSetSize <- 1000L
maxCrossCorrelationDistance <- 500L

outputDir <- file.path("inst", "extdata")
scriptDir <- file.path("inst", "scripts")

set.seed(randomSeed)


# ---- Dependency check -------------------------------------------------------

# Everything here runs offline, so none of these belong in DESCRIPTION.
requiredPackages <- c(
  "chipseqDBData",
  "AnnotationHub",
  "ATACseqTFEA",
  "BSgenome.Mmusculus.UCSC.mm10",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "GenomicFeatures",
  "GenomicRanges",
  "GenomeInfoDb",
  "IRanges",
  "Rsamtools",
  "BiocGenerics",
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


# ---- BAM files and sample sheet ---------------------------------------------

# NFYAData() downloads to the ExperimentHub cache and returns Name, Description
# and a List of BamFile objects.
nfyaData <- chipseqDBData::NFYAData()

if (nrow(nfyaData) != 5) {
  stop("NFYAData() returned ", nrow(nfyaData), " libraries; 5 were expected.")
}

bamPaths <- vapply(
  as.list(nfyaData$Path),
  function(bamFile) BiocGenerics::path(bamFile),
  character(1)
)

# Conditions and replicate numbers come straight out of the Description strings
# ("NF-YA ESC (1)", "NF-YA TN (2)", "Input").
sampleSheet <-
  data.frame(
    sampleName = as.character(nfyaData$Name),
    description = as.character(nfyaData$Description),
    bamPath = unname(bamPaths),
    stringsAsFactors = FALSE
  ) |>
  dplyr::mutate(
    condition = dplyr::case_when(
      grepl("ESC", description) ~ "ESC",
      grepl("TN", description) ~ "TN",
      TRUE ~ "input"
    ),
    replicate = as.integer(
      dplyr::if_else(
        grepl("\\([0-9]+\\)", description),
        sub(".*\\(([0-9]+)\\).*", "\\1", description),
        "1"
      )
    ),
    isInput = condition == "input"
  ) |>
  dplyr::arrange(condition, replicate)

if (!all(c("ESC", "TN", "input") %in% sampleSheet$condition)) {
  stop("Could not assign ESC, TN and input conditions from the NFYAData descriptions.")
}

chipSamples <- sampleSheet |> dplyr::filter(!isInput)
inputSample <- sampleSheet |> dplyr::filter(isInput)


# ---- Exclusion list ---------------------------------------------------------

# mm10 ENCODE blacklist v2 (Amemiya et al. 2019), read straight from the
# Boyle Lab release rather than through AnnotationHub.
blacklistUrl <- paste0(
  "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/mm10-blacklist.v2.bed.gz"
)

blacklist <-
  rtracklayer::import(blacklistUrl, format = "BED") |>
  GenomeInfoDb::keepSeqlevels(targetChromosomes, pruning.mode = "coarse")

GenomeInfoDb::genome(blacklist) <- "mm10"


# ---- Fragment extension length ----------------------------------------------

# Single-end reads, so the extension length is estimated once here and hard-coded
# into the vignette rather than recomputed at build time.
crossCorrelation <- csaw::correlateReads(
  chipSamples$bamPath[1],
  max.dist = maxCrossCorrelationDistance,
  param = csaw::readParam(dedup = TRUE, restrict = targetChromosomes)
)

extensionLength <- csaw::maximizeCcf(crossCorrelation)
message("Estimated fragment extension length: ", extensionLength, " bp")


# ---- NF-Y binding sites -----------------------------------------------------

# ATACseqTFEA ships a merged PWMatrixList; keep the NF-Y entries (CCAAT box).
motifList <- readRDS(
  system.file("extdata", "PWMatrixList.rds", package = "ATACseqTFEA", mustWork = TRUE)
)

nfyMotifs <- motifList[grepl("NFY", names(motifList), ignore.case = TRUE)]

if (length(nfyMotifs) == 0) {
  stop("No NF-Y motif found in the ATACseqTFEA PWMatrixList. Fall back to JASPAR via TFBSTools.")
}

message("Scanning ", length(nfyMotifs), " NF-Y motif(s) over ",
        paste(targetChromosomes, collapse = " and "), ".")

# Arguments are passed positionally to match the prepareBindingSites() signature.
# If memory becomes an issue, loop over chromosomes and pass the grange argument.
motifSites <- ATACseqTFEA::prepareBindingSites(
  nfyMotifs,
  BSgenome.Mmusculus.UCSC.mm10::Mmusculus,
  targetChromosomes,
  p.cutoff = motifPvalueCutoff
)

motifSites <- GenomicRanges::granges(motifSites)
GenomicRanges::strand(motifSites) <- "*"


# ---- Promoters --------------------------------------------------------------

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

promoterRegions <-
  GenomicFeatures::promoters(
    txdb,
    upstream = promoterUpstream,
    downstream = promoterDownstream
  ) |>
  GenomeInfoDb::keepSeqlevels(targetChromosomes, pruning.mode = "coarse") |>
  GenomicRanges::reduce(ignore.strand = TRUE)

transcriptStarts <-
  GenomicFeatures::transcripts(txdb) |>
  GenomeInfoDb::keepSeqlevels(targetChromosomes, pruning.mode = "coarse") |>
  GenomicRanges::resize(width = 1L, fix = "start")

GenomicRanges::strand(transcriptStarts) <- "*"
transcriptStarts <- unique(GenomicRanges::granges(transcriptStarts))


# ---- Region universe --------------------------------------------------------

# Candidate anchors are motif midpoints and transcription start sites. Each one
# becomes a fixed-width window so that counts stay comparable across sets.
candidateAnchors <- c(
  GenomicRanges::resize(motifSites, width = 1L, fix = "center"),
  transcriptStarts
)

GenomeInfoDb::seqlevels(candidateAnchors, pruning.mode = "coarse") <- targetChromosomes
GenomeInfoDb::seqinfo(candidateAnchors) <-
  GenomeInfoDb::seqinfo(BSgenome.Mmusculus.UCSC.mm10::Mmusculus)[targetChromosomes]

candidateRegions <-
  candidateAnchors |>
  GenomicRanges::resize(width = regionWidth, fix = "center") |>
  GenomicRanges::trim() |>
  sort()

# Windows clipped at a chromosome edge are dropped rather than kept ragged.
candidateRegions <- candidateRegions[GenomicRanges::width(candidateRegions) == regionWidth]

# disjointBins assigns overlapping ranges to different bins, so the first bin is
# a non-overlapping subset of the candidates.
candidateRegions <- candidateRegions[GenomicRanges::disjointBins(candidateRegions) == 1L]

# Anything touching the exclusion list goes.
candidateRegions <- candidateRegions[
  !IRanges::overlapsAny(candidateRegions, blacklist, ignore.strand = TRUE)
]


# ---- Null control regions ---------------------------------------------------

# Random windows matched in width, kept clear of motifs, promoters, the exclusion
# list and the rest of the universe. These should come out non-significant and
# are what shows the method is not calling hits everywhere.
chromosomeLengths <-
  GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10::Mmusculus)[targetChromosomes]

sampledStarts <- lapply(targetChromosomes, function(chromosome) {
  nDraws <- nullSetSize * 10L
  GenomicRanges::GRanges(
    seqnames = chromosome,
    ranges = IRanges::IRanges(
      start = sample.int(chromosomeLengths[[chromosome]] - regionWidth, nDraws),
      width = regionWidth
    )
  )
})

nullCandidates <- do.call(c, sampledStarts)
GenomeInfoDb::seqinfo(nullCandidates) <-
  GenomeInfoDb::seqinfo(BSgenome.Mmusculus.UCSC.mm10::Mmusculus)[targetChromosomes]

occupiedRanges <- c(
  GenomicRanges::granges(candidateRegions),
  GenomicRanges::granges(motifSites),
  GenomicRanges::granges(promoterRegions),
  GenomicRanges::granges(blacklist)
)

nullRegions <- nullCandidates[
  !IRanges::overlapsAny(nullCandidates, occupiedRanges, ignore.strand = TRUE)
]

nullRegions <- nullRegions[GenomicRanges::disjointBins(sort(nullRegions)) == 1L]

if (length(nullRegions) < nullSetSize) {
  stop("Only ", length(nullRegions), " null regions survived filtering; ",
       nullSetSize, " were requested.")
}

nullRegions <- sort(sample(nullRegions, nullSetSize))


# ---- Region table and region sets -------------------------------------------

allRegions <- sort(c(GenomicRanges::granges(candidateRegions), nullRegions))

regions <-
  data.frame(
    seqnames = as.character(GenomicRanges::seqnames(allRegions)),
    start = GenomicRanges::start(allRegions),
    end = GenomicRanges::end(allRegions),
    stringsAsFactors = FALSE
  ) |>
  dplyr::mutate(
    regionId = sprintf("region_%05d", dplyr::row_number()),
    hasMotif = IRanges::overlapsAny(allRegions, motifSites, ignore.strand = TRUE),
    isPromoter = IRanges::overlapsAny(allRegions, promoterRegions, ignore.strand = TRUE),
    isNull = IRanges::overlapsAny(allRegions, nullRegions, type = "equal")
  ) |>
  dplyr::select(regionId, seqnames, start, end, hasMotif, isPromoter, isNull)

# Four sets: the expected responders, a distal comparison, a motif-free contrast
# and the null control.
regionSets <-
  dplyr::bind_rows(
    regions |>
      dplyr::filter(hasMotif, isPromoter) |>
      dplyr::mutate(setName = "promoter_CCAAT"),
    regions |>
      dplyr::filter(hasMotif, !isPromoter) |>
      dplyr::mutate(setName = "distal_CCAAT"),
    regions |>
      dplyr::filter(!hasMotif, isPromoter) |>
      dplyr::mutate(setName = "promoter_noMotif"),
    regions |>
      dplyr::filter(isNull) |>
      dplyr::mutate(setName = "random_null")
  ) |>
  dplyr::select(setName, regionId) |>
  dplyr::arrange(setName, regionId)

setSizes <- regionSets |> dplyr::count(setName, name = "nRegions")
print(setSizes)

if (any(setSizes$nRegions < 100)) {
  stop("At least one region set holds fewer than 100 regions. Loosen the motif ",
       "p-value cutoff or widen the promoter window.")
}


# ---- Counting ---------------------------------------------------------------

# NOTE: adjust the argument names below once countReads() and countBackground()
# are frozen. Reads are single-end and duplicates are marked but not removed in
# the source BAM files.
regionsGRanges <- GenomicRanges::makeGRangesFromDataFrame(
  regions,
  keep.extra.columns = TRUE
)
names(regionsGRanges) <- regions$regionId

readCounts <- RegionSetDE::countReads(
  bamFiles = chipSamples$bamPath,
  sampleNames = chipSamples$sampleName,
  regions = regionsGRanges,
  pairedEnd = FALSE,
  extensionLength = extensionLength,
  ignoreDuplicates = TRUE,
  blacklist = blacklist
)

backgroundCounts <- RegionSetDE::countBackground(
  bamFiles = c(chipSamples$bamPath, inputSample$bamPath),
  sampleNames = c(chipSamples$sampleName, inputSample$sampleName),
  regions = regionsGRanges,
  pairedEnd = FALSE,
  extensionLength = extensionLength,
  ignoreDuplicates = TRUE,
  blacklist = blacklist
)


# ---- Save -------------------------------------------------------------------

objectsToSave <- list(
  "nfya_sampleSheet.rds" = sampleSheet |> dplyr::select(-bamPath),
  "nfya_regions.rds" = regions,
  "nfya_regionSets.rds" = regionSets,
  "nfya_blacklist.rds" = blacklist,
  "nfya_readCounts.rds" = readCounts,
  "nfya_backgroundCounts.rds" = backgroundCounts
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
  targetChromosomes = targetChromosomes,
  genomeBuild = "mm10",
  regionWidth = regionWidth,
  promoterWindow = c(upstream = promoterUpstream, downstream = promoterDownstream),
  motifPvalueCutoff = motifPvalueCutoff,
  extensionLength = extensionLength,
  blacklistRecord = blacklistId,
  motifNames = names(nfyMotifs),
  sourceAccession = "GSE25532",
  sourcePackage = paste0(
    "chipseqDBData ",
    as.character(utils::packageVersion("chipseqDBData"))
  )
)

saveRDS(buildMetadata, file.path(outputDir, "nfya_buildMetadata.rds"), compress = "xz")

writeLines(
  utils::capture.output(utils::sessionInfo()),
  file.path(scriptDir, "make-data-sessionInfo.txt")
)


# ---- Size report ------------------------------------------------------------

sizeReport <-
  data.frame(
    file = list.files(outputDir, pattern = "^nfya_"),
    stringsAsFactors = FALSE
  ) |>
  dplyr::mutate(sizeKb = round(file.size(file.path(outputDir, file)) / 1024, 1)) |>
  dplyr::arrange(dplyr::desc(sizeKb))

print(sizeReport)
message("Total: ", round(sum(sizeReport$sizeKb) / 1024, 2), " MB in ", outputDir)


