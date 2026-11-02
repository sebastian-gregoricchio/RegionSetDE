# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

## Builds the annotation the vignettes run on, so that neither of them has to depend on a
## full TxDb at build time:
##
##   inst/extdata/euratrans_genes.rds          gene start sites of the rat window (main vignette)
##   inst/extdata/euratrans_geneSets.rds       GO biological process sets over those genes
##   inst/extdata/lncapAR/chr19_genes.rds      gene start sites of the human window (peak vignette)
##   inst/extdata/lncapAR/chr19_geneSets.rds   GO biological process sets over those genes
##
## Each file is a slice of the annotation the vignette names in its own chunks, so what a
## vignette shows and what it runs on come from the same source. Run from the package root:
##
##   Rscript inst/scripts/make-example-annotation.R
##
## Needs GenomicFeatures, AnnotationDbi, GO.db, the two TxDb packages and the two OrgDb ones.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(AnnotationDbi)
})

extdataDirectory <- file.path("inst", "extdata")
if (!dir.exists(extdataDirectory)) {
  stop("Run this from the root of the package, '", extdataDirectory, "' does not exist.", call. = FALSE)
}


#' Gene start sites of one window, with their symbols
geneSitesOfWindow <- function(txdb, orgdb, chromosome, windowStart, windowEnd, keyType) {

  startSites <- GenomicRanges::resize(GenomicFeatures::genes(txdb), width = 1, fix = "start")

  windowRange <- GenomicRanges::GRanges(chromosome, IRanges::IRanges(windowStart, windowEnd))
  startSites <- IRanges::subsetByOverlaps(startSites, windowRange, ignore.strand = TRUE)

  if (length(startSites) == 0) {
    stop("No gene start site falls in ", chromosome, ":", windowStart, "-", windowEnd, ".", call. = FALSE)
  }

  # The symbols are for reading only, nothing is ever matched on them
  geneSymbols <- tryCatch(suppressMessages(AnnotationDbi::mapIds(orgdb,
                                                                 keys = names(startSites),
                                                                 column = "SYMBOL",
                                                                 keytype = keyType,
                                                                 multiVals = "first")),
                          error = function(e) {
                            message("  symbols unavailable (", conditionMessage(e), "), leaving them out")
                            return(NULL)
                          })

  S4Vectors::mcols(startSites)$gene.id <- names(startSites)
  S4Vectors::mcols(startSites)$gene.symbol <-
    if (is.null(geneSymbols)) {NA_character_} else {unname(geneSymbols[names(startSites)])}

  names(startSites) <- NULL
  startSites <- GenomeInfoDb::keepSeqlevels(startSites, chromosome, pruning.mode = "coarse")

  return(startSites)
}


#' GO biological process terms over a set of genes, kept at a size a set level test can use
geneSetsOfGenes <- function(orgdb, geneIds, keyType, minSize = 5L, maxSize = 40L, setNumber = 25L) {

  goMapping <- suppressMessages(AnnotationDbi::select(orgdb,
                                                      keys = geneIds,
                                                      columns = c("GOALL", "ONTOLOGYALL"),
                                                      keytype = keyType))

  goColumn <- if ("GOALL" %in% colnames(goMapping)) {"GOALL"} else {"GO"}
  ontologyColumn <- if ("ONTOLOGYALL" %in% colnames(goMapping)) {"ONTOLOGYALL"} else {"ONTOLOGY"}

  goMapping <- goMapping[!is.na(goMapping[[goColumn]]) & goMapping[[ontologyColumn]] == "BP", , drop = FALSE]
  goMapping <- unique(goMapping[, c(goColumn, keyType)])

  setsByTerm <- split(goMapping[[keyType]], goMapping[[goColumn]])

  setSizes <- lengths(setsByTerm)
  setsByTerm <- setsByTerm[setSizes >= minSize & setSizes <= maxSize]

  if (length(setsByTerm) == 0) {
    stop("No GO biological process term of a usable size covers these genes.", call. = FALSE)
  }

  # The most specific terms first, those being the informative ones, then a readable handful
  setsByTerm <- utils::head(setsByTerm[order(lengths(setsByTerm))], setNumber)

  goTerms <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,
                                                    keys = names(setsByTerm),
                                                    columns = "TERM",
                                                    keytype = "GOID"))

  names(setsByTerm) <- goTerms$TERM[match(names(setsByTerm), goTerms$GOID)]

  return(setsByTerm[!is.na(names(setsByTerm))])
}


#' Write one pair of files and print the lines SOURCES.md wants
writeAnnotation <- function(label, txdbName, orgdbName, keyType, chromosome,
                            windowStart, windowEnd, genesFile, setsFile) {

  message("\n== ", label, " ==")
  for (packageName in c(txdbName, orgdbName, "GO.db")) {
    if (!requireNamespace(packageName, quietly = TRUE)) {
      stop("The '", packageName, "' package is needed and is not installed.", call. = FALSE)
    }
  }

  txdb <- getExportedValue(txdbName, txdbName)
  orgdb <- getExportedValue(orgdbName, orgdbName)

  message("Window: ", chromosome, ":", format(windowStart, big.mark = ","),
          "-", format(windowEnd, big.mark = ","))

  startSites <- geneSitesOfWindow(txdb = txdb, orgdb = orgdb, chromosome = chromosome,
                                  windowStart = windowStart, windowEnd = windowEnd, keyType = keyType)
  message(length(startSites), " gene start sites kept.")
  saveRDS(startSites, genesFile, compress = "xz")

  geneSets <- geneSetsOfGenes(orgdb = orgdb,
                              geneIds = unique(S4Vectors::mcols(startSites)$gene.id),
                              keyType = keyType)
  message(length(geneSets), " gene sets kept, holding ",
          min(lengths(geneSets)), " to ", max(lengths(geneSets)), " genes.")
  saveRDS(geneSets, setsFile, compress = "xz")

  message("\nFor SOURCES.md:")
  message("  ", basename(genesFile), ": gene start sites of ", chromosome, ":",
          format(windowStart, big.mark = ","), "-", format(windowEnd, big.mark = ","),
          ", a slice of ", txdbName, " ", utils::packageVersion(txdbName),
          " with symbols from ", orgdbName, " ", utils::packageVersion(orgdbName),
          ". Built by inst/scripts/make-example-annotation.R on ", format(Sys.Date()), ".")
  message("  ", basename(setsFile), ": GO biological process terms over those genes, from ",
          orgdbName, " ", utils::packageVersion(orgdbName), " and GO.db ",
          utils::packageVersion("GO.db"), ", kept at 5 to 40 genes, the 25 most specific. ",
          "Same script, same date.")
}


#-------------------------------#
# The rat example                #
#-------------------------------#
regionTable <- readRDS(file.path(extdataDirectory, "euratrans_regions.rds"))
ratChromosome <- unique(as.character(regionTable$seqnames))

if (length(ratChromosome) != 1) {
  stop("The rat example spans more than one chromosome: ", paste(ratChromosome, collapse = ", "), ".", call. = FALSE)
}

windowFlank <- 1e6L

writeAnnotation(label = "Rat, main vignette",
                txdbName = "TxDb.Rnorvegicus.UCSC.rn4.ensGene",
                orgdbName = "org.Rn.eg.db",
                keyType = "ENSEMBL",
                chromosome = ratChromosome,
                windowStart = max(1L, min(regionTable$start) - windowFlank),
                windowEnd = max(regionTable$end) + windowFlank,
                genesFile = file.path(extdataDirectory, "euratrans_genes.rds"),
                setsFile = file.path(extdataDirectory, "euratrans_geneSets.rds"))


#-------------------------------#
# The human example              #
#-------------------------------#
peakDirectory <- file.path(extdataDirectory, "lncapAR")
peakFiles <- list.files(peakDirectory, pattern = "_peaks\\.narrowPeak$", full.names = TRUE)

if (length(peakFiles) == 0) {
  stop("No narrowPeak file found in '", peakDirectory, "'.", call. = FALSE)
}

peakTable <- do.call(rbind, lapply(peakFiles, utils::read.delim, header = FALSE, comment.char = "#"))
peakChromosome <- unique(as.character(peakTable$V1))

if (length(peakChromosome) != 1) {
  stop("The example peaks span more than one chromosome: ", paste(peakChromosome, collapse = ", "), ".", call. = FALSE)
}

# The peaks are named the Ensembl way and the TxDb the UCSC way, the table follows the TxDb
# so that the vignette shows annotateRegions() reconciling the two
ucscChromosome <- if (grepl("^chr", peakChromosome)) {peakChromosome} else {paste0("chr", peakChromosome)}

writeAnnotation(label = "Human, peak vignette",
                txdbName = "TxDb.Hsapiens.UCSC.hg38.knownGene",
                orgdbName = "org.Hs.eg.db",
                keyType = "ENTREZID",
                chromosome = ucscChromosome,
                windowStart = max(1L, min(peakTable$V2) - windowFlank),
                windowEnd = max(peakTable$V3) + windowFlank,
                genesFile = file.path(peakDirectory, "chr19_genes.rds"),
                setsFile = file.path(peakDirectory, "chr19_geneSets.rds"))

message("\nDone.")
