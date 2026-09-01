#' @title loadCounts
#'
#' @description Imports a count matrix computed outside R, for instance by featureCounts, bedtools multicov or deeptools multiBigwigSummary, and attaches it to the regions of a \code{RegionSetDE} object. The rows of the matrix are matched to the regions either by coordinates or by identifier.
#'
#' @param regionSet \code{RegionSetDE} object returned by \code{\link{loadRegions}}, or a named \code{GRangesList}.
#' @param counts Matrix, data.frame or path to a tab separated file with the counts. Lines starting with \code{#}, such as the featureCounts header, are skipped.
#' @param sampleNames Character vector with the sample names. Default: \code{NULL}, the names of the count columns are used.
#' @param sampleMetadata Data.frame with the sample annotation, stored in the \code{colData}. When it contains a \code{sample} column the rows are matched by name, otherwise they must follow the order of the count columns. Default: \code{NULL}.
#' @param countColumns Character vector or numeric positions of the columns holding the counts. Default: \code{NULL}, every numeric column that is not a coordinate or a standard annotation column.
#' @param coordinateColumns Character vector or numeric positions of the three columns holding chromosome, start and end. Used only when \code{matchBy = "coordinates"}. Default: \code{NULL}, detected from the column names.
#' @param idColumn String or numeric position of the column holding the region identifiers. Used only when \code{matchBy = "id"}. Default: \code{NULL}, the row names are used.
#' @param matchBy String indicating how the rows of \code{counts} are assigned to the regions, either \code{"coordinates"} or \code{"id"}. Default: \code{"coordinates"}.
#' @param startsAt Numeric value with the coordinate system of the count table, \code{1} for 1-based tables such as featureCounts, \code{0} for BED-like tables. Default: \code{1}.
#' @param tileWidth Numeric value with the width of the tiles, to be set when the external matrix has been computed on tiles rather than on whole regions. Default: \code{NULL}.
#' @param partialTiles Logical value indicating whether the trailing shorter tile of each region has been kept. Default: \code{TRUE}.
#' @param missingRegions String indicating what to do with the regions absent from the count table, one among \code{"stop"}, \code{"zero"} or \code{"drop"}. Default: \code{"stop"}.
#' @param librarySizes Numeric vector with the library size of each sample, in the same order as the count columns. Default: \code{NULL}, the column sums of the imported table are used.
#' @param header Logical value indicating whether the file carries a header line. Ignored when \code{counts} is not a file path. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.counts} object with one row per region, or per tile, and one column per sample.
#'
#' @details The column sums of the imported table are a poor substitute for the real library sizes, since they only cover the regions present in the file. When the sequencing depth is known it should be passed through \code{librarySizes}, otherwise the normalisation should rely on factors estimated elsewhere. Rows of the count table that match no region are ignored, which makes it safe to import a genome wide matrix and keep only the sets of interest.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' regionRanges <- SummarizedExperiment::rowRanges(counts)
#'
#' # A count table as it would come out of featureCounts or a coverage tool
#' countTable <- data.frame(
#'   seqnames = as.character(GenomicRanges::seqnames(regionRanges)),
#'   start = GenomicRanges::start(regionRanges),
#'   end = GenomicRanges::end(regionRanges),
#'   SummarizedExperiment::assay(counts, "counts"),
#'   check.names = FALSE)
#'
#' regions <- splitLoadRegions(regionRanges, splitBy = "region.set",
#'                             genomeAssembly = "rn4", verbose = FALSE)
#'
#' reloaded <- loadCounts(regions, counts = countTable, verbose = FALSE)
#' reloaded
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countReads}}, \code{\link{countBigwig}}
#'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
#' @importFrom dplyr select all_of mutate
#' @importFrom utils read.delim
#'
#' @export loadCounts

loadCounts <-
  function(regionSet,
           counts,
           sampleNames = NULL,
           sampleMetadata = NULL,
           countColumns = NULL,
           coordinateColumns = NULL,
           idColumn = NULL,
           matchBy = "coordinates",
           startsAt = 1,
           tileWidth = NULL,
           partialTiles = TRUE,
           missingRegions = "stop",
           librarySizes = NULL,
           header = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    matchBy <- tolower(matchBy[1])
    if (!(matchBy %in% c("coordinates", "id"))) {
      stop("The 'matchBy' parameter must be either 'coordinates' or 'id'.", call. = FALSE)
    }

    missingRegions <- tolower(missingRegions[1])
    if (!(missingRegions %in% c("stop", "zero", "drop"))) {
      stop("The 'missingRegions' parameter must be one among 'stop', 'zero' or 'drop'.", call. = FALSE)
    }

    startsAt <- as.integer(startsAt[1])
    if (!(startsAt %in% c(0L, 1L))) {
      stop("The 'startsAt' parameter must be either 0 or 1.", call. = FALSE)
    }

    #------------------------#
    # Read the count table   #
    #------------------------#
    if (is.character(counts) & length(counts) == 1) {
      if (!file.exists(counts)) {
        stop("The file '", counts, "' does not exist.", call. = FALSE)
      }
      countTable <- utils::read.delim(file = counts, header = header, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
    } else if (is.matrix(counts)) {
      countTable <- as.data.frame(counts, stringsAsFactors = FALSE)
      countTable$region.id <- rownames(counts)
    } else if (is.data.frame(counts)) {
      countTable <- as.data.frame(counts, stringsAsFactors = FALSE)
    } else {
      stop("The 'counts' parameter must be a matrix, a data.frame or the path to a count table.", call. = FALSE)
    }

    #-------------------------------#
    # Identify the columns to use   #
    #-------------------------------#
    # Positions are accepted everywhere, working on names afterwards keeps the rest of the code readable
    if (is.numeric(countColumns)) {countColumns <- colnames(countTable)[countColumns]}
    if (is.numeric(coordinateColumns)) {coordinateColumns <- colnames(countTable)[coordinateColumns]}
    if (is.numeric(idColumn)) {idColumn <- colnames(countTable)[idColumn]}

    annotationColumns <- c("geneid", "chr", "chrom", "chromosome", "seqnames", "start", "end", "stop",
                           "strand", "length", "width", "name", "score", "id", "region.id", "region.set", "tile.id")

    if (is.null(countColumns)) {
      isCountColumn <- vapply(countTable, is.numeric, logical(1)) & !(tolower(colnames(countTable)) %in% annotationColumns)
      countColumns <- colnames(countTable)[isCountColumn]
    }

    if (length(countColumns) == 0) {
      stop("No count column could be identified, provide them through the 'countColumns' parameter.", call. = FALSE)
    }

    absentColumns <- setdiff(countColumns, colnames(countTable))
    if (length(absentColumns) > 0) {
      stop("The following columns are missing from the count table: ", paste(absentColumns, collapse = ", "), ".", call. = FALSE)
    }

    #-------------------------------#
    # Build the matching keys       #
    #-------------------------------#
    allRegions <- .flattenRegionSets(regionSet = regionSet,
                                     tileWidth = tileWidth,
                                     partialTiles = partialTiles,
                                     verbose = verbose)

    if (matchBy == "coordinates") {
      if (is.null(coordinateColumns)) {
        chromosomeColumn <- colnames(countTable)[tolower(colnames(countTable)) %in% c("chr", "chrom", "chromosome", "seqnames")][1]
        startColumn <- colnames(countTable)[tolower(colnames(countTable)) %in% c("start")][1]
        endColumn <- colnames(countTable)[tolower(colnames(countTable)) %in% c("end", "stop")][1]
        coordinateColumns <- c(chromosomeColumn, startColumn, endColumn)
      }

      if (any(is.na(coordinateColumns)) | length(coordinateColumns) != 3) {
        stop("The coordinate columns could not be identified, provide them through the 'coordinateColumns' parameter.", call. = FALSE)
      }

      # BED tables count from zero, the shift aligns them to the 1-based GRanges coordinates
      tableKey <- paste0(as.character(countTable[[coordinateColumns[1]]]), ":",
                         as.numeric(countTable[[coordinateColumns[2]]]) + (1L - startsAt), "-",
                         as.numeric(countTable[[coordinateColumns[3]]]))

      regionMatchKey <- paste0(as.character(GenomeInfoDb::seqnames(allRegions)), ":",
                               BiocGenerics::start(allRegions), "-",
                               BiocGenerics::end(allRegions))
    } else {
      if (is.null(idColumn)) {
        if (is.null(rownames(countTable))) {
          stop("The count table has no row names, provide the identifiers through the 'idColumn' parameter.", call. = FALSE)
        }
        tableKey <- rownames(countTable)
      } else {
        tableKey <- as.character(countTable[[idColumn]])
      }

      regionMatchKey <- as.character(S4Vectors::mcols(allRegions)$region.id)
    }

    if (any(duplicated(tableKey))) {
      stop("The count table contains duplicated regions, they must be unique to be assigned unambiguously.", call. = FALSE)
    }

    #-------------------------------#
    # Match the rows to the regions #
    #-------------------------------#
    matchIndex <- match(regionMatchKey, tableKey)
    unmatchedRegions <- sum(is.na(matchIndex))

    if (unmatchedRegions > 0) {
      if (missingRegions == "stop") {
        stop(unmatchedRegions, " regions out of ", length(allRegions), " have no matching row in the count table. ",
             "Check the 'startsAt' parameter, or set 'missingRegions' to 'zero' or 'drop'.", call. = FALSE)
      }
      if (isTRUE(verbose)) {
        message(unmatchedRegions, " regions have no matching row in the count table and have been ",
                ifelse(missingRegions == "drop", "removed.", "set to zero."))
      }
    }

    countMatrix <- as.matrix(dplyr::select(countTable, dplyr::all_of(countColumns)))
    fullColumnSums <- colSums(countMatrix, na.rm = TRUE)

    countMatrix <- countMatrix[matchIndex, , drop = FALSE]

    if (unmatchedRegions > 0) {
      if (missingRegions == "zero") {
        countMatrix[is.na(matchIndex), ] <- 0
      } else {
        countMatrix <- countMatrix[!is.na(matchIndex), , drop = FALSE]
        allRegions <- allRegions[!is.na(matchIndex)]
      }
    }

    #---------------------#
    # Assemble the object #
    #---------------------#
    sampleTable <- .buildSampleTable(files = countColumns,
                                     sampleNames = sampleNames,
                                     sampleMetadata = sampleMetadata,
                                     fileColumn = "source.column",
                                     extensionPattern = "\\.bam$")

    if (!is.null(librarySizes)) {
      if (length(librarySizes) != length(countColumns)) {
        stop("The 'librarySizes' parameter must have one value per count column.", call. = FALSE)
      }
      sampleTable <- dplyr::mutate(sampleTable, library.size = as.numeric(librarySizes))
    } else {
      # The column sums only cover the imported regions, they are a stand-in until real depths are supplied
      sampleTable <- dplyr::mutate(sampleTable, library.size = as.numeric(fullColumnSums))
    }

    newParameters <- list(loadCounts = list(source = ifelse(is.character(counts) & length(counts) == 1, counts, "in-memory table"),
                                            countColumns = countColumns,
                                            matchBy = matchBy,
                                            startsAt = startsAt,
                                            tileWidth = tileWidth,
                                            partialTiles = partialTiles,
                                            missingRegions = missingRegions,
                                            librarySizes.supplied = !is.null(librarySizes)))

    countsObject <- .newCountsObject(countMatrix = countMatrix,
                                     regions = allRegions,
                                     sampleTable = sampleTable,
                                     provenance = .provenanceSlots(regionSet),
                                     countingLevel = "region",
                                     newParameters = newParameters,
                                     metadataList = list(signal.type = "external"))

    if (isTRUE(verbose)) {
      message("Imported ", nrow(countMatrix), " regions and ", ncol(countMatrix), " samples.")
    }

    return(countsObject)
  } # END function
