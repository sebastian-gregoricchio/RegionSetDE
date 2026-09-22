# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title countTable
#'
#' @description Returns the values of a counts object as a table, one row per region or per tile, with the coordinates and the annotation of each row next to the values of the samples. The values can be raw or normalised, and come in wide or long format, or as a plain matrix. On a tiled object the tiles can be combined back into their region.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results} or \code{RegionSetDE.setResults} object, or either of the two list classes holding several contrasts. For the results the counts are the ones carried inside them, which requires the test to have been run with \code{carryCounts = TRUE}.
#' @param level String indicating whether the table must have one row per region (\code{"region"}) or one row per tile (\code{"tile"}). On a tiled object \code{"region"} combines the tiles of each region into a single row. Default: \code{"region"}.
#' @param normalized Logical value to indicate whether the normalised values must be returned instead of the raw ones. Default: \code{FALSE}.
#' @param format String with the shape of the output: \code{"wide"} gives one column per sample, \code{"long"} one row per row of the object and sample, with the \code{colData} of the samples attached, and \code{"matrix"} a numeric matrix with the rows named after the regions. Default: \code{"wide"}.
#' @param set Character vector with the names of the region sets to keep. Default: \code{NULL}, all of them.
#' @param tileSummary String indicating how the tiles are combined into their region when \code{level = "region"}, one among \code{"sum"}, \code{"mean"}, \code{"max"} and \code{"min"}. Default: \code{NULL}, read from the way the object was counted.
#' @param extraColumns Annotation carried by the regions that must be added to the table, after the coordinates. Either \code{TRUE} for every column of the \code{rowData} beyond the ones the package writes itself, \code{FALSE} for none, or a character vector naming the ones wanted. Ignored when \code{format = "matrix"}. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#' @param ... Arguments passed on to the method for \code{RegionSetDE.counts} objects.
#'
#' @return With \code{format = "wide"}, a data.frame with one row per region, or per tile, described by the \code{region.set}, \code{region.id}, \code{tile.id} (tiles only), \code{seqnames}, \code{start}, \code{end} and \code{width} columns, followed by the annotation of the regions and by one column per sample. When the tiles have been combined, \code{n.tiles} reports how many of them each region was built from. With \code{format = "long"}, a data.frame with the same description repeated for every sample, a \code{sample} column, the values in a column named after the assay they were read from (\code{counts}, or \code{norm.counts} for the normalised ones), and the \code{colData} of the samples. With \code{format = "matrix"}, a numeric matrix with one column per sample and the rows named \code{"set|id"}, or \code{"set|id|tileN"} for the tiles.
#'
#' @details With \code{level = "region"} on a tiled object the tiles of each region are combined into one row, spanning from the first to the last tile present. Read counts are summed. Signal read from bigWig files follows the \code{summaryFunction} used by \code{\link{countBigwig}}: a sum stays a sum, a mean is averaged with the width of each tile as its weight, so that a shorter trailing tile counts for the bases it covers, and maxima and minima stay maxima and minima. \code{tileSummary} overrides the rule, for instance for a matrix imported through \code{\link{loadCounts}} that holds a mean signal rather than counts.
#'
#' For bigWig signal the combination is exact, since the values are integrated base by base. The one exception is a mean computed with \code{missingAsZero = FALSE}, which is taken over the covered bases only, while the tiles are weighted by their full width. For reads the combination is not exact. \code{\link{countReads}} counts a fragment in every row it overlaps, so a fragment lying across the border between two tiles is counted in both, and the sum over the tiles is higher than the count of the same region taken whole. The excess grows with the fragment length relative to the tile width: a few percent with tiles of several kb, about double when the tiles are as narrow as the fragments. It is similar across libraries with similar fragment lengths, so the values still compare between samples, but the number of fragments falling in a region can only come from counting the regions without tiles.
#'
#' \code{\link{testRegionSets}} combines the tiles differently. It averages them (\code{tileHandling = "collapse"}) so that every region weighs the same within its set, whatever its width. A table of counts describes the regions themselves, and the sum is the value closest to a region counted in one piece.
#'
#' The normalised assay holds the raw counts divided by the scaling factor of each sample, so the region values built from it equal the combined raw counts divided by the same factor. The exception is \code{method = "loess"}, where every row carries its own offset: each tile is corrected at its own abundance before the tiles are combined, which is not what a loess fit on the region counts would return.
#'
#' After \code{\link{filterRegions}} a region holds only the tiles that passed the filter, and its value covers those tiles alone.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' # Raw counts, one column per sample
#' head(countTable(counts), 3)
#'
#' # Normalised values in long format, with the sample annotation attached for ggplot2
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#' longTable <- countTable(counts, normalized = TRUE, format = "long")
#' head(longTable, 3)
#'
#' # A plain matrix of one set, for ComplexHeatmap or any other tool
#' countMatrix <- countTable(counts, normalized = TRUE, format = "matrix", set = "promoterCpG")
#' dim(countMatrix)
#'
#' # On a tiled object the tiles can be kept, or summed back into their region
#' bamFile <- system.file("extdata", "ex1.bam", package = "Rsamtools")
#'
#' exampleRegions <- GenomicRanges::GRanges(
#'   seqnames = rep(c("seq1", "seq2"), each = 2),
#'   ranges = IRanges::IRanges(start = rep(c(1, 800), 2), width = 400))
#'
#' exampleRegions$setName <- rep(c("firstSet", "secondSet"), each = 2)
#'
#' exampleSets <- splitLoadRegions(exampleRegions, splitBy = "setName",
#'                                 seqlevelsStyle = NULL, verbose = FALSE)
#'
#' tiledCounts <- countReads(exampleSets, bamFiles = bamFile,
#'                           sampleNames = "example", tileWidth = 100,
#'                           verbose = FALSE)
#'
#' countTable(tiledCounts, level = "tile")
#' countTable(tiledCounts, level = "region")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{normalizeCounts}}, \code{\link{resultCounts}}, \code{\link{fitCounts}}
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export

setGeneric(name = "countTable", def = function(object, ...) {standardGeneric("countTable")})


#' @rdname countTable
#' @export

setMethod(f = "countTable",
          signature = "RegionSetDE.counts",
          definition = function(object,
                                level = "region",
                                normalized = FALSE,
                                format = "wide",
                                set = NULL,
                                tileSummary = NULL,
                                extraColumns = TRUE,
                                verbose = TRUE) {
            return(.buildCountTable(counts = object,
                                    level = level,
                                    normalized = normalized,
                                    format = format,
                                    set = set,
                                    tileSummary = tileSummary,
                                    extraColumns = extraColumns,
                                    verbose = verbose))
          })


#' @rdname countTable
#' @export

setMethod(f = "countTable",
          signature = "RegionSetDE.fit",
          definition = function(object, ...) {
            # The counts the model was fitted on, filtered and restricted to the samples used
            return(countTable(fitCounts(object), ...))
          })


#' @rdname countTable
#' @export

setMethod(f = "countTable",
          signature = "RegionSetDE.results",
          definition = function(object, ...) {
            return(countTable(.carriedCounts(object), ...))
          })


#' @rdname countTable
#' @export

setMethod(f = "countTable",
          signature = "RegionSetDE.setResults",
          definition = function(object, ...) {
            return(countTable(.carriedCounts(object), ...))
          })


#' @rdname countTable
#' @export

setMethod(f = "countTable",
          signature = "RegionSetDE.resultsList",
          definition = function(object, ...) {
            # Every contrast comes from the same fit, so any of them holds the same counts
            return(countTable(.carriedCounts(object), ...))
          })


#' @rdname countTable
#' @export

setMethod(f = "countTable",
          signature = "RegionSetDE.setResultsList",
          definition = function(object, ...) {
            return(countTable(.carriedCounts(object), ...))
          })




#' @title .buildCountTable
#'
#' @description Builds the table returned by \code{countTable} from a \code{RegionSetDE.counts} object.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param level String, either \code{"region"} or \code{"tile"}.
#' @param normalized Logical value indicating whether the normalised assay must be read.
#' @param format String, one among \code{"wide"}, \code{"long"} and \code{"matrix"}.
#' @param set Character vector with the region sets to keep, or \code{NULL}.
#' @param tileSummary String with the rule combining the tiles, or \code{NULL}.
#' @param extraColumns \code{TRUE}, \code{FALSE} or a character vector with the annotation columns wanted.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A data.frame or a numeric matrix, see \code{countTable}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment assay assayNames colData rowData rowRanges
#' @importFrom S4Vectors metadata mcols
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end width
#' @importFrom dplyr filter mutate group_by summarise first n row_number relocate left_join
#' @importFrom rlang .data
#'
#' @keywords internal

.buildCountTable <-
  function(counts,
           level = "region",
           normalized = FALSE,
           format = "wide",
           set = NULL,
           tileSummary = NULL,
           extraColumns = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!(level[1] %in% c("region", "tile"))) {
      stop("The 'level' parameter must be either 'region' or 'tile'.", call. = FALSE)
    }
    level <- level[1]

    if (!(format[1] %in% c("wide", "long", "matrix"))) {
      stop("The 'format' parameter must be one of 'wide', 'long', 'matrix'.", call. = FALSE)
    }
    format <- format[1]

    if (!is.logical(normalized) | length(normalized) != 1 | anyNA(normalized)) {
      stop("The 'normalized' parameter must be a single logical value.", call. = FALSE)
    }

    isTiled <- identical(counts@counting.level, "tile")

    if (level == "tile" & !isTiled) {
      stop("The counts were not tiled, no tile level table is available.", call. = FALSE)
    }

    #-------------------------------#
    # Assay to read                 #
    #-------------------------------#
    assayName <- "counts"

    if (isTRUE(normalized)) {
      normalizationInfo <- S4Vectors::metadata(counts)$normalization

      if (is.null(normalizationInfo)) {
        stop("The object carries no normalisation, run normalizeCounts() first.", call. = FALSE)
      }

      assayName <- normalizationInfo$normalized.assay

      if (!(assayName %in% SummarizedExperiment::assayNames(counts))) {
        stop("The normalised assay '", assayName, "' is absent from the object.", call. = FALSE)
      }
    }

    #-------------------------------#
    # Rows of the requested sets    #
    #-------------------------------#
    rowRangesObject <- SummarizedExperiment::rowRanges(counts)

    # The counting functions name every row, the fallback only covers objects assembled by hand
    rowKeys <- rownames(counts)
    if (is.null(rowKeys)) {
      rowKeys <- paste(S4Vectors::mcols(rowRangesObject)$region.set, S4Vectors::mcols(rowRangesObject)$region.id, sep = "|")
      if (isTiled) {
        rowKeys <- paste(rowKeys, paste0("tile", S4Vectors::mcols(rowRangesObject)$tile.id), sep = "|")
      }
    }

    rowTable <- data.frame(row.index = seq_along(rowRangesObject),
                           row.key = rowKeys,
                           region.set = as.character(S4Vectors::mcols(rowRangesObject)$region.set),
                           region.id = as.character(S4Vectors::mcols(rowRangesObject)$region.id),
                           tile.id = S4Vectors::mcols(rowRangesObject)$tile.id,
                           seqnames = as.character(GenomeInfoDb::seqnames(rowRangesObject)),
                           start = BiocGenerics::start(rowRangesObject),
                           end = BiocGenerics::end(rowRangesObject),
                           width = BiocGenerics::width(rowRangesObject),
                           stringsAsFactors = FALSE)

    if (!is.null(set)) {
      absentSets <- setdiff(set, unique(rowTable$region.set))
      if (length(absentSets) > 0) {
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
      }
      rowTable <- dplyr::filter(rowTable, .data$region.set %in% set)
    }

    valueMatrix <- as.matrix(SummarizedExperiment::assay(counts, assayName))[rowTable$row.index, , drop = FALSE]

    # A sample named like a coordinate would give the wide table two columns with the same name
    identifierColumns <- c("region.set", "region.id", "tile.id", "seqnames", "start", "end", "width", "n.tiles")
    sampleNames <- colnames(valueMatrix)

    if (format == "wide" & any(sampleNames %in% identifierColumns)) {
      stop("The following sample names are also the names of identifier columns, use format = 'long': ",
           paste(intersect(sampleNames, identifierColumns), collapse = ", "), ".", call. = FALSE)
    }

    #-------------------------------#
    # Annotation of the rows        #
    #-------------------------------#
    sampleTable <- as.data.frame(SummarizedExperiment::colData(counts))
    sampleTable$sample <- sampleNames
    sampleTable <- dplyr::relocate(sampleTable, "sample")

    annotationTable <- data.frame(row.names = seq_len(nrow(rowTable)))

    if (format != "matrix") {
      reservedColumns <- if (format == "wide") {c(identifierColumns, sampleNames)} else {c(identifierColumns, assayName, colnames(sampleTable))}

      annotationTable <- .extraRowColumns(rowTable = as.data.frame(SummarizedExperiment::rowData(counts))[rowTable$row.index, , drop = FALSE],
                                          extraColumns = extraColumns,
                                          reserved = reservedColumns,
                                          verbose = verbose)
    }

    #-------------------------------#
    # Tiles back into their region  #
    #-------------------------------#
    if (level == "region" & isTiled) {
      summaryRule <- .tileSummaryRule(counts = counts, tileSummary = tileSummary)

      # Matching the keys rather than relying on the row order keeps the tiles together whatever the sorting
      regionKey <- paste(rowTable$region.set, rowTable$region.id, sep = "|")
      rowTable <- dplyr::mutate(rowTable,
                                table.row = dplyr::row_number(),
                                region.index = match(regionKey, unique(regionKey)))

      valueMatrix <- .combineTileValues(valueMatrix = valueMatrix,
                                        regionIndex = rowTable$region.index,
                                        tileWidths = rowTable$width,
                                        summaryRule = summaryRule)

      rowTable <- dplyr::summarise(dplyr::group_by(rowTable, .data$region.index),
                                   first.row = dplyr::first(.data$table.row),
                                   region.set = dplyr::first(.data$region.set),
                                   region.id = dplyr::first(.data$region.id),
                                   seqnames = dplyr::first(.data$seqnames),
                                   start = min(.data$start),
                                   end = max(.data$end),
                                   n.tiles = dplyr::n(),
                                   .groups = "drop")

      rowTable <- as.data.frame(rowTable, stringsAsFactors = FALSE)
      rowTable$width <- rowTable$end - rowTable$start + 1
      rowTable$row.key <- paste(rowTable$region.set, rowTable$region.id, sep = "|")

      # The tiles inherit the annotation of their region, so the first one speaks for all of them
      annotationTable <- annotationTable[rowTable$first.row, , drop = FALSE]

      if (isTRUE(verbose)) {
        signalType <- S4Vectors::metadata(counts)$signal.type
        if (summaryRule == "sum" & (is.null(signalType) || signalType == "reads")) {
          message("The tiles have been summed into their regions. A fragment overlapping two tiles is counted in both, ",
                  "so these values are higher than the counts of the same regions taken whole.")
        } else if (identical(signalType, "external") & is.null(tileSummary)) {
          message("The imported tiles have been summed into their regions, set 'tileSummary' if they hold a mean signal rather than counts.")
        }
      }
    }

    rownames(valueMatrix) <- rowTable$row.key
    colnames(valueMatrix) <- sampleNames

    #-------------------------------#
    # Shape of the output           #
    #-------------------------------#
    if (format == "matrix") {
      return(valueMatrix)
    }

    # tile.id means nothing once the tiles are combined, and n.tiles nothing before
    descriptionColumns <- if (level == "tile") {
      c("region.set", "region.id", "tile.id", "seqnames", "start", "end", "width")
    } else if (isTiled) {
      c("region.set", "region.id", "seqnames", "start", "end", "width", "n.tiles")
    } else {
      c("region.set", "region.id", "seqnames", "start", "end", "width")
    }

    descriptionTable <- rowTable[, descriptionColumns, drop = FALSE]
    rownames(annotationTable) <- NULL

    if (ncol(annotationTable) > 0) {
      descriptionTable <- cbind(descriptionTable, annotationTable)
    }

    if (format == "wide") {
      wideTable <- cbind(descriptionTable, as.data.frame(valueMatrix, optional = TRUE))
      colnames(wideTable) <- c(colnames(descriptionTable), sampleNames)
      rownames(wideTable) <- rowTable$row.key
      return(wideTable)
    }

    # One row per row and sample, the samples in the order of the object
    longTable <- descriptionTable[rep(seq_len(nrow(descriptionTable)), times = ncol(valueMatrix)), , drop = FALSE]
    longTable$sample <- rep(sampleNames, each = nrow(valueMatrix))
    longTable[[assayName]] <- as.vector(valueMatrix)
    rownames(longTable) <- NULL

    # A colData column named like a row column would otherwise be duplicated by the join
    collidingColumns <- setdiff(intersect(colnames(sampleTable), colnames(longTable)), "sample")
    if (length(collidingColumns) > 0) {
      colnames(sampleTable)[colnames(sampleTable) %in% collidingColumns] <- paste0(collidingColumns, ".sample")

      if (isTRUE(verbose)) {
        message("The following colData columns share a name with a row column and carry the suffix '.sample': ",
                paste(collidingColumns, collapse = ", "), ".")
      }
    }

    longTable <- dplyr::left_join(longTable, sampleTable, by = "sample")

    return(longTable)
  } # END function




#' @title .tileSummaryRule
#'
#' @description Decides how the tiles of a region are combined, from the argument when it is given and from the way the object was counted otherwise.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param tileSummary String with the rule requested, or \code{NULL}.
#'
#' @return A string, one among \code{"sum"}, \code{"mean"}, \code{"max"} and \code{"min"}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom S4Vectors metadata
#'
#' @keywords internal

.tileSummaryRule <-
  function(counts,
           tileSummary = NULL) {

    if (!is.null(tileSummary)) {
      if (!(tileSummary[1] %in% c("sum", "mean", "max", "min"))) {
        stop("The 'tileSummary' parameter must be one of 'sum', 'mean', 'max', 'min'.", call. = FALSE)
      }
      return(tileSummary[1])
    }

    # bigWig signal was summarised per tile, and the region needs the same summary
    if (identical(S4Vectors::metadata(counts)$signal.type, "bigwig")) {
      bigwigSummary <- counts@parameters$countBigwig$summaryFunction
      if (!is.null(bigwigSummary)) {
        return(bigwigSummary)
      }
    }

    return("sum")
  } # END function




#' @title .combineTileValues
#'
#' @description Combines the rows of a value matrix over the tiles of each region.
#'
#' @param valueMatrix Numeric matrix with one row per tile.
#' @param regionIndex Integer vector giving the region of every tile, numbered in order of first appearance.
#' @param tileWidths Numeric vector with the width of every tile.
#' @param summaryRule String, one among \code{"sum"}, \code{"mean"}, \code{"max"} and \code{"min"}.
#'
#' @return A numeric matrix with one row per region, in the order of \code{regionIndex}.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.combineTileValues <-
  function(valueMatrix,
           regionIndex,
           tileWidths,
           summaryRule = "sum") {

    if (summaryRule == "sum") {
      return(rowsum(x = valueMatrix, group = regionIndex, reorder = TRUE))
    }

    # Weighting by width lets the shorter trailing tile count for the bases it covers, not for a whole tile
    if (summaryRule == "mean") {
      weightedSums <- rowsum(x = valueMatrix * tileWidths, group = regionIndex, reorder = TRUE)
      return(weightedSums / as.numeric(rowsum(x = tileWidths, group = regionIndex, reorder = TRUE)))
    }

    # Maxima and minima have no rowsum counterpart, so they are taken one sample at a time
    summaryFunction <- if (summaryRule == "max") {max} else {min}
    regionFactor <- factor(regionIndex, levels = sort(unique(regionIndex)))

    combinedList <- lapply(seq_len(ncol(valueMatrix)),
                           function(j) {as.numeric(tapply(valueMatrix[, j], regionFactor, summaryFunction))})

    return(matrix(unlist(combinedList), nrow = nlevels(regionFactor), ncol = ncol(valueMatrix)))
  } # END function




#' @title .carriedCounts
#'
#' @description Returns the counts carried by a results object, and stops when the test left them behind.
#'
#' @param results \code{RegionSetDE.results}, \code{RegionSetDE.setResults} or either of the two list classes.
#'
#' @return A \code{RegionSetDE.counts} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.carriedCounts <-
  function(results) {

    carriedCounts <- resultCounts(results)

    if (ncol(carriedCounts) == 0) {
      stop("The result carries no counts, run the test with carryCounts = TRUE or call countTable() on the fit.", call. = FALSE)
    }

    return(carriedCounts)
  } # END function
