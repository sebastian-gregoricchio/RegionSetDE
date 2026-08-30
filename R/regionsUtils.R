#' @title renameBedColumns
#'
#' @description Assigns the standard BED column names to the first columns of a data.frame, leaving the remaining ones untouched. Useful to convert a headerless table, read with generic \code{V1}, \code{V2}, \code{V3} names, into a data.frame ready for \code{GenomicRanges::makeGRangesFromDataFrame}.
#'
#' @param table A data.frame whose first columns hold the genomic coordinates.
#' @param bedFormat Numeric value indicating how many columns must be renamed, one among \code{3} (seqnames, start, end), \code{6} (adding name, score, strand), \code{9} (adding thickStart, thickEnd, itemRgb) or \code{12} (adding blockCount, blockSizes, blockStarts). Default: \code{3}.
#'
#' @return The input data.frame with the first columns renamed.
#'
#' @examples
#' \dontrun{
#' peaks <- utils::read.delim("peaks.bed", header = FALSE)
#' peaks <- renameBedColumns(peaks, bedFormat = 6)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @export renameBedColumns

renameBedColumns <-
  function(table,
           bedFormat = 3) {

    ### Check of the arguments
    if (!is.data.frame(table)) {
      stop("The 'table' parameter must be a data.frame.", call. = FALSE)
    }

    bedFormat <- bedFormat[1]
    if (!(bedFormat %in% c(3, 6, 9, 12))) {
      stop("The 'bedFormat' parameter must be one among 3, 6, 9 or 12.", call. = FALSE)
    }

    if (ncol(table) < bedFormat) {
      stop(paste0("The table contains ", ncol(table), " columns, less than the ", bedFormat, " required by the requested BED format."), call. = FALSE)
    }

    ### Renaming of the coordinate columns, the extra ones keep their original names
    bedColumns <- c("seqnames", "start", "end", "name", "score", "strand",
                    "thickStart", "thickEnd", "itemRgb",
                    "blockCount", "blockSizes", "blockStarts")

    colnames(table)[seq_len(bedFormat)] <- bedColumns[seq_len(bedFormat)]

    # Renaming a column onto a name already used further right would break any subsequent selection
    if (any(duplicated(colnames(table)))) {
      stop(paste0("The renaming produces duplicated column names: ",
                  paste(unique(colnames(table)[duplicated(colnames(table))]), collapse = ", "), "."), call. = FALSE)
    }

    return(table)
  } # END function




#' @title .readRegionFile
#'
#' @description Internal reader of BED-like and tabular files, shared by \code{loadRegions} and \code{splitLoadRegions}. Headerless files are treated as BED and their coordinates are converted from the 0-based to the 1-based convention.
#'
#' @param filePath String with the path to the file.
#' @param header Logical value to indicate whether the file carries a header line.
#' @param asGRanges Logical value to indicate whether the output must be a \code{GRanges} rather than a data.frame.
#'
#' @return A \code{GRanges} or a data.frame with the regions, or \code{NULL} when the file cannot be parsed.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter mutate %>%
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom utils read.delim
#'
#' @keywords internal
#' @noRd

.readRegionFile <-
  function(filePath,
           header = FALSE,
           asGRanges = TRUE) {

    if (!file.exists(filePath)) {
      stop(paste0("The file '", filePath, "' does not exist."), call. = FALSE)
    }

    bed <- utils::read.delim(filePath, header = header, comment.char = "#", stringsAsFactors = FALSE, quote = "")

    if (ncol(bed) < 3) {
      stop(paste0("The file '", filePath, "' contains less than 3 columns."), call. = FALSE)
    }

    ### Headerless files follow the BED specification, only the first 6 columns have a defined meaning
    if (header == FALSE) {
      bed <- renameBedColumns(bed, bedFormat = ifelse(ncol(bed) >= 6, 6, 3))

      # Drop the track/browser headers and any row without usable coordinates
      bed <-
        bed %>%
        dplyr::filter(!grepl("^(track|browser|#)", seqnames)) %>%
        dplyr::filter(grepl("^[0-9]+$", trimws(start)), grepl("^[0-9]+$", trimws(end))) %>%
        dplyr::mutate(seqnames = as.character(seqnames),
                      start = as.numeric(start) + 1, # BED files are 0-based, GRanges are 1-based
                      end = as.numeric(end))

      # narrowPeak/broadPeak often store a dot in the strand column
      if ("strand" %in% colnames(bed)) {
        bed <- dplyr::mutate(bed, strand = ifelse(strand %in% c("+", "-", "*"), strand, "*"))
      }
    }

    if (nrow(bed) == 0) {
      stop(paste0("No valid region found in '", filePath, "'."), call. = FALSE)
    }

    if (asGRanges == FALSE) {return(bed)}

    return(GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE))
  } # END function
