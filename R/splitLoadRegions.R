#' @title splitLoadRegions
#'
#' @description Imports a collection of genomic regions stored in a single object or file and splits them into individual region sets according to the values of one of its columns. The regions can be provided as a \code{GRanges}, a data.frame or the path to a BED-like/tabular file. All the arguments controlling the import are passed to \code{\link{loadRegions}}.
#'
#' @param regions A \code{GRanges} object, a data.frame with chromosome/start/end columns, or a string indicating the path to a BED-like or tabular file containing all the regions.
#' @param splitBy String with the name of the column/metadata field containing the region set labels, or numeric value with its position in the input. The chromosome column (\code{"seqnames"}) is accepted as well, to split the regions by chromosome. Default: \code{"name"}, the fourth column of a BED file.
#' @param header Logical value to indicate whether the file contains a header line. Used only when \code{regions} is a file path. Notice that headerless files are read as BED and their coordinates are converted from 0-based to 1-based, while data.frames and \code{GRanges} are assumed to be 1-based already. Default: \code{FALSE}.
#' @param selectedSets Character vector with the labels of the sets to import, the others being discarded. Default: \code{NULL}, all the sets are imported.
#' @param minRegionsPerSet Numeric value indicating the minimum number of regions that a set must contain to be imported. Sets falling below this threshold are dropped with a warning. Default: \code{1}.
#' @param maxSets Numeric value indicating the maximum number of sets tolerated. Splitting on a continuous column, such as the score, would generate thousands of single-region sets, therefore the loading is interrupted above this threshold. Default: \code{100}.
#' @param keepSplitColumn Logical value to indicate whether the splitting column must be kept in the metadata. Since it is constant within each set, it is redundant with the set name. Default: \code{FALSE}.
#' @param keepMetadata Logical value to indicate whether the metadata columns must be kept. Default: \code{TRUE}.
#' @param sortRegions Logical value to indicate whether the regions must be sorted by coordinate. Default: \code{TRUE}.
#' @param reduceRegions Logical value to indicate whether the overlapping regions within the same set must be collapsed (\code{IRanges::reduce}). Notice that metadata are lost upon reduction. Default: \code{FALSE}.
#' @param removeDuplicatedRegions Logical value to indicate whether the regions with identical coordinates within the same set must be collapsed to a single entry. Default: \code{TRUE}.
#' @param duplicatedSets String indicating how to handle two or more sets containing exactly the same regions: \code{"stop"}, \code{"remove"} or \code{"keep"}. Default: \code{"stop"}.
#' @param seqlevelsStyle String indicating the chromosome naming style to apply, one among \code{"UCSC"} (chr1), \code{"Ensembl"} (1) or \code{"NCBI"}. Default: \code{"UCSC"}.
#' @param genomeAssembly String indicating the genome assembly to store in the ranges seqinfo, e.g. \code{"hg38"}. Default: \code{NULL}, no assembly is assigned.
#' @param outputFormat String indicating the class of the returned object, one among \code{"RegionSetDE"}, \code{"GRangesList"} or \code{"list"}. Default: \code{"RegionSetDE"}.
#' @param verbose Logical value to indicate whether the loading messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE} object or, depending on \code{outputFormat}, a named \code{GRangesList} or a named list of \code{GRanges}.
#'
#' @examples
#' \dontrun{
#' regions <- splitLoadRegions("peaks/all_peaks_annotated.bed",
#'                             splitBy = "name",
#'                             genomeAssembly = "hg38")
#'
#' regions <- splitLoadRegions(grAllPeaks,
#'                             splitBy = "cluster",
#'                             minRegionsPerSet = 50)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadRegions}}
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom methods is
#'
#' @export splitLoadRegions

splitLoadRegions <-
  function(regions,
           splitBy = "name",
           header = FALSE,
           selectedSets = NULL,
           minRegionsPerSet = 1,
           maxSets = 100,
           keepSplitColumn = FALSE,
           keepMetadata = TRUE,
           sortRegions = TRUE,
           reduceRegions = FALSE,
           removeDuplicatedRegions = TRUE,
           duplicatedSets = "stop",
           seqlevelsStyle = "UCSC",
           genomeAssembly = NULL,
           outputFormat = "RegionSetDE",
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (length(splitBy) != 1) {
      stop("The 'splitBy' parameter must be a single column name or a single column position.", call. = FALSE)
    }

    if (minRegionsPerSet < 1) {
      stop("The 'minRegionsPerSet' parameter must be at least 1.", call. = FALSE)
    }

    #-------------------------------------------#
    # Collect the regions and their set labels  #
    #-------------------------------------------#
    # The labels are read before any coercion, so that the splitting column can also be a coordinate column
    if (is.character(regions) & length(regions) == 1) {
      ### File input, kept as a table so that the column positions match the ones of the file
      regionTable <- .readRegionFile(filePath = regions, header = header, asGRanges = FALSE)

      if (is.numeric(splitBy)) {
        if (splitBy > ncol(regionTable)) {
          stop(paste0("The 'splitBy' position (", splitBy, ") exceeds the number of columns of the file (", ncol(regionTable), ")."), call. = FALSE)
        }
        splitBy <- colnames(regionTable)[splitBy]
      }

      if (!(splitBy %in% colnames(regionTable))) {
        stop(paste0("The column '", splitBy, "' is not available in the file. The columns found are: ",
                    paste(colnames(regionTable), collapse = ", "), "."), call. = FALSE)
      }

      setLabels <- as.character(regionTable[[splitBy]])
      allRegions <- GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE)

    } else if (is.data.frame(regions)) {
      ### Data.frame input
      if (is.numeric(splitBy)) {
        if (splitBy > ncol(regions)) {
          stop(paste0("The 'splitBy' position (", splitBy, ") exceeds the number of columns of the data.frame (", ncol(regions), ")."), call. = FALSE)
        }
        splitBy <- colnames(regions)[splitBy]
      }

      if (!(splitBy %in% colnames(regions))) {
        stop(paste0("The column '", splitBy, "' is not available in the data.frame. The columns found are: ",
                    paste(colnames(regions), collapse = ", "), "."), call. = FALSE)
      }

      setLabels <- as.character(regions[[splitBy]])
      allRegions <- GenomicRanges::makeGRangesFromDataFrame(regions, keep.extra.columns = TRUE)

    } else if (methods::is(regions, "GRanges")) {
      ### GRanges input, the labels are searched among the metadata columns
      metadataColumns <- colnames(S4Vectors::mcols(regions))

      if (is.numeric(splitBy)) {
        if (splitBy > length(metadataColumns)) {
          stop(paste0("The 'splitBy' position (", splitBy, ") exceeds the number of metadata columns of the GRanges (", length(metadataColumns), ")."), call. = FALSE)
        }
        splitBy <- metadataColumns[splitBy]
      }

      if (splitBy %in% metadataColumns) {
        setLabels <- as.character(S4Vectors::mcols(regions)[[splitBy]])
      } else if (splitBy == "seqnames") {
        setLabels <- as.character(GenomeInfoDb::seqnames(regions))
      } else {
        stop(paste0("The column '", splitBy, "' is not available in the GRanges. The metadata columns found are: ",
                    ifelse(length(metadataColumns) == 0, "none", paste(metadataColumns, collapse = ", ")), "."), call. = FALSE)
      }

      allRegions <- regions

    } else {
      stop("The 'regions' parameter must be a single GRanges, a data.frame or the path to a single file. Use 'loadRegions' to import a collection of separate sets.", call. = FALSE)
    }

    #-----------------------------#
    # Discard the unlabelled rows #
    #-----------------------------#
    # Regions without a label belong to no set, dropping them silently would hide a truncated import
    unlabelled <- which(is.na(setLabels) | setLabels == "")
    if (length(unlabelled) > 0) {
      warning(paste0(length(unlabelled), " regions carry no value in the column '", splitBy, "' and have been discarded."), call. = FALSE)
      allRegions <- allRegions[-unlabelled]
      setLabels <- setLabels[-unlabelled]
    }

    if (length(allRegions) == 0) {
      stop(paste0("No region left after discarding the entries without a value in the column '", splitBy, "'."), call. = FALSE)
    }

    # Splitting on a continuous column would generate one set per region
    if (length(unique(setLabels)) > maxSets) {
      stop(paste0("The column '", splitBy, "' contains ", length(unique(setLabels)), " distinct values, above the 'maxSets' threshold (", maxSets,
                  "). Please check that the splitting column holds the region set names."), call. = FALSE)
    }

    #-------------------#
    # Split the regions #
    #-------------------#
    # The splitting column is constant within each set, therefore redundant with the set name
    if (keepSplitColumn == FALSE & splitBy %in% colnames(S4Vectors::mcols(allRegions))) {
      S4Vectors::mcols(allRegions)[[splitBy]] <- NULL
    }

    # Index subsetting keeps the metadata untouched and follows the order of appearance of the labels
    setLevels <- unique(setLabels)
    regionSets <- lapply(X = setLevels, FUN = function(x) {return(allRegions[setLabels == x])})
    names(regionSets) <- setLevels

    ### Selection of the sets requested by the user
    if (!is.null(selectedSets)) {
      missingSets <- setdiff(selectedSets, names(regionSets))
      if (length(missingSets) > 0) {
        stop(paste0("The following sets are not present in the column '", splitBy, "': ", paste(missingSets, collapse = ", "),
                    ". The sets available are: ", paste(names(regionSets), collapse = ", "), "."), call. = FALSE)
      }
      regionSets <- regionSets[selectedSets]
    }

    ### Removal of the underpowered sets
    setSizes <- vapply(regionSets, length, numeric(1))
    if (any(setSizes < minRegionsPerSet)) {
      smallSets <- names(regionSets)[setSizes < minRegionsPerSet]
      warning(paste0("Region sets discarded because containing less than ", minRegionsPerSet, " regions: ",
                     paste0(smallSets, " (", setSizes[smallSets], ")", collapse = ", "), "."), call. = FALSE)
      regionSets <- regionSets[setSizes >= minRegionsPerSet]
    }

    if (length(regionSets) == 0) {
      stop("No region set left after the filtering, please lower the 'minRegionsPerSet' threshold.", call. = FALSE)
    }

    if (verbose == TRUE) {
      message(paste0("Regions split by '", splitBy, "' into ", length(regionSets), " sets."))
    }

    #-----------------------#
    # Import the split sets #
    #-----------------------#
    # The import is delegated to loadRegions, which handles the chromosome names, the duplicates and the output class
    return(loadRegions(regions = regionSets,
                       keepMetadata = keepMetadata,
                       sortRegions = sortRegions,
                       reduceRegions = reduceRegions,
                       removeDuplicatedRegions = removeDuplicatedRegions,
                       duplicatedSets = duplicatedSets,
                       seqlevelsStyle = seqlevelsStyle,
                       genomeAssembly = genomeAssembly,
                       skipInvalid = FALSE,
                       outputFormat = outputFormat,
                       verbose = verbose))
  } # END function
