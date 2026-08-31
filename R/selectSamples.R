#' @title selectSamples
#'
#' @description Restricts a \code{RegionSetDE.counts} object to a subset of its samples. The conditions are written as they would be in \code{dplyr::filter} and are evaluated on the \code{colData}, so any column of the sample metadata can be used. The regions are left untouched, only the columns change.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param ... Conditions evaluated on the \code{colData}, in the syntax of \code{dplyr::filter}, e.g. \code{mark == "H3K27ac"}. Several conditions are combined with AND.
#' @param samples Character vector with the names of the samples to keep, or a numeric vector of column positions. Applied before the conditions. Default: \code{NULL}, all the samples.
#' @param dropNormalization Logical value to indicate whether the normalisation stored in the object must be discarded. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.counts} object with the selected columns.
#'
#' @details Marks, and more generally experiments run on different antibodies or different assays, should not share a model. The dispersion, the dynamic range and the meaning of the scaling factors all differ between them, so a fit that pools them borrows information across rows that have nothing to say about each other. The selection therefore belongs upstream of \code{\link{normalizeCounts}} rather than at test time, and this is why \code{dropNormalization} defaults to \code{TRUE}: factors estimated over a set of samples that no longer exists describe a library composition that no longer exists either. The raw counts are never modified, so re-normalising costs one call.
#'
#' The background bins stored in the metadata, when present, are subset along with the regions.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' brownNorway <- selectSamples(counts, condition == "BN", verbose = FALSE)
#' SummarizedExperiment::colData(brownNorway)
#'
#' firstReplicates <- selectSamples(counts, biologicalReplicate == "bio2",
#'                                  verbose = FALSE)
#' SummarizedExperiment::colData(firstReplicates)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{splitSamples}}, \code{\link{normalizeCounts}}, \code{\link{fitRegions}}
#'
#' @importFrom SummarizedExperiment colData assayNames assays assays<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom dplyr filter mutate
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export selectSamples

selectSamples <-
  function(counts,
           ...,
           samples = NULL,
           dropNormalization = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    if (ncol(counts) == 0) {
      stop("The 'counts' object contains no sample.", call. = FALSE)
    }

    #-------------------------------#
    # Resolve the columns to keep   #
    #-------------------------------#
    colTable <- as.data.frame(SummarizedExperiment::colData(counts))
    colTable <- dplyr::mutate(colTable,
                              sample.name = colnames(counts),
                              column.index = seq_len(ncol(counts)))

    if (!is.null(samples)) {
      if (is.character(samples)) {
        absentSamples <- setdiff(samples, colnames(counts))
        if (length(absentSamples) > 0) {
          stop(paste0("The following samples are absent from the object: ", paste(absentSamples, collapse = ", "), "."), call. = FALSE)
        }
        colTable <- dplyr::filter(colTable, .data$sample.name %in% samples)
      } else {
        samples <- as.integer(samples)
        if (any(is.na(samples)) | any(samples < 1) | any(samples > ncol(counts))) {
          stop("The 'samples' parameter contains positions outside the range of the columns.", call. = FALSE)
        }
        colTable <- dplyr::filter(colTable, .data$column.index %in% samples)
      }
    }

    # The conditions are evaluated in the caller environment, so external objects stay visible
    colTable <- dplyr::filter(colTable, ...)

    if (nrow(colTable) == 0) {
      stop("No sample satisfies the selection.", call. = FALSE)
    }

    #-------------------------------#
    # Subset the object             #
    #-------------------------------#
    selectedCounts <- counts[, colTable$column.index]

    # The bins carry the same columns as the regions and are used by the background normalisation
    backgroundBins <- S4Vectors::metadata(selectedCounts)$background
    if (!is.null(backgroundBins)) {
      S4Vectors::metadata(selectedCounts)$background <- backgroundBins[, colTable$column.index]
    }

    #-------------------------------#
    # Drop the stale normalisation  #
    #-------------------------------#
    normalizationInfo <- S4Vectors::metadata(selectedCounts)$normalization

    if (isTRUE(dropNormalization) & !is.null(normalizationInfo)) {
      remainingAssays <- setdiff(SummarizedExperiment::assayNames(selectedCounts),
                                 c(normalizationInfo$normalized.assay, "offset"))
      SummarizedExperiment::assays(selectedCounts) <- SummarizedExperiment::assays(selectedCounts)[remainingAssays]

      SummarizedExperiment::colData(selectedCounts)$norm.factor <- NULL
      SummarizedExperiment::colData(selectedCounts)$scaling.factor <- NULL
      S4Vectors::metadata(selectedCounts)$normalization <- NULL

      if (isTRUE(verbose)) {
        message(paste0("The normalisation computed with '", normalizationInfo$method,
                       "' has been removed, run normalizeCounts on the selected samples."))
      }
    } else if (!is.null(normalizationInfo) & isTRUE(verbose)) {
      warning("The scaling factors have been kept but were estimated on the whole set of samples.", call. = FALSE)
    }

    selectedCounts@parameters <- c(selectedCounts@parameters,
                                   list(selectSamples = list(samples = colTable$sample.name,
                                                             dropNormalization = dropNormalization)))

    if (isTRUE(verbose)) {
      message(paste0("Kept ", nrow(colTable), " samples out of ", ncol(counts), ": ",
                     paste(colTable$sample.name, collapse = ", "), "."))
    }

    return(selectedCounts)
  } # END function




#' @title splitSamples
#'
#' @description Splits a \code{RegionSetDE.counts} object into a list of objects, one per level of a column of the \code{colData}. Convenient when several marks or assays have been counted together over the same regions and each of them needs its own normalisation and its own fit.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param by String, or character vector, with the names of the \code{colData} columns defining the groups. When more than one is given the levels are combined.
#' @param dropNormalization Logical value to indicate whether the normalisation stored in the object must be discarded in each piece. Default: \code{TRUE}.
#' @param minSamples Numeric value with the minimum number of samples a group must contain to be returned. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A named list of \code{RegionSetDE.counts} objects.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' byStrain <- splitSamples(counts, by = "condition", verbose = FALSE)
#' names(byStrain)
#' byStrain[[1]]
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{selectSamples}}, \code{\link{fitRegions}}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter mutate
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export splitSamples

splitSamples <-
  function(counts,
           by,
           dropNormalization = TRUE,
           minSamples = 1,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    colTable <- as.data.frame(SummarizedExperiment::colData(counts))
    absentColumns <- setdiff(by, colnames(colTable))
    if (length(absentColumns) > 0) {
      stop(paste0("The following columns are absent from the colData: ", paste(absentColumns, collapse = ", "), "."), call. = FALSE)
    }

    #-------------------------------#
    # Build the groups              #
    #-------------------------------#
    colTable <- dplyr::mutate(colTable,
                              column.index = seq_len(ncol(counts)),
                              group.label = apply(colTable[, by, drop = FALSE], MARGIN = 1, FUN = paste, collapse = "_"))

    groupLabels <- unique(colTable$group.label)

    countsList <-
      lapply(groupLabels,
             function(groupLabel) {
               groupTable <- dplyr::filter(colTable, .data$group.label == groupLabel)

               if (nrow(groupTable) < minSamples) {
                 return(NULL)
               }

               return(selectSamples(counts = counts,
                                    samples = groupTable$column.index,
                                    dropNormalization = dropNormalization,
                                    verbose = FALSE))
             })

    names(countsList) <- groupLabels
    countsList <- countsList[!vapply(countsList, is.null, logical(1))]

    if (length(countsList) == 0) {
      stop("No group reaches 'minSamples'.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      groupSizes <- vapply(countsList, ncol, integer(1))
      message(paste0("Split into ", length(countsList), " groups: ",
                     paste(paste0(names(groupSizes), " (n=", groupSizes, ")"), collapse = ", "), "."))
    }

    return(countsList)
  } # END function
