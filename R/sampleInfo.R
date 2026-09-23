# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title sampleInfo
#'
#' @description Returns the sample table of an object as a data.frame: one row per sample, with the annotation of the sample sheet and what the package added along the way, such as the library sizes, the layout of each library and, once the counts are normalised, the scaling factors. It is \code{SummarizedExperiment::colData()} in a form \code{dplyr} and \code{ggplot2} take directly.
#'
#' @param object \code{RegionSetDE.counts} or \code{RegionSetDE.fit} object, or any of the results classes, whose carried counts are used.
#' @param columns Character vector with the columns to keep, in that order. Default: \code{NULL}, all of them.
#'
#' @return A data.frame with one row per sample.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' sampleInfo(counts)
#'
#' counts <- normalizeCounts(counts, method = "TMM", verbose = FALSE)
#' sampleInfo(counts, columns = c("sample", "condition", "library.size", "scaling.factor"))
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{libInfo}}, \code{\link{contrastInfo}}, \code{\link{countTable}}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#'
#' @export sampleInfo

sampleInfo <-
  function(object,
           columns = NULL) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    # The results reach their samples through the counts they carry
    if (methods::is(object, "RegionSetDE.fit")) {
      object <- fitCounts(object)
    } else if (any(vapply(c("RegionSetDE.results", "RegionSetDE.resultsList", "RegionSetDE.setResults", "RegionSetDE.setResultsList"),
                          function(className) {methods::is(object, className)}, logical(1)))) {
      object <- .carriedCounts(object)
    }

    if (!methods::is(object, "RegionSetDE.counts")) {
      stop("The 'object' parameter must be a RegionSetDE.counts, a RegionSetDE.fit or a results object.", call. = FALSE)
    }

    sampleTable <- as.data.frame(SummarizedExperiment::colData(object), optional = TRUE)

    #------------------------#
    # Columns to keep        #
    #------------------------#
    if (!is.null(columns)) {
      missingColumns <- setdiff(columns, colnames(sampleTable))

      if (length(missingColumns) > 0) {
        stop("The following columns are not in the sample table: ", paste(missingColumns, collapse = ", "),
             ". Available: ", paste(colnames(sampleTable), collapse = ", "), ".", call. = FALSE)
      }

      sampleTable <- sampleTable[, columns, drop = FALSE]
    }

    # The sample names are a column already, as row names they would only be dropped by dplyr
    if ("sample" %in% colnames(SummarizedExperiment::colData(object))) {
      rownames(sampleTable) <- NULL
    }

    return(sampleTable)
  } # END function
