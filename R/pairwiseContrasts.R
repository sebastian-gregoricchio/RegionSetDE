# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title pairwiseContrasts
#'
#' @description Writes out the contrasts between the levels of a column of the sample table, ready to be handed to \code{\link{testRegions}} or \code{\link{testRegionSets}}: every pair of levels, or every level against a reference.
#'
#' @param object \code{RegionSetDE.counts} or \code{RegionSetDE.fit} object, whose \code{colData} holds the column. A data.frame with the sample annotation is accepted as well.
#' @param column String with the name of the column.
#' @param levels Character vector with the levels to use, in the order that sets the direction of the contrasts: each level is compared against the ones before it. Default: \code{NULL}, the levels of the column when it is a factor, and the order in which they first appear otherwise.
#' @param reference String with a level to compare every other level against. When given, only those contrasts are written, instead of every pair. Default: \code{NULL}.
#'
#' @return A named list of character vectors \code{c(column, levelA, levelB)}, one per contrast, each read as \code{levelA} against \code{levelB}. The names are \code{levelA_vs_levelB}, and they become the names of the contrasts in the results.
#'
#' @details With \emph{k} levels there are \emph{k(k - 1) / 2} pairs, which grows quickly, and every pair is a separate family of tests. The correction for multiple testing is applied within each contrast, never across them, so the more contrasts are tested the more of the reported regions are false positives overall. A reference level keeps the number at \emph{k - 1}, and is usually what a treatment design is asking.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' pairwiseContrasts(counts, column = "condition")
#'
#' \dontrun{
#' sampleSheet <- loadExampleData("peakSheet")
#'
#' # Every pair of the three conditions
#' pairwiseContrasts(sampleSheet, column = "condition")
#'
#' # Every treatment against the vehicle only
#' pairwiseContrasts(sampleSheet, column = "condition", reference = "DMSO")
#'
#' results <- testRegions(fit, contrast = pairwiseContrasts(fit, column = "condition"))
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegions}}, \code{\link{testRegionSets}}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#' @importFrom utils combn
#'
#' @export pairwiseContrasts

pairwiseContrasts <-
  function(object,
           column,
           levels = NULL,
           reference = NULL) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (methods::is(object, "RegionSetDE.fit")) {
      object <- fitCounts(object)
    }

    sampleTable <-
      if (methods::is(object, "RegionSetDE.counts")) {
        as.data.frame(SummarizedExperiment::colData(object), optional = TRUE)
      } else if (is.data.frame(object)) {
        object
      } else {
        stop("The 'object' parameter must be a RegionSetDE.counts or a RegionSetDE.fit object, or a data.frame.", call. = FALSE)
      }

    if (!is.character(column) | length(column) != 1) {
      stop("The 'column' parameter must be a single column name.", call. = FALSE)
    }

    if (!(column %in% colnames(sampleTable))) {
      stop("The column '", column, "' is not in the sample table.", call. = FALSE)
    }

    #------------------------#
    # Levels and their order #
    #------------------------#
    # The order decides the sign of every fold change, so it is taken from the data rather than sorted
    columnValues <- sampleTable[[column]]
    presentLevels <- if (is.factor(columnValues)) {levels(droplevels(columnValues))} else {unique(as.character(columnValues[!is.na(columnValues)]))}

    if (is.null(levels)) {
      levels <- presentLevels
    }

    unknownLevels <- setdiff(levels, presentLevels)
    if (length(unknownLevels) > 0) {
      stop("The following levels are not in the column '", column, "': ", paste(unknownLevels, collapse = ", "), ".", call. = FALSE)
    }

    if (length(levels) < 2) {
      stop("At least two levels are needed to write a contrast.", call. = FALSE)
    }

    #------------------------#
    # Write the contrasts    #
    #------------------------#
    if (!is.null(reference)) {
      if (!(reference %in% levels)) {
        stop("The reference '", reference, "' is not among the levels: ", paste(levels, collapse = ", "), ".", call. = FALSE)
      }

      contrastList <- lapply(setdiff(levels, reference), function(level) {c(column, level, reference)})
    } else {
      # Each level against the ones before it, so the first one plays the reference in every pair it enters
      levelPairs <- utils::combn(seq_along(levels), 2)
      contrastList <- lapply(seq_len(ncol(levelPairs)),
                             function(i) {c(column, levels[levelPairs[2, i]], levels[levelPairs[1, i]])})
    }

    names(contrastList) <- vapply(contrastList, function(contrast) {paste(contrast[2], "vs", contrast[3], sep = "_")}, character(1))

    return(contrastList)
  } # END function
