# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title contrastInfo
#'
#' @description Describes the contrasts of a results object in one table: the engine that tested them, the two groups compared and how many samples each holds, and the distribution the test statistic follows with its degrees of freedom. Together with the \code{stat} column of \code{\link{resultsTable}} it is what a power or sample size analysis needs, for instance with power4peaks.
#'
#' @param results \code{RegionSetDE.results} object, or a \code{RegionSetDE.resultsList} holding several contrasts.
#'
#' @return A data.frame with one row per contrast:
#' \itemize{
#'   \item \code{contrast}: the name of the contrast in the list, or its description for a single one.
#'   \item \code{contrast.description}: what the contrast compares.
#'   \item \code{engine}: \code{"edgeR"}, \code{"voom"}, \code{"limma"}, \code{"dream"} or \code{"deseq2"}.
#'   \item \code{column}, \code{group1}, \code{group2}: the column of the sample table the contrast separates and its two levels, the fold change being \code{group1} over \code{group2}.
#'   \item \code{n.group1}, \code{n.group2}: the number of samples in each group.
#'   \item \code{stat.distribution}: the distribution of \code{stat} under the null hypothesis, \code{"f"}, \code{"chisq"}, \code{"t"} or \code{"norm"}.
#'   \item \code{df1}, \code{df2}: its degrees of freedom, the median over the regions when they differ between regions, as the prior degrees of freedom of a robust fit make them do.
#'   \item \code{n.regions}: the number of rows tested.
#' }
#' The group columns are \code{NA} when the contrast is not a difference between two levels of one column, a numeric vector over the coefficients for instance, and \code{stat.distribution} is \code{NA} for the threshold tests run with \code{lfcThreshold > 0}.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' contrastInfo(results)
#'
#' # The statistics themselves, one per region
#' head(resultsTable(results)$stat)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegions}}, \code{\link{resultsTable}}, \code{\link{contrastName}}
#'
#' @importFrom stats median
#' @importFrom dplyr bind_rows
#' @importFrom methods is
#'
#' @export contrastInfo

contrastInfo <-
  function(results) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (methods::is(results, "RegionSetDE.resultsList")) {
      infoList <- lapply(names(results@results),
                         function(contrastKey) {.singleContrastInfo(results = results@results[[contrastKey]], contrastKey = contrastKey)})
      return(dplyr::bind_rows(infoList))
    }

    if (!methods::is(results, "RegionSetDE.results")) {
      stop("The 'results' parameter must be a RegionSetDE.results or a RegionSetDE.resultsList object.", call. = FALSE)
    }

    return(.singleContrastInfo(results = results, contrastKey = results@contrast))
  } # END function




#' @title .singleContrastInfo
#'
#' @description Builds the row of \code{\link{contrastInfo}} for one contrast.
#'
#' @param results \code{RegionSetDE.results} object.
#' @param contrastKey String with the name the contrast goes under.
#'
#' @return A data.frame with one row.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats median na.omit
#'
#' @keywords internal

.singleContrastInfo <-
  function(results,
           contrastKey) {

    resultTable <- results@results
    contrastGroups <- results@contrast.groups

    # A value per region collapses to one, the way a power analysis takes it
    summariseColumn <- function(columnName) {
      if (!(columnName %in% colnames(resultTable)) || all(is.na(resultTable[[columnName]]))) {
        return(NA_real_)
      }
      stats::median(resultTable[[columnName]], na.rm = TRUE)
    }

    statDistribution <-
      if ("stat.distribution" %in% colnames(resultTable)) {
        unique(stats::na.omit(as.character(resultTable$stat.distribution)))
      } else {
        character(0)
      }

    hasGroups <- length(contrastGroups$groups) == 2
    groupSizes <- if (length(contrastGroups$n.samples) == 2) {as.integer(contrastGroups$n.samples)} else {c(NA_integer_, NA_integer_)}

    return(data.frame(contrast = contrastKey,
                      contrast.description = results@contrast,
                      engine = results@engine,
                      column = if (hasGroups) {contrastGroups$column} else {NA_character_},
                      group1 = if (hasGroups) {contrastGroups$groups[1]} else {NA_character_},
                      group2 = if (hasGroups) {contrastGroups$groups[2]} else {NA_character_},
                      n.group1 = groupSizes[1],
                      n.group2 = groupSizes[2],
                      stat.distribution = if (length(statDistribution) == 1) {statDistribution} else {NA_character_},
                      df1 = summariseColumn("df1"),
                      df2 = summariseColumn("df2"),
                      n.regions = nrow(resultTable),
                      stringsAsFactors = FALSE))
  } # END function
