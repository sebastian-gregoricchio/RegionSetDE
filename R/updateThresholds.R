# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title updateThresholds
#'
#' @description Labels the regions of a result again with other cut-offs, without running the test a second time. The \code{diff.status} column is filled anew from the \code{FDR} and \code{log2FC} given, and the thresholds stored in the object are replaced, so that the tables, the plots and the export all follow the new ones.
#'
#' @param results \code{RegionSetDE.results} or \code{RegionSetDE.resultsList} object returned by \code{\link{testRegions}}, or \code{RegionSetDE.setResults} or \code{RegionSetDE.setResultsList} object returned by \code{\link{testRegionSets}} or \code{\link{testSetContrast}}.
#' @param FDR Numeric value with the new adjusted p-value cut-off. Default: \code{NULL}, the one stored in the object.
#' @param log2FC Numeric value with the new absolute log2 fold change cut-off. Not available for the results of the region sets, which carry no label per region. Default: \code{NULL}, the one stored in the object.
#' @param contrast Character vector with the names or the positions of the contrasts to update, when \code{results} holds several of them. Default: \code{NULL}, all of them.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return An object of the same class as \code{results}, with the new labels in \code{diff.status} and the new cut-offs in the \code{thresholds} slot and in the parameters recorded for the test. The previous cut-offs are kept in \code{parameters$updateThresholds}.
#'
#' @details The cut-offs only decide the labels. The p-values, the adjusted p-values and the fold changes come from the model and do not move, which is why nothing has to be fitted or tested again. What does need a new test is \code{lfcThreshold} of \code{\link{testRegions}}, which moves the fold change inside the test itself and changes the p-values; it stays as it was and is reported next to the new cut-offs.
#'
#' Everything reading the labels afterwards sees the new ones: \code{\link{resultsTable}}, \code{\link{topRegions}}, \code{\link{plotVolcano}} and \code{\link{plotResultsMA}}, which draw the dashed lines at the stored cut-offs, \code{\link{peakOccupancyTable}}, \code{\link{plotProfile}} and \code{\link{exportResults}}, which writes the cut-offs down with the rest of the parameters. On tiled results the region labels are updated, while \code{n.tiles.up} and \code{n.tiles.down} keep counting the tiles changing within the region at the fixed rate csaw uses, since they come from the combination of the tiles and not from the labels.
#'
#' For the results of the region sets only \code{FDR} applies: it decides which sets \code{\link{plotSetEffect}} and the printed summary call significant.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' table(resultsTable(results)$diff.status)
#'
#' # Stricter cut-offs, with no new test
#' strictResults <- updateThresholds(results, FDR = 0.01, log2FC = 1)
#' table(resultsTable(strictResults)$diff.status)
#'
#' # The cut-offs travel with the object, the plots draw them
#' plotVolcano(strictResults)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegions}}, \code{\link{testRegionSets}}, \code{\link{resultsTable}}
#'
#' @importFrom methods is
#'
#' @export updateThresholds

updateThresholds <-
  function(results,
           FDR = NULL,
           log2FC = NULL,
           contrast = NULL,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!is.null(FDR) && (!is.numeric(FDR) | length(FDR) != 1 || is.na(FDR) || FDR <= 0 || FDR > 1)) {
      stop("The 'FDR' parameter must be a single number above 0 and not above 1.", call. = FALSE)
    }

    if (!is.null(log2FC) && (!is.numeric(log2FC) | length(log2FC) != 1 || is.na(log2FC) || log2FC < 0)) {
      stop("The 'log2FC' parameter must be a single number of at least 0.", call. = FALSE)
    }

    if (is.null(FDR) & is.null(log2FC)) {
      stop("Give at least one new cut-off, 'FDR' or 'log2FC'.", call. = FALSE)
    }

    #------------------------#
    # One contrast or many   #
    #------------------------#
    if (methods::is(results, "RegionSetDE.results") | methods::is(results, "RegionSetDE.setResults")) {
      if (!is.null(contrast)) {
        stop("The object holds a single contrast, leave 'contrast' empty.", call. = FALSE)
      }
      return(.relabelContrast(results = results, FDR = FDR, log2FC = log2FC, verbose = verbose))
    }

    if (methods::is(results, "RegionSetDE.resultsList") | methods::is(results, "RegionSetDE.setResultsList")) {
      contrastNames <- names(results@results)

      # The contrasts are named as testRegions named them, or taken by position
      contrastIndex <- if (is.null(contrast)) {
        seq_along(contrastNames)
      } else if (is.numeric(contrast)) {
        if (any(is.na(contrast) | contrast < 1 | contrast > length(contrastNames))) {
          stop("The 'contrast' parameter contains positions outside the range of the contrasts.", call. = FALSE)
        }
        as.integer(contrast)
      } else {
        absentContrasts <- setdiff(contrast, contrastNames)
        if (length(absentContrasts) > 0) {
          stop("The following contrasts are absent from the object: ", paste(absentContrasts, collapse = ", "),
               ". The contrasts are: ", paste(contrastNames, collapse = ", "), ".", call. = FALSE)
        }
        match(contrast, contrastNames)
      }

      for (i in contrastIndex) {
        if (isTRUE(verbose)) {message(contrastNames[i], ":")}
        results@results[[i]] <- .relabelContrast(results = results@results[[i]], FDR = FDR, log2FC = log2FC, verbose = verbose)
      }

      return(results)
    }

    stop("The 'results' parameter must be an object returned by testRegions(), testRegionSets() or testSetContrast().", call. = FALSE)
  } # END function




#' @title .relabelContrast
#'
#' @description Applies new cut-offs to the result of a single contrast: new labels for the regions, new thresholds in the object and in its parameters.
#'
#' @param results \code{RegionSetDE.results} or \code{RegionSetDE.setResults} object.
#' @param FDR Numeric value with the new adjusted p-value cut-off, or \code{NULL} to keep the stored one.
#' @param log2FC Numeric value with the new absolute log2 fold change cut-off, or \code{NULL} to keep the stored one.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return The object with the new labels and thresholds.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr mutate case_when
#' @importFrom rlang .data
#' @importFrom methods is validObject
#'
#' @keywords internal

.relabelContrast <-
  function(results,
           FDR,
           log2FC,
           verbose) {

    previousThresholds <- results@thresholds
    newFDR <- if (is.null(FDR)) {previousThresholds$FDR} else {FDR}

    #-------------------------------#
    # Region sets: the FDR alone    #
    #-------------------------------#
    if (methods::is(results, "RegionSetDE.setResults")) {
      if (!is.null(log2FC)) {
        stop("The results of the region sets carry no label per region, only 'FDR' can be updated.", call. = FALSE)
      }

      results@thresholds$FDR <- newFDR
      testStep <- intersect(c("testRegionSets", "testSetContrast"), names(results@parameters))
      for (stepName in testStep) {
        results@parameters[[stepName]]$FDR <- newFDR
      }

      results@parameters$updateThresholds <- list(previous.FDR = previousThresholds$FDR, FDR = newFDR)

      if (isTRUE(verbose)) {
        significanceColumn <- intersect(c("camera.FDR", "fry.FDR"), colnames(results@results))[1]
        message("FDR cut-off moved from ", previousThresholds$FDR, " to ", newFDR, ": ",
                sum(results@results[[significanceColumn]] < newFDR, na.rm = TRUE), " of ", nrow(results@results),
                " sets below it (", significanceColumn, ").")
      }

      methods::validObject(results)
      return(results)
    }

    #-------------------------------#
    # Regions: new labels           #
    #-------------------------------#
    newLog2FC <- if (is.null(log2FC)) {previousThresholds$log2FC} else {log2FC}

    # The arguments share their names with two columns of the table, which dplyr would read first
    FDRthreshold <- newFDR
    log2FCthreshold <- newLog2FC

    results@results <- dplyr::mutate(results@results,
                                     diff.status = dplyr::case_when(.data$FDR < FDRthreshold & .data$log2FC > log2FCthreshold ~ "up",
                                                                    .data$FDR < FDRthreshold & .data$log2FC < (-log2FCthreshold) ~ "down",
                                                                    TRUE ~ "null"))
    results@results$diff.status <- factor(results@results$diff.status, levels = c("down", "null", "up"))

    #-------------------------------#
    # Thresholds and parameters     #
    #-------------------------------#
    results@thresholds$FDR <- newFDR
    results@thresholds$log2FC <- newLog2FC

    if (!is.null(results@parameters$testRegions)) {
      results@parameters$testRegions$FDR <- newFDR
      results@parameters$testRegions$log2FC <- newLog2FC
    }

    results@parameters$updateThresholds <- list(previous.FDR = previousThresholds$FDR,
                                                previous.log2FC = previousThresholds$log2FC,
                                                FDR = newFDR,
                                                log2FC = newLog2FC,
                                                lfcThreshold = previousThresholds$lfcThreshold)

    if (isTRUE(verbose)) {
      statusTable <- table(results@results$diff.status)
      message("Relabelled with FDR < ", newFDR, " and |log2FC| > ", newLog2FC, ": ",
              statusTable[["up"]], " up and ", statusTable[["down"]], " down out of ", nrow(results@results), " regions.")

      # The threshold inside the test shaped the p-values, relabelling cannot move it
      if (isTRUE(previousThresholds$lfcThreshold > 0)) {
        message("The p-values were computed against lfcThreshold = ", previousThresholds$lfcThreshold, ", which only a new test can change.")
      }
    }

    methods::validObject(results)
    return(results)
  } # END function
