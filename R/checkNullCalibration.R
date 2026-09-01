#' @title checkNullCalibration
#'
#' @description Runs the contrast of a fit on rows that should not respond to it, and reports how many of them come out significant anyway. On a design with no replicates this is the only empirical check there is on the dispersion that was supplied, and it should be run before any of the results are believed.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrast Contrast to test, in the syntax accepted by \code{\link{testRegions}}.
#' @param source String with where the null rows come from, one of \code{"background"}, \code{"regionSet"} and \code{"supplied"}. Default: \code{NULL}, the source the dispersion of the fit came from.
#' @param regionSets Character vector with the names of the sets used as null rows. Only for \code{source = "regionSet"}. Default: \code{NULL}.
#' @param index Integer vector with the positions of the null rows. Only for \code{source = "supplied"}. Default: \code{NULL}.
#' @param subset Integer vector restricting the check to a subset of the rows \code{source} selects. Default: \code{NULL}, the rows the fit held out of its estimate.
#' @param minCount Numeric value with the average count a null row must carry to be used. Default: \code{10}.
#' @param maxRows Numeric value with the number of null rows kept. Default: \code{50000}.
#' @param FDR Numeric value with the adjusted p-value cut-off the count is reported at. Default: \code{0.05}.
#' @param pValue Numeric value with the raw p-value cut-off the count is reported at. Default: \code{0.05}.
#' @param strata Numeric value with the number of abundance strata the calibration is broken down over. Default: \code{5}.
#' @param suggest Logical value to indicate whether a dispersion that would calibrate the null must be searched for when the current one does not. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with the p-values of the null rows, the proportions passing the two cut-offs, the same broken down by abundance, the dispersion the fit used, a suggested one when it is off, and the verdict.
#'
#' @details The logic is the same one that made the dispersion estimable in the first place. Rows assumed not to respond are fitted with the same design, the same offsets, the same dispersion and the same test as the regions, and their p-values are read against what they should look like if nothing is happening: uniform between zero and one, with a fraction \code{pValue} of them below the raw cut-off and essentially none surviving the correction.
#' Too many small p-values means the dispersion is too low, the fit is treating library-to-library variation as treatment effect, and the regions are inheriting the same optimism at the same rate. Too few means the opposite and costs power rather than credibility. The ratio of what was seen to what was expected is reported as the inflation, and one is what a calibrated fit looks like.
#' Checking a dispersion on the rows it was fitted to would be calibrated by construction, so \code{\link{estimateNullDispersion}} holds half of them back and this function uses that half by default. Passing \code{subset} yourself overrides it, and a check run against a region set believed to be invariant, with the dispersion estimated on the genome, is stronger still.
#' The breakdown by abundance is where the answer usually is when the check fails. Overdispersion concentrated in the weakest rows means the filter was too loose and \code{\link{filterRegions}} is the fix; raising the dispersion globally would instead cost real signal in the abundant rows to pay for noise in the faint ones.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' # The background bins carry no strain effect, so the p-values should be flat
#' calibration <- checkNullCalibration(fit,
#'                                     contrast = c("condition", "SHR", "BN"),
#'                                     source = "background",
#'                                     verbose = FALSE)
#' calibration
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{estimateNullDispersion}}, \code{\link{fitRegions}}, \code{\link{plotNullCalibration}}, \code{\link{filterRegions}}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom edgeR DGEList scaleOffset
#' @importFrom stats p.adjust median quantile
#' @importFrom methods is
#'
#' @export checkNullCalibration

checkNullCalibration <-
  function(fit,
           contrast,
           source = NULL,
           regionSets = NULL,
           index = NULL,
           subset = NULL,
           minCount = 10,
           maxRows = 50000,
           FDR = 0.05,
           pValue = 0.05,
           strata = 5,
           suggest = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(fit, "RegionSetDE.fit")) {
      stop("The 'fit' parameter must be a RegionSetDE.fit object.", call. = FALSE)
    }

    if (is.null(fit@dispersion$common)) {
      stop("The fit carries no common dispersion, so there is nothing to calibrate.", call. = FALSE)
    }

    #-------------------------------#
    # Rows the fit has not seen     #
    #-------------------------------#
    fitSource <- if (is.null(fit@dispersion$source)) {"background"} else {fit@dispersion$source}
    source <- if (is.null(source)) {fitSource} else {source}

    usedHoldout <- FALSE
    if (is.null(subset) & identical(source, fitSource) & length(fit@dispersion$holdout.index) > 0) {
      subset <- fit@dispersion$holdout.index
      usedHoldout <- TRUE
    }

    contrastObject <- .resolveContrast(contrast = contrast,
                                       design = fit@design,
                                       colData = SummarizedExperiment::colData(fit@counts))

    nullObject <- .nullMatrix(counts = fit@counts, source = source, regionSets = regionSets,
                              index = index, subset = subset, minCount = minCount, maxRows = maxRows)

    if (nrow(nullObject$counts) < 100) {
      stop("Only ", nrow(nullObject$counts), " null rows reach ", minCount,
           " counts, which is too few to calibrate against.", call. = FALSE)
    }

    #-------------------------------#
    # The same model on the null    #
    #-------------------------------#
    # Same design, same offsets, same dispersion and the same test: only the rows change
    dgeList <- edgeR::DGEList(counts = nullObject$counts, lib.size = nullObject$library.size)
    dgeList <- edgeR::scaleOffset(y = dgeList, offset = nullObject$offset)

    useQuasiLikelihood <- !identical(fit@fit$test, "lrt")

    nullTable <- .nullTest(dgeList = dgeList,
                           design = fit@design,
                           contrastVector = contrastObject$vector,
                           dispersion = fit@dispersion$common,
                           useQuasiLikelihood = useQuasiLikelihood)

    nullPvalues <- nullTable$PValue
    nullFDR <- stats::p.adjust(nullPvalues, method = "BH")

    #-------------------------------#
    # Where the miscalibration is   #
    #-------------------------------#
    rawProportion <- mean(nullPvalues < pValue, na.rm = TRUE)

    abundanceStratum <- .widthStratum(rowWidths = nullObject$abundance, strataNumber = strata)

    strataTable <-
      do.call(what = rbind,
              args = lapply(sort(unique(abundanceStratum)),
                            function(stratumName) {
                              stratumRows <- abundanceStratum == stratumName
                              return(data.frame(stratum = stratumName,
                                                n.rows = sum(stratumRows),
                                                median.abundance = stats::median(nullObject$abundance[stratumRows]),
                                                proportion = mean(nullPvalues[stratumRows] < pValue, na.rm = TRUE),
                                                stringsAsFactors = FALSE))
                            }))
    strataTable$inflation <- strataTable$proportion / pValue

    verdict <- if (rawProportion > 2 * pValue) {
      "too liberal"
    } else if (rawProportion < 0.5 * pValue) {
      "too conservative"
    } else {
      "calibrated"
    }

    #-------------------------------#
    # What would calibrate it       #
    #-------------------------------#
    suggestedDispersion <- NA_real_
    if (isTRUE(suggest) & verdict != "calibrated") {
      suggestedDispersion <- .suggestDispersion(dgeList = dgeList,
                                                design = fit@design,
                                                contrastVector = contrastObject$vector,
                                                dispersion = fit@dispersion$common,
                                                useQuasiLikelihood = useQuasiLikelihood,
                                                target = pValue)
    }

    calibration <- list(p.value = nullPvalues,
                        FDR = nullFDR,
                        log2FC = nullTable$logFC,
                        abundance = nullObject$abundance,
                        n.rows = length(nullPvalues),
                        proportion.raw = rawProportion,
                        inflation = rawProportion / pValue,
                        proportion.FDR = mean(nullFDR < FDR, na.rm = TRUE),
                        n.FDR = sum(nullFDR < FDR, na.rm = TRUE),
                        median.log2FC = stats::median(nullTable$logFC, na.rm = TRUE),
                        strata = strataTable,
                        dispersion = fit@dispersion$common,
                        suggested.dispersion = suggestedDispersion,
                        source = source,
                        independent = usedHoldout | !identical(source, fitSource),
                        contrast = contrastObject$label,
                        thresholds = list(FDR = FDR, p.value = pValue),
                        verdict = verdict)

    #-------------------------------#
    # Report                        #
    #-------------------------------#
    if (isTRUE(verbose)) {
      message(calibration$n.rows, " null rows from '", source, "'",
                     if (isTRUE(usedHoldout)) {" (held out of the estimate)"} else {""}, ": ",
                     round(100 * rawProportion, 2), "% below p = ", pValue,
                     " (expected ", round(100 * pValue, 2), "%), inflation ", signif(calibration$inflation, 3),
                     ", ", calibration$n.FDR, " below FDR = ", FDR, ".")

      if (verdict == "too liberal") {
        message("The dispersion of ", signif(fit@dispersion$common, 3), " is too low.",
                if (is.finite(suggestedDispersion)) {
                  paste0(" Around ", signif(suggestedDispersion, 3), " (BCV ", signif(sqrt(suggestedDispersion), 2),
                         ") would calibrate these rows.")
                } else {""})
      } else if (verdict == "too conservative") {
        message("The dispersion is higher than the null rows need, which costs power but not credibility.")
      } else {
        message("The null p-values are close to uniform, which is what a usable dispersion looks like.")
      }

      # Overdispersion piling up in the faint rows is a filtering problem, not a dispersion problem
      if (nrow(strataTable) > 1 && strataTable$inflation[1] > 2 * strataTable$inflation[nrow(strataTable)]) {
        message("The excess sits in the weakest rows, so raising 'minCount' in filterRegions is a better answer than raising the dispersion.")
      }

      if (isFALSE(calibration$independent)) {
        message("These are the rows the dispersion was fitted to, so this check is not independent of it.")
      }
    }

    return(calibration)
  } # END function




#' @title .nullTest
#'
#' @description Fits and tests the null rows the same way the fit tested the regions.
#'
#' @param dgeList \code{DGEList} with the null rows.
#' @param design Design matrix.
#' @param contrastVector Numeric vector with the contrast.
#' @param dispersion Numeric value with the dispersion of the fit.
#' @param useQuasiLikelihood Logical value indicating whether the fit used the quasi-likelihood F test.
#'
#' @return A data.frame with the \code{logFC} and \code{PValue} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom edgeR glmFit glmLRT glmQLFit glmQLFTest
#'
#' @keywords internal

.nullTest <-
  function(dgeList,
           design,
           contrastVector,
           dispersion,
           useQuasiLikelihood = FALSE) {

    # Checking a likelihood ratio against a quasi-likelihood F would calibrate a test that was never run
    if (isTRUE(useQuasiLikelihood)) {
      nullFit <- edgeR::glmQLFit(y = dgeList, design = design, dispersion = dispersion, robust = TRUE)
      return(edgeR::glmQLFTest(glmfit = nullFit, contrast = contrastVector)$table)
    }

    nullFit <- edgeR::glmFit(y = dgeList, design = design, dispersion = dispersion)
    return(edgeR::glmLRT(glmfit = nullFit, contrast = contrastVector)$table)
  } # END function




#' @title .suggestDispersion
#'
#' @description Searches for the dispersion that would leave the expected fraction of the null rows below the p-value cut-off.
#'
#' @param dgeList \code{DGEList} with the null rows.
#' @param design Design matrix.
#' @param contrastVector Numeric vector with the contrast.
#' @param dispersion Numeric value with the dispersion the fit used, which anchors the search.
#' @param useQuasiLikelihood Logical value indicating whether the fit used the quasi-likelihood F test.
#' @param target Numeric value with the proportion aimed at.
#' @param steps Numeric value with the number of dispersions tried. Default: \code{12}.
#'
#' @return A numeric value, \code{NA} when the target lies outside the range searched.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.suggestDispersion <-
  function(dgeList,
           design,
           contrastVector,
           dispersion,
           useQuasiLikelihood = FALSE,
           target = 0.05,
           steps = 12) {

    # A grid around the current value, since the proportion moves monotonically with the dispersion but not analytically
    candidateValues <- exp(seq(from = log(dispersion / 8), to = log(dispersion * 8), length.out = steps))

    candidateProportions <-
      vapply(candidateValues,
             function(candidateDispersion) {
               nullTable <- try(.nullTest(dgeList = dgeList, design = design, contrastVector = contrastVector,
                                          dispersion = candidateDispersion, useQuasiLikelihood = useQuasiLikelihood),
                                silent = TRUE)

               if (inherits(nullTable, "try-error")) {
                 return(NA_real_)
               }
               return(mean(nullTable$PValue < target, na.rm = TRUE))
             },
             numeric(1))

    if (all(is.na(candidateProportions))) {
      return(NA_real_)
    }

    if (min(candidateProportions, na.rm = TRUE) > target | max(candidateProportions, na.rm = TRUE) < target) {
      return(NA_real_)
    }

    return(candidateValues[which.min(abs(candidateProportions - target))])
  } # END function




#' @title plotNullCalibration
#'
#' @description Draws the p-values of the rows that should not respond to a contrast, against the flat distribution they would follow if the dispersion were right.
#'
#' @param calibration List returned by \code{\link{checkNullCalibration}}.
#' @param style String with the kind of plot, one of \code{"histogram"}, \code{"qq"} and \code{"abundance"}. Default: \code{"histogram"}.
#' @param bins Numeric value with the number of bins of the histogram. Default: \code{50}.
#' @param colours Character vector of length two with the colours of the data and of the expected line. Default: \code{NULL}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}, the contrast.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}, the dispersion and the verdict.
#' @param legendPosition String with the position of the legend. Default: \code{"none"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object.
#'
#' @details A flat histogram sitting on the expected line is the target. Bars piling up near zero mean the dispersion is too low and the fit is calling library-to-library variation a treatment effect; the same optimism is in the regions, at the same rate. A histogram sagging near zero means the opposite, which costs sensitivity and nothing else.
#'
#' The quantile plot says the same thing with more resolution in the tail, which is where the regions that end up in a figure come from. The abundance plot says where the trouble is: a flat line across the strata means the dispersion is wrong everywhere, while a line rising towards the weak rows means the filter is too loose and \code{\link{filterRegions}} is the fix.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' calibration <- checkNullCalibration(fit,
#'                                     contrast = c("condition", "SHR", "BN"),
#'                                     source = "background",
#'                                     verbose = FALSE)
#'
#' # A flat histogram is what a well calibrated test looks like
#' plotNullCalibration(calibration)
#'
#' plotNullCalibration(calibration, style = "qq")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{checkNullCalibration}}, \code{\link{estimateNullDispersion}}
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_hline geom_point geom_line geom_abline labs theme element_blank element_rect element_line
#' @importFrom rlang .data
#' @importFrom stats ppoints
#'
#' @export plotNullCalibration

plotNullCalibration <-
  function(calibration,
           style = "histogram",
           bins = 50,
           colours = NULL,
           title = NULL,
           subtitle = NULL,
           legendPosition = "none",
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!is.list(calibration) | is.null(calibration$p.value)) {
      stop("The 'calibration' parameter must be the list returned by checkNullCalibration.", call. = FALSE)
    }

    if (!(style %in% c("histogram", "qq", "abundance"))) {
      stop("The 'style' parameter must be one of 'histogram', 'qq', 'abundance'.", call. = FALSE)
    }

    colourValues <- if (is.null(colours)) {c("grey55", "#B2182B")} else {colours}

    if (is.null(subtitle)) {
      subtitle <- sprintf("dispersion %.3g (BCV %.2f), inflation %.2f, *%s*%s",
                          calibration$dispersion, sqrt(calibration$dispersion),
                          calibration$inflation, calibration$verdict,
                          if (isFALSE(calibration$independent)) {" (not independent)"} else {""})
    }

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    if (style == "histogram") {
      calibrationPlot <-
        ggplot2::ggplot(data = data.frame(p.value = calibration$p.value),
                        mapping = ggplot2::aes(x = .data$p.value)) +
        ggplot2::geom_histogram(bins = bins, fill = colourValues[1], colour = NA) +
        ggplot2::geom_hline(yintercept = calibration$n.rows / bins,
                            linetype = "dashed", linewidth = 0.5, colour = colourValues[2]) +
        ggplot2::labs(x = "p-value of the null rows", y = "Rows")

    } else if (style == "qq") {
      # The tail is where the regions that reach a figure come from, and a quantile plot shows it properly
      observedValues <- sort(calibration$p.value)
      plotTable <- data.frame(expected = -log10(stats::ppoints(length(observedValues))),
                              observed = -log10(observedValues))
      plotTable <- plotTable[.thinIndex(n = nrow(plotTable), maxPoints = 20000), , drop = FALSE]

      calibrationPlot <-
        ggplot2::ggplot(data = plotTable, mapping = ggplot2::aes(x = .data$expected, y = .data$observed)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.5, colour = colourValues[2]) +
        ggplot2::geom_point(size = 0.7, colour = colourValues[1], stroke = NA) +
        ggplot2::labs(x = "Expected -log<sub>10</sub>(p)", y = "Observed -log<sub>10</sub>(p)")

    } else {
      calibrationPlot <-
        ggplot2::ggplot(data = calibration$strata,
                        mapping = ggplot2::aes(x = .data$median.abundance, y = .data$inflation)) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, colour = colourValues[2]) +
        ggplot2::geom_line(linewidth = 0.5, colour = colourValues[1]) +
        ggplot2::geom_point(size = 2.5, colour = colourValues[1], stroke = NA) +
        ggplot2::labs(x = "Median abundance of the stratum (log<sub>2</sub> CPM)",
                      y = "Observed over expected")
    }

    calibrationPlot <- calibrationPlot +
      ggplot2::labs(title = if (is.null(title)) {calibration$contrast} else {title},
                    subtitle = subtitle) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize) +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5, colour = "black"),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "gray"))

    return(calibrationPlot)
  } # END function
