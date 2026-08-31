#' @title plotNormComparison
#'
#' @description Compares the scaling factors that different normalisation methods give for the same object, without modifying it. The factors already stored in the object, whether estimated or supplied by hand, are shown alongside the others so that a manual set of factors can be checked against the automatic ones.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param methods Character vector with the methods to compare, among those accepted by \code{\link{normalizeCounts}} that need no extra input: \code{"librarySize"}, \code{"TMM"}, \code{"TMMwsp"}, \code{"RLE"}, \code{"upperQuartile"} and \code{"background"}. Default: \code{c("librarySize", "TMM", "RLE", "background")}.
#' @param plotType String with the type of plot, either \code{"factors"} to show one point per sample and method, or \code{"ma"} to show the counts of each sample against a reference with the factors drawn as horizontal lines. Default: \code{"factors"}.
#' @param referenceSample String or numeric position of the sample used as reference by the MA plot. Default: \code{NULL}, the sample with the median depth.
#' @param useBackground Logical value indicating whether the MA plot must be drawn on the background bins rather than on the regions. Default: \code{FALSE}.
#' @param useRegionSets Character vector with the names of the region sets used to estimate the factors. Default: \code{NULL}, all of them.
#' @param minCount Numeric value with the minimum total count required for a region to take part in the estimation. Default: \code{1}.
#' @param priorCount Numeric value added to the counts before the log transformation of the MA plot. Default: \code{1}.
#' @param maxRegions Numeric value with the maximum number of regions drawn in the MA plot, thinned at regular intervals when exceeded. Default: \code{20000}.
#' @param pointSize Numeric value with the size of the points of the factor plot. Default: \code{3}.
#' @param facetScales String passed to the facets of the MA plot, one among \code{"fixed"}, \code{"free"}, \code{"free_x"} or \code{"free_y"}. Default: \code{"fixed"}.
#' @param title String with the title of the plot. Default: \code{NULL}, a title describing the plot type.
#' @param returnData Logical value indicating whether the table behind the plot must be returned instead of the plot. Default: \code{FALSE}.
#'
#' @return A \code{ggplot} object, or a data.frame when \code{returnData} is \code{TRUE}.
#'
#' @details The MA plot is drawn on the raw counts, so what it shows is the bias before any correction: a cloud sitting away from zero means that the sample and the reference differ by more than a common factor. The line of each method marks the shift that method would subtract, which makes the comparison direct, a line running through the middle of the cloud describes the data while one sitting off to the side does not. Methods estimated on the regions follow the regions by construction, so the interesting comparison is against \code{"background"} or against a manual set of factors, which are free to disagree.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' # Scaling factors from four methods side by side, before committing to one
#' plotNormComparison(counts)
#'
#' # The numbers behind the panel
#' head(plotNormComparison(counts, returnData = TRUE))
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{normalizeCounts}}, \code{\link{plotSetMA}}
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom S4Vectors metadata
#' @importFrom dplyr bind_rows filter mutate select left_join
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_hline geom_bin2d facet_wrap labs theme element_text element_blank scale_y_continuous scale_fill_gradient
#' @importFrom stats median
#' @importFrom methods is
#'
#' @export plotNormComparison

plotNormComparison <-
  function(counts,
           methods = c("librarySize", "TMM", "RLE", "background"),
           plotType = "factors",
           referenceSample = NULL,
           useBackground = FALSE,
           useRegionSets = NULL,
           minCount = 1,
           priorCount = 1,
           maxRegions = 20000,
           pointSize = 3,
           facetScales = "fixed",
           title = NULL,
           returnData = FALSE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    plotType <- tolower(plotType[1])
    if (!(plotType %in% c("factors", "ma"))) {
      stop("The 'plotType' parameter must be either 'factors' or 'ma'.", call. = FALSE)
    }

    facetScales <- tolower(facetScales[1])
    if (!(facetScales %in% c("fixed", "free", "free_x", "free_y"))) {
      stop("The 'facetScales' parameter must be one among 'fixed', 'free', 'free_x' or 'free_y'.", call. = FALSE)
    }

    availableMethods <- c("librarySize", "TMM", "TMMwsp", "RLE", "upperQuartile", "background")
    unknownMethods <- setdiff(methods, availableMethods)
    if (length(unknownMethods) > 0) {
      stop(paste0("The following methods cannot be compared here: ", paste(unknownMethods, collapse = ", "),
                  ". Only the ones needing no extra input are accepted: ", paste(availableMethods, collapse = ", "), "."), call. = FALSE)
    }

    #-----------------------------------#
    # Collect the factors of each method #
    #-----------------------------------#
    factorList <-
      lapply(methods,
             function(oneMethod) {
               # A method can legitimately fail, background without bins for instance, and that must not stop the comparison
               estimated <- tryCatch(normalizeCounts(counts = counts,
                                                     method = oneMethod,
                                                     useRegionSets = useRegionSets,
                                                     minCount = minCount,
                                                     verbose = FALSE),
                                     error = function(e) {
                                       warning(paste0("The '", oneMethod, "' method could not be applied: ", conditionMessage(e)), call. = FALSE)
                                       return(NULL)
                                     })

               if (is.null(estimated)) {return(NULL)}

               return(data.frame(sample = colnames(counts),
                                 method = oneMethod,
                                 scaling.factor = as.numeric(SummarizedExperiment::colData(estimated)$scaling.factor),
                                 stringsAsFactors = FALSE))
             })

    # The factors already in the object are the ones the analysis will actually use
    appliedFactors <- SummarizedExperiment::colData(counts)$scaling.factor

    if (!is.null(appliedFactors) & !all(is.na(appliedFactors))) {
      appliedMethod <- S4Vectors::metadata(counts)$normalization$method
      factorList <- c(factorList,
                      list(data.frame(sample = colnames(counts),
                                      method = paste0("applied (", ifelse(is.null(appliedMethod), "unknown", appliedMethod), ")"),
                                      scaling.factor = as.numeric(appliedFactors),
                                      stringsAsFactors = FALSE)))
    }

    factorTable <- dplyr::bind_rows(factorList)

    if (nrow(factorTable) == 0) {
      stop("None of the requested methods could be applied to this object.", call. = FALSE)
    }

    factorTable <- dplyr::mutate(factorTable, sample = factor(.data$sample, levels = colnames(counts)))

    #---------------------------#
    # One point per sample      #
    #---------------------------#
    if (plotType == "factors") {
      if (isTRUE(returnData)) {return(factorTable)}

      normalizationPlot <-
        ggplot2::ggplot(data = factorTable,
                        mapping = ggplot2::aes(x = .data$sample, y = .data$scaling.factor, color = .data$method, group = .data$method)) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
        ggplot2::geom_line(alpha = 0.6, show.legend = FALSE) +
        ggplot2::geom_point(size = pointSize, stroke = NA) +
        ggplot2::scale_y_continuous(trans = "log2") +
        ggplot2::labs(title = ifelse(is.null(title), "Scaling factors across methods", title),
                      x = NULL,
                      y = "Scaling factor (divisor, log<sub>2</sub> scale)",
                      color = "Method") +
        .regionSetTheme(legendPosition = "right") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.line = ggplot2::element_blank(),
                       panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.5),
                       panel.grid.major.x = ggplot2::element_line(colour = "gray", linewidth = 0.2),
                       panel.grid.minor.y = ggplot2::element_line(colour = "gray", linewidth = 0.2))

      return(normalizationPlot)
    }

    #---------------------------#
    # Counts against a reference #
    #---------------------------#
    if (isTRUE(useBackground)) {
      backgroundBins <- S4Vectors::metadata(counts)$background
      if (is.null(backgroundBins)) {
        stop("No background bin is stored in the object, run 'countBackground' before using 'useBackground'.", call. = FALSE)
      }
      valueMatrix <- SummarizedExperiment::assay(backgroundBins, "counts")
    } else {
      valueMatrix <- SummarizedExperiment::assay(counts, "counts")
    }

    referenceIndex <- .resolveSampleIndex(sample = referenceSample, sampleNames = colnames(counts))
    if (is.null(referenceIndex)) {
      referenceIndex <- which.min(abs(colSums(valueMatrix) - stats::median(colSums(valueMatrix))))[1]
    }

    # Thinning at regular intervals keeps the picture and avoids drawing millions of points
    rowIndex <- unique(round(seq(from = 1, to = nrow(valueMatrix), length.out = min(maxRegions, nrow(valueMatrix)))))
    logMatrix <- log2(valueMatrix[rowIndex, , drop = FALSE] + priorCount)

    maTable <-
      dplyr::bind_rows(lapply(setdiff(seq_len(ncol(logMatrix)), referenceIndex),
                              function(sampleIndex) {
                                data.frame(sample = colnames(counts)[sampleIndex],
                                           A = (logMatrix[, sampleIndex] + logMatrix[, referenceIndex]) / 2,
                                           M = logMatrix[, sampleIndex] - logMatrix[, referenceIndex],
                                           stringsAsFactors = FALSE)
                              }))

    # Regions empty in both samples pile up on the origin and hide the rest of the cloud
    maTable <- dplyr::filter(maTable, .data$A > log2(priorCount + 1))
    maTable <- dplyr::mutate(maTable, sample = factor(.data$sample, levels = colnames(counts)))

    referenceFactors <- dplyr::filter(factorTable, .data$sample == colnames(counts)[referenceIndex])
    referenceFactors <- dplyr::select(referenceFactors, method = "method", reference.factor = "scaling.factor")

    lineTable <- dplyr::filter(factorTable, .data$sample != colnames(counts)[referenceIndex])
    lineTable <- dplyr::left_join(lineTable, referenceFactors, by = "method")

    # The line marks how far the method would move the sample relative to the reference
    lineTable <- dplyr::mutate(lineTable, shift = log2(.data$scaling.factor / .data$reference.factor))

    if (isTRUE(returnData)) {return(list(ma = maTable, lines = lineTable))}

    maPlot <-
      ggplot2::ggplot(data = maTable, mapping = ggplot2::aes(x = .data$A, y = .data$M)) +
      ggplot2::geom_bin2d(bins = 80) +
      ggplot2::scale_fill_gradient(low = "gray85", high = "gray15", trans = "log10") +
      ggplot2::geom_hline(yintercept = 0, color = "gray40") +
      ggplot2::geom_hline(data = lineTable, mapping = ggplot2::aes(yintercept = .data$shift, color = .data$method), linewidth = 0.7) +
      ggplot2::facet_wrap(facets = ~ sample, scales = facetScales) +
      ggplot2::labs(title = ifelse(is.null(title), paste0("Counts against ", colnames(counts)[referenceIndex]), title),
                    x = "Average log<sub>2</sub> count",
                    y = "log<sub>2</sub> count ratio",
                    color = "Method",
                    fill = "Regions") +
      .regionSetTheme(legendPosition = "right")

    return(maPlot)
  } # END function
