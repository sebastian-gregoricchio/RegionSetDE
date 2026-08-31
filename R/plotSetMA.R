#' @title plotSetMA
#'
#' @description Draws the log ratio of two groups of samples against their average abundance, with the regions of a set highlighted over the others. A set whose cloud sits away from zero while the rest stays on it is the picture the package is built to produce, and seeing it before any test tells whether the normalisation has already decided the answer.
#'
#' @param counts \code{RegionSetDE.counts} object, normalised by \code{\link{normalizeCounts}}.
#' @param set Character vector with the names of the region sets to highlight. Default: \code{NULL}, one panel per set with the others greyed behind.
#' @param groupBy String with the name of the \code{colData} column defining the groups. Default: \code{NULL}, accepted only when the object holds two samples.
#' @param contrast Character vector of length two with the levels of \code{groupBy} to compare, given as \code{c(numerator, denominator)}. Default: \code{NULL}, the two levels found in the column, in alphabetical order.
#' @param assayName String with the assay to plot. Default: \code{NULL}, the normalised assay recorded by \code{\link{normalizeCounts}}, or the raw counts when the object has not been normalised.
#' @param priorCount Numeric value added to the values before the log transformation. Default: \code{1}.
#' @param minCount Numeric value with the minimum total count required to draw a region. Default: \code{1}.
#' @param showTrend Logical value indicating whether a loess trend must be drawn over the highlighted regions. Default: \code{TRUE}.
#' @param showMedian Logical value indicating whether the median log ratio of the highlighted regions must be drawn as a dashed line. Default: \code{TRUE}.
#' @param highlightColor String with the colour of the highlighted regions. Default: \code{"#B22222"}.
#' @param backgroundColor String with the colour of the other regions. Default: \code{"gray75"}.
#' @param pointSize Numeric value with the size of the points. Default: \code{1.2}.
#' @param maxPoints Numeric value with the maximum number of points used for the grey backdrop and for the trend fit, thinned at regular intervals when exceeded. The highlighted regions are always drawn in full. Default: \code{10000}.
#' @param facetScales String passed to the facets, one among \code{"fixed"}, \code{"free"}, \code{"free_x"} or \code{"free_y"}. Default: \code{"fixed"}.
#' @param title String with the title of the plot. Default: \code{NULL}, no title.
#' @param returnData Logical value indicating whether the table behind the plot must be returned instead of the plot. Default: \code{FALSE}.
#'
#' @return A \code{ggplot} object, or a data.frame when \code{returnData} is \code{TRUE}.
#'
#' @details The ratio is computed on the mean of each group, so it carries no dispersion estimate and no test, it only shows where the regions sit. When the object has not been normalised the plot runs on the raw counts and a global offset between the groups is expected, that offset is what the normalisation removes.
#'
#' The median line is the number to read: a set whose median sits at zero behaves like the rest of the genome under the chosen normalisation, whichever way the individual regions scatter. When every set shares the same non zero median, the shift belongs to the normalisation rather than to the biology.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#' # One panel per region set, each against the rest of the data
#' plotSetMA(counts, groupBy = "condition")
#'
#' plotSetMA(counts, set = "promoterCpG", groupBy = "condition")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{normalizeCounts}}, \code{\link{plotNormComparison}}
#'
#' @importFrom SummarizedExperiment assay assayNames colData rowData
#' @importFrom S4Vectors metadata
#' @importFrom dplyr filter mutate group_by summarise select slice n
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_smooth facet_wrap labs
#' @importFrom stats median
#' @importFrom methods is
#'
#' @export plotSetMA

plotSetMA <-
  function(counts,
           set = NULL,
           groupBy = NULL,
           contrast = NULL,
           assayName = NULL,
           priorCount = 1,
           minCount = 1,
           showTrend = TRUE,
           showMedian = TRUE,
           highlightColor = "#B22222",
           backgroundColor = "gray75",
           pointSize = 1.2,
           maxPoints = 10000,
           facetScales = "fixed",
           title = NULL,
           returnData = FALSE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    facetScales <- tolower(facetScales[1])
    if (!(facetScales %in% c("fixed", "free", "free_x", "free_y"))) {
      stop("The 'facetScales' parameter must be one among 'fixed', 'free', 'free_x' or 'free_y'.", call. = FALSE)
    }

    # The normalised values are the ones the plot is meant to judge, the raw counts are only a fallback
    if (is.null(assayName)) {
      assayName <- S4Vectors::metadata(counts)$normalization$normalized.assay
      if (is.null(assayName)) {assayName <- "counts"}
    }

    if (!(assayName %in% SummarizedExperiment::assayNames(counts))) {
      stop(paste0("The assay '", assayName, "' is absent from the object."), call. = FALSE)
    }

    valueMatrix <- SummarizedExperiment::assay(counts, assayName)

    #---------------------------#
    # Groups to be compared     #
    #---------------------------#
    if (is.null(groupBy)) {
      if (ncol(counts) != 2) {
        stop("The 'groupBy' parameter is required whenever the object holds more than two samples.", call. = FALSE)
      }
      groupVector <- colnames(counts)
      contrast <- if (is.null(contrast)) {colnames(counts)} else {contrast}
    } else {
      if (!(groupBy %in% colnames(SummarizedExperiment::colData(counts)))) {
        stop(paste0("The column '", groupBy, "' is absent from the sample metadata."), call. = FALSE)
      }
      groupVector <- as.character(SummarizedExperiment::colData(counts)[[groupBy]])
      if (is.null(contrast)) {contrast <- sort(unique(groupVector))}
    }

    if (length(contrast) != 2) {
      stop("The 'contrast' parameter must name exactly two groups, as c(numerator, denominator).", call. = FALSE)
    }

    absentGroups <- setdiff(contrast, groupVector)
    if (length(absentGroups) > 0) {
      stop(paste0("The following groups are absent from the samples: ", paste(absentGroups, collapse = ", "), "."), call. = FALSE)
    }

    numeratorColumns <- which(groupVector == contrast[1])
    denominatorColumns <- which(groupVector == contrast[2])

    #---------------------------#
    # Ratio and abundance       #
    #---------------------------#
    numeratorMean <- rowMeans(valueMatrix[, numeratorColumns, drop = FALSE])
    denominatorMean <- rowMeans(valueMatrix[, denominatorColumns, drop = FALSE])

    maTable <- data.frame(region.set = as.character(SummarizedExperiment::rowData(counts)$region.set),
                          region.id = as.character(SummarizedExperiment::rowData(counts)$region.id),
                          total.count = as.numeric(rowSums(SummarizedExperiment::assay(counts, "counts"))),
                          A = (log2(numeratorMean + priorCount) + log2(denominatorMean + priorCount)) / 2,
                          M = log2(numeratorMean + priorCount) - log2(denominatorMean + priorCount),
                          stringsAsFactors = FALSE)

    maTable <- dplyr::filter(maTable, .data$total.count >= minCount)

    if (nrow(maTable) == 0) {
      stop("No region passes 'minCount', there is nothing to draw.", call. = FALSE)
    }

    if (!is.null(set)) {
      absentSets <- setdiff(set, unique(maTable$region.set))
      if (length(absentSets) > 0) {
        stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }
    }

    #---------------------------#
    # Highlight and panels      #
    #---------------------------#
    # With no set named every set gets its own panel, with the whole object greyed underneath as a reference
    if (is.null(set)) {
      highlightTable <- dplyr::mutate(maTable, panel = .data$region.set)
      backdropTable <- dplyr::select(maTable, .data$A, .data$M)
    } else {
      highlightTable <- dplyr::mutate(dplyr::filter(maTable, .data$region.set %in% set), panel = paste(set, collapse = ", "))
      backdropTable <- dplyr::select(dplyr::filter(maTable, !(.data$region.set %in% set)), .data$A, .data$M)
    }

    # The backdrop carries no panel column, so it is redrawn in every facet: thinning it keeps that cost flat
    backdropTable <- backdropTable[.thinIndex(n = nrow(backdropTable), maxPoints = maxPoints), , drop = FALSE]

    if (isTRUE(returnData)) {
      return(dplyr::mutate(maTable, highlighted = .data$region.set %in% if (is.null(set)) {unique(maTable$region.set)} else {set}))
    }

    medianTable <- dplyr::summarise(dplyr::group_by(highlightTable, .data$panel), median.M = stats::median(.data$M), .groups = "drop")

    #---------------------------#
    # Assemble the plot         #
    #---------------------------#
    setPlot <-
      ggplot2::ggplot(mapping = ggplot2::aes(x = .data$A, y = .data$M)) +
      ggplot2::geom_point(data = backdropTable, color = backgroundColor, size = pointSize, alpha = 0.4, stroke = NA) +
      ggplot2::geom_point(data = highlightTable, color = highlightColor, size = pointSize, alpha = 0.6, stroke = NA) +
      ggplot2::geom_hline(yintercept = 0, color = "gray30")

    if (isTRUE(showMedian)) {
      setPlot <- setPlot + ggplot2::geom_hline(data = medianTable,
                                               mapping = ggplot2::aes(yintercept = .data$median.M),
                                               color = highlightColor, linetype = "dashed")
    }

    if (isTRUE(showTrend)) {
      # The loess cost grows with the square of the points and one fit is run per panel, a thinned sample draws the same line
      trendTable <- dplyr::slice(dplyr::group_by(highlightTable, .data$panel),
                                 .thinIndex(n = dplyr::n(), maxPoints = maxPoints))

      setPlot <- setPlot + ggplot2::geom_smooth(data = trendTable, method = "loess", formula = y ~ x,
                                                se = FALSE, color = "black", linewidth = 0.6)
    }

    setPlot <-
      setPlot +
      ggplot2::facet_wrap(facets = ~ panel, scales = facetScales) +
      ggplot2::labs(title = title,
                    x = paste0("Average log<sub>2</sub> signal (", assayName, ")"),
                    y = paste0("log<sub>2</sub> ", contrast[1], " / ", contrast[2]),
                    caption = "Dashed line: median of the highlighted regions") +
      .regionSetTheme(legendPosition = "none")

    return(setPlot)
  } # END function
