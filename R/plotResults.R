#' @title plotVolcano
#'
#' @description Draws the log2 fold change of a contrast against the significance, one panel per region set, with the number of changing regions written in the top corners of each panel.
#'
#' @param results \code{RegionSetDE.results} or \code{RegionSetDE.resultsList} object.
#' @param set Character vector with the names of the region sets to draw. Default: \code{NULL}, all of them.
#' @param contrast String with the name of the contrast to draw, or its position, when \code{results} holds several of them. Default: \code{NULL}.
#' @param colourBy String with the variable driving the colour, either \code{"diff.status"} or \code{"region.set"}. Default: \code{"diff.status"}.
#' @param facetBySet Logical value to indicate whether each region set must get its own panel. Default: \code{TRUE}.
#' @param facetScales String with the scales of the panels, one among \code{"fixed"}, \code{"free"}, \code{"free_x"} and \code{"free_y"}. Default: \code{"fixed"}.
#' @param yValue String with the quantity on the y axis, either \code{"FDR"} or \code{"p.value"}. Default: \code{"FDR"}.
#' @param FDR Numeric value with the adjusted p-value cut-off drawn as a line. Default: \code{NULL}, the threshold stored in the object.
#' @param log2FC Numeric value with the absolute log2 fold change cut-off drawn as a line. Default: \code{NULL}, the threshold stored in the object.
#' @param showCounts Logical value to indicate whether the number of changing regions must be written in the top corners of each panel. Default: \code{TRUE}.
#' @param labelTop Numeric value with the number of top regions to label, per panel. Default: \code{0}.
#' @param labelColumn String with the column holding the labels. Default: \code{"region.id"}.
#' @param colours Named character vector with the colours. Default: \code{NULL}, a grey, blue and red palette for \code{"diff.status"}.
#' @param pointSize Numeric value with the size of the points. Default: \code{0.8}.
#' @param maxPoints Numeric value with the number of non-changing points drawn per panel. Default: \code{20000}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}, the contrast.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}.
#' @param legendPosition String with the position of the legend. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object.
#'
#' @details The y axis is the adjusted p-value by default. Showing the raw p-value while drawing the cut-off line at the FDR puts two different quantities on the same figure, which is where most misread volcano plots come from.
#'
#' Only the points labelled \code{"null"} are thinned by \code{maxPoints}, and the thinning happens inside each panel so that a small set keeps all of its points. Everything passing the thresholds is drawn. The thinning is deterministic, so the figure does not change between calls.
#'
#' The counts in the corners come from the full table, before any thinning, and they are the number of regions labelled \code{"down"} on the left and \code{"up"} on the right.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' plotVolcano(results)
#'
#' # One set, with the strongest regions labelled
#' plotVolcano(results, set = "promoterCpG", labelTop = 5)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegions}}, \code{\link{plotResultsMA}}, \code{\link{topRegions}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline geom_text scale_colour_manual facet_wrap labs guides guide_legend
#' @importFrom dplyr filter arrange group_by slice_head ungroup
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export plotVolcano

plotVolcano <-
  function(results,
           set = NULL,
           contrast = NULL,
           colourBy = "diff.status",
           facetBySet = TRUE,
           facetScales = "fixed",
           yValue = "FDR",
           FDR = NULL,
           log2FC = NULL,
           showCounts = TRUE,
           labelTop = 0,
           labelColumn = "region.id",
           colours = NULL,
           pointSize = 0.8,
           maxPoints = 20000,
           title = NULL,
           subtitle = NULL,
           legendPosition = "right",
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    results <- .pickResults(results = results, contrast = contrast)

    if (!(colourBy %in% c("diff.status", "region.set"))) {
      stop("The 'colourBy' parameter must be either 'diff.status' or 'region.set'.", call. = FALSE)
    }

    if (!(yValue %in% c("FDR", "p.value"))) {
      stop("The 'yValue' parameter must be either 'FDR' or 'p.value'.", call. = FALSE)
    }

    if (!(facetScales %in% c("fixed", "free", "free_x", "free_y"))) {
      stop("The 'facetScales' parameter must be one among 'fixed', 'free', 'free_x' or 'free_y'.", call. = FALSE)
    }

    plotTable <- .prepareResultTable(results = results, set = set, FDR = FDR, log2FC = log2FC)

    FDRthreshold <- if (is.null(FDR)) {results@thresholds$FDR} else {FDR}
    log2FCthreshold <- if (is.null(log2FC)) {results@thresholds$log2FC} else {log2FC}

    plotTable$y.value <- -log10(plotTable[[yValue]])

    statusColours <- .diffStatusColours(colours = colours)

    #-------------------------------#
    # Thin the unchanging cloud     #
    #-------------------------------#
    drawnTable <- rbind(.thinBySet(regionTable = dplyr::filter(plotTable, .data$diff.status == "null"),
                                   maxPoints = maxPoints,
                                   bySet = facetBySet),
                        dplyr::filter(plotTable, .data$diff.status != "null"))

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    volcanoPlot <-
      ggplot2::ggplot(data = drawnTable,
                      mapping = ggplot2::aes(x = .data$log2FC, y = .data$y.value, colour = .data[[colourBy]])) +
      ggplot2::geom_point(size = pointSize, alpha = 0.6, stroke = NA) +
      ggplot2::geom_hline(yintercept = -log10(FDRthreshold), linetype = "dashed", linewidth = 0.3, colour = "black") +
      ggplot2::labs(x = paste0("log<sub>2</sub> fold change (", results@contrast, ")"),
                    y = paste0("-log<sub>10</sub>(", if (yValue == "FDR") {"FDR"} else {"p-value"}, ")"),
                    colour = if (colourBy == "diff.status") {"Status"} else {"Region set"},
                    title = if (is.null(title)) {results@contrast} else {title},
                    subtitle = subtitle) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = max(c(pointSize, 3)), alpha = 1))) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize)

    if (log2FCthreshold > 0) {
      volcanoPlot <- volcanoPlot +
        ggplot2::geom_vline(xintercept = c(-log2FCthreshold, log2FCthreshold),
                            linetype = "dashed", linewidth = 0.3, colour = "black")
    }

    if (colourBy == "diff.status") {
      volcanoPlot <- volcanoPlot + ggplot2::scale_colour_manual(values = statusColours, drop = FALSE)
    } else if (!is.null(colours)) {
      volcanoPlot <- volcanoPlot + ggplot2::scale_colour_manual(values = colours)
    }

    if (isTRUE(facetBySet)) {
      volcanoPlot <- volcanoPlot + ggplot2::facet_wrap(facets = ~ region.set, scales = facetScales)
    }

    #-------------------------------#
    # Counts in the corners         #
    #-------------------------------#
    if (isTRUE(showCounts)) {
      countTable <- .diffCounts(regionTable = plotTable, bySet = facetBySet)

      volcanoPlot <- volcanoPlot +
        ggplot2::geom_text(data = countTable,
                           mapping = ggplot2::aes(x = -Inf, y = Inf, label = paste0("n = ", .data$n.down)),
                           inherit.aes = FALSE, hjust = -0.15, vjust = 1.4,
                           size = baseSize / 4.5, colour = statusColours[["down"]]) +
        ggplot2::geom_text(data = countTable,
                           mapping = ggplot2::aes(x = Inf, y = Inf, label = paste0("n = ", .data$n.up)),
                           inherit.aes = FALSE, hjust = 1.15, vjust = 1.4,
                           size = baseSize / 4.5, colour = statusColours[["up"]])
    }

    #-------------------------------#
    # Labels on the top regions     #
    #-------------------------------#
    if (labelTop > 0) {
      if (!(labelColumn %in% colnames(plotTable))) {
        stop(paste0("The column '", labelColumn, "' is absent from the results table."), call. = FALSE)
      }
      if (!requireNamespace("ggrepel", quietly = TRUE)) {
        stop("The 'ggrepel' package is needed to label the points.", call. = FALSE)
      }

      labelTable <- dplyr::arrange(dplyr::filter(plotTable, .data$diff.status != "null"), .data$FDR, .data$p.value)
      labelTable <- if (isTRUE(facetBySet)) {
        dplyr::ungroup(dplyr::slice_head(dplyr::group_by(labelTable, .data$region.set), n = labelTop))
      } else {
        dplyr::slice_head(labelTable, n = labelTop)
      }

      volcanoPlot <- volcanoPlot +
        ggrepel::geom_text_repel(data = labelTable,
                                 mapping = ggplot2::aes(label = .data[[labelColumn]]),
                                 size = baseSize / 4.5, colour = "black", max.overlaps = Inf, show.legend = FALSE)
    }

    return(volcanoPlot)
  } # END function




#' @title plotResultsMA
#'
#' @description Draws the log2 fold change of a contrast against the average signal, one panel per region set, which shows whether the response depends on how much signal a region carried to begin with.
#'
#' @param results \code{RegionSetDE.results} or \code{RegionSetDE.resultsList} object.
#' @param set Character vector with the names of the region sets to draw. Default: \code{NULL}, all of them.
#' @param contrast String with the name of the contrast to draw, or its position, when \code{results} holds several of them. Default: \code{NULL}.
#' @param colourBy String with the variable driving the colour, either \code{"diff.status"} or \code{"region.set"}. Default: \code{"diff.status"}.
#' @param facetBySet Logical value to indicate whether each region set must get its own panel. Default: \code{TRUE}.
#' @param facetScales String with the scales of the panels, one among \code{"fixed"}, \code{"free"}, \code{"free_x"} and \code{"free_y"}. Default: \code{"fixed"}.
#' @param FDR Numeric value with the adjusted p-value cut-off used to label the points. Default: \code{NULL}, the threshold stored in the object.
#' @param log2FC Numeric value with the absolute log2 fold change cut-off used to label the points. Default: \code{NULL}, the threshold stored in the object.
#' @param showCounts Logical value to indicate whether the number of changing regions must be written in the corners of each panel, on the right hand side and on the same side of zero as the regions they count. Default: \code{TRUE}.
#' @param showTrend Logical value to indicate whether a loess trend must be drawn. Default: \code{TRUE}.
#' @param colours Named character vector with the colours. Default: \code{NULL}.
#' @param pointSize Numeric value with the size of the points. Default: \code{0.8}.
#' @param maxPoints Numeric value with the number of non-changing points drawn per panel. Default: \code{20000}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}, the contrast.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}.
#' @param legendPosition String with the position of the legend. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object.
#'
#' @details A trend that leaves zero at one end of the abundance range is the usual sign that the normalisation has not done its job, and it is easier to see here than on any summary statistic. This plot answers a different question from \code{\link{plotSetMA}}, which compares samples before any model is fitted: here the y axis is a fitted coefficient rather than a difference between two libraries.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' plotResultsMA(results)
#'
#' plotResultsMA(results, set = "promoterCpG", facetBySet = FALSE)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotVolcano}}, \code{\link{plotSetMA}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_smooth geom_text scale_colour_manual facet_wrap labs guides guide_legend
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export plotResultsMA

plotResultsMA <-
  function(results,
           set = NULL,
           contrast = NULL,
           colourBy = "diff.status",
           facetBySet = TRUE,
           facetScales = "fixed",
           FDR = NULL,
           log2FC = NULL,
           showCounts = TRUE,
           showTrend = TRUE,
           colours = NULL,
           pointSize = 0.8,
           maxPoints = 20000,
           title = NULL,
           subtitle = NULL,
           legendPosition = "right",
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    results <- .pickResults(results = results, contrast = contrast)

    if (!(colourBy %in% c("diff.status", "region.set"))) {
      stop("The 'colourBy' parameter must be either 'diff.status' or 'region.set'.", call. = FALSE)
    }

    if (!(facetScales %in% c("fixed", "free", "free_x", "free_y"))) {
      stop("The 'facetScales' parameter must be one among 'fixed', 'free', 'free_x' or 'free_y'.", call. = FALSE)
    }

    plotTable <- .prepareResultTable(results = results, set = set, FDR = FDR, log2FC = log2FC)

    if (!("average.signal" %in% colnames(plotTable))) {
      stop("The results table carries no 'average.signal' column.", call. = FALSE)
    }

    statusColours <- .diffStatusColours(colours = colours)

    #-------------------------------#
    # Thin the unchanging cloud     #
    #-------------------------------#
    drawnTable <- rbind(.thinBySet(regionTable = dplyr::filter(plotTable, .data$diff.status == "null"),
                                   maxPoints = maxPoints,
                                   bySet = facetBySet),
                        dplyr::filter(plotTable, .data$diff.status != "null"))

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    maPlot <-
      ggplot2::ggplot(data = drawnTable,
                      mapping = ggplot2::aes(x = .data$average.signal, y = .data$log2FC, colour = .data[[colourBy]])) +
      ggplot2::geom_point(size = pointSize, alpha = 0.6, stroke = NA) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
      ggplot2::labs(x = "Average signal (log<sub>2</sub> CPM)",
                    y = paste0("log<sub>2</sub> fold change (", results@contrast, ")"),
                    colour = if (colourBy == "diff.status") {"Status"} else {"Region set"},
                    title = if (is.null(title)) {results@contrast} else {title},
                    subtitle = subtitle) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = max(c(pointSize, 3)), alpha = 1))) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize)

    if (colourBy == "diff.status") {
      maPlot <- maPlot + ggplot2::scale_colour_manual(values = statusColours, drop = FALSE)
    } else if (!is.null(colours)) {
      maPlot <- maPlot + ggplot2::scale_colour_manual(values = colours)
    }

    if (isTRUE(showTrend)) {
      # One loess per set, on a thinned sample: the fit costs the square of the points and draws the same line either way
      trendTable <- .thinBySet(regionTable = plotTable, maxPoints = maxPoints, bySet = TRUE)

      maPlot <- maPlot +
        ggplot2::geom_smooth(data = trendTable,
                             mapping = ggplot2::aes(group = .data$region.set),
                             method = "loess", formula = y ~ x, se = FALSE,
                             colour = "black", linewidth = 0.5)
    }

    if (isTRUE(facetBySet)) {
      maPlot <- maPlot + ggplot2::facet_wrap(facets = ~ region.set, scales = facetScales)
    }

    #-------------------------------#
    # Counts in the corners         #
    #-------------------------------#
    if (isTRUE(showCounts)) {
      countTable <- .diffCounts(regionTable = plotTable, bySet = facetBySet)

      # The y axis carries the fold change here, so the labels sit on the side they describe
      maPlot <- maPlot +
        ggplot2::geom_text(data = countTable,
                           mapping = ggplot2::aes(x = Inf, y = Inf, label = paste0("n = ", .data$n.up)),
                           inherit.aes = FALSE, hjust = 1.15, vjust = 1.4,
                           size = baseSize / 4.5, colour = statusColours[["up"]]) +
        ggplot2::geom_text(data = countTable,
                           mapping = ggplot2::aes(x = Inf, y = -Inf, label = paste0("n = ", .data$n.down)),
                           inherit.aes = FALSE, hjust = 1.15, vjust = -0.6,
                           size = baseSize / 4.5, colour = statusColours[["down"]])
    }

    return(maPlot)
  } # END function




#' @title plotRegion
#'
#' @description Draws the signal of a single region across the samples. On a tiled object the values are drawn along the coordinates, one line per sample, which shows whether the whole region moved or only part of it. On a non-tiled object the region carries one value per sample and those are drawn as points.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results} or \code{RegionSetDE.resultsList} object. A result carries both the values and the statistics, so nothing else has to be passed.
#' @param region String identifying the region, written as \code{"set|id"} or as the region identifier alone when it is unique across the sets. A \code{GRanges} of length one is accepted as well, in which case the overlapping rows are drawn.
#' @param counts \code{RegionSetDE.counts} object holding the values, when \code{object} carries none. Default: \code{NULL}.
#' @param contrast String with the name of the contrast to annotate with, or its position, when \code{object} holds several of them. Default: \code{NULL}.
#' @param groupBy String with the name of a \code{colData} column driving the colour, e.g. \code{"condition"}. Default: \code{NULL}, one colour per sample.
#' @param assay String with the name of the assay to draw. Default: \code{NULL}, the normalised assay when present, the raw counts otherwise.
#' @param log2Scale Logical value to indicate whether the values must be drawn on a log2 scale. Default: \code{TRUE}.
#' @param summarise Logical value to indicate whether the replicates of a group must be summarised rather than drawn one by one: a mean line with a ribbon along a tiled region, a mean with its spread next to the individual points on a region counted as a single row. Requires \code{groupBy}. Default: \code{FALSE}.
#' @param colours Named character vector with the colours. Default: \code{NULL}.
#' @param pointSize Numeric value with the size of the points, on a non-tiled region. Default: \code{3}.
#' @param rotateX Logical value to indicate whether the labels of the x axis must be angled. Default: \code{TRUE}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}, the region identifier.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}, the statistics of the region when they are available.
#' @param legendPosition String with the position of the legend. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object.
#'
#' @details The values come from the object and nothing is re-read from the BAM or bigWig files, so the resolution of the plot is the resolution of the counting. A region counted as a single row gives a single point per sample, which is the honest picture of what the model saw.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' topRegion <- topRegions(results, n = 1, FDR = 1)$region.id
#'
#' plotRegion(results, region = topRegion, groupBy = "condition")
#'
#' # Summarised to one point per group rather than one per sample
#' plotRegion(results, region = topRegion, groupBy = "condition", summarise = TRUE)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{topRegions}}, \code{\link{testRegions}}, \code{\link{plotTopHeatmap}}
#'
#' @importFrom SummarizedExperiment assay assayNames colData rowData rowRanges
#' @importFrom S4Vectors metadata queryHits
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon geom_errorbar geom_crossbar labs scale_colour_manual scale_fill_manual theme element_blank element_rect element_line
#' @importFrom dplyr filter mutate group_by summarise ungroup
#' @importFrom rlang .data
#' @importFrom stats sd
#' @importFrom methods is
#'
#' @export plotRegion

plotRegion <-
  function(object,
           region,
           counts = NULL,
           contrast = NULL,
           groupBy = NULL,
           assay = NULL,
           log2Scale = TRUE,
           summarise = FALSE,
           colours = NULL,
           pointSize = 3,
           rotateX = TRUE,
           title = NULL,
           subtitle = NULL,
           legendPosition = "right",
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    resolvedObject <- .resolveCounts(object = object, counts = counts, contrast = contrast)
    counts <- resolvedObject$counts
    results <- resolvedObject$results

    if (isTRUE(summarise) & is.null(groupBy)) {
      stop("The 'summarise' parameter needs a 'groupBy' column.", call. = FALSE)
    }

    #-------------------------------#
    # Locate the region             #
    #-------------------------------#
    regionRows <- .locateRegion(counts = counts, region = region)

    #-------------------------------#
    # Extract the values            #
    #-------------------------------#
    assay <- .defaultAssay(counts = counts, assay = assay)

    valueMatrix <- as.matrix(SummarizedExperiment::assay(counts, assay))[regionRows$row.index, , drop = FALSE]
    rowRangesObject <- SummarizedExperiment::rowRanges(counts)[regionRows$row.index]

    plotTable <- data.frame(sample = rep(colnames(counts), each = nrow(valueMatrix)),
                            position = rep((BiocGenerics::start(rowRangesObject) + BiocGenerics::end(rowRangesObject)) / 2,
                                           times = ncol(valueMatrix)),
                            value = as.numeric(valueMatrix),
                            stringsAsFactors = FALSE)

    if (isTRUE(log2Scale)) {
      plotTable <- dplyr::mutate(plotTable, value = log2(.data$value + 1))
    }

    #-------------------------------#
    # Attach the sample metadata    #
    #-------------------------------#
    colTable <- as.data.frame(SummarizedExperiment::colData(counts))
    colTable$sample <- colnames(counts)

    if (!is.null(groupBy)) {
      if (!(groupBy %in% colnames(colTable))) {
        stop(paste0("The column '", groupBy, "' is absent from the colData."), call. = FALSE)
      }
      plotTable$group <- as.character(colTable[[groupBy]])[match(plotTable$sample, colTable$sample)]
    } else {
      plotTable$group <- plotTable$sample
    }

    #-------------------------------#
    # Statistics of the region      #
    #-------------------------------#
    isTiled <- nrow(valueMatrix) > 1
    yLabel <- paste0(if (isTRUE(log2Scale)) {"log<sub>2</sub> "} else {""}, assay)

    if (is.null(subtitle) & !is.null(results)) {
      statisticsRow <- dplyr::filter(resultsTable(results),
                                     paste(.data$region.set, .data$region.id, sep = "|") == unique(regionRows$region.key))
      if (nrow(statisticsRow) == 1) {
        subtitle <- sprintf("*%s*: log<sub>2</sub>FC %.2f, FDR %.1e",
                            results@contrast, statisticsRow$log2FC[1], statisticsRow$FDR[1])
      }
    }

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    summaryTable <- NULL
    if (isTRUE(summarise)) {
      summaryTable <- dplyr::ungroup(dplyr::summarise(dplyr::group_by(plotTable, .data$group, .data$position),
                                                      mean.value = mean(.data$value),
                                                      sd.value = stats::sd(.data$value),
                                                      .groups = "drop"))
      summaryTable$sd.value[is.na(summaryTable$sd.value)] <- 0
    }

    # A ribbon needs a coordinate to run along, which only a tiled region has
    hasRibbon <- isTiled & isTRUE(summarise)

    if (isTRUE(hasRibbon)) {
      # One line per group, with the spread of the replicates around it, when the individual samples clutter the panel
      regionPlot <-
        ggplot2::ggplot(data = summaryTable,
                        mapping = ggplot2::aes(x = .data$position, y = .data$mean.value,
                                               colour = .data$group, fill = .data$group)) +
        ggplot2::geom_ribbon(mapping = ggplot2::aes(ymin = .data$mean.value - .data$sd.value,
                                                    ymax = .data$mean.value + .data$sd.value),
                             alpha = 0.2, colour = NA) +
        ggplot2::geom_line(linewidth = 0.6)

    } else if (isTRUE(isTiled)) {
      regionPlot <-
        ggplot2::ggplot(data = plotTable,
                        mapping = ggplot2::aes(x = .data$position, y = .data$value,
                                               colour = .data$group, group = .data$sample)) +
        ggplot2::geom_line(linewidth = 0.5)

    } else if (isTRUE(summarise)) {
      # The region holds one value per sample, so the summary is a mean with its spread rather than a profile
      regionPlot <-
        ggplot2::ggplot(data = plotTable,
                        mapping = ggplot2::aes(x = .data$group, y = .data$value, colour = .data$group)) +
        ggplot2::geom_errorbar(data = summaryTable,
                               mapping = ggplot2::aes(x = .data$group,
                                                      ymin = .data$mean.value - .data$sd.value,
                                                      ymax = .data$mean.value + .data$sd.value),
                               inherit.aes = FALSE, width = 0.15, linewidth = 0.4, colour = "black") +
        ggplot2::geom_crossbar(data = summaryTable,
                               mapping = ggplot2::aes(x = .data$group, y = .data$mean.value,
                                                      ymin = .data$mean.value, ymax = .data$mean.value),
                               inherit.aes = FALSE, width = 0.35, linewidth = 0.4, colour = "black") +
        ggplot2::geom_point(size = pointSize, stroke = NA)

    } else {
      regionPlot <-
        ggplot2::ggplot(data = plotTable,
                        mapping = ggplot2::aes(x = .data$group, y = .data$value, colour = .data$group)) +
        ggplot2::geom_point(size = pointSize, stroke = NA)
    }

    # Only the summarised plot maps a fill, naming the others would raise an unknown label warning
    labelList <- list(x = if (isTiled) {paste0(as.character(GenomeInfoDb::seqnames(rowRangesObject))[1], " (bp)")} else {""},
                      y = yLabel,
                      colour = if (is.null(groupBy)) {"Sample"} else {groupBy},
                      title = if (is.null(title)) {unique(regionRows$region.key)} else {title},
                      subtitle = subtitle)

    if (isTRUE(hasRibbon)) {
      labelList$fill <- if (is.null(groupBy)) {"Sample"} else {groupBy}
    }

    regionPlot <- regionPlot +
      do.call(what = ggplot2::labs, args = labelList) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize, rotateX = rotateX) +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5, colour = "black"),
                     axis.ticks.x = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "gray"))

    if (!is.null(colours)) {
      regionPlot <- regionPlot + ggplot2::scale_colour_manual(values = colours)

      if (isTRUE(hasRibbon)) {
        regionPlot <- regionPlot + ggplot2::scale_fill_manual(values = colours)
      }
    }

    return(regionPlot)
  } # END function




#' @title .resultsTheme
#'
#' @description Extends the theme shared by the package with the pieces the result plots need: a centred subtitle rendered as markdown, axis labels in black at full size, and optionally angled labels on the x axis. The weight of the subtitle is written out rather than left to the inheritance, which otherwise picks up the bold of the title.
#'
#' @param legendPosition String with the position of the legend. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#' @param rotateX Logical value to indicate whether the labels of the x axis must be angled. Default: \code{FALSE}.
#'
#' @return A \code{ggplot2} theme.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom ggplot2 theme element_text rel
#' @importFrom ggtext element_markdown
#'
#' @keywords internal

.resultsTheme <-
  function(legendPosition = "right",
           baseSize = 12,
           rotateX = FALSE) {

    # The axis labels come out at rel(0.8) by default, which is smaller than the panels of a figure usually want
    axisTextX <- ggplot2::element_text(colour = "black", size = ggplot2::rel(1),
                                       angle = if (isTRUE(rotateX)) {45} else {0},
                                       hjust = if (isTRUE(rotateX)) {1} else {0.5},
                                       vjust = if (isTRUE(rotateX)) {1} else {0.5})

    return(.regionSetTheme(legendPosition = legendPosition, baseSize = baseSize) +
             ggplot2::theme(plot.subtitle = ggtext::element_markdown(face = "plain", hjust = 0.5),
                            axis.text.x = axisTextX,
                            axis.text.y = ggplot2::element_text(colour = "black", size = ggplot2::rel(1))))
  } # END function




#' @title .diffStatusColours
#'
#' @description Returns the palette of the \code{diff.status} column, completed with the defaults when only part of it is given.
#'
#' @param colours Named character vector, or \code{NULL}.
#'
#' @return A named character vector with the \code{down}, \code{null} and \code{up} entries.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.diffStatusColours <-
  function(colours = NULL) {

    defaultColours <- c("down" = "#2166AC", "null" = "grey70", "up" = "#B2182B")

    if (is.null(colours)) {
      return(defaultColours)
    }

    if (is.null(names(colours))) {
      if (length(colours) != 3) {
        stop("The 'colours' parameter must be named, or hold three values for down, null and up.", call. = FALSE)
      }
      names(colours) <- c("down", "null", "up")
    }

    defaultColours[names(colours)] <- colours
    return(defaultColours)
  } # END function




#' @title .diffCounts
#'
#' @description Counts the changing regions of every set, for the annotation written in the corners of the panels.
#'
#' @param regionTable Data.frame with the \code{region.set} and \code{diff.status} columns.
#' @param bySet Logical value indicating whether the counts must be split by region set.
#'
#' @return A data.frame with the \code{n.up} and \code{n.down} columns, carrying \code{region.set} when \code{bySet} is \code{TRUE}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr group_by summarise ungroup
#' @importFrom rlang .data
#'
#' @keywords internal

.diffCounts <-
  function(regionTable,
           bySet = TRUE) {

    # The counts describe the whole table, not the thinned one, otherwise they would depend on maxPoints
    if (isFALSE(bySet)) {
      return(data.frame(n.up = sum(regionTable$diff.status == "up"),
                        n.down = sum(regionTable$diff.status == "down"),
                        stringsAsFactors = FALSE))
    }

    return(as.data.frame(dplyr::ungroup(dplyr::summarise(dplyr::group_by(regionTable, .data$region.set),
                                                         n.up = sum(.data$diff.status == "up"),
                                                         n.down = sum(.data$diff.status == "down"),
                                                         .groups = "drop"))))
  } # END function




#' @title .thinBySet
#'
#' @description Thins a table of regions down to a number of points that a panel can hold, inside each region set when the plot is faceted.
#'
#' @param regionTable Data.frame with a \code{region.set} column.
#' @param maxPoints Numeric value with the number of points kept per panel.
#' @param bySet Logical value indicating whether the thinning must happen inside each set.
#'
#' @return A data.frame with at most \code{maxPoints} rows per panel.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.thinBySet <-
  function(regionTable,
           maxPoints = 20000,
           bySet = TRUE) {

    if (nrow(regionTable) == 0) {
      return(regionTable)
    }

    if (isFALSE(bySet)) {
      return(regionTable[.thinIndex(n = nrow(regionTable), maxPoints = maxPoints), , drop = FALSE])
    }

    # Thinning the pool would empty a small set to make room for a large one, each panel gets its own budget
    thinnedList <-
      lapply(unique(regionTable$region.set),
             function(setName) {
               setTable <- regionTable[regionTable$region.set == setName, , drop = FALSE]
               return(setTable[.thinIndex(n = nrow(setTable), maxPoints = maxPoints), , drop = FALSE])
             })

    return(do.call(what = rbind, args = thinnedList))
  } # END function




#' @title .prepareResultTable
#'
#' @description Pulls the table out of a results object and restricts it to a subset of region sets, relabelling \code{diff.status} when the thresholds differ from the stored ones.
#'
#' @param results \code{RegionSetDE.results} object.
#' @param set Character vector with the region sets to keep, or \code{NULL}.
#' @param FDR Numeric value with the adjusted p-value cut-off, or \code{NULL}.
#' @param log2FC Numeric value with the log2 fold change cut-off, or \code{NULL}.
#'
#' @return A data.frame ready to be plotted.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter mutate case_when
#' @importFrom rlang .data
#'
#' @keywords internal

.prepareResultTable <-
  function(results,
           set = NULL,
           FDR = NULL,
           log2FC = NULL) {

    plotTable <- results@results

    if (!is.null(set)) {
      absentSets <- setdiff(set, unique(plotTable$region.set))
      if (length(absentSets) > 0) {
        stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }
      plotTable <- dplyr::filter(plotTable, .data$region.set %in% set)
    }

    # A threshold given here only changes the labels, the p-values and the correction stay the ones of the test
    if (!is.null(FDR) | !is.null(log2FC)) {
      FDRthreshold <- if (is.null(FDR)) {results@thresholds$FDR} else {FDR}
      log2FCthreshold <- if (is.null(log2FC)) {results@thresholds$log2FC} else {log2FC}

      plotTable <- dplyr::mutate(plotTable,
                                 diff.status = dplyr::case_when(.data$FDR < FDRthreshold & .data$log2FC > log2FCthreshold ~ "up",
                                                                .data$FDR < FDRthreshold & .data$log2FC < (-log2FCthreshold) ~ "down",
                                                                TRUE ~ "null"))
      plotTable$diff.status <- factor(plotTable$diff.status, levels = c("down", "null", "up"))
    }

    if (nrow(plotTable) == 0) {
      stop("The selection leaves no region to draw.", call. = FALSE)
    }

    return(plotTable)
  } # END function




#' @title .resolveCounts
#'
#' @description Returns the counts and, when there are any, the statistics a plotting function has to work with, whatever kind of object it was handed.
#'
#' @param object Any object of the package holding counts: \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results}, \code{RegionSetDE.setResults}, or one of the two list classes.
#' @param counts \code{RegionSetDE.counts} object overriding the one carried by \code{object}. Default: \code{NULL}.
#' @param contrast String with the name of a contrast, or its position. Default: \code{NULL}.
#'
#' @return A list with the \code{counts} and the \code{results} elements, the second one \code{NULL} when no statistics are available.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods is
#'
#' @keywords internal

.resolveCounts <-
  function(object,
           counts = NULL,
           contrast = NULL) {

    if (methods::is(object, "RegionSetDE.counts")) {
      return(list(counts = object, results = NULL))
    }

    if (methods::is(object, "RegionSetDE.fit")) {
      return(list(counts = object@counts, results = NULL))
    }

    resultClasses <- c("RegionSetDE.results", "RegionSetDE.resultsList",
                       "RegionSetDE.setResults", "RegionSetDE.setResultsList")

    if (any(vapply(resultClasses, function(x) {methods::is(object, x)}, logical(1)))) {
      results <- .pickResults(results = object, contrast = contrast)

      # A result built with carryCounts = FALSE has the statistics but not the values behind them
      resolvedCounts <- if (!is.null(counts)) {counts} else {results@counts}

      if (ncol(resolvedCounts) == 0) {
        stop("The result carries no counts, pass them through 'counts' or run the test with carryCounts = TRUE.", call. = FALSE)
      }

      return(list(counts = resolvedCounts, results = results))
    }

    stop("The object must be a counts, a fit or a results object of the package.", call. = FALSE)
  } # END function




#' @title .locateRegion
#'
#' @description Finds the rows of a counts object belonging to one region.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param region String with the region identifier, written as \code{"set|id"} or as the identifier alone, or a \code{GRanges} of length one.
#'
#' @return A data.frame with the rows of the region and their positions in the object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment rowData rowRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors queryHits
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @keywords internal

.locateRegion <-
  function(counts,
           region) {

    rowTable <- as.data.frame(SummarizedExperiment::rowData(counts))
    rowTable$row.index <- seq_len(nrow(counts))
    rowTable$region.key <- paste(rowTable$region.set, rowTable$region.id, sep = "|")

    if (methods::is(region, "GRanges")) {
      if (length(region) != 1) {
        stop("The 'region' parameter must be a GRanges of length one.", call. = FALSE)
      }
      overlapIndex <- S4Vectors::queryHits(GenomicRanges::findOverlaps(query = SummarizedExperiment::rowRanges(counts),
                                                                       subject = region))
      regionRows <- rowTable[overlapIndex, , drop = FALSE]
      regionLabel <- paste0(as.character(GenomeInfoDb::seqnames(region)), ":",
                            BiocGenerics::start(region), "-", BiocGenerics::end(region))

    } else {
      regionRows <- dplyr::filter(rowTable, .data$region.key == region | .data$region.id == region)
      regionLabel <- region
    }

    if (nrow(regionRows) == 0) {
      stop(paste0("No row matches '", regionLabel, "'."), call. = FALSE)
    }

    # An identifier shared by several sets points at two different rows, the set has to be named
    if (length(unique(regionRows$region.key)) > 1) {
      stop(paste0("The identifier matches several region sets (",
                  paste(unique(regionRows$region.set), collapse = ", "), "), write it as 'set|id'."), call. = FALSE)
    }

    return(regionRows)
  } # END function




#' @title .defaultAssay
#'
#' @description Picks the assay a plot should draw, preferring the normalised values when the object holds them.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param assay String with the name of an assay, or \code{NULL}.
#'
#' @return A string with the name of the assay.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment assayNames
#' @importFrom S4Vectors metadata
#'
#' @keywords internal

.defaultAssay <-
  function(counts,
           assay = NULL) {

    if (is.null(assay)) {
      normalizationInfo <- S4Vectors::metadata(counts)$normalization
      assay <- if (!is.null(normalizationInfo)) {normalizationInfo$normalized.assay} else {"counts"}
    }

    if (!(assay %in% SummarizedExperiment::assayNames(counts))) {
      stop(paste0("The assay '", assay, "' is absent from the object."), call. = FALSE)
    }

    return(assay)
  } # END function
