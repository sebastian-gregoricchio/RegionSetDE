# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

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
#' @param log2FC Numeric value with the absolute log2 fold change cut-off used to label the points, drawn as two dashed vertical lines at \code{-log2FC} and \code{log2FC} when above zero. Default: \code{NULL}, the threshold stored in the object, which is the \code{log2FC} given to \code{\link{testRegions}}.
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
        stop("The column '", labelColumn, "' is absent from the results table.", call. = FALSE)
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
#' @param log2FC Numeric value with the absolute log2 fold change cut-off used to label the points, drawn as two dashed horizontal lines when above zero. Default: \code{NULL}, the threshold stored in the object.
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
      # DESeq2 reports the mean of the normalised counts, the other engines log2 counts per million
      ggplot2::labs(x = if (identical(results@engine, "deseq2")) {"Average signal (log<sub>2</sub> normalised count + 1)"} else {"Average signal (log<sub>2</sub> CPM)"},
                    y = paste0("log<sub>2</sub> fold change (", results@contrast, ")"),
                    colour = if (colourBy == "diff.status") {"Status"} else {"Region set"},
                    title = if (is.null(title)) {results@contrast} else {title},
                    subtitle = subtitle) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = max(c(pointSize, 3)), alpha = 1))) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize)

    # The fold change cut-off of the up and down labels, the same one the volcano draws as vertical lines
    log2FCthreshold <- if (is.null(log2FC)) {results@thresholds$log2FC} else {log2FC}

    if (isTRUE(log2FCthreshold > 0)) {
      maPlot <- maPlot +
        ggplot2::geom_hline(yintercept = c(-log2FCthreshold, log2FCthreshold),
                            linetype = "dashed", linewidth = 0.3, colour = "black")
    }

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
#' @description Draws the signal of a single region across the samples. On a tiled object the values are drawn along the coordinates, one line per sample, which shows whether the whole region moved or only part of it. On a non-tiled object the region carries one value per sample and those are drawn as points, with brackets between the groups when \code{pairwiseTest} is set.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results} or \code{RegionSetDE.resultsList} object. A result carries both the values and the statistics, so nothing else has to be passed.
#' @param region String identifying the region, written as \code{"set|id"} or as the region identifier alone when it is unique across the sets. A \code{GRanges} of length one is accepted as well, in which case the overlapping rows are drawn.
#' @param counts \code{RegionSetDE.counts} object holding the values, when \code{object} carries none. Default: \code{NULL}.
#' @param contrast String with the name of the contrast to annotate with, or its position, when \code{object} holds several of them. With \code{pairwiseTest = "model"} it keeps the bracket of that contrast only. Default: \code{NULL}, which with \code{pairwiseTest = "model"} draws one bracket per contrast of the object.
#' @param groupBy String with the name of a \code{colData} column driving the colour, e.g. \code{"condition"}. Its levels are also the groups compared by \code{pairwiseTest}. Default: \code{NULL}, one colour per sample.
#' @param assay String with the name of the assay to draw. Default: \code{NULL}, the normalised assay when present, the raw counts otherwise.
#' @param log2Scale Logical value to indicate whether the values must be drawn on a log2 scale. Default: \code{TRUE}.
#' @param summarise Logical value to indicate whether the replicates of a group must be summarised rather than drawn one by one: a mean line with a ribbon along a tiled region, a mean with its spread next to the individual points on a region counted as a single row. Requires \code{groupBy}. Default: \code{FALSE}.
#' @param pairwiseTest String with what the brackets between the groups of \code{groupBy} report, one among \code{"t.test"}, \code{"wilcox.test"} and \code{"model"}. The first two test the plotted values, with the Welch t-test or the Wilcoxon rank-sum test, or with their paired versions when \code{pairBy} is given. \code{"model"} writes instead the log2 fold change and the FDR of the fitted contrasts, and needs a results object. Only for a region counted as a single row. Default: \code{NULL}, no brackets.
#' @param pairBy String with the name of a \code{colData} column matching the samples of two groups, e.g. \code{"replicate"}, for a paired test. A value can appear once per group, and a sample whose value has no partner in the other group is left out of that comparison. Default: \code{NULL}, unpaired tests.
#' @param comparisons List of character vectors of length two, naming the pairs of groups joined by a bracket. Default: \code{NULL}, every pair of groups, or with \code{pairwiseTest = "model"} every pair tested by a contrast.
#' @param pAdjustMethod String with the correction applied to the p-values of the brackets, one of \code{stats::p.adjust.methods}. It covers the comparisons drawn in the plot, not the regions tested, and it is not used with \code{pairwiseTest = "model"}. Default: \code{"none"}.
#' @param pLabel String indicating how the values are written on the brackets: \code{"p.value"} for the number, \code{"stars"} for the significance symbols (\code{ns} above 0.05, then \code{*}, \code{**}, \code{***} and \code{****} up to 0.05, 0.01, 0.001 and 0.0001). Default: \code{"p.value"}.
#' @param pDecimals Numeric value with the number of decimals of the values written on the brackets. Values below 0.1 are written in scientific notation, with the exponent as a superscript. Default: \code{2}.
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
#' The brackets drawn by \code{"t.test"} and \code{"wilcox.test"} answer a smaller question than the model does. The test runs on the values as they are drawn, on the log2 scale when \code{log2Scale = TRUE}, and it knows nothing of the offsets, of the moderated dispersion, of the covariates in the design or of the thousands of regions tested next to this one. Its p-value and the FDR in the subtitle will often disagree. The region is also usually picked from \code{\link{topRegions}} on the same samples, so the test describes the region and does not confirm the result. Where it helps is on counts and fit objects, which carry no statistics of their own, and between groups that no contrast compared. On raw counts the function warns, since a difference in sequencing depth between the groups is read as a difference in signal.
#'
#' Few replicates limit these tests more than they limit the model. With three samples per group the exact Wilcoxon test cannot go below p = 0.1, and with three pairs the signed-rank test cannot go below 0.25; a message says so when the groups are that small. With \code{pairBy} the samples are matched on the value of that column, whatever their order in the object. A test paired by replicate asks what \code{design = ~ replicate + condition} asks in \code{\link{fitRegions}}, and the model is the place for that pairing when the numbers are meant to be reported.
#'
#' \code{"model"} writes on each bracket the log2 fold change and the FDR that \code{\link{testRegions}} gave to the region, one bracket per contrast comparing two levels of \code{groupBy}. A contrast qualifies when it was written as \code{c("column", "groupA", "groupB")}, or when it reduces to that difference, which the result records in its \code{contrast.groups} slot. The caption says which test, or which fit, the brackets come from.
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
#' # Fold change and FDR of the fitted contrast on the bracket
#' plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "model")
#'
#' # Welch t-test on the plotted values
#' plotRegion(results, region = topRegion, groupBy = "condition", pairwiseTest = "t.test")
#'
#' # Paired t-test; the example has no matched replicates, so the pairs are made up to show the call
#' pairedCounts <- resultCounts(results)
#' pairedCounts$pair <- c("p1", "p2", "p1", "p2")
#'
#' plotRegion(pairedCounts, region = topRegion, groupBy = "condition",
#'            pairwiseTest = "t.test", pairBy = "pair")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{topRegions}}, \code{\link{testRegions}}, \code{\link{plotTopHeatmap}}, \code{\link{plotSetSignal}}
#'
#' @importFrom SummarizedExperiment assay assayNames colData rowData rowRanges
#' @importFrom S4Vectors metadata queryHits
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon geom_errorbar geom_crossbar geom_segment expand_limits labs scale_colour_manual scale_fill_manual theme element_blank element_rect element_line
#' @importFrom ggtext geom_richtext element_markdown
#' @importFrom grid unit
#' @importFrom dplyr filter mutate group_by summarise ungroup
#' @importFrom rlang .data
#' @importFrom stats sd p.adjust.methods
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
           pairwiseTest = NULL,
           pairBy = NULL,
           comparisons = NULL,
           pAdjustMethod = "none",
           pLabel = "p.value",
           pDecimals = 2,
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
    if (!is.null(pairwiseTest)) {
      if (length(pairwiseTest) != 1 || !(pairwiseTest %in% c("t.test", "wilcox.test", "model"))) {
        stop("The 'pairwiseTest' parameter must be one of 't.test', 'wilcox.test', 'model'.", call. = FALSE)
      }

      if (is.null(groupBy)) {
        stop("The 'pairwiseTest' parameter needs a 'groupBy' column.", call. = FALSE)
      }

      if (length(pAdjustMethod) != 1 || !(pAdjustMethod %in% stats::p.adjust.methods)) {
        stop("The 'pAdjustMethod' parameter must be one of ", paste(stats::p.adjust.methods, collapse = ", "), ".", call. = FALSE)
      }

      if (length(pLabel) != 1 || !(pLabel %in% c("p.value", "stars"))) {
        stop("The 'pLabel' parameter must be either 'p.value' or 'stars'.", call. = FALSE)
      }

      if (!is.numeric(pDecimals) || length(pDecimals) != 1 || is.na(pDecimals) || pDecimals < 0) {
        stop("The 'pDecimals' parameter must be a single number, zero or above.", call. = FALSE)
      }

    } else if (!is.null(pairBy) | !is.null(comparisons)) {
      warning("No 'pairwiseTest' was given, the 'pairBy' and 'comparisons' parameters have been ignored.", call. = FALSE)
    }

    # The model brackets read every contrast of a list, so a single one is picked only when it is named
    modelResults <- NULL
    if (identical(pairwiseTest, "model")) {
      modelResults <- .modelResultsList(object = object, contrast = contrast)
    }

    resolvedObject <- .resolveCounts(object = object,
                                     counts = counts,
                                     contrast = if (length(modelResults) > 1) {1} else {contrast})
    counts <- resolvedObject$counts

    # With several contrasts on the brackets, none of them can speak for the subtitle alone
    results <- if (length(modelResults) > 1) {NULL} else {resolvedObject$results}

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
        stop("The column '", groupBy, "' is absent from the colData.", call. = FALSE)
      }
      plotTable$group <- as.character(colTable[[groupBy]])[match(plotTable$sample, colTable$sample)]
    } else {
      plotTable$group <- plotTable$sample
    }

    # The order is fixed here, so that the brackets land where the axis puts the groups
    plotTable$group <- factor(plotTable$group, levels = sort(unique(plotTable$group)))

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
    # Pairwise comparisons          #
    #-------------------------------#
    comparisonTable <- NULL
    captionText <- NULL

    if (!is.null(pairwiseTest)) {
      # A bracket joins two positions of the axis, and along a tiled region the axis holds the coordinates
      if (isTRUE(isTiled)) {
        stop("The pairwise comparisons need a region counted as a single row, on a tiled region the x axis holds the coordinates.", call. = FALSE)
      }

      if (pairwiseTest == "model") {
        if (!is.null(pairBy)) {
          warning("The model brackets follow the design of the fit, the 'pairBy' parameter has been ignored.", call. = FALSE)
        }

        if (pAdjustMethod != "none") {
          warning("The FDR of the model is already corrected over the regions, the 'pAdjustMethod' parameter has been ignored.", call. = FALSE)
        }

        comparisonTable <- .modelBrackets(resultsList = modelResults,
                                          regionKey = unique(regionRows$region.key),
                                          groupBy = groupBy,
                                          groupLevels = levels(plotTable$group),
                                          comparisons = comparisons)

        comparisonTable$label <- .bracketLabels(pValues = comparisonTable$FDR,
                                                prefix = "FDR",
                                                pLabel = pLabel,
                                                pDecimals = pDecimals,
                                                log2FC = comparisonTable$log2FC)

        engineNames <- unique(vapply(modelResults, function(x) {x@engine}, character(1)))
        captionText <- paste0("log<sub>2</sub>FC and FDR from the ", paste(engineNames, collapse = " and "), " fit")

      } else {
        # Coverage declared as not count-like is taken as already scaled, while counts need the normalisation first
        objectMetadata <- S4Vectors::metadata(counts)
        if (!isFALSE(objectMetadata$count.like) & !identical(assay, objectMetadata$normalization$normalized.assay)) {
          warning("The test runs on the '", assay, "' assay, which is not normalised: a difference in sequencing depth between the groups ",
                  "is read as a difference in signal. Run normalizeCounts() first, or pick the normalised values through 'assay'.", call. = FALSE)
        }

        comparisonTable <- .pairwiseTests(plotTable = plotTable,
                                          colTable = colTable,
                                          test = pairwiseTest,
                                          pairBy = pairBy,
                                          comparisons = comparisons,
                                          groupLevels = levels(plotTable$group),
                                          pAdjustMethod = pAdjustMethod)

        # One comparison leaves nothing to correct for, the label keeps the plain p-value then
        isAdjusted <- pAdjustMethod != "none" & nrow(comparisonTable) > 1

        comparisonTable$label <- .bracketLabels(pValues = if (isTRUE(isAdjusted)) {comparisonTable$p.adjusted} else {comparisonTable$p.value},
                                                prefix = if (isTRUE(isAdjusted)) {"p<sub>adj</sub>"} else {"p"},
                                                pLabel = pLabel,
                                                pDecimals = pDecimals)

        testName <- if (pairwiseTest == "t.test") {
          if (is.null(pairBy)) {"Welch t-test"} else {"Paired t-test"}
        } else {
          if (is.null(pairBy)) {"Wilcoxon rank-sum test"} else {"Wilcoxon signed-rank test"}
        }

        captionText <- paste0(testName, " on ", yLabel,
                              if (!is.null(pairBy)) {paste0(", paired by ", pairBy)} else {""},
                              if (isTRUE(isAdjusted)) {paste0(", ", pAdjustMethod, " adjustment over ", nrow(comparisonTable), " comparisons")} else {""},
                              "<br>Not corrected for the number of regions tested")
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

    #-------------------------------#
    # Brackets over the groups      #
    #-------------------------------#
    if (!is.null(comparisonTable)) {
      # The error bars of a summarised plot can reach above the points, the brackets start above both
      valueRange <- range(c(plotTable$value,
                            if (!is.null(summaryTable)) {summaryTable$mean.value + summaryTable$sd.value} else {NULL}),
                          na.rm = TRUE)

      comparisonTable <- .placeBrackets(bracketTable = comparisonTable,
                                        groupLevels = levels(plotTable$group),
                                        valueRange = valueRange,
                                        labelLines = if (pairwiseTest == "model") {2} else {1})

      regionPlot <- regionPlot +
        ggplot2::geom_segment(data = comparisonTable,
                              mapping = ggplot2::aes(x = .data$x.start, xend = .data$x.end,
                                                     y = .data$y.bracket, yend = .data$y.bracket),
                              inherit.aes = FALSE, linewidth = 0.4, colour = "grey20") +
        ggplot2::geom_segment(data = comparisonTable,
                              mapping = ggplot2::aes(x = .data$x.start, xend = .data$x.start,
                                                     y = .data$y.bracket, yend = .data$y.bracket - .data$y.tick),
                              inherit.aes = FALSE, linewidth = 0.4, colour = "grey20") +
        ggplot2::geom_segment(data = comparisonTable,
                              mapping = ggplot2::aes(x = .data$x.end, xend = .data$x.end,
                                                     y = .data$y.bracket, yend = .data$y.bracket - .data$y.tick),
                              inherit.aes = FALSE, linewidth = 0.4, colour = "grey20") +
        ggtext::geom_richtext(data = comparisonTable,
                              mapping = ggplot2::aes(x = (.data$x.start + .data$x.end) / 2,
                                                     y = .data$y.bracket, label = .data$label),
                              inherit.aes = FALSE, vjust = 0, size = baseSize / 4.5,
                              colour = "grey20", fill = NA, label.color = NA,
                              label.padding = grid::unit(rep(1, 4), "pt")) +
        ggplot2::expand_limits(y = comparisonTable$y.top[1])
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

    if (!is.null(captionText)) {
      labelList$caption <- captionText
    }

    regionPlot <- regionPlot +
      do.call(what = ggplot2::labs, args = labelList) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize, rotateX = rotateX) +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5, colour = "black"),
                     axis.ticks.x = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "gray"),
                     plot.caption = ggtext::element_markdown(colour = "gray30"))

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
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
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
      stop("No row matches '", regionLabel, "'.", call. = FALSE)
    }

    # An identifier shared by several sets points at two different rows, the set has to be named
    if (length(unique(regionRows$region.key)) > 1) {
      stop("The identifier matches several region sets (",
           paste(unique(regionRows$region.set), collapse = ", "), "), write it as 'set|id'.", call. = FALSE)
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
      stop("The assay '", assay, "' is absent from the object.", call. = FALSE)
    }

    return(assay)
  } # END function




#' @title .modelResultsList
#'
#' @description Collects the contrasts whose statistics the model brackets of \code{\link{plotRegion}} are written from: the single contrast of a result, every contrast of a list, or the one named in \code{contrast}.
#'
#' @param object Object handed to \code{\link{plotRegion}}.
#' @param contrast String with the name of a contrast, or its position, or \code{NULL}. Default: \code{NULL}.
#'
#' @return A named list of \code{RegionSetDE.results} objects.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats setNames
#' @importFrom methods is
#'
#' @keywords internal

.modelResultsList <-
  function(object,
           contrast = NULL) {

    if (methods::is(object, "RegionSetDE.results")) {
      return(stats::setNames(list(object), object@contrast))
    }

    if (!methods::is(object, "RegionSetDE.resultsList")) {
      stop("The model brackets need the statistics of a region level test, pass the object returned by testRegions().", call. = FALSE)
    }

    if (is.null(contrast)) {
      return(object@results)
    }

    pickedResults <- .pickResults(results = object, contrast = contrast)
    contrastName <- if (is.numeric(contrast)) {names(object@results)[as.integer(contrast[1])]} else {contrast[1]}

    return(stats::setNames(list(pickedResults), contrastName))
  } # END function




#' @title .modelBrackets
#'
#' @description Reads, for every contrast comparing two levels of the grouping column, the log2 fold change and the FDR the fit gave to one region, which is what the model brackets of \code{\link{plotRegion}} write.
#'
#' @param resultsList Named list of \code{RegionSetDE.results} objects.
#' @param regionKey String with the region, written as \code{"set|id"}.
#' @param groupBy String with the column grouping the samples on the axis.
#' @param groupLevels Character vector with the groups, in the order of the axis.
#' @param comparisons List of character vectors of length two, or \code{NULL}. Default: \code{NULL}.
#'
#' @return A data.frame with one row per bracket and the \code{contrast}, \code{group1}, \code{group2}, \code{log2FC} and \code{FDR} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom methods .hasSlot
#'
#' @keywords internal

.modelBrackets <-
  function(resultsList,
           regionKey,
           groupBy,
           groupLevels,
           comparisons = NULL) {

    #-------------------------------#
    # One bracket per contrast      #
    #-------------------------------#
    bracketList <-
      lapply(names(resultsList),
             function(contrastName) {
               contrastResults <- resultsList[[contrastName]]

               # An object saved before the slot existed carries no groups, and is treated as such
               contrastGroups <- if (methods::.hasSlot(contrastResults, "contrast.groups")) {contrastResults@contrast.groups} else {list()}

               # Only a difference between two levels of the grouping column describes two positions of the axis
               if (!identical(contrastGroups$column, groupBy) || !all(contrastGroups$groups %in% groupLevels)) {
                 return(NULL)
               }

               statisticsRow <- dplyr::filter(resultsTable(contrastResults),
                                              paste(.data$region.set, .data$region.id, sep = "|") == regionKey)

               if (nrow(statisticsRow) != 1) {
                 return(NULL)
               }

               return(data.frame(contrast = contrastName,
                                 group1 = contrastGroups$groups[1],
                                 group2 = contrastGroups$groups[2],
                                 log2FC = statisticsRow$log2FC[1],
                                 FDR = statisticsRow$FDR[1],
                                 stringsAsFactors = FALSE))
             })

    missedContrasts <- names(resultsList)[vapply(bracketList, is.null, logical(1))]
    bracketTable <- do.call(what = rbind, args = bracketList)

    if (is.null(bracketTable)) {
      stop("None of the contrasts compares two levels of '", groupBy, "' on this region. The model brackets need a contrast written as c('",
           groupBy, "', 'groupA', 'groupB'), or one that reduces to that difference.", call. = FALSE)
    }

    if (length(missedContrasts) > 0) {
      message("No bracket for ", paste(missedContrasts, collapse = ", "), ": the contrast does not compare two levels of '",
              groupBy, "', or holds no statistics for this region.")
    }

    #-------------------------------#
    # Pairs asked for               #
    #-------------------------------#
    if (!is.null(comparisons)) {
      comparisons <- .checkComparisons(comparisons = comparisons, groupLevels = groupLevels)

      keptRows <- vapply(seq_len(nrow(bracketTable)),
                         function(i) {
                           any(vapply(comparisons,
                                      function(x) {setequal(x, c(bracketTable$group1[i], bracketTable$group2[i]))},
                                      logical(1)))
                         },
                         logical(1))

      bracketTable <- bracketTable[keptRows, , drop = FALSE]

      if (nrow(bracketTable) == 0) {
        stop("None of the requested comparisons matches a pair of groups tested by the contrasts.", call. = FALSE)
      }
    }

    rownames(bracketTable) <- NULL
    return(bracketTable)
  } # END function




#' @title .pairwiseTests
#'
#' @description Runs a two-sample test between each pair of groups on the values of one region, as they are drawn, pairing the samples through a column of the metadata when asked to.
#'
#' @param plotTable Data.frame with one row per sample and the \code{sample}, \code{value} and \code{group} columns.
#' @param colTable Data.frame with the sample metadata and a \code{sample} column.
#' @param test String with the test, either \code{"t.test"} or \code{"wilcox.test"}. Default: \code{"t.test"}.
#' @param pairBy String with the column pairing the samples, or \code{NULL} for unpaired tests. Default: \code{NULL}.
#' @param comparisons List of character vectors of length two, or \code{NULL} for every pair of groups. Default: \code{NULL}.
#' @param groupLevels Character vector with the groups, in the order of the axis.
#' @param pAdjustMethod String with the correction applied across the comparisons. Default: \code{"none"}.
#'
#' @return A data.frame with one row per comparison and the \code{group1}, \code{group2}, \code{n.first}, \code{n.second}, \code{p.value} and \code{p.adjusted} columns. With \code{pairBy} the two sizes are the number of pairs.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom stats t.test wilcox.test p.adjust
#'
#' @keywords internal

.pairwiseTests <-
  function(plotTable,
           colTable,
           test = "t.test",
           pairBy = NULL,
           comparisons = NULL,
           groupLevels,
           pAdjustMethod = "none") {

    comparisons <- .checkComparisons(comparisons = comparisons, groupLevels = groupLevels)

    #-------------------------------#
    # Values and pairing column     #
    #-------------------------------#
    sampleTable <- data.frame(sample = plotTable$sample,
                              value = plotTable$value,
                              group = as.character(plotTable$group),
                              stringsAsFactors = FALSE)

    if (!is.null(pairBy)) {
      if (!(pairBy %in% colnames(colTable))) {
        stop("The column '", pairBy, "' is absent from the colData.", call. = FALSE)
      }
      sampleTable$pair <- as.character(colTable[[pairBy]])[match(sampleTable$sample, colTable$sample)]
    }

    #-------------------------------#
    # One test per pair of groups   #
    #-------------------------------#
    testList <-
      lapply(comparisons,
             function(groupPair) {
               firstTable <- dplyr::filter(sampleTable, .data$group == groupPair[1])
               secondTable <- dplyr::filter(sampleTable, .data$group == groupPair[2])

               if (is.null(pairBy)) {
                 firstValues <- firstTable$value
                 secondValues <- secondTable$value

                 if (test == "t.test" & (length(firstValues) < 2 | length(secondValues) < 2)) {
                   stop("The t-test between '", groupPair[1], "' and '", groupPair[2],
                        "' needs at least two samples per group, leave the pair out through 'comparisons'.", call. = FALSE)
                 }

               } else {
                 # A value repeated inside a group offers two partners for one sample, and nothing to choose between them
                 for (groupTable in list(firstTable, secondTable)) {
                   repeatedValues <- unique(groupTable$pair[duplicated(groupTable$pair) & !is.na(groupTable$pair)])
                   if (length(repeatedValues) > 0) {
                     stop("The value '", repeatedValues[1], "' of '", pairBy, "' appears more than once in the group '",
                          groupTable$group[1], "', the samples cannot be paired.", call. = FALSE)
                   }
                 }

                 sharedPairs <- sort(intersect(firstTable$pair[!is.na(firstTable$pair)],
                                               secondTable$pair[!is.na(secondTable$pair)]))

                 if (length(sharedPairs) < 2) {
                   stop("Fewer than two pairs of samples are matched by '", pairBy, "' between '", groupPair[1], "' and '", groupPair[2],
                        "'. Leave the pair out through 'comparisons', or run the test unpaired.", call. = FALSE)
                 }

                 unmatchedSamples <- c(firstTable$sample[!(firstTable$pair %in% sharedPairs)],
                                       secondTable$sample[!(secondTable$pair %in% sharedPairs)])

                 if (length(unmatchedSamples) > 0) {
                   warning("No partner in '", pairBy, "' for ", paste(unmatchedSamples, collapse = ", "),
                           ", left out of the comparison between '", groupPair[1], "' and '", groupPair[2], "'.", call. = FALSE)
                 }

                 firstValues <- firstTable$value[match(sharedPairs, firstTable$pair)]
                 secondValues <- secondTable$value[match(sharedPairs, secondTable$pair)]
               }

               # The exact Wilcoxon distribution needs untied values, with ties the normal approximation is used
               exactTest <- NULL

               if (test == "wilcox.test") {
                 if (is.null(pairBy)) {
                   hasTies <- anyDuplicated(c(firstValues, secondValues)) > 0
                 } else {
                   valueDifferences <- firstValues - secondValues
                   hasTies <- any(valueDifferences == 0) | anyDuplicated(abs(valueDifferences)) > 0
                 }

                 if (isTRUE(hasTies)) {
                   exactTest <- FALSE
                   message("Tied values between '", groupPair[1], "' and '", groupPair[2],
                           "', the Wilcoxon p-value comes from the normal approximation.")

                 } else if (is.null(pairBy)) {
                   smallestP <- min(1, 2 / choose(length(firstValues) + length(secondValues), length(firstValues)))
                   if (smallestP >= 0.05) {
                     message("With ", length(firstValues), " and ", length(secondValues), " samples, the exact Wilcoxon test between '",
                             groupPair[1], "' and '", groupPair[2], "' cannot go below p = ", signif(smallestP, 2), ".")
                   }

                 } else {
                   smallestP <- min(1, 2 / 2^length(firstValues))
                   if (smallestP >= 0.05) {
                     message("With ", length(firstValues), " pairs, the exact Wilcoxon test between '",
                             groupPair[1], "' and '", groupPair[2], "' cannot go below p = ", signif(smallestP, 2), ".")
                   }
                 }
               }

               # Constant values leave the test without a variance, the bracket then reads NA
               testOutcome <-
                 tryCatch(expr = {
                   testResult <- if (test == "t.test") {
                     stats::t.test(x = firstValues, y = secondValues, paired = !is.null(pairBy))
                   } else {
                     stats::wilcox.test(x = firstValues, y = secondValues, paired = !is.null(pairBy), exact = exactTest)
                   }
                   list(p.value = testResult$p.value, error = NULL)
                 },
                 error = function(e) {
                   return(list(p.value = NA_real_, error = conditionMessage(e)))
                 })

               if (is.na(testOutcome$p.value)) {
                 warning("No p-value between '", groupPair[1], "' and '", groupPair[2], "' (",
                         if (is.null(testOutcome$error)) {"the test returned NA"} else {testOutcome$error},
                         "), the bracket reads NA.", call. = FALSE)
               }

               return(data.frame(group1 = groupPair[1],
                                 group2 = groupPair[2],
                                 n.first = length(firstValues),
                                 n.second = length(secondValues),
                                 p.value = as.numeric(testOutcome$p.value),
                                 stringsAsFactors = FALSE))
             })

    testTable <- do.call(what = rbind, args = testList)

    # The correction covers the brackets of this plot, the regions tested elsewhere do not enter it
    testTable$p.adjusted <- stats::p.adjust(testTable$p.value, method = pAdjustMethod)

    return(testTable)
  } # END function




#' @title .checkComparisons
#'
#' @description Checks the pairs of groups asked for through \code{comparisons}, or builds every pair when none is given.
#'
#' @param comparisons List of character vectors of length two, or \code{NULL}. Default: \code{NULL}.
#' @param groupLevels Character vector with the groups, in the order of the axis.
#'
#' @return A list of character vectors of length two, each pair appearing once.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom utils combn
#'
#' @keywords internal

.checkComparisons <-
  function(comparisons = NULL,
           groupLevels) {

    if (length(groupLevels) < 2) {
      stop("The 'groupBy' column holds a single group, there is nothing to compare.", call. = FALSE)
    }

    if (is.null(comparisons)) {
      return(utils::combn(groupLevels, 2, simplify = FALSE))
    }

    if (!is.list(comparisons) || any(vapply(comparisons, length, integer(1)) != 2)) {
      stop("The 'comparisons' parameter must be a list of character vectors of length two.", call. = FALSE)
    }

    comparisons <- lapply(comparisons, as.character)

    absentGroups <- setdiff(unlist(comparisons), groupLevels)
    if (length(absentGroups) > 0) {
      stop("The following groups are absent from the object: ", paste(absentGroups, collapse = ", "), ".", call. = FALSE)
    }

    if (any(vapply(comparisons, function(x) {x[1] == x[2]}, logical(1)))) {
      stop("Each element of 'comparisons' must name two different groups.", call. = FALSE)
    }

    # The same pair written twice, in either order, would draw two brackets carrying the same number
    pairKeys <- vapply(comparisons, function(x) {paste(sort(x), collapse = "\t")}, character(1))

    return(comparisons[!duplicated(pairKeys)])
  } # END function




#' @title .placeBrackets
#'
#' @description Gives each bracket its horizontal span on the axis and a height above the data. A bracket shares a level with the others unless the two would overlap, in which case it moves up one level.
#'
#' @param bracketTable Data.frame with the \code{group1} and \code{group2} columns.
#' @param groupLevels Character vector with the groups, in the order of the axis.
#' @param valueRange Numeric vector with the lowest and the highest value drawn.
#' @param labelLines Numeric value with the number of lines of each label. Default: \code{1}.
#'
#' @return The input data.frame with the \code{x.start}, \code{x.end}, \code{y.bracket}, \code{y.tick} and \code{y.top} columns added.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.placeBrackets <-
  function(bracketTable,
           groupLevels,
           valueRange,
           labelLines = 1) {

    firstPositions <- match(bracketTable$group1, groupLevels)
    secondPositions <- match(bracketTable$group2, groupLevels)

    bracketTable$x.start <- pmin(firstPositions, secondPositions)
    bracketTable$x.end <- pmax(firstPositions, secondPositions)

    # Short brackets sit low, and a bracket goes up a level only when it would overlap one already placed
    bracketLevels <- rep(0, nrow(bracketTable))

    for (i in order(bracketTable$x.end - bracketTable$x.start, bracketTable$x.start)) {
      candidateLevel <- 1
      while (any(bracketLevels == candidateLevel &
                 bracketTable$x.start <= bracketTable$x.end[i] &
                 bracketTable$x.end >= bracketTable$x.start[i])) {
        candidateLevel <- candidateLevel + 1
      }
      bracketLevels[i] <- candidateLevel
    }

    valueSpan <- diff(valueRange)
    if (!is.finite(valueSpan) | valueSpan == 0) {
      valueSpan <- 1
    }

    # Each level leaves room for the label written over the bracket below it
    levelStep <- valueSpan * (0.06 + 0.1 * labelLines)

    bracketTable$y.bracket <- valueRange[2] + valueSpan * 0.08 + levelStep * (bracketLevels - 1)
    bracketTable$y.tick <- valueSpan * 0.02
    bracketTable$y.top <- max(bracketTable$y.bracket) + levelStep

    return(bracketTable)
  } # END function




#' @title .bracketLabels
#'
#' @description Writes the labels of the brackets as markdown, either as numbers or as significance symbols, with the log2 fold change above them when it is given.
#'
#' @param pValues Numeric vector with the values to write.
#' @param prefix String naming the value, e.g. \code{"p"} or \code{"FDR"}. Default: \code{"p"}.
#' @param pLabel String, either \code{"p.value"} or \code{"stars"}. Default: \code{"p.value"}.
#' @param pDecimals Numeric value with the number of decimals. Default: \code{2}.
#' @param log2FC Numeric vector with the fold changes written on a first line, or \code{NULL}. Default: \code{NULL}.
#'
#' @return A character vector with one label per value.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.bracketLabels <-
  function(pValues,
           prefix = "p",
           pLabel = "p.value",
           pDecimals = 2,
           log2FC = NULL) {

    if (pLabel == "stars") {
      # Escaped, since an asterisk is markdown emphasis and three of them make a horizontal rule
      valueLabels <- as.character(cut(pValues,
                                      breaks = c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf),
                                      labels = c("\\*\\*\\*\\*", "\\*\\*\\*", "\\*\\*", "\\*", "ns")))
      valueLabels[is.na(valueLabels)] <- "NA"

    } else {
      formattedValues <- .formatPvalue(p = pValues, decimals = pDecimals)
      valueLabels <- paste0(prefix, ifelse(startsWith(formattedValues, "&lt;"), " ", " = "), formattedValues)
    }

    if (!is.null(log2FC)) {
      valueLabels <- paste0("log<sub>2</sub>FC = ", sprintf("%.2f", log2FC), "<br>", valueLabels)
    }

    return(valueLabels)
  } # END function




#' @title .formatPvalue
#'
#' @description Formats p-values for a markdown label, with a fixed number of decimals from 0.1 upwards and in scientific notation below, e.g. 3.20 x 10 to the -2 written with a superscript exponent. A value that underflowed to zero is written as a bound.
#'
#' @param p Numeric vector with the p-values.
#' @param decimals Numeric value with the number of decimals. Default: \code{2}.
#'
#' @return A character vector with one formatted value per p-value, \code{"NA"} for the missing ones.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.formatPvalue <-
  function(p,
           decimals = 2) {

    decimals <- max(0, as.integer(round(decimals)))

    return(vapply(p,
                  function(value) {
                    if (is.na(value)) {
                      return("NA")
                    }

                    # Zero is what the arithmetic ran out at, not what the test found
                    prefix <- ""
                    if (value <= 0) {
                      value <- .Machine$double.xmin
                      prefix <- "&lt; "
                    }

                    exponent <- floor(log10(value))
                    mantissa <- round(value / 10^exponent, decimals)

                    # Rounding can carry the mantissa to 10, which belongs to the next power
                    if (mantissa >= 10) {
                      mantissa <- mantissa / 10
                      exponent <- exponent + 1
                    }

                    if (exponent > -2) {
                      return(paste0(prefix, formatC(round(value, decimals), format = "f", digits = decimals)))
                    }

                    # The entity keeps the label in ASCII, a literal sign breaks the markdown parser under a non UTF-8 locale
                    return(paste0(prefix, formatC(mantissa, format = "f", digits = decimals), "&times;10<sup>", exponent, "</sup>"))
                  },
                  character(1),
                  USE.NAMES = FALSE))
  } # END function
