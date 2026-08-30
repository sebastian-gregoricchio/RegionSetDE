#' @title plotSetEffect
#'
#' @description Draws the effect size of every region set with its confidence interval, which is the figure the set level conclusion should rest on. The interval carries the inflation for the correlation between the regions, so a set of thirty thousand promoters does not come out looking thirty thousand times more certain than a set of two hundred enhancers.
#'
#' @param setResults \code{RegionSetDE.setResults} or \code{RegionSetDE.setResultsList} object.
#' @param value String with the quantity on the axis, either \code{"delta.log2FC"} (the difference with the background) or \code{"mean.log2FC"} (the shift away from zero). Default: \code{"delta.log2FC"}.
#' @param contrast String with the name of the contrast to draw, or its position, when \code{setResults} holds several of them. Default: \code{NULL}.
#' @param colourBy String with the variable driving the colour, either \code{"FDR"}, \code{"direction"} or \code{"none"}. Default: \code{"FDR"}.
#' @param FDR Numeric value with the adjusted p-value cut-off used by the colouring. Default: \code{NULL}, the threshold stored in the object.
#' @param showN Logical value to indicate whether the number of regions must be written next to each set. Default: \code{TRUE}.
#' @param orderBy String with the ordering of the sets, either \code{"effect"} or \code{"name"}. Default: \code{"effect"}.
#' @param colours Character vector of length two with the colours of the significant and non-significant sets. Default: \code{NULL}.
#' @param pointSize Numeric value with the size of the points. Default: \code{2.5}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}, the contrast.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}.
#' @param legendPosition String with the position of the legend. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object.
#'
#' @details \code{"mean.log2FC"} is a self-contained quantity and moves with any global shift of the mark, including one left behind by an imperfect normalisation. \code{"delta.log2FC"} is the difference between the set and what it was compared against, and a scaling error common to both cancels out of it. When the two tell different stories, the second is the one that survives a reviewer.
#'
#' Both levels of the colour scale are kept in the legend even when only one of them occurs, so that a figure in which nothing reaches the cut-off still says what the cut-off was.
#'
#' @examples
#' \dontrun{
#' plotSetEffect(setRes)
#'
#' plotSetEffect(setRes, value = "mean.log2FC", colourBy = "direction")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegionSets}}, \code{\link{plotSetDistribution}}, \code{\link{plotSetSignal}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbarh geom_vline geom_text labs theme element_blank element_rect element_line scale_colour_manual guides guide_legend
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom stats reorder
#' @importFrom methods is
#'
#' @export plotSetEffect

plotSetEffect <-
  function(setResults,
           value = "delta.log2FC",
           contrast = NULL,
           colourBy = "FDR",
           FDR = NULL,
           showN = TRUE,
           orderBy = "effect",
           colours = NULL,
           pointSize = 2.5,
           title = NULL,
           subtitle = NULL,
           legendPosition = "right",
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    setResults <- .pickResults(results = setResults, contrast = contrast)

    if (!methods::is(setResults, "RegionSetDE.setResults")) {
      stop("The 'setResults' parameter must be a RegionSetDE.setResults object.", call. = FALSE)
    }

    if (!(value %in% c("delta.log2FC", "mean.log2FC"))) {
      stop("The 'value' parameter must be either 'delta.log2FC' or 'mean.log2FC'.", call. = FALSE)
    }

    if (!(colourBy %in% c("FDR", "direction", "none"))) {
      stop("The 'colourBy' parameter must be one of 'FDR', 'direction', 'none'.", call. = FALSE)
    }

    plotTable <- setResults@results
    FDRthreshold <- if (is.null(FDR)) {setResults@thresholds$FDR} else {FDR}

    #-------------------------------#
    # Label of every row            #
    #-------------------------------#
    plotTable$set.label <- if (setResults@test == "setContrast") {
      paste(plotTable$set.1, "vs", plotTable$set.2)
    } else {
      plotTable$region.set
    }

    plotTable$n.label <- if (setResults@test == "setContrast") {
      paste0("n=", plotTable$n.regions.1, "/", plotTable$n.regions.2)
    } else {
      paste0("n=", plotTable$n.regions)
    }

    #-------------------------------#
    # The interval follows delta    #
    #-------------------------------#
    # The bounds were computed on the difference, drawing them around the self-contained mean would misplace them
    plotTable <- dplyr::mutate(plotTable,
                               effect = .data[[value]],
                               lower = if (value == "delta.log2FC") {.data$CI.lower} else {NA_real_},
                               upper = if (value == "delta.log2FC") {.data$CI.upper} else {NA_real_})

    if (orderBy == "effect") {
      plotTable$set.label <- stats::reorder(plotTable$set.label, plotTable$effect)
    } else {
      plotTable$set.label <- factor(plotTable$set.label, levels = sort(unique(as.character(plotTable$set.label)), decreasing = TRUE))
    }

    #-------------------------------#
    # Colour                        #
    #-------------------------------#
    significanceColumn <- if ("camera.FDR" %in% colnames(plotTable)) {"camera.FDR"} else {"fry.FDR"}

    if (colourBy == "FDR") {
      significanceLevels <- c(paste0("FDR < ", FDRthreshold), paste0("FDR \u2265 ", FDRthreshold))

      # Both levels stay in the scale, so a panel where nothing passes still declares the cut-off it used
      plotTable$colour.group <- factor(ifelse(plotTable[[significanceColumn]] < FDRthreshold,
                                              significanceLevels[1], significanceLevels[2]),
                                       levels = significanceLevels)

      colourValues <- if (is.null(colours)) {c("#B2182B", "grey60")} else {colours}
      names(colourValues) <- significanceLevels

    } else if (colourBy == "direction") {
      plotTable$colour.group <- factor(ifelse(plotTable$effect > 0, "up", "down"), levels = c("down", "up"))
      colourValues <- if (is.null(colours)) {c("down" = "#2166AC", "up" = "#B2182B")} else {colours}

    } else {
      plotTable$colour.group <- factor("all")
      colourValues <- c("all" = "black")
    }

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    effectPlot <-
      ggplot2::ggplot(data = plotTable,
                      mapping = ggplot2::aes(x = .data$effect, y = .data$set.label, colour = .data$colour.group)) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "black") +
      ggplot2::geom_point(size = pointSize, stroke = NA) +
      ggplot2::scale_colour_manual(values = colourValues, drop = FALSE) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = max(c(pointSize, 3))))) +
      ggplot2::labs(x = paste0(sub("log2FC", "log<sub>2</sub>FC", value), " (", setResults@contrast, ")"),
                    y = "",
                    colour = "",
                    title = if (is.null(title)) {setResults@contrast} else {title},
                    subtitle = subtitle,
                    caption = paste0(paste(setResults@methods, collapse = " and "),
                                     ", intervals inflated for the correlation between regions")) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize) +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.5, colour = "black"),
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.2, colour = "gray"))

    if (value == "delta.log2FC") {
      effectPlot <- effectPlot +
        ggplot2::geom_errorbarh(mapping = ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                                height = 0.15, linewidth = 0.5)
    }

    if (isTRUE(showN)) {
      effectPlot <- effectPlot +
        ggplot2::geom_text(mapping = ggplot2::aes(label = .data$n.label),
                           hjust = -0.2, vjust = -1, size = baseSize / 4.5, colour = "grey30", show.legend = FALSE)
    }

    return(effectPlot)
  } # END function




#' @title plotSetDistribution
#'
#' @description Draws the whole distribution of the per-region log2 fold changes of every set, rather than its summary. A set whose mean has moved because a small group of regions collapsed looks nothing like one whose regions all shifted a little, and only the distribution tells the two apart.
#'
#' @param setResults \code{RegionSetDE.setResults} or \code{RegionSetDE.setResultsList} object.
#' @param style String with the kind of plot, one of \code{"violin"}, \code{"boxplot"} and \code{"ecdf"}. Default: \code{"violin"}.
#' @param set Character vector with the names of the region sets to draw. Default: \code{NULL}, all of them.
#' @param contrast String with the name of the contrast to draw, or its position, when \code{setResults} holds several of them. Default: \code{NULL}.
#' @param annotate Logical value to indicate whether the effect size and the adjusted p-value must be written above each set. Ignored for \code{style = "ecdf"}. Default: \code{TRUE}.
#' @param colours Named character vector with one colour per region set. Default: \code{NULL}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}, the contrast.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}.
#' @param legendPosition String with the position of the legend. Default: \code{NULL}, no legend for the violins and the boxplots, where the sets are already named on the axis, and a legend for the cumulative curves, where they are not.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' plotSetDistribution(setRes)
#'
#' plotSetDistribution(setRes, style = "ecdf")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegionSets}}, \code{\link{plotSetEffect}}, \code{\link{plotSetSignal}}
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot stat_ecdf geom_hline geom_vline geom_text labs theme element_blank scale_colour_manual scale_fill_manual
#' @importFrom dplyr filter mutate
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export plotSetDistribution

plotSetDistribution <-
  function(setResults,
           style = "violin",
           set = NULL,
           contrast = NULL,
           annotate = TRUE,
           colours = NULL,
           title = NULL,
           subtitle = NULL,
           legendPosition = NULL,
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    setResults <- .pickResults(results = setResults, contrast = contrast)

    if (!methods::is(setResults, "RegionSetDE.setResults")) {
      stop("The 'setResults' parameter must be a RegionSetDE.setResults object.", call. = FALSE)
    }

    if (!(style %in% c("violin", "boxplot", "ecdf"))) {
      stop("The 'style' parameter must be one of 'violin', 'boxplot', 'ecdf'.", call. = FALSE)
    }

    # The sets are written on the x axis of a violin, a legend would repeat them
    if (is.null(legendPosition)) {
      legendPosition <- if (style == "ecdf") {"right"} else {"none"}
    }

    plotTable <- setResults@regionStats

    if (!is.null(set)) {
      absentSets <- setdiff(set, unique(plotTable$region.set))
      if (length(absentSets) > 0) {
        stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }
      plotTable <- dplyr::filter(plotTable, .data$region.set %in% set)
    }

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    if (style == "ecdf") {
      distributionPlot <-
        ggplot2::ggplot(data = plotTable,
                        mapping = ggplot2::aes(x = .data$log2FC, colour = .data$region.set)) +
        ggplot2::stat_ecdf(geom = "step", linewidth = 0.7) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "black") +
        ggplot2::labs(x = paste0("log<sub>2</sub> fold change (", setResults@contrast, ")"),
                      y = "Cumulative fraction of regions",
                      colour = "") +
        .resultsTheme(legendPosition = legendPosition, baseSize = baseSize)

    } else {
      distributionPlot <-
        ggplot2::ggplot(data = plotTable,
                        mapping = ggplot2::aes(x = .data$region.set, y = .data$log2FC,
                                               colour = .data$region.set, fill = .data$region.set)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "black") +
        ggplot2::labs(x = "",
                      y = paste0("log<sub>2</sub> fold change (", setResults@contrast, ")"),
                      colour = "", fill = "") +
        .resultsTheme(legendPosition = legendPosition, baseSize = baseSize, rotateX = TRUE) +
        ggplot2::theme(axis.ticks.x = ggplot2::element_blank())

      if (style == "violin") {
        distributionPlot <- distributionPlot +
          ggplot2::geom_violin(alpha = 0.3, linewidth = 0.4, scale = "width") +
          ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.4)
      } else {
        distributionPlot <- distributionPlot +
          ggplot2::geom_boxplot(alpha = 0.3, outlier.size = 0.3, outlier.stroke = NA, linewidth = 0.4)
      }

      #-------------------------------#
      # Annotation over each set      #
      #-------------------------------#
      if (isTRUE(annotate) & setResults@test == "set") {
        distributionPlot <- distributionPlot +
          ggplot2::geom_text(data = .setAnnotation(setResults = setResults,
                                                   setNames = unique(plotTable$region.set),
                                                   yPosition = max(plotTable$log2FC, na.rm = TRUE) * 1.05),
                             mapping = ggplot2::aes(x = .data$region.set, y = .data$y.position, label = .data$label),
                             inherit.aes = FALSE, size = baseSize / 4.5, colour = "grey20", vjust = 0)
      }
    }

    distributionPlot <- distributionPlot +
      ggplot2::labs(title = if (is.null(title)) {setResults@contrast} else {title},
                    subtitle = subtitle)

    if (!is.null(colours)) {
      distributionPlot <- distributionPlot +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours)
    }

    return(distributionPlot)
  } # END function




#' @title plotSetSignal
#'
#' @description Draws the signal of every sample over the regions of each set, one violin per sample, with a bracket joining the groups being compared and the set level fold change and p-value written on it. Where \code{\link{plotSetDistribution}} shows the fold changes the model estimated, this one shows the values those estimates came from, so a set that moved because one replicate is out of line is visible rather than hidden inside a mean.
#'
#' @param setResults \code{RegionSetDE.setResults} or \code{RegionSetDE.setResultsList} object. A \code{RegionSetDE.fit} or a \code{RegionSetDE.counts} object is accepted as well, in which case no annotation is written.
#' @param counts \code{RegionSetDE.counts} object holding the values, when \code{setResults} carries none. Default: \code{NULL}.
#' @param contrast String with the name of the contrast to draw, or its position, when \code{setResults} holds several of them. Default: \code{NULL}.
#' @param set Character vector with the names of the region sets to draw. Default: \code{NULL}, all of them.
#' @param groupBy String with the name of a \code{colData} column driving the colour and the ordering of the samples. Default: \code{NULL}, the column the contrast separates.
#' @param groupOrder Character vector with the levels of \code{groupBy} in the order they must appear, which is also the order the colour families are handed out in. Default: \code{NULL}, the reference level of the contrast first, alphabetical when the contrast does not name two levels.
#' @param comparisons List of character vectors of length two, naming the groups joined by a bracket. Only the pair the contrast actually compares is labelled. Default: \code{NULL}, that pair.
#' @param valueColumn String with the fold change written on the bracket, either \code{"mean.log2FC"} (the shift of the set away from zero, which is what the two sides of the bracket differ by) or \code{"delta.log2FC"} (the difference with the background). Default: \code{"mean.log2FC"}.
#' @param style String with the kind of plot, either \code{"violin"} or \code{"boxplot"}. Default: \code{"violin"}.
#' @param assay String with the name of the assay to draw. Default: \code{NULL}, the normalised assay when present, the raw counts otherwise.
#' @param log2Scale Logical value to indicate whether the values must be drawn on a log2 scale. Default: \code{TRUE}.
#' @param annotate Logical value to indicate whether the bracket and its label must be drawn. Default: \code{TRUE}.
#' @param facetScales String with the scales of the panels, one among \code{"fixed"}, \code{"free"}, \code{"free_x"} and \code{"free_y"}. Default: \code{"fixed"}.
#' @param maxRegions Numeric value with the number of regions drawn per set. Default: \code{20000}.
#' @param colours Named character vector with one colour per sample, overriding the shades built from the groups. Default: \code{NULL}.
#' @param baseColours Character vector with one base colour per group, from which the shades of the replicates are built. Default: \code{NULL}, blue then red, matching the down and up colours of \code{\link{plotVolcano}}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}, the contrast.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}.
#' @param legendPosition String with the position of the legend. Default: \code{"none"}, the samples are already named on the axis.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object.
#'
#' @details The violins are drawn on the values as they are, without centring each region on its own mean, so what separates two groups here is the difference in absolute signal rather than the per-region contrast the model tested. Those two disagree whenever a set is heterogeneous: a handful of very strong regions dominate the shape of a violin while contributing one row each to the fold change. The number on the bracket comes from the set level test, not from the violins under it, which is why a small visible shift can carry a decisive p-value and a large one need not.
#'
#' \code{valueColumn} decides which fold change that number is. \code{"mean.log2FC"} is the shift of the set away from zero under the contrast, which is the quantity the two sides of the bracket differ by and the reason it is the default. \code{"delta.log2FC"} is the difference between the set and the background it was compared against, a different comparison that the bracket does not draw.
#'
#' The fold change and the p-value belong to one comparison, the one the contrast declared, so they are written only when the axis is grouped by the variable that contrast separates. Grouping by anything else, a replicate or a batch, still colours and orders the violins but drops the bracket: the numbers would describe a comparison the figure is no longer showing. The default for \code{groupBy} is therefore the column the contrast came from, read from the object rather than guessed.
#'
#' Each group is given a colour family and the replicates inside it a shade, running from light to the base colour in the order the samples appear. The first group gets the blue of the down regions of \code{\link{plotVolcano}} and the second the red of the up ones, so a figure made of both reads consistently; the reference level of the contrast is the one that comes out blue.
#'
#' Sets larger than \code{maxRegions} are thinned, deterministically, before the violin is computed. The annotation is unaffected, since it is read from the test rather than recomputed here.
#'
#' @examples
#' \dontrun{
#' # groupBy defaults to the column the contrast separates
#' plotSetSignal(setRes)
#'
#' plotSetSignal(setRes, style = "boxplot", facetScales = "free_y")
#'
#' # Ordered by replicate: the violins are recoloured, the bracket is dropped
#' plotSetSignal(setRes, groupBy = "replicate")
#'
#' # Three groups, two brackets
#' plotSetSignal(setRes, groupBy = "treatment",
#'               comparisons = list(c("DMSO", "EPZ"), c("DMSO", "COMBO")))
#'
#' # Before any test has been run
#' plotSetSignal(counts, groupBy = "treatment", annotate = FALSE)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotSetDistribution}}, \code{\link{plotSetEffect}}, \code{\link{plotRegion}}
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_segment geom_text facet_wrap labs theme element_blank scale_colour_manual scale_fill_manual scale_y_continuous expansion
#' @importFrom ggtext geom_richtext
#' @importFrom dplyr filter mutate
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export plotSetSignal

plotSetSignal <-
  function(setResults,
           counts = NULL,
           contrast = NULL,
           set = NULL,
           groupBy = NULL,
           groupOrder = NULL,
           comparisons = NULL,
           valueColumn = "mean.log2FC",
           style = "violin",
           assay = NULL,
           log2Scale = TRUE,
           annotate = TRUE,
           facetScales = "fixed",
           maxRegions = 20000,
           colours = NULL,
           baseColours = NULL,
           title = NULL,
           subtitle = NULL,
           legendPosition = "none",
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!(style %in% c("violin", "boxplot"))) {
      stop("The 'style' parameter must be either 'violin' or 'boxplot'.", call. = FALSE)
    }

    if (!(facetScales %in% c("fixed", "free", "free_x", "free_y"))) {
      stop("The 'facetScales' parameter must be one among 'fixed', 'free', 'free_x' or 'free_y'.", call. = FALSE)
    }

    if (!(valueColumn %in% c("mean.log2FC", "delta.log2FC"))) {
      stop("The 'valueColumn' parameter must be either 'mean.log2FC' or 'delta.log2FC'.", call. = FALSE)
    }

    resolvedObject <- .resolveCounts(object = setResults, counts = counts, contrast = contrast)
    counts <- resolvedObject$counts
    setResults <- resolvedObject$results

    #-------------------------------#
    # Long table of the values      #
    #-------------------------------#
    assay <- .defaultAssay(counts = counts, assay = assay)

    rowTable <- as.data.frame(SummarizedExperiment::rowData(counts))
    rowTable$row.index <- seq_len(nrow(counts))

    if (!is.null(set)) {
      absentSets <- setdiff(set, unique(as.character(rowTable$region.set)))
      if (length(absentSets) > 0) {
        stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }
      rowTable <- dplyr::filter(rowTable, .data$region.set %in% set)
    }

    # A violin over a hundred thousand regions draws the same shape as one over twenty thousand, at a fraction of the cost
    rowTable <- .thinBySet(regionTable = rowTable, maxPoints = maxRegions, bySet = TRUE)

    valueMatrix <- as.matrix(SummarizedExperiment::assay(counts, assay))[rowTable$row.index, , drop = FALSE]

    plotTable <- data.frame(sample = rep(colnames(counts), each = nrow(valueMatrix)),
                            region.set = rep(as.character(rowTable$region.set), times = ncol(valueMatrix)),
                            value = as.numeric(valueMatrix),
                            stringsAsFactors = FALSE)

    if (isTRUE(log2Scale)) {
      plotTable <- dplyr::mutate(plotTable, value = log2(.data$value + 1))
    }

    #-------------------------------#
    # Order the samples by group    #
    #-------------------------------#
    colTable <- as.data.frame(SummarizedExperiment::colData(counts))
    colTable$sample <- colnames(counts)

    # The contrast knows which variable it separates, and that is the only grouping its numbers describe
    contrastInfo <- if (is.null(setResults)) {list()} else {setResults@contrast.groups}
    contrastColumn <- contrastInfo$column
    contrastGroups <- contrastInfo$groups

    if (is.null(groupBy)) {
      groupBy <- contrastColumn
    }

    if (!is.null(groupBy)) {
      if (!(groupBy %in% colnames(colTable))) {
        stop(paste0("The column '", groupBy, "' is absent from the colData."), call. = FALSE)
      }
      colTable$group <- as.character(colTable[[groupBy]])
    } else {
      colTable$group <- colTable$sample
    }

    groupsMatchContrast <- !is.null(contrastColumn) && identical(groupBy, contrastColumn)

    groupOrder <- .resolveGroupOrder(groupOrder = groupOrder,
                                     groupLevels = unique(colTable$group),
                                     contrastGroups = if (isTRUE(groupsMatchContrast)) {contrastGroups} else {NULL})

    # Ordering the axis by group puts the replicates of a condition next to each other
    colTable$group <- factor(colTable$group, levels = groupOrder)
    colTable <- colTable[order(colTable$group, colTable$sample), , drop = FALSE]

    plotTable$group <- as.character(colTable$group)[match(plotTable$sample, colTable$sample)]
    plotTable$sample <- factor(plotTable$sample, levels = colTable$sample)

    sampleColours <- if (is.null(colours)) {
      .groupShades(colTable = colTable, baseColours = baseColours)
    } else {
      colours
    }

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    contrastLabel <- if (is.null(setResults)) {NULL} else {setResults@contrast}

    signalPlot <-
      ggplot2::ggplot(data = plotTable,
                      mapping = ggplot2::aes(x = .data$sample, y = .data$value,
                                             colour = .data$sample, fill = .data$sample)) +
      ggplot2::facet_wrap(facets = ~ region.set, scales = facetScales) +
      ggplot2::scale_colour_manual(values = sampleColours) +
      ggplot2::scale_fill_manual(values = sampleColours) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
      ggplot2::labs(x = "",
                    y = paste0(if (isTRUE(log2Scale)) {"log<sub>2</sub> "} else {""}, assay),
                    colour = "", fill = "",
                    title = if (is.null(title)) {contrastLabel} else {title},
                    subtitle = subtitle) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize, rotateX = TRUE) +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank())

    if (style == "violin") {
      signalPlot <- signalPlot +
        ggplot2::geom_violin(alpha = 0.35, linewidth = 0.4, scale = "width") +
        ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.4)
    } else {
      signalPlot <- signalPlot +
        ggplot2::geom_boxplot(alpha = 0.35, outlier.size = 0.3, outlier.stroke = NA, linewidth = 0.4)
    }

    #-------------------------------#
    # Bracket over the groups       #
    #-------------------------------#
    if (isTRUE(annotate) & !is.null(setResults)) {
      if (setResults@test != "set") {
        warning("The annotation is written from a set level test, and this object holds a contrast between sets.", call. = FALSE)

      } else if (isFALSE(groupsMatchContrast)) {
        # Labelling a bracket over the replicates with the treatment effect would attach the number to the wrong comparison
        message(paste0("The samples are grouped by '", if (is.null(groupBy)) {"sample"} else {groupBy},
                       "' while the contrast compares ",
                       if (is.null(contrastColumn)) {"coefficients of the design"} else {paste0("'", contrastColumn, "'")},
                       ", so no fold change is written."))

      } else {
        bracketTable <- .signalBrackets(plotTable = plotTable,
                                        colTable = colTable,
                                        setResults = setResults,
                                        comparisons = comparisons,
                                        contrastGroups = contrastGroups,
                                        groupOrder = groupOrder,
                                        valueColumn = valueColumn,
                                        perFacet = grepl("free", facetScales) & facetScales != "free_x")

        if (!is.null(bracketTable)) {
          signalPlot <- signalPlot +
            ggplot2::geom_segment(data = bracketTable,
                                  mapping = ggplot2::aes(x = .data$x.start, xend = .data$x.end,
                                                         y = .data$y.bracket, yend = .data$y.bracket),
                                  inherit.aes = FALSE, linewidth = 0.4, colour = "grey20") +
            ggplot2::geom_segment(data = bracketTable,
                                  mapping = ggplot2::aes(x = .data$x.start, xend = .data$x.start,
                                                         y = .data$y.bracket, yend = .data$y.bracket - .data$y.tick),
                                  inherit.aes = FALSE, linewidth = 0.4, colour = "grey20") +
            ggplot2::geom_segment(data = bracketTable,
                                  mapping = ggplot2::aes(x = .data$x.end, xend = .data$x.end,
                                                         y = .data$y.bracket, yend = .data$y.bracket - .data$y.tick),
                                  inherit.aes = FALSE, linewidth = 0.4, colour = "grey20") +
            ggtext::geom_richtext(data = bracketTable,
                                  mapping = ggplot2::aes(x = (.data$x.start + .data$x.end) / 2,
                                                         y = .data$y.bracket, label = .data$label),
                                  inherit.aes = FALSE, vjust = 0, size = baseSize / 4.5,
                                  colour = "grey20", fill = NA, label.color = NA,
                                  label.padding = grid::unit(rep(1, 4), "pt"))
        }
      }
    }

    return(signalPlot)
  } # END function




#' @title .resolveGroupOrder
#'
#' @description Decides the order the groups appear in, putting the reference level of the contrast on the left so that it gets the blue family.
#'
#' @param groupOrder Character vector given by the user, or \code{NULL}.
#' @param groupLevels Character vector with the groups present in the object.
#' @param contrastGroups Character vector of length two with the levels the contrast compares, the first being the one the fold change is positive for. Default: \code{NULL}.
#'
#' @return A character vector with every group, ordered.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resolveGroupOrder <-
  function(groupOrder = NULL,
           groupLevels,
           contrastGroups = NULL) {

    if (!is.null(groupOrder)) {
      absentGroups <- setdiff(groupOrder, groupLevels)
      if (length(absentGroups) > 0) {
        stop(paste0("The following groups are absent from the object: ", paste(absentGroups, collapse = ", "), "."), call. = FALSE)
      }
      return(c(groupOrder, setdiff(groupLevels, groupOrder)))
    }

    #-------------------------------#
    # Read it from the contrast     #
    #-------------------------------#
    # The fold change is positive for the first level, so the second one is the reference and belongs on the left
    if (!is.null(contrastGroups) && length(contrastGroups) == 2 && all(contrastGroups %in% groupLevels)) {
      return(c(contrastGroups[2], contrastGroups[1], setdiff(groupLevels, contrastGroups)))
    }

    return(sort(groupLevels))
  } # END function




#' @title .groupShades
#'
#' @description Builds one colour per sample, a family per group and a shade per replicate inside it.
#'
#' @param colTable Data.frame with the \code{sample} and \code{group} columns, already ordered.
#' @param baseColours Character vector with one base colour per group, or \code{NULL}.
#'
#' @return A named character vector with one colour per sample.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom grDevices colorRampPalette col2rgb rgb
#'
#' @keywords internal

.groupShades <-
  function(colTable,
           baseColours = NULL) {

    # Blue then red, the same pair the volcano and the MA plot use for the down and up regions
    if (is.null(baseColours)) {
      baseColours <- c("#2166AC", "#B2182B", "#1B7837", "#E08214", "#762A83", "#01665E")
    }

    groupLevels <- levels(colTable$group)
    if (is.null(groupLevels)) {
      groupLevels <- unique(as.character(colTable$group))
    }

    shadeList <-
      lapply(seq_along(groupLevels),
             function(i) {
               groupSamples <- colTable$sample[as.character(colTable$group) == groupLevels[i]]
               baseColour <- baseColours[((i - 1) %% length(baseColours)) + 1]

               # A single replicate has no shading to do, and colorRampPalette of length one returns the light end
               groupShades <- if (length(groupSamples) < 2) {
                 baseColour
               } else {
                 grDevices::colorRampPalette(c(.lightenColour(colour = baseColour, amount = 0.55), baseColour))(length(groupSamples))
               }

               names(groupShades) <- groupSamples
               return(groupShades)
             })

    return(unlist(shadeList))
  } # END function




#' @title .lightenColour
#'
#' @description Mixes a colour with white, used to build the light end of the shades of a group.
#'
#' @param colour String with the colour.
#' @param amount Numeric value between 0 and 1 with the fraction of white.
#'
#' @return A string with the mixed colour.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom grDevices col2rgb rgb
#'
#' @keywords internal

.lightenColour <-
  function(colour,
           amount = 0.5) {

    colourChannels <- as.numeric(grDevices::col2rgb(colour)) / 255
    mixedChannels <- colourChannels + (1 - colourChannels) * amount

    return(grDevices::rgb(red = mixedChannels[1], green = mixedChannels[2], blue = mixedChannels[3]))
  } # END function




#' @title .signalBrackets
#'
#' @description Places the brackets joining the groups compared in a signal plot, and writes the fold change and the adjusted p-value of the set on each of them.
#'
#' @param plotTable Data.frame with the values being drawn.
#' @param colTable Data.frame with the samples and their groups, ordered.
#' @param setResults \code{RegionSetDE.setResults} object.
#' @param comparisons List of character vectors of length two, or \code{NULL}.
#' @param contrastGroups Character vector of length two with the levels the contrast compares, or \code{NULL}.
#' @param groupOrder Character vector with the groups, ordered.
#' @param valueColumn String with the fold change written on the bracket.
#' @param perFacet Logical value indicating whether the height must be computed inside each panel.
#'
#' @return A data.frame with one row per bracket, or \code{NULL} when no comparison can be drawn.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom stats setNames
#'
#' @keywords internal

.signalBrackets <-
  function(plotTable,
           colTable,
           setResults,
           comparisons = NULL,
           contrastGroups = NULL,
           groupOrder,
           valueColumn = "mean.log2FC",
           perFacet = FALSE) {

    #-------------------------------#
    # Pairs to join                 #
    #-------------------------------#
    if (is.null(contrastGroups) | length(contrastGroups) != 2) {
      return(NULL)
    }

    if (is.null(comparisons)) {
      comparisons <- list(contrastGroups)
    }

    if (!is.list(comparisons) | any(vapply(comparisons, length, integer(1)) != 2)) {
      stop("The 'comparisons' parameter must be a list of character vectors of length two.", call. = FALSE)
    }

    absentGroups <- setdiff(unlist(comparisons), groupOrder)
    if (length(absentGroups) > 0) {
      stop(paste0("The following groups are absent from the object: ", paste(absentGroups, collapse = ", "), "."), call. = FALSE)
    }

    # A pair the contrast did not compare has no fold change of its own, so it gets no bracket rather than a borrowed one
    comparisons <- comparisons[vapply(comparisons, function(x) {setequal(x, contrastGroups)}, logical(1))]

    if (length(comparisons) == 0) {
      message("None of the requested comparisons matches the pair the contrast tested, no bracket has been drawn.")
      return(NULL)
    }

    significanceColumn <- if ("camera.FDR" %in% colnames(setResults@results)) {"camera.FDR"} else {"fry.FDR"}
    testName <- if (significanceColumn == "camera.FDR") {"camera"} else {"fry"}

    samplePositions <- stats::setNames(seq_along(levels(plotTable$sample)), levels(plotTable$sample))
    setNames <- unique(plotTable$region.set)

    #-------------------------------#
    # One bracket per set and pair  #
    #-------------------------------#
    bracketList <-
      lapply(setNames,
             function(setName) {
               setRow <- dplyr::filter(setResults@results, .data$region.set == setName)
               if (nrow(setRow) != 1) {
                 return(NULL)
               }

               # With free scales every panel has its own ceiling, otherwise they share the one of the whole plot
               panelValues <- if (isTRUE(perFacet)) {
                 plotTable$value[plotTable$region.set == setName]
               } else {
                 plotTable$value
               }

               valueRange <- range(panelValues, na.rm = TRUE)
               valueSpan <- diff(valueRange)
               if (!is.finite(valueSpan) | valueSpan == 0) {
                 valueSpan <- 1
               }

               pairList <-
                 lapply(seq_along(comparisons),
                        function(i) {
                          firstSamples <- colTable$sample[as.character(colTable$group) == comparisons[[i]][1]]
                          secondSamples <- colTable$sample[as.character(colTable$group) == comparisons[[i]][2]]

                          return(data.frame(region.set = setName,
                                            x.start = mean(samplePositions[firstSamples]),
                                            x.end = mean(samplePositions[secondSamples]),
                                            y.bracket = valueRange[2] + valueSpan * (0.04 + 0.12 * (i - 1)),
                                            y.tick = valueSpan * 0.02,
                                            label = sprintf("log<sub>2</sub>FC = %.2f<br>%s FDR = %.1e",
                                                            setRow[[valueColumn]][1], testName, setRow[[significanceColumn]][1]),
                                            stringsAsFactors = FALSE))
                        })

               return(do.call(what = rbind, args = pairList))
             })

    bracketTable <- do.call(what = rbind, args = bracketList)

    if (is.null(bracketTable)) {
      return(NULL)
    }

    if (nrow(bracketTable) == 0) {
      return(NULL)
    }

    return(bracketTable)
  } # END function


#' @title plotUniverseMatching
#'
#' @description Compares a region set with the rows it is compared against, on width and on abundance. The two distributions should sit on top of each other; where they do not, the competitive test is partly reading the difference between the intervals rather than the difference between the biologies.
#'
#' @param object \code{RegionSetDE.fit} or \code{RegionSetDE.setResults} object, both of which carry the universe they used.
#' @param universe \code{RegionSetDE.universe} object. Default: \code{NULL}, the one stored in \code{object}.
#' @param contrast String with the name of the contrast, or its position, when \code{object} holds several of them. Default: \code{NULL}.
#' @param set Character vector with the names of the region sets to draw. Default: \code{NULL}, all of them.
#' @param covariate String with the covariate on the axis, either \code{"width"} or \code{"abundance"}. Default: \code{"abundance"}.
#' @param colours Character vector of length two with the colours of the set and of the rows it is compared against. Default: \code{NULL}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}.
#' @param legendPosition String with the position of the legend. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object.
#'
#' @details One panel per region set, holding two density curves: the regions of the set and the rows it is compared against. Read it as a check on the matching rather than as a result. Curves lying on top of each other mean the two groups are comparable on that covariate, so a difference found by \code{\link{testRegionSets}} cannot be attributed to it. Curves offset from each other mean the matching did not find enough eligible rows in some strata, which happens when a set occupies a corner of the width or abundance range that nothing else reaches, and the numbers behind it are in the \code{diagnostics} slot of the universe.
#'
#' @examples
#' \dontrun{
#' # Straight off the fit, before any test has been run
#' plotUniverseMatching(fit, covariate = "abundance")
#'
#' setRes <- testRegionSets(fit, contrast = "conditionCOMBO")
#' plotUniverseMatching(setRes, covariate = "width")
#'
#' # On a universe built by hand
#' plotUniverseMatching(fit, universe = makeSetUniverse(fit, match = "width"))
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{makeSetUniverse}}, \code{\link{testRegionSets}}
#'
#' @importFrom SummarizedExperiment assay colData rowData rowRanges
#' @importFrom BiocGenerics width
#' @importFrom ggplot2 ggplot aes geom_density facet_wrap labs scale_colour_manual scale_fill_manual scale_x_log10
#' @importFrom dplyr filter bind_rows mutate
#' @importFrom rlang .data
#' @importFrom stats median
#' @importFrom methods is
#'
#' @export plotUniverseMatching

plotUniverseMatching <-
  function(object,
           universe = NULL,
           contrast = NULL,
           set = NULL,
           covariate = "abundance",
           colours = NULL,
           title = NULL,
           legendPosition = "right",
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    resolvedObject <- .resolveCounts(object = object, counts = NULL, contrast = contrast)
    counts <- resolvedObject$counts

    # The fit and the set level result both carry the universe they used
    if (is.null(universe)) {
      universe <- if (methods::is(object, "RegionSetDE.fit")) {
        object@universe
      } else if (methods::is(resolvedObject$results, "RegionSetDE.setResults")) {
        resolvedObject$results@universe
      } else {
        NULL
      }
    }

    if (!methods::is(universe, "RegionSetDE.universe")) {
      stop("No universe is available, pass one built with makeSetUniverse through 'universe'.", call. = FALSE)
    }

    if (length(universe@index) == 0) {
      stop("The universe is empty, which happens when the object holds a single region set.", call. = FALSE)
    }

    if (!(covariate %in% c("width", "abundance"))) {
      stop("The 'covariate' parameter must be either 'width' or 'abundance'.", call. = FALSE)
    }

    #-------------------------------#
    # Row level covariates          #
    #-------------------------------#
    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, 1))
    rowWidths <- BiocGenerics::width(SummarizedExperiment::rowRanges(counts))

    librarySizes <- SummarizedExperiment::colData(counts)$library.size
    if (is.null(librarySizes) | any(is.na(librarySizes))) {
      librarySizes <- colSums(countMatrix)
    }

    rowTable <- data.frame(region.set = as.character(SummarizedExperiment::rowData(counts)$region.set),
                           row.index = seq_len(nrow(counts)),
                           width = rowWidths,
                           abundance = .regionAbundance(countMatrix = countMatrix,
                                                        librarySizes = librarySizes,
                                                        rowWidths = rowWidths,
                                                        referenceWidth = stats::median(rowWidths)),
                           stringsAsFactors = FALSE)

    setNames <- if (is.null(set)) {names(universe@index)} else {set}
    absentSets <- setdiff(setNames, names(universe@index))
    if (length(absentSets) > 0) {
      stop(paste0("The universe holds no entry for: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
    }

    #-------------------------------#
    # Set against its background    #
    #-------------------------------#
    plotTable <-
      dplyr::bind_rows(lapply(setNames,
                              function(setName) {
                                setIndex <- rowTable$row.index[rowTable$region.set == setName]

                                # The universe holds the set as well, the curve to compare against is what is left
                                comparisonIndex <- setdiff(universe@index[[setName]], setIndex)

                                setRows <- dplyr::mutate(rowTable[setIndex, , drop = FALSE],
                                                         panel = setName, origin = "set")
                                comparisonRows <- dplyr::mutate(rowTable[comparisonIndex, , drop = FALSE],
                                                                panel = setName, origin = "compared against")
                                return(rbind(setRows, comparisonRows))
                              }))

    plotTable$origin <- factor(plotTable$origin, levels = c("set", "compared against"))

    colourValues <- if (is.null(colours)) {c("set" = "#B2182B", "compared against" = "grey55")} else {colours}

    matchingPlot <-
      ggplot2::ggplot(data = plotTable,
                      mapping = ggplot2::aes(x = .data[[covariate]], colour = .data$origin, fill = .data$origin)) +
      ggplot2::geom_density(alpha = 0.25, linewidth = 0.5) +
      ggplot2::facet_wrap(facets = ~ panel, scales = "free_y") +
      ggplot2::scale_colour_manual(values = colourValues, drop = FALSE) +
      ggplot2::scale_fill_manual(values = colourValues, drop = FALSE) +
      ggplot2::labs(x = if (covariate == "width") {"Region width (bp)"} else {"Width-adjusted abundance (log<sub>2</sub> CPM)"},
                    y = "Density", colour = "", fill = "",
                    title = if (is.null(title)) {paste0("Universe matching on ", covariate)} else {title}) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize)

    if (covariate == "width") {
      matchingPlot <- matchingPlot + ggplot2::scale_x_log10()
    }

    return(matchingPlot)
  } # END function




#' @title .setAnnotation
#'
#' @description Builds the label written above each region set, carrying the effect size and the adjusted p-value of the set level test.
#'
#' @param setResults \code{RegionSetDE.setResults} object.
#' @param setNames Character vector with the sets present in the plot.
#' @param yPosition Numeric value with the height the labels are drawn at.
#'
#' @return A data.frame with the \code{region.set}, \code{label} and \code{y.position} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter mutate
#' @importFrom rlang .data
#'
#' @keywords internal

.setAnnotation <-
  function(setResults,
           setNames,
           yPosition) {

    significanceColumn <- if ("camera.FDR" %in% colnames(setResults@results)) {"camera.FDR"} else {"fry.FDR"}
    testName <- if (significanceColumn == "camera.FDR") {"camera"} else {"fry"}

    # The label names the test, since a self-contained and a competitive FDR are not the same claim
    annotationTable <- dplyr::mutate(setResults@results,
                                     label = sprintf("delta %.2f\n%s FDR %.1e",
                                                     .data$delta.log2FC, testName, .data[[significanceColumn]]),
                                     y.position = yPosition)

    return(as.data.frame(dplyr::filter(annotationTable, .data$region.set %in% setNames)))
  } # END function
