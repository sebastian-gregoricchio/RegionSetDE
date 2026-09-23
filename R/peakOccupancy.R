# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title peakOccupancyTable
#'
#' @description Crosses the consensus a region came from with what happened to it in the test: how many of the regions called in both groups changed, how many of the ones called in a single group did, and in which direction. It works on the results of regions built by \code{\link{loadConsensusPeaks}}, which carry the occupancy of every group.
#'
#' @param results \code{RegionSetDE.results} object, or a \code{RegionSetDE.resultsList} together with \code{contrast}.
#' @param contrast String with the name of a contrast, or its position, when \code{results} holds several of them. Default: \code{NULL}.
#' @param groups Character vector with the consensus groups the regions are classified by. Default: \code{NULL}, the two groups the contrast compares when they are among them, every group otherwise.
#' @param by String indicating what the regions are grouped by, either \code{"group"}, the combination of consensus groups covering the region, or \code{"samples"}, the number of samples that carried a peak on it. Default: \code{"group"}.
#' @param set Character vector with the names of the region sets kept. Default: \code{NULL}, all of them.
#' @param FDR Numeric value with the adjusted p-value below which a region counts as changed. Default: \code{NULL}, the threshold the test was run with.
#' @param log2FC Numeric value with the log2 fold change a region must reach. Default: \code{NULL}, the threshold the test was run with.
#'
#' @return A data.frame with one row per occupancy class: the class itself, the number of regions in it, how many of them came out \code{down}, \code{null} and \code{up}, and the percentage that changed in either direction.
#'
#' @details The table says where the differences sit, which is the question a peak-based experiment usually starts from: whether the condition rearranged the shared peaks or whether it mostly gained and lost whole sites. A region called in one group only and coming out unchanged is worth as much attention as the reverse, since it means the peak caller drew a line the counts do not support.
#'
#' What the table is not is a test of occupancy. Whether a peak is called in a sample depends on the depth of that library and on the threshold of the caller as much as on the chromatin, so the number of group-specific regions is not an effect size and comparing those numbers between groups is not evidence of anything. The evidence is in the columns beside them, which come from the counts.
#'
#' @examples
#' if (requireNamespace("consensusRegions", quietly = TRUE)) {
#'   # AR binding without ligand and after 24 hours of R1881
#'   sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#'   sampleSheet <- dplyr::filter(sampleSheet, condition %in% c("DMSO", "R1881_24h"))
#'
#'   consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)
#'
#'   counts <- countReads(consensus, sampleSheet = sampleSheet, verbose = FALSE)
#'   counts <- countBackground(counts, binSize = 10000, verbose = FALSE)
#'   counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#'   fit <- fitRegions(counts, design = ~ condition, verbose = FALSE)
#'   results <- testRegions(fit, contrast = c("condition", "R1881_24h", "DMSO"), verbose = FALSE)
#'
#'   peakOccupancyTable(results)
#'
#'   # By number of samples calling a peak, whatever the group
#'   peakOccupancyTable(results, by = "samples")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotPeakOccupancy}}, \code{\link{loadConsensusPeaks}}, \code{\link{plotPeakUpset}}
#'
#' @importFrom dplyr group_by summarise arrange n
#' @importFrom rlang .data
#'
#' @export peakOccupancyTable

peakOccupancyTable <-
  function(results,
           contrast = NULL,
           groups = NULL,
           by = "group",
           set = NULL,
           FDR = NULL,
           log2FC = NULL) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    results <- .pickResults(results = results, contrast = contrast)

    by <- tolower(by[1])
    if (!(by %in% c("group", "samples"))) {
      stop("The 'by' parameter must be either 'group' or 'samples'.", call. = FALSE)
    }

    resultTable <- .prepareResultTable(results = results, set = set, FDR = FDR, log2FC = log2FC)

    #-------------------------------#
    # Class of every region         #
    #-------------------------------#
    resultTable$occupancy <- .peakOccupancyClass(resultTable = resultTable,
                                                 groups = groups,
                                                 by = by,
                                                 contrastGroups = results@contrast.groups)

    #-------------------------------#
    # Count the directions          #
    #-------------------------------#
    occupancyCounts <-
      dplyr::summarise(dplyr::group_by(resultTable, .data$occupancy),
                       n.regions = dplyr::n(),
                       down = sum(.data$diff.status == "down"),
                       null = sum(.data$diff.status == "null"),
                       up = sum(.data$diff.status == "up"),
                       .groups = "drop")

    occupancyCounts <- dplyr::arrange(occupancyCounts, .data$occupancy)
    occupancyCounts$percent.changed <- round(100 * (occupancyCounts$down + occupancyCounts$up) / occupancyCounts$n.regions, 1)

    return(as.data.frame(occupancyCounts, stringsAsFactors = FALSE))
  } # END function




#' @title plotPeakOccupancy
#'
#' @description Draws the table of \code{\link{peakOccupancyTable}} as stacked bars, one per occupancy class, split by the direction of the change.
#'
#' @param results \code{RegionSetDE.results} object, or a \code{RegionSetDE.resultsList} together with \code{contrast}.
#' @param contrast String with the name of a contrast, or its position, when \code{results} holds several of them. Default: \code{NULL}.
#' @param groups Character vector with the consensus groups the regions are classified by. Default: \code{NULL}, the two groups the contrast compares when they are among them, every group otherwise.
#' @param by String indicating what the regions are grouped by, either \code{"group"} or \code{"samples"}. Default: \code{"group"}.
#' @param set Character vector with the names of the region sets kept. Default: \code{NULL}, all of them.
#' @param FDR Numeric value with the adjusted p-value below which a region counts as changed. Default: \code{NULL}, the threshold the test was run with.
#' @param log2FC Numeric value with the log2 fold change a region must reach. Default: \code{NULL}, the threshold the test was run with.
#' @param proportion Logical value to indicate whether the bars must be scaled to the same height, showing the composition of each class rather than its size. Default: \code{FALSE}.
#' @param showCounts Logical value to indicate whether the number of regions must be written on the bars. Default: \code{TRUE}.
#' @param colours Named character vector with the colours of \code{down}, \code{null} and \code{up}. Default: \code{NULL}, the palette of the other plots of the package.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}, the contrast.
#' @param legendPosition String with the position of the legend. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#' @param returnData Logical value indicating whether the table behind the plot must be returned instead of the plot. Default: \code{FALSE}.
#'
#' @return A \code{ggplot} object, or a data.frame when \code{returnData} is \code{TRUE}.
#'
#' @details With \code{proportion = TRUE} the bars all reach the same height and what they show is the composition of each class. That is the version to read when the shared regions outnumber the group-specific ones by an order of magnitude, which they usually do, and the absolute bars leave the small classes invisible.
#'
#' @examples
#' if (requireNamespace("consensusRegions", quietly = TRUE)) {
#'   # AR binding without ligand and after 24 hours of R1881
#'   sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#'   sampleSheet <- dplyr::filter(sampleSheet, condition %in% c("DMSO", "R1881_24h"))
#'
#'   consensus <- loadConsensusPeaks(sampleSheet, groupBy = "condition", seqlevelsStyle = "Ensembl", verbose = FALSE)
#'
#'   counts <- countReads(consensus, sampleSheet = sampleSheet, verbose = FALSE)
#'   counts <- countBackground(counts, binSize = 10000, verbose = FALSE)
#'   counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#'   fit <- fitRegions(counts, design = ~ condition, verbose = FALSE)
#'   results <- testRegions(fit, contrast = c("condition", "R1881_24h", "DMSO"), verbose = FALSE)
#'
#'   plotPeakOccupancy(results)
#'   plotPeakOccupancy(results, proportion = TRUE)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{peakOccupancyTable}}, \code{\link{plotPeakUpset}}, \code{\link{plotVolcano}}
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text labs scale_fill_manual scale_y_continuous position_stack theme element_blank
#' @importFrom scales percent
#' @importFrom rlang .data
#'
#' @export plotPeakOccupancy

plotPeakOccupancy <-
  function(results,
           contrast = NULL,
           groups = NULL,
           by = "group",
           set = NULL,
           FDR = NULL,
           log2FC = NULL,
           proportion = FALSE,
           showCounts = TRUE,
           colours = NULL,
           title = NULL,
           subtitle = NULL,
           legendPosition = "right",
           baseSize = 12,
           returnData = FALSE) {

    #------------------------#
    # Table behind the plot  #
    #------------------------#
    results <- .pickResults(results = results, contrast = contrast)

    occupancyCounts <- peakOccupancyTable(results = results, groups = groups, by = by,
                                          set = set, FDR = FDR, log2FC = log2FC)

    if (isTRUE(returnData)) {
      return(occupancyCounts)
    }

    # One row per class and direction, which is the shape the stacked bars need
    plotTable <- data.frame(occupancy = rep(occupancyCounts$occupancy, times = 3),
                            diff.status = factor(rep(c("down", "null", "up"), each = nrow(occupancyCounts)),
                                                 levels = c("down", "null", "up")),
                            n.regions = c(occupancyCounts$down, occupancyCounts$null, occupancyCounts$up),
                            stringsAsFactors = FALSE)

    plotTable$occupancy <- factor(plotTable$occupancy, levels = occupancyCounts$occupancy)

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    statusColours <- .diffStatusColours(colours = colours)
    barPosition <- if (isTRUE(proportion)) {"fill"} else {"stack"}

    occupancyPlot <-
      ggplot2::ggplot(data = plotTable,
                      mapping = ggplot2::aes(x = .data$occupancy, y = .data$n.regions, fill = .data$diff.status)) +
      ggplot2::geom_col(position = barPosition, width = 0.7) +
      ggplot2::scale_fill_manual(values = statusColours, drop = FALSE) +
      ggplot2::labs(x = if (by == "group") {"consensus"} else {"samples with a peak"},
                    y = if (isTRUE(proportion)) {"share of the regions"} else {"regions"},
                    fill = "",
                    title = title,
                    subtitle = if (is.null(subtitle)) {results@contrast} else {subtitle}) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize)

    if (isTRUE(proportion)) {
      occupancyPlot <- occupancyPlot + ggplot2::scale_y_continuous(labels = scales::percent)
    } else {
      # A count of regions has no halves, so neither do the breaks
      occupancyPlot <- occupancyPlot +
        ggplot2::scale_y_continuous(breaks = function(axisLimits) {unique(floor(pretty(axisLimits)))},
                                    labels = function(breakValues) {format(breakValues, big.mark = ",", trim = TRUE)})
    }

    # The counts are written on the pieces they belong to, in whichever of black and white reads on them
    if (isTRUE(showCounts)) {
      plotTable$label.colour <- .contrastColour(colour = statusColours[as.character(plotTable$diff.status)])

      occupancyPlot <- occupancyPlot +
        ggplot2::geom_text(data = plotTable,
                           mapping = ggplot2::aes(label = ifelse(.data$n.regions > 0, format(.data$n.regions, big.mark = ",", trim = TRUE), "")),
                           position = ggplot2::position_stack(vjust = 0.5),
                           colour = plotTable$label.colour,
                           size = baseSize / 4, show.legend = FALSE) +
        ggplot2::theme(axis.ticks.x = ggplot2::element_blank())
    }

    return(occupancyPlot)
  } # END function




#' @title .peakOccupancyClass
#'
#' @description Gives every region the class it is summarised under, either the combination of consensus groups covering it or the number of samples that carried a peak on it.
#'
#' @param resultTable Data.frame with the results, carrying the \code{peak.*} columns.
#' @param groups Character vector with the consensus groups used, or \code{NULL}.
#' @param by String, either \code{"group"} or \code{"samples"}.
#' @param contrastGroups List with the \code{column} and the two \code{groups} of the contrast, possibly empty.
#'
#' @return A factor with one value per row of the table.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.peakOccupancyClass <-
  function(resultTable,
           groups = NULL,
           by = "group",
           contrastGroups = list()) {

    #-------------------------------#
    # Number of samples             #
    #-------------------------------#
    if (by == "samples") {
      if (!("peak.samples" %in% colnames(resultTable))) {
        stop("The results carry no 'peak.samples' column, the regions were not built by loadConsensusPeaks().", call. = FALSE)
      }

      sampleCounts <- as.integer(resultTable$peak.samples)
      return(factor(sampleCounts, levels = sort(unique(sampleCounts))))
    }

    #-------------------------------#
    # Groups of the consensus       #
    #-------------------------------#
    occupancyColumns <- setdiff(grep("^peak\\.", colnames(resultTable), value = TRUE), c("peak.groups", "peak.samples"))

    if (length(occupancyColumns) == 0) {
      stop("The results carry no occupancy column, the regions were not built by loadConsensusPeaks().", call. = FALSE)
    }

    availableGroups <- sub("^peak\\.", "", occupancyColumns)

    # The two groups being contrasted are the ones the question is about, when the consensus knows them
    if (is.null(groups) && length(contrastGroups$groups) == 2 && all(make.names(contrastGroups$groups) %in% availableGroups)) {
      groups <- make.names(contrastGroups$groups)
    }

    if (is.null(groups)) {
      groups <- availableGroups
    }

    groups <- make.names(groups)
    absentGroups <- setdiff(groups, availableGroups)

    if (length(absentGroups) > 0) {
      stop("The following groups have no occupancy column in the results: ", paste(absentGroups, collapse = ", "),
           ". Available: ", paste(availableGroups, collapse = ", "), ".", call. = FALSE)
    }

    occupancyMatrix <- as.matrix(resultTable[, paste0("peak.", groups), drop = FALSE])
    mode(occupancyMatrix) <- "logical"

    # Shared and group specific are the two readings people want, the combinations in between are named as they are
    classLabels <- apply(occupancyMatrix, MARGIN = 1,
                         FUN = function(occupied) {
                           if (all(occupied)) {return(if (length(groups) > 1) {"shared"} else {groups[1]})}
                           if (!any(occupied)) {return("none")}
                           if (sum(occupied) == 1) {return(paste(groups[occupied], "only"))}
                           return(paste(groups[occupied], collapse = " + "))
                         })

    # Shared first, then the single groups in the order they were given, then the mixtures, then the regions without a peak
    sharedLabel <- if (length(groups) > 1) {"shared"} else {groups[1]}
    singleLabels <- paste(groups, "only")
    mixedLabels <- setdiff(unique(classLabels), c(sharedLabel, singleLabels, "none"))
    classLevels <- intersect(c(sharedLabel, singleLabels, sort(mixedLabels), "none"), unique(classLabels))

    return(factor(classLabels, levels = classLevels))
  } # END function
