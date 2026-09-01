#' @title plotTopHeatmap
#'
#' @description Draws the signal of the regions responding most strongly to a contrast, as a heatmap with one block of rows per region set. The values come from the counts carried by the result, the samples are annotated from the \code{colData}, and the log2 fold change of every region is drawn next to it.
#'
#' @param results \code{RegionSetDE.results} or \code{RegionSetDE.resultsList} object.
#' @param n Numeric value with the number of regions taken from each region set. Default: \code{25}.
#' @param counts \code{RegionSetDE.counts} object holding the values, when the result carries none. Default: \code{NULL}.
#' @param contrast String with the name of the contrast to draw, or its position, when \code{results} holds several of them. Default: \code{NULL}.
#' @param set Character vector with the names of the region sets to draw. Default: \code{NULL}, all of them.
#' @param sortBy String with the ranking used to pick the regions, one of \code{"log2FC"}, \code{"FDR"} and \code{"stat"}. Default: \code{"log2FC"}.
#' @param direction String restricting the regions to one direction of change, one of \code{"both"}, \code{"up"} and \code{"down"}. Default: \code{"both"}.
#' @param FDR Numeric value with the adjusted p-value cut-off applied before the ranking. Default: \code{NULL}, the threshold stored in the object.
#' @param log2FC Numeric value with the absolute log2 fold change cut-off applied before the ranking. Default: \code{NULL}, the threshold stored in the object.
#' @param assay String with the name of the assay to draw. Default: \code{NULL}, the normalised assay when present.
#' @param scaleRows Logical value to indicate whether every row must be centred and scaled, so that the colour shows the shift between samples rather than how much signal the region carries. Default: \code{TRUE}.
#' @param annotationColumns Character vector with the \code{colData} columns drawn above the heatmap. Default: \code{NULL}, none.
#' @param showLog2FC Logical value to indicate whether the log2 fold change of every region must be drawn as a bar next to the rows. Default: \code{TRUE}.
#' @param border Logical value to indicate whether a frame must be drawn around each block of rows and around the barplot. Default: \code{TRUE}.
#' @param showRegionNames Logical value to indicate whether the region identifiers must be written. Default: \code{FALSE}.
#' @param clusterColumns Logical value to indicate whether the samples must be clustered rather than kept in the order of the object. Default: \code{FALSE}.
#' @param clusterRowsWithinSet Logical value to indicate whether the rows must be clustered inside each set block. Default: \code{TRUE}.
#' @param colours Character vector of three colours for the low, middle and high ends of the scale. Default: \code{NULL}.
#' @param limits Numeric vector of length two with the range of the colour scale, either value possibly \code{NA} to take that end from the data. Values outside are drawn at the nearest end rather than dropped, and how many were is reported. Default: \code{NULL}, plus and minus two for scaled rows and the first and last percentiles otherwise.
#' @param title String with the title of the heatmap. Default: \code{NULL}, the contrast.
#' @param ... Further arguments passed to \code{ComplexHeatmap::Heatmap}.
#'
#' @return A \code{Heatmap} object, drawn when printed.
#'
#' @details Ranking by \code{"log2FC"} sorts on the effect size among the regions that already passed the FDR cut-off, which is usually what a figure of this kind is meant to show. Ranking by \code{"FDR"} gives the regions the model is most certain about, and on a well powered object those are often the ones with the smallest effects.
#' \code{limits} fixes the range the colour scale covers, with \code{NA} on either end meaning that end follows the data. Anything past a limit is drawn at the end colour rather than left blank, which keeps the cell visible at the cost of hiding how far past it went, so the number of cells it happened to is reported.
#' \code{scaleRows} decides what the colour means, and the two answers are not interchangeable. Scaled rows show the pattern between the samples and put a region moving from 2 to 4 reads next to one moving from 200 to 400; unscaled rows keep the abundance visible and let a broad domain dominate the palette. The default is the scaled one because the figure is about a contrast, and the abundance is available in the \code{average.signal} column and in \code{\link{plotResultsMA}}.
#' The bars carry the colour of the direction they point to and no outline of their own, so at forty rows per block they stay readable; the frame around them is the annotation box, drawn to match the blocks of the heatmap.
#' On a tiled object the rows of the counts are tiles, so a region contributes several of them. Only the tile carrying the region level statistic is drawn, which is the representative tile chosen by the Simes combination in \code{\link{testRegions}}.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' # FDR = 1 ranks the regions instead of filtering them, which this small
#' # example dataset needs to fill a heatmap
#' plotTopHeatmap(results, n = 25, FDR = 1)
#'
#' plotTopHeatmap(results, n = 15, set = "promoterCpG", FDR = 1,
#'                annotationColumns = "condition")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{topRegions}}, \code{\link{plotRegion}}, \code{\link{plotVolcano}}
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom dplyr filter arrange group_by slice_head ungroup desc mutate
#' @importFrom rlang .data
#' @importFrom stats median
#' @importFrom methods is
#'
#' @export plotTopHeatmap

plotTopHeatmap <-
  function(results,
           n = 25,
           counts = NULL,
           contrast = NULL,
           set = NULL,
           sortBy = "log2FC",
           direction = "both",
           FDR = NULL,
           log2FC = NULL,
           assay = NULL,
           scaleRows = TRUE,
           annotationColumns = NULL,
           showLog2FC = TRUE,
           border = TRUE,
           showRegionNames = FALSE,
           clusterColumns = FALSE,
           clusterRowsWithinSet = TRUE,
           colours = NULL,
           limits = NULL,
           title = NULL,
           ...) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE) | !requireNamespace("circlize", quietly = TRUE)) {
      stop("The 'ComplexHeatmap' and 'circlize' packages are needed to draw the heatmap.", call. = FALSE)
    }

    if (!(sortBy %in% c("log2FC", "FDR", "stat"))) {
      stop("The 'sortBy' parameter must be one of 'log2FC', 'FDR', 'stat'.", call. = FALSE)
    }

    resolvedObject <- .resolveCounts(object = results, counts = counts, contrast = contrast)
    counts <- resolvedObject$counts
    results <- resolvedObject$results

    #-------------------------------#
    # Pick the regions              #
    #-------------------------------#
    setNames <- if (is.null(set)) {unique(as.character(results@results$region.set))} else {set}

    topTable <-
      do.call(what = rbind,
              args = lapply(setNames,
                            function(setName) {
                              setTop <- try(topRegions(results = results, n = n, set = setName,
                                                       sortBy = sortBy, direction = direction,
                                                       FDR = FDR, log2FC = log2FC),
                                            silent = TRUE)

                              # A set with nothing passing the thresholds drops out rather than stopping the figure
                              if (inherits(setTop, "try-error") | nrow(setTop) == 0) {
                                return(NULL)
                              }
                              return(setTop)
                            }))

    # rbind over an empty list returns NULL, which has no nrow to compare
    if (is.null(topTable)) {
      stop("No region passes the thresholds in any of the sets.", call. = FALSE)
    }

    if (nrow(topTable) == 0) {
      stop("No region passes the thresholds in any of the sets.", call. = FALSE)
    }

    topTable$region.key <- paste(topTable$region.set, topTable$region.id, sep = "|")

    #-------------------------------#
    # Rows of the counts to draw    #
    #-------------------------------#
    rowTable <- as.data.frame(SummarizedExperiment::rowData(counts))
    rowTable$row.index <- seq_len(nrow(counts))
    rowTable$region.key <- paste(rowTable$region.set, rowTable$region.id, sep = "|")

    matrixIndex <- .representativeRow(rowTable = rowTable, topTable = topTable, results = results)

    absentRegions <- sum(is.na(matrixIndex))
    if (absentRegions > 0) {
      warning(absentRegions, " regions of the result are absent from the counts and have been dropped.", call. = FALSE)
      topTable <- topTable[!is.na(matrixIndex), , drop = FALSE]
      matrixIndex <- matrixIndex[!is.na(matrixIndex)]
    }

    #-------------------------------#
    # Values                        #
    #-------------------------------#
    assay <- .defaultAssay(counts = counts, assay = assay)
    valueMatrix <- log2(as.matrix(SummarizedExperiment::assay(counts, assay))[matrixIndex, , drop = FALSE] + 1)
    rownames(valueMatrix) <- topTable$region.id
    colnames(valueMatrix) <- colnames(counts)

    legendName <- paste0("log2 ", assay)
    legendTitle <- bquote(log[2] ~ .(assay))

    if (isTRUE(scaleRows)) {
      # A row with the same value everywhere has no scale, dividing by it would fill the heatmap with NaN
      rowCentres <- rowMeans(valueMatrix)
      rowScales <- apply(valueMatrix, MARGIN = 1, FUN = stats::sd)
      rowScales[!is.finite(rowScales) | rowScales == 0] <- 1

      valueMatrix <- (valueMatrix - rowCentres) / rowScales
      legendName <- "z-score"
      legendTitle <- "z-score"
    }

    #-------------------------------#
    # Colours                       #
    #-------------------------------#
    paletteColours <- if (is.null(colours)) {c("#2166AC", "white", "#B2182B")} else {colours}
    if (length(paletteColours) != 3) {
      stop("The 'colours' parameter must hold three values, for the low, middle and high ends of the scale.", call. = FALSE)
    }

    defaultLimits <- if (isTRUE(scaleRows)) {
      c(-2, 2)
    } else {
      as.numeric(stats::quantile(valueMatrix, probs = c(0.01, 0.99), na.rm = TRUE))
    }

    scaleLimits <- .resolveScaleLimits(values = defaultLimits,
                                       limits = limits,
                                       drawnValues = as.numeric(valueMatrix),
                                       label = if (isTRUE(scaleRows)) {"z-score"} else {"value"})

    # colorRamp2 already draws anything past the ends at the end colour, which is the behaviour wanted here
    colourFunction <- circlize::colorRamp2(breaks = c(scaleLimits[1], mean(scaleLimits), scaleLimits[2]),
                                           colors = paletteColours)

    #-------------------------------#
    # Annotations                   #
    #-------------------------------#
    topAnnotation <- NULL
    if (!is.null(annotationColumns)) {
      colTable <- as.data.frame(SummarizedExperiment::colData(counts))
      absentColumns <- setdiff(annotationColumns, colnames(colTable))
      if (length(absentColumns) > 0) {
        stop("The following columns are absent from the colData: ", paste(absentColumns, collapse = ", "), ".", call. = FALSE)
      }

      topAnnotation <- ComplexHeatmap::HeatmapAnnotation(df = colTable[, annotationColumns, drop = FALSE],
                                                         annotation_name_side = "left")
    }

    rightAnnotation <- NULL
    if (isTRUE(showLog2FC)) {
      statusColours <- .diffStatusColours(colours = NULL)

      # col = NA drops the outline of every bar, the frame around the panel comes from 'border' instead
      barGraphics <- grid::gpar(fill = ifelse(topTable$log2FC > 0, statusColours[["up"]], statusColours[["down"]]),
                                col = NA)

      rightAnnotation <- ComplexHeatmap::rowAnnotation(log2FC = ComplexHeatmap::anno_barplot(x = topTable$log2FC,
                                                                                             baseline = 0,
                                                                                             gp = barGraphics,
                                                                                             border = border,
                                                                                             bar_width = 1,
                                                                                             axis_param = list(gp = grid::gpar(fontsize = 7))),
                                                       annotation_label = list(log2FC = expression(log[2] * "FC")),
                                                       annotation_name_rot = 0,
                                                       annotation_name_gp = grid::gpar(fontsize = 9))
    }

    #-------------------------------#
    # Build the heatmap             #
    #-------------------------------#
    heatmapObject <- ComplexHeatmap::Heatmap(matrix = valueMatrix,
                                             name = legendName,
                                             col = colourFunction,
                                             heatmap_legend_param = list(title = legendTitle),
                                             column_title = if (is.null(title)) {results@contrast} else {title},
                                             row_split = factor(topTable$region.set, levels = setNames),
                                             row_title_rot = 0,
                                             cluster_rows = clusterRowsWithinSet,
                                             cluster_columns = clusterColumns,
                                             show_row_names = showRegionNames,
                                             row_names_gp = grid::gpar(fontsize = 7),
                                             column_names_gp = grid::gpar(fontsize = 9),
                                             top_annotation = topAnnotation,
                                             right_annotation = rightAnnotation,
                                             border = border,
                                             ...)

    return(heatmapObject)
  } # END function




#' @title .representativeRow
#'
#' @description Finds, for every region of a result table, the row of the counts that carries its statistic. On a tiled object that is the tile the combination picked, and on a non-tiled one it is the region itself.
#'
#' @param rowTable Data.frame with the rows of the counts and their positions.
#' @param topTable Data.frame with the regions to draw.
#' @param results \code{RegionSetDE.results} object.
#'
#' @return An integer vector with one position per region, \code{NA} when the region is absent from the counts.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#' @keywords internal

.representativeRow <-
  function(rowTable,
           topTable,
           results) {

    isCombined <- isTRUE(results@combination$applied)

    # A region counted as a single row is its own representative, nothing has to be chosen
    if (!isCombined) {
      return(match(topTable$region.key, rowTable$region.key))
    }

    tileTable <- results@tiles
    tileTable$region.key <- paste(tileTable$region.set, tileTable$region.id, sep = "|")

    matchedIndex <-
      vapply(seq_len(nrow(topTable)),
             function(i) {
               candidateRows <- rowTable[rowTable$region.key == topTable$region.key[i], , drop = FALSE]

               if (nrow(candidateRows) == 0) {
                 return(NA_integer_)
               }

               # rep.tile.start comes from the tile the Simes combination reported, so the row matches the statistic
               if ("rep.tile.start" %in% colnames(topTable)) {
                 representativeRow <- candidateRows$row.index[
                   BiocGenerics::start(SummarizedExperiment::rowRanges(results@counts))[candidateRows$row.index] == topTable$rep.tile.start[i]]

                 if (length(representativeRow) == 1) {
                   return(as.integer(representativeRow))
                 }
               }

               return(as.integer(candidateRows$row.index[1]))
             },
             integer(1))

    return(matchedIndex)
  } # END function
