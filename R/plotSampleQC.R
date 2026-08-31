#' @title plotRegionPCA
#'
#' @description Places the samples on the first principal components of the region signal, which is the fastest way to see whether the conditions separate, whether the replicates pair, and whether either of those is really the sequencing depth in disguise.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit} or any result object of the package.
#' @param set Character vector with the names of the region sets used. Default: \code{NULL}, all of them.
#' @param contrast String with the name of a contrast, or its position, when \code{object} holds several of them. Default: \code{NULL}.
#' @param colourBy String with the name of a \code{colData} column driving the colour. Default: \code{NULL}.
#' @param shapeBy String with the name of a \code{colData} column driving the shape. Default: \code{NULL}.
#' @param labelBy String with the name of a \code{colData} column written next to the points, or \code{"sample"}. Default: \code{"sample"}.
#' @param useOffsets Logical value to indicate whether the normalisation stored in the object must be applied. Default: \code{TRUE}.
#' @param compareOffsets Logical value to indicate whether the same ordination must be drawn twice, once with the normalisation and once on the library sizes alone. Default: \code{FALSE}.
#' @param facetBySet Logical value to indicate whether each region set must get its own ordination. Default: \code{FALSE}.
#' @param dimensions Numeric vector of length two with the components drawn. Default: \code{c(1, 2)}.
#' @param topRegions Numeric value with the number of most variable regions the ordination is computed on. Default: \code{2000}.
#' @param pointSize Numeric value with the size of the points. Default: \code{3}.
#' @param colours Named character vector with the colours. Default: \code{NULL}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}.
#' @param legendPosition String with the position of the legend. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object, carrying the coordinates and the variance explained as the \code{pca} attribute, in the \code{x.variance} and \code{y.variance} columns. The variance explained is written on the axis titles when a single panel is drawn, and in the panel labels when there are several, since every panel recomputes its own components.
#'
#' @details \code{compareOffsets} is the argument worth using. Scaling factors estimated outside the object, from a spike-in or a greenlist, impose a grouping of their own, and when that grouping happens to match the replicates it is indistinguishable from a batch effect until the two ordinations are put side by side. A separation that survives the normalisation being removed is in the data; one that appears only with it is the factors writing themselves into the ordination, and blocking on it would be blocking on an artefact.
#'
#' The regions are the same in both panels, chosen once by variance on the normalised values, so what differs between them is the transformation and not the selection. Restricting to the most variable rows is what makes an ordination read the structure rather than the depth, and \code{topRegions} controls how aggressively.
#'
#' Marks and assays should not share an ordination any more than they share a model. Split with \code{\link{splitSamples}} first, or use \code{facetBySet} when the sets themselves are the question.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#' plotRegionPCA(counts, colourBy = "condition", shapeBy = "sex")
#'
#' # Restricted to one set, which is where the strain effect should show
#' plotRegionPCA(counts, set = "promoterCpG", colourBy = "condition")

#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotSampleCorrelation}}, \code{\link{normalizeCounts}}, \code{\link{splitSamples}}
#'
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline facet_wrap facet_grid labs theme element_blank element_rect element_line scale_colour_manual guides guide_legend
#' @importFrom stats prcomp
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export plotRegionPCA

plotRegionPCA <-
  function(object,
           set = NULL,
           contrast = NULL,
           colourBy = NULL,
           shapeBy = NULL,
           labelBy = "sample",
           useOffsets = TRUE,
           compareOffsets = FALSE,
           facetBySet = FALSE,
           dimensions = c(1, 2),
           topRegions = 2000,
           pointSize = 3,
           colours = NULL,
           title = NULL,
           subtitle = NULL,
           legendPosition = "right",
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    counts <- .resolveCounts(object = object, counts = NULL, contrast = contrast)$counts

    if (ncol(counts) < 3) {
      stop("At least three samples are needed for an ordination.", call. = FALSE)
    }

    if (length(dimensions) != 2 | any(dimensions < 1)) {
      stop("The 'dimensions' parameter must hold two positive component numbers.", call. = FALSE)
    }

    colTable <- .sampleTable(counts = counts, colourBy = colourBy, shapeBy = shapeBy, labelBy = labelBy)
    panelSets <- .panelSets(counts = counts, set = set, facetBySet = facetBySet)
    offsetPanels <- if (isTRUE(compareOffsets)) {c(TRUE, FALSE)} else {useOffsets}

    #-------------------------------#
    # One ordination per panel      #
    #-------------------------------#
    coordinateList <- list()

    for (setName in names(panelSets)) {
      # The rows are chosen once, on the normalised values, so the panels differ by transformation and not by selection
      rowIndex <- .topVariableRows(counts = counts, rowIndex = panelSets[[setName]], topRegions = topRegions)

      for (offsetFlag in offsetPanels) {
        logMatrix <- .countsLogMatrix(counts = counts, useOffsets = offsetFlag)[rowIndex, , drop = FALSE]

        pcaObject <- stats::prcomp(x = t(logMatrix), center = TRUE, scale. = FALSE)
        varianceShare <- round(100 * pcaObject$sdev^2 / sum(pcaObject$sdev^2), 1)

        if (max(dimensions) > ncol(pcaObject$x)) {
          stop(paste0("Only ", ncol(pcaObject$x), " components exist, which is fewer than requested."), call. = FALSE)
        }

        offsetLabel <- if (isTRUE(offsetFlag)) {"normalised"} else {"library size only"}

        coordinateList[[paste(setName, offsetLabel)]] <-
          data.frame(sample = colnames(counts),
                     x.value = pcaObject$x[, dimensions[1]],
                     y.value = pcaObject$x[, dimensions[2]],
                     x.variance = varianceShare[dimensions[1]],
                     y.variance = varianceShare[dimensions[2]],
                     region.set = setName,
                     offset.label = offsetLabel,
                     panel = sprintf("%s%s — PC%d %.1f%%, PC%d %.1f%%",
                                     if (setName == "all") {""} else {paste0(setName, ", ")},
                                     offsetLabel, dimensions[1], varianceShare[dimensions[1]],
                                     dimensions[2], varianceShare[dimensions[2]]),
                     n.regions = length(rowIndex),
                     stringsAsFactors = FALSE)
      }
    }

    plotTable <- merge(x = do.call(what = rbind, args = coordinateList), y = colTable, by = "sample", sort = FALSE)
    plotTable$offset.label <- factor(plotTable$offset.label, levels = c("normalised", "library size only"))

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    # Each panel recomputes its own components, so the variance explained can go on
    # the axes only when there is one of them. With several it stays in the strip.
    if (length(unique(plotTable$panel)) > 1) {
      xAxisLabel <- paste0("PC", dimensions[1])
      yAxisLabel <- paste0("PC", dimensions[2])
    } else {
      xAxisLabel <- sprintf("PC%d (%.1f%%)", dimensions[1], plotTable$x.variance[1])
      yAxisLabel <- sprintf("PC%d (%.1f%%)", dimensions[2], plotTable$y.variance[1])
    }

    pcaPlot <-
      ggplot2::ggplot(data = plotTable,
                      mapping = ggplot2::aes(x = .data$x.value, y = .data$y.value)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "gray60") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "gray60") +
      ggplot2::geom_point(mapping = ggplot2::aes(colour = .data$colour.group, shape = .data$shape.group),
                          size = pointSize, stroke = NA) +
      ggplot2::labs(x = xAxisLabel,
                    y = yAxisLabel,
                    colour = if (is.null(colourBy)) {""} else {colourBy},
                    shape = if (is.null(shapeBy)) {""} else {shapeBy},
                    title = title,
                    subtitle = subtitle,
                    caption = paste0("computed on the ", plotTable$n.regions[1], " most variable regions")) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = max(c(pointSize, 3))))) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize) +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.75, colour = "black"),
                     aspect.ratio = 1)

    if (is.null(colourBy)) {
      pcaPlot <- pcaPlot + ggplot2::guides(colour = "none")
    } else if (!is.null(colours)) {
      pcaPlot <- pcaPlot + ggplot2::scale_colour_manual(values = colours)
    }

    if (is.null(shapeBy)) {
      pcaPlot <- pcaPlot + ggplot2::guides(shape = "none")
    }

    if (length(unique(plotTable$panel)) > 1) {
      pcaPlot <- pcaPlot + ggplot2::facet_wrap(facets = ~ panel, scales = "free")
    }

    #-------------------------------#
    # Labels on the points          #
    #-------------------------------#
    if (!is.null(labelBy) && requireNamespace("ggrepel", quietly = TRUE)) {
      pcaPlot <- pcaPlot +
        ggrepel::geom_text_repel(mapping = ggplot2::aes(label = .data$point.label),
                                 size = baseSize / 4.5, colour = "grey20", max.overlaps = Inf, show.legend = FALSE)
    }

    attr(pcaPlot, "pca") <- plotTable
    return(pcaPlot)
  } # END function




#' @title plotSampleCorrelation
#'
#' @description Draws the pairwise correlation between the samples over the region signal, which with a handful of libraries often reads more clearly than an ordination: replicates of a condition should sit closer to each other than to anything else.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit} or any result object of the package.
#' @param set Character vector with the names of the region sets used. Default: \code{NULL}, all of them.
#' @param contrast String with the name of a contrast, or its position, when \code{object} holds several of them. Default: \code{NULL}.
#' @param method String with the correlation, one of \code{"spearman"}, \code{"pearson"} and \code{"kendall"}. Default: \code{"spearman"}.
#' @param groupBy String with the name of a \code{colData} column defining the groups whose within and between correlations are summarised in the panel label. Default: \code{NULL}.
#' @param useOffsets Logical value to indicate whether the normalisation stored in the object must be applied. Default: \code{TRUE}.
#' @param compareOffsets Logical value to indicate whether the same matrix must be drawn twice, once with the normalisation and once on the library sizes alone. Default: \code{FALSE}.
#' @param facetBySet Logical value to indicate whether each region set must get its own matrix. Default: \code{FALSE}.
#' @param cluster Logical value to indicate whether the samples must be ordered by hierarchical clustering rather than kept in the order of the object. Default: \code{TRUE}.
#' @param clusteringMethod String with the agglomeration passed to \code{stats::hclust}. Default: \code{"complete"}.
#' @param topRegions Numeric value with the number of most variable regions the correlation is computed on. Default: \code{NULL}, all of them.
#' @param excludeDiagonal Logical value to indicate whether the diagonal must be left empty. Default: \code{FALSE}.
#' @param palette Character vector with the colours of the scale. Default: \code{NULL}, \code{viridisLite::mako(100, direction = -1)}.
#' @param limits Numeric vector of length two with the range of the colour scale, either value possibly \code{NA} to take that end from the data. Values outside are drawn at the nearest end rather than dropped, and how many were is reported. Default: \code{NULL}, the range of the values off the diagonal.
#' @param showValues Logical value to indicate whether the correlations must be written in the cells. Default: \code{TRUE}.
#' @param valuesColour String with the colour of the written values. Default: \code{NULL}, black or white on each cell depending on how dark it is.
#' @param valuesSize Numeric value with the font size of the written values. Default: \code{2.5}.
#' @param digits Numeric value with the number of decimals written. Default: \code{2}.
#' @param title String with the title of the plot, rendered as markdown. Default: \code{NULL}.
#' @param subtitle String with the subtitle of the plot, rendered as markdown. Default: \code{NULL}.
#' @param legendPosition String with the position of the legend. Default: \code{"right"}.
#' @param baseSize Numeric value with the base font size. Default: \code{12}.
#'
#' @return A \code{ggplot} object, carrying the correlation matrices as the \code{correlation} attribute.
#'
#' @details The scale runs over the values off the diagonal rather than from zero to one, because every sample correlates with itself perfectly and every pair of libraries from the same assay correlates highly. A scale anchored at zero turns the whole matrix one shade and hides the differences that matter. \code{limits} takes that decision back, and either end can be left as \code{NA} to be read from the data: \code{c(NA, 1)} fixes the top at one and lets the bottom follow the values.
#'
#' Values outside \code{limits} are drawn at the nearest end of the scale rather than left blank, so a cell that falls below the floor still shows as the extreme colour. That hides how far below it went, which is why the number of cells it happened to is reported.
#'
#' The palette is sequential, since a correlation has a low end and a high end and nothing meaningful in the middle. A diverging scale with white at the centre reads that midpoint as an absence, which on a matrix where everything sits between 0.9 and 1 is exactly wrong.
#'
#' With \code{groupBy}, the mean correlation within a group and between groups is written in the panel label. Within above between is what a usable experiment looks like; the two being equal says the condition effect is small next to the replicate noise, and that is the answer about whether to block, regardless of what an ordination suggests.
#'
#' The clustering order is taken from the first panel and reused in the others, so that a comparison across \code{compareOffsets} shows the values changing rather than the rows moving. No dendrogram is drawn for the same reason: a single dendrogram cannot describe several panels, and one per panel would defeat the comparison.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#' plotSampleCorrelation(counts, groupBy = "condition")
#'
#' # Pearson on the CpG island promoters only
#' plotSampleCorrelation(counts, set = "promoterCpG", method = "pearson")

#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotRegionPCA}}, \code{\link{normalizeCounts}}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes geom_tile geom_text facet_wrap labs coord_fixed theme element_blank element_rect element_text scale_x_discrete scale_y_discrete scale_fill_gradientn after_scale
#' @importFrom stats cor hclust as.dist
#' @importFrom scales squish
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export plotSampleCorrelation

plotSampleCorrelation <-
  function(object,
           set = NULL,
           contrast = NULL,
           method = "spearman",
           groupBy = NULL,
           useOffsets = TRUE,
           compareOffsets = FALSE,
           facetBySet = FALSE,
           cluster = TRUE,
           clusteringMethod = "complete",
           topRegions = NULL,
           excludeDiagonal = FALSE,
           palette = NULL,
           limits = NULL,
           showValues = TRUE,
           valuesColour = NULL,
           valuesSize = 2.5,
           digits = 2,
           title = NULL,
           subtitle = NULL,
           legendPosition = "right",
           baseSize = 12) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    counts <- .resolveCounts(object = object, counts = NULL, contrast = contrast)$counts

    if (!(method %in% c("spearman", "pearson", "kendall"))) {
      stop("The 'method' parameter must be one of 'spearman', 'pearson', 'kendall'.", call. = FALSE)
    }

    if (ncol(counts) < 2) {
      stop("At least two samples are needed for a correlation.", call. = FALSE)
    }

    if (is.null(palette)) {
      if (!requireNamespace("viridisLite", quietly = TRUE)) {
        stop("The 'viridisLite' package is needed for the default palette, or pass one through 'palette'.", call. = FALSE)
      }
      palette <- viridisLite::mako(n = 100, direction = -1)
    }

    panelSets <- .panelSets(counts = counts, set = set, facetBySet = facetBySet)
    offsetPanels <- if (isTRUE(compareOffsets)) {c(TRUE, FALSE)} else {useOffsets}

    groupVector <- NULL
    if (!is.null(groupBy)) {
      colTable <- as.data.frame(SummarizedExperiment::colData(counts))
      if (!(groupBy %in% colnames(colTable))) {
        stop(paste0("The column '", groupBy, "' is absent from the colData."), call. = FALSE)
      }
      groupVector <- as.character(colTable[[groupBy]])
    }

    #-------------------------------#
    # One matrix per panel          #
    #-------------------------------#
    correlationList <- list()
    matrixList <- list()
    sampleOrder <- NULL

    for (setName in names(panelSets)) {
      rowIndex <- .topVariableRows(counts = counts, rowIndex = panelSets[[setName]], topRegions = topRegions)

      for (offsetFlag in offsetPanels) {
        logMatrix <- .countsLogMatrix(counts = counts, useOffsets = offsetFlag)[rowIndex, , drop = FALSE]
        correlationMatrix <- stats::cor(logMatrix, method = method)

        # The order comes from the first panel and is reused, so a comparison shows the values moving and not the rows
        if (is.null(sampleOrder)) {
          sampleOrder <- if (isTRUE(cluster) & ncol(correlationMatrix) > 2) {
            colnames(correlationMatrix)[stats::hclust(d = stats::as.dist(1 - correlationMatrix),
                                                      method = clusteringMethod)$order]
          } else {
            colnames(correlationMatrix)
          }
        }

        panelLabel <- paste0(if (setName == "all") {""} else {paste0(setName, ", ")},
                             if (isTRUE(offsetFlag)) {"normalised"} else {"library size only"},
                             .correlationSummary(correlationMatrix = correlationMatrix,
                                                 groupVector = groupVector, digits = digits))

        matrixList[[panelLabel]] <- correlationMatrix

        correlationList[[panelLabel]] <-
          data.frame(sample.x = rep(sampleOrder, times = length(sampleOrder)),
                     sample.y = rep(sampleOrder, each = length(sampleOrder)),
                     correlation = as.numeric(correlationMatrix[sampleOrder, sampleOrder]),
                     panel = panelLabel,
                     stringsAsFactors = FALSE)
      }
    }

    plotTable <- do.call(what = rbind, args = correlationList)

    if (isTRUE(excludeDiagonal)) {
      plotTable <- plotTable[plotTable$sample.x != plotTable$sample.y, , drop = FALSE]
    }

    plotTable$sample.x <- factor(plotTable$sample.x, levels = sampleOrder)
    plotTable$sample.y <- factor(plotTable$sample.y, levels = rev(sampleOrder))

    #-------------------------------#
    # A scale over the real range   #
    #-------------------------------#
    offDiagonal <- plotTable$correlation[as.character(plotTable$sample.x) != as.character(plotTable$sample.y)]

    scaleLimits <- .resolveScaleLimits(values = offDiagonal,
                                       limits = limits,
                                       drawnValues = plotTable$correlation,
                                       label = "correlation")

    #-------------------------------#
    # Build the plot                #
    #-------------------------------#
    correlationPlot <-
      ggplot2::ggplot(data = plotTable,
                      mapping = ggplot2::aes(x = .data$sample.x, y = .data$sample.y, fill = .data$correlation)) +
      ggplot2::geom_tile() +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::coord_fixed() +
      ggplot2::labs(x = NULL, y = NULL,
                    fill = paste0(toupper(substring(method, 1, 1)), substring(method, 2), "\ncoefficient"),
                    title = title,
                    subtitle = subtitle) +
      .resultsTheme(legendPosition = legendPosition, baseSize = baseSize, rotateX = TRUE) +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     axis.ticks = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.75, colour = "black"),
                     aspect.ratio = 1)

    #-------------------------------#
    # Values, readable on any cell  #
    #-------------------------------#
    if (isTRUE(showValues)) {
      correlationPlot <- correlationPlot +
        if (is.null(valuesColour)) {
          # after_scale reads the colour the cell was actually filled with, so the text follows the palette
          ggplot2::geom_text(mapping = ggplot2::aes(label = format(round(.data$correlation, digits), nsmall = digits),
                                                    colour = ggplot2::after_scale(.contrastColour(.data$fill))),
                             size = valuesSize)
        } else {
          ggplot2::geom_text(mapping = ggplot2::aes(label = format(round(.data$correlation, digits), nsmall = digits)),
                             colour = valuesColour, size = valuesSize)
        }
    }

    # The scale is added last so that after_scale has a filled cell to read
    correlationPlot <- correlationPlot +
      ggplot2::scale_fill_gradientn(colours = palette,
                                    limits = scaleLimits,
                                    oob = scales::squish,
                                    na.value = "grey90")

    if (length(unique(plotTable$panel)) > 1) {
      correlationPlot <- correlationPlot + ggplot2::facet_wrap(facets = ~ panel)
    }

    attr(correlationPlot, "correlation") <- matrixList
    return(correlationPlot)
  } # END function




#' @title .contrastColour
#'
#' @description Returns black or white for every colour given, whichever of the two reads against it.
#'
#' @param colour Character vector with the colours of the cells.
#'
#' @return A character vector of the same length, holding \code{"black"} and \code{"white"}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom farver get_channel
#'
#' @keywords internal

.contrastColour <-
  function(colour) {

    # The luminance in HCL is what the eye reads, so it is what decides between black and white text
    luminance <- farver::get_channel(colour = colour, channel = "l", space = "hcl")

    textColour <- rep("black", length(colour))
    textColour[!is.na(luminance) & luminance < 50] <- "white"

    return(textColour)
  } # END function




#' @title .resolveScaleLimits
#'
#' @description Works out the range of a colour scale, filling in either end from the data when it was left open, and saying how many values will be drawn at the ends rather than at their own position.
#'
#' @param values Numeric vector the default range is read from.
#' @param limits Numeric vector of length two, either element possibly \code{NA}, or \code{NULL}.
#' @param drawnValues Numeric vector with every value that will be drawn, used to count the ones falling outside. Default: \code{NULL}, \code{values}.
#' @param label String naming the quantity, used in the message. Default: \code{"value"}.
#'
#' @return A numeric vector of length two.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resolveScaleLimits <-
  function(values,
           limits = NULL,
           drawnValues = NULL,
           label = "value") {

    dataRange <- range(values, na.rm = TRUE)

    if (is.null(limits)) {
      return(dataRange)
    }

    if (length(limits) != 2) {
      stop("The 'limits' parameter must hold two values, either of which may be NA.", call. = FALSE)
    }

    # NA on one end means that end follows the data, which is what makes c(NA, 1) useful
    limits <- as.numeric(limits)
    limits[is.na(limits)] <- dataRange[is.na(limits)]

    if (limits[1] >= limits[2]) {
      stop("The lower limit must sit below the upper one.", call. = FALSE)
    }

    #-------------------------------#
    # Say what is being squashed    #
    #-------------------------------#
    drawnValues <- if (is.null(drawnValues)) {values} else {drawnValues}
    outsideCount <- sum(drawnValues < limits[1] | drawnValues > limits[2], na.rm = TRUE)

    if (outsideCount > 0) {
      message(paste0(outsideCount, " ", label, "s fall outside the scale and are drawn at its ends, ",
                     "so how far past they went is not visible."))
    }

    return(limits)
  } # END function


#' @title .countsLogMatrix
#'
#' @description Turns a counts object into log2 counts per million, with or without the normalisation it carries.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param useOffsets Logical value to indicate whether the stored normalisation must be applied.
#' @param priorCount Numeric value with the prior count added before taking the logarithm. Default: \code{2}.
#'
#' @return A numeric matrix with one row per region and one column per sample.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment assay colData
#'
#' @keywords internal

.countsLogMatrix <-
  function(counts,
           useOffsets = TRUE,
           priorCount = 2) {

    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, "counts"))

    librarySizes <- SummarizedExperiment::colData(counts)$library.size
    if (is.null(librarySizes) | any(is.na(librarySizes))) {
      librarySizes <- colSums(countMatrix)
    }

    offsetMatrix <- if (isTRUE(useOffsets)) {.fitOffsets(counts = counts, useOffsets = TRUE, verbose = FALSE)} else {NULL}

    if (is.null(offsetMatrix)) {
      offsetMatrix <- matrix(data = rep(log(librarySizes), each = nrow(countMatrix)),
                             nrow = nrow(countMatrix), ncol = ncol(countMatrix))
    } else {
      # The offsets are centred on the library sizes so that the values come out on a CPM scale
      offsetMatrix <- offsetMatrix - rowMeans(offsetMatrix) + mean(log(librarySizes))
    }

    logMatrix <- log2(countMatrix + priorCount) - offsetMatrix / log(2) + log2(1e6)
    dimnames(logMatrix) <- dimnames(countMatrix)

    return(logMatrix)
  } # END function




#' @title .panelSets
#'
#' @description Splits the rows of a counts object into the panels a sample level plot will draw.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param set Character vector with the region sets to keep, or \code{NULL}.
#' @param facetBySet Logical value indicating whether every set gets a panel of its own.
#'
#' @return A named list of integer vectors, one per panel.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment rowData
#'
#' @keywords internal

.panelSets <-
  function(counts,
           set = NULL,
           facetBySet = FALSE) {

    rowSets <- as.character(SummarizedExperiment::rowData(counts)$region.set)

    if (!is.null(set)) {
      absentSets <- setdiff(set, unique(rowSets))
      if (length(absentSets) > 0) {
        stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }
      keptRows <- which(rowSets %in% set)
    } else {
      keptRows <- seq_along(rowSets)
    }

    if (isFALSE(facetBySet)) {
      return(list(all = keptRows))
    }

    panelList <- lapply(unique(rowSets[keptRows]), function(setName) {keptRows[rowSets[keptRows] == setName]})
    names(panelList) <- unique(rowSets[keptRows])

    return(panelList)
  } # END function




#' @title .topVariableRows
#'
#' @description Picks the most variable rows of a panel, which is what makes an ordination read the structure between samples rather than the differences in depth.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param rowIndex Integer vector with the rows of the panel.
#' @param topRegions Numeric value with the number of rows kept, or \code{NULL} for all of them.
#'
#' @return An integer vector with the rows kept.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats var
#'
#' @keywords internal

.topVariableRows <-
  function(counts,
           rowIndex,
           topRegions = 2000) {

    # A NULL here means every row, and testing it against a length would leave the branch with nothing to read
    if (is.null(topRegions)) {
      return(rowIndex)
    }

    if (length(rowIndex) <= topRegions) {
      return(rowIndex)
    }

    # The variance is read off the normalised values, so both panels of a comparison work on the same rows
    logMatrix <- .countsLogMatrix(counts = counts, useOffsets = TRUE)[rowIndex, , drop = FALSE]
    rowVariance <- apply(logMatrix, MARGIN = 1, FUN = stats::var)

    return(rowIndex[order(rowVariance, decreasing = TRUE)[seq_len(topRegions)]])
  } # END function




#' @title .sampleTable
#'
#' @description Assembles the sample metadata a sample level plot maps onto colour, shape and labels.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param colourBy String with a \code{colData} column, or \code{NULL}.
#' @param shapeBy String with a \code{colData} column, or \code{NULL}.
#' @param labelBy String with a \code{colData} column, \code{"sample"}, or \code{NULL}.
#'
#' @return A data.frame with one row per sample.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment colData
#'
#' @keywords internal

.sampleTable <-
  function(counts,
           colourBy = NULL,
           shapeBy = NULL,
           labelBy = "sample") {

    colTable <- as.data.frame(SummarizedExperiment::colData(counts))
    colTable$sample <- colnames(counts)

    for (columnName in c(colourBy, shapeBy, labelBy)) {
      if (!is.null(columnName) && columnName != "sample" && !(columnName %in% colnames(colTable))) {
        stop(paste0("The column '", columnName, "' is absent from the colData."), call. = FALSE)
      }
    }

    return(data.frame(sample = colTable$sample,
                      colour.group = if (is.null(colourBy)) {"all"} else {as.character(colTable[[colourBy]])},
                      shape.group = if (is.null(shapeBy)) {"all"} else {as.character(colTable[[shapeBy]])},
                      point.label = if (is.null(labelBy)) {""} else if (labelBy == "sample") {colTable$sample} else {as.character(colTable[[labelBy]])},
                      stringsAsFactors = FALSE))
  } # END function




#' @title .correlationSummary
#'
#' @description Summarises a correlation matrix as the mean within a group against the mean between groups, for the label of a panel.
#'
#' @param correlationMatrix Numeric matrix of correlations.
#' @param groupVector Character vector with the group of every sample, or \code{NULL}.
#' @param digits Numeric value with the number of decimals written.
#'
#' @return A string, empty when no grouping was given.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.correlationSummary <-
  function(correlationMatrix,
           groupVector = NULL,
           digits = 3) {

    if (is.null(groupVector)) {
      return("")
    }

    samePair <- outer(groupVector, groupVector, FUN = "==")
    offDiagonal <- upper.tri(correlationMatrix)

    withinValues <- correlationMatrix[samePair & offDiagonal]
    betweenValues <- correlationMatrix[!samePair & offDiagonal]

    # A single sample per group leaves no within-group pair to average
    if (length(withinValues) == 0 | length(betweenValues) == 0) {
      return("")
    }

    return(sprintf("\nwithin %s, between %s",
                   format(round(mean(withinValues), digits), nsmall = digits),
                   format(round(mean(betweenValues), digits), nsmall = digits)))
  } # END function
