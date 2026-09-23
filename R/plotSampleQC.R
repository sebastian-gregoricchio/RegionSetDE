# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title plotRegionPCA
#'
#' @description Places the samples on the first principal components of the region signal, which is the fastest way to see whether the conditions separate, whether the replicates pair, and whether either of those is really the sequencing depth in disguise.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit} or any result object of the package, or the list returned by \code{\link{computeSamplePCA}}, in which case \code{set}, \code{useOffsets}, \code{compareOffsets}, \code{facetBySet} and \code{topRegions} are those of the computation.
#' @param set Character vector with the names of the region sets used. Default: \code{NULL}, all of them.
#' @param contrast String with the name of a contrast, or its position, when \code{object} holds several of them. Default: \code{NULL}.
#' @param colourBy String with the name of a \code{colData} column driving the colour. Default: \code{NULL}.
#' @param shapeBy String with the name of a \code{colData} column driving the shape. Default: \code{NULL}.
#' @param labelBy String with the name of a \code{colData} column written next to the points, or \code{"sample"}. Default: \code{"sample"}.
#' @param useOffsets Logical value to indicate whether the normalisation stored in the object must be applied, \code{FALSE} scaling the samples by their library sizes alone. Default: \code{TRUE}.
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
#' @seealso \code{\link{computeSamplePCA}}, \code{\link{plotSampleCorrelation}}, \code{\link{normalizeCounts}}, \code{\link{splitSamples}}
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
    if (length(dimensions) != 2 | any(dimensions < 1)) {
      stop("The 'dimensions' parameter must hold two positive component numbers.", call. = FALSE)
    }

    # A computation done beforehand is drawn as it is, the panels only make sense on an object with counts
    precomputed <- is.list(object) && all(c("scores", "variance") %in% names(object))

    #-------------------------------#
    # One ordination per panel      #
    #-------------------------------#
    pcaList <- list()

    # The transformation is named when the two are being compared, or when the values are not the normalised ones
    nameOffsets <- isTRUE(compareOffsets) | isFALSE(useOffsets)

    if (isTRUE(precomputed)) {
      pcaList[["computed"]] <- object
      panelNames <- "all"
      offsetLabels <- if (isFALSE(object$parameters$useOffsets)) {"library size only"} else {"normalised"}
      nameOffsets <- isFALSE(object$parameters$useOffsets)
      colTable <- .sampleTable(sampleAnnotation = object$scores, colourBy = colourBy, shapeBy = shapeBy, labelBy = labelBy)
    } else {
      counts <- .resolveCounts(object = object, counts = NULL, contrast = contrast)$counts

      if (ncol(counts) < 3) {
        stop("At least three samples are needed for an ordination.", call. = FALSE)
      }

      colTable <- .sampleTable(sampleAnnotation = .sampleAnnotation(counts = counts), colourBy = colourBy, shapeBy = shapeBy, labelBy = labelBy)
      panelSets <- .panelSets(counts = counts, set = set, facetBySet = facetBySet)
      offsetPanels <- if (isTRUE(compareOffsets)) {c(TRUE, FALSE)} else {useOffsets}

      .normalisationNotice(counts = counts, useOffsets = useOffsets | isTRUE(compareOffsets), verbose = TRUE)

      panelNames <- character(0)
      offsetLabels <- character(0)

      # The rows are chosen on the normalised values in every panel, so the panels differ by transformation and not by selection
      for (setName in names(panelSets)) {
        for (offsetFlag in offsetPanels) {
          pcaList[[length(pcaList) + 1]] <- computeSamplePCA(object = counts,
                                                             set = if (setName == "all") {set} else {setName},
                                                             useOffsets = offsetFlag,
                                                             topRegions = topRegions,
                                                             verbose = FALSE)
          panelNames <- c(panelNames, setName)
          offsetLabels <- c(offsetLabels, if (isTRUE(offsetFlag)) {"normalised"} else {"library size only"})
        }
      }
    }

    coordinateList <- list()

    for (i in seq_along(pcaList)) {
      scoreTable <- pcaList[[i]]$scores
      varianceShare <- round(pcaList[[i]]$variance$variance.percent, 1)
      componentNames <- paste0("PC", dimensions)

      if (!all(componentNames %in% colnames(scoreTable))) {
        stop("Only ", nrow(pcaList[[i]]$variance), " components exist, which is fewer than requested.", call. = FALSE)
      }

      panelLabel <- paste(c(if (panelNames[i] == "all") {NULL} else {panelNames[i]},
                            if (isTRUE(nameOffsets)) {offsetLabels[i]} else {NULL}), collapse = ", ")

      coordinateList[[i]] <-
        data.frame(sample = scoreTable$sample,
                   x.value = scoreTable[[componentNames[1]]],
                   y.value = scoreTable[[componentNames[2]]],
                   x.variance = varianceShare[dimensions[1]],
                   y.variance = varianceShare[dimensions[2]],
                   region.set = panelNames[i],
                   offset.label = offsetLabels[i],
                   panel = sprintf("%sPC%d %.1f%%, PC%d %.1f%%",
                                   if (panelLabel == "") {""} else {paste0(panelLabel, " - ")},
                                   dimensions[1], varianceShare[dimensions[1]],
                                   dimensions[2], varianceShare[dimensions[2]]),
                   n.regions = pcaList[[i]]$parameters$n.regions,
                   stringsAsFactors = FALSE)
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




#' @title computeSamplePCA
#'
#' @description Computes the principal components of the samples of an object on the log2 signal of its most variable regions, normalised or raw, and returns the coordinates, the variance explained and the loadings, ready for \code{\link{plotRegionPCA}} or for any other use.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit} or any result object of the package.
#' @param set Character vector with the names of the region sets used. Default: \code{NULL}, all of them.
#' @param contrast String with the name of a contrast, or its position, when \code{object} holds several of them. Default: \code{NULL}.
#' @param useOffsets Logical value to indicate whether the normalisation stored in the object must be applied, \code{FALSE} scaling the samples by their library sizes alone. Default: \code{TRUE}.
#' @param topRegions Numeric value with the number of most variable regions the ordination is computed on. Default: \code{2000}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with \code{scores}, a data.frame with one row per sample, its coordinates on every component and its \code{colData}; \code{variance}, a data.frame with the standard deviation, the percentage of variance and the cumulative percentage of every component; \code{loadings}, the matrix of the weights of the regions on the components; and \code{parameters}, with the region sets, the normalisation and the number of regions used.
#'
#' @details The values are log2 counts per million, with a prior count of 2 added to every count and, when \code{useOffsets = TRUE}, the scaling factors or the offsets stored by \code{\link{normalizeCounts}} applied. The regions are centred and not scaled, so a region weighs as much as it varies. An object that was never normalised is scaled by the library sizes whatever \code{useOffsets} says, and a message reports it. The most variable regions are chosen on the normalised values, so the raw and the normalised ordinations of one object are computed on the same rows.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#' samplePCA <- computeSamplePCA(counts, topRegions = 1000)
#' samplePCA$variance
#' head(samplePCA$scores[, c("sample", "PC1", "PC2", "condition")])
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotRegionPCA}}, \code{\link{computeSampleCorrelation}}
#'
#' @importFrom stats prcomp
#' @importFrom dplyr left_join
#'
#' @export computeSamplePCA

computeSamplePCA <-
  function(object,
           set = NULL,
           contrast = NULL,
           useOffsets = TRUE,
           topRegions = 2000,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    counts <- .resolveCounts(object = object, counts = NULL, contrast = contrast)$counts

    if (ncol(counts) < 3) {
      stop("At least three samples are needed for an ordination.", call. = FALSE)
    }

    .normalisationNotice(counts = counts, useOffsets = useOffsets, verbose = verbose)

    #-------------------------------#
    # Components of the samples     #
    #-------------------------------#
    rowIndex <- .topVariableRows(counts = counts,
                                 rowIndex = .panelSets(counts = counts, set = set, facetBySet = FALSE)$all,
                                 topRegions = topRegions)

    logMatrix <- .countsLogMatrix(counts = counts, useOffsets = useOffsets)[rowIndex, , drop = FALSE]

    # Centring costs one degree of freedom, so the last component of a full rank matrix is numerical noise and is left out
    pcaObject <- stats::prcomp(x = t(logMatrix), center = TRUE, scale. = FALSE,
                               rank. = min(ncol(counts) - 1, nrow(logMatrix)))
    varianceShare <- 100 * pcaObject$sdev^2 / sum(pcaObject$sdev^2)

    #-------------------------------#
    # Coordinates and annotation    #
    #-------------------------------#
    scoreTable <- data.frame(sample = colnames(counts), pcaObject$x, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    sampleTable <- .sampleAnnotation(counts = counts)

    # A colData column named like a component would otherwise be duplicated by the join
    collidingColumns <- setdiff(intersect(colnames(sampleTable), colnames(scoreTable)), "sample")
    if (length(collidingColumns) > 0) {
      colnames(sampleTable)[colnames(sampleTable) %in% collidingColumns] <- paste0(collidingColumns, ".sample")
    }

    scoreTable <- dplyr::left_join(scoreTable, sampleTable, by = "sample")

    varianceTable <- data.frame(component = colnames(pcaObject$x),
                                standard.deviation = pcaObject$sdev[seq_len(ncol(pcaObject$x))],
                                variance.percent = varianceShare[seq_len(ncol(pcaObject$x))],
                                cumulative.percent = cumsum(varianceShare)[seq_len(ncol(pcaObject$x))],
                                stringsAsFactors = FALSE)

    return(list(scores = scoreTable,
                variance = varianceTable,
                loadings = pcaObject$rotation,
                parameters = list(set = set,
                                  useOffsets = useOffsets,
                                  topRegions = topRegions,
                                  n.regions = length(rowIndex))))
  } # END function




#' @title computeSampleCorrelation
#'
#' @description Computes the correlation between the samples of an object on the log2 signal of its regions, normalised or raw, and returns it together with the annotation of the samples, ready for \code{\link{plotSampleCorrelation}} or for any other use.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit} or any result object of the package.
#' @param set Character vector with the names of the region sets used. Default: \code{NULL}, all of them.
#' @param contrast String with the name of a contrast, or its position, when \code{object} holds several of them. Default: \code{NULL}.
#' @param method String with the correlation, one of \code{"spearman"}, \code{"pearson"} and \code{"kendall"}. Default: \code{"spearman"}.
#' @param useOffsets Logical value to indicate whether the normalisation stored in the object must be applied, \code{FALSE} scaling the samples by their library sizes alone. Default: \code{TRUE}.
#' @param topRegions Numeric value with the number of most variable regions the correlation is computed on. Default: \code{NULL}, all of them.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with \code{correlation}, the matrix of the correlations between the samples; \code{samples}, the \code{colData} of the object with a \code{sample} column in the order of the matrix; and \code{parameters}, with the method, the normalisation, the region sets and the number of regions used.
#'
#' @details The values are log2 counts per million, with a prior count of 2 added to every count and, when \code{useOffsets = TRUE}, the scaling factors or the offsets stored by \code{\link{normalizeCounts}} applied. An object that was never normalised is scaled by the library sizes whatever \code{useOffsets} says, and a message reports it. The most variable regions are chosen on the normalised values, so the raw and the normalised correlations of one object are computed on the same rows.
#'
#' A correlation does not see a single factor per sample. Scaling a library moves its log values by a constant, which leaves the three correlations exactly where they were, so \code{useOffsets} changes the matrix only when the normalisation holds one offset per region, as \code{method = "loess"} and offsets supplied from outside do. An ordination is a different matter, and \code{\link{plotRegionPCA}} does move with the factors.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#' sampleCorrelation <- computeSampleCorrelation(counts, method = "pearson")
#' round(sampleCorrelation$correlation, 3)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{plotSampleCorrelation}}, \code{\link{computeSamplePCA}}
#'
#' @importFrom stats cor
#'
#' @export computeSampleCorrelation

computeSampleCorrelation <-
  function(object,
           set = NULL,
           contrast = NULL,
           method = "spearman",
           useOffsets = TRUE,
           topRegions = NULL,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    counts <- .resolveCounts(object = object, counts = NULL, contrast = contrast)$counts

    if (!(method[1] %in% c("spearman", "pearson", "kendall"))) {
      stop("The 'method' parameter must be one of 'spearman', 'pearson', 'kendall'.", call. = FALSE)
    }
    method <- method[1]

    if (ncol(counts) < 2) {
      stop("At least two samples are needed for a correlation.", call. = FALSE)
    }

    .normalisationNotice(counts = counts, useOffsets = useOffsets, verbose = verbose)

    #-------------------------------#
    # Correlation of the samples    #
    #-------------------------------#
    rowIndex <- .topVariableRows(counts = counts,
                                 rowIndex = .panelSets(counts = counts, set = set, facetBySet = FALSE)$all,
                                 topRegions = topRegions)

    logMatrix <- .countsLogMatrix(counts = counts, useOffsets = useOffsets)[rowIndex, , drop = FALSE]
    correlationMatrix <- stats::cor(logMatrix, method = method)

    return(list(correlation = correlationMatrix,
                samples = .sampleAnnotation(counts = counts),
                parameters = list(method = method,
                                  useOffsets = useOffsets,
                                  set = set,
                                  topRegions = topRegions,
                                  n.regions = length(rowIndex))))
  } # END function




#' @title plotSampleCorrelation
#'
#' @description Draws the correlation between the samples as a heatmap, clustered and annotated with any column of the sample table, such as the condition, the treatment or the replicate. The correlation can be computed on the normalised or on the raw signal, on every region or on the most variable ones.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit} or any result object of the package, or the list returned by \code{\link{computeSampleCorrelation}}.
#' @param set Character vector with the names of the region sets used. Default: \code{NULL}, all of them.
#' @param contrast String with the name of a contrast, or its position, when \code{object} holds several of them. Default: \code{NULL}.
#' @param method String with the correlation, one of \code{"spearman"}, \code{"pearson"} and \code{"kendall"}. Default: \code{"spearman"}.
#' @param groupBy String with the name of a \code{colData} column defining the groups whose within and between correlations are summarised in the title. Default: \code{NULL}.
#' @param annotationColumns Character vector with the \code{colData} columns drawn as annotation bars above and beside the heatmap, for instance \code{c("condition", "replicate")}. Default: \code{NULL}, the \code{groupBy} column when there is one.
#' @param annotationColours Named list with the colours of the annotation columns, as \code{ComplexHeatmap} takes them: a named vector for a categorical column, a \code{circlize::colorRamp2} function for a numeric one. The columns left out take the palette of the package. Default: \code{NULL}.
#' @param useOffsets Logical value to indicate whether the normalisation stored in the object must be applied, \code{FALSE} scaling the samples by their library sizes alone. Default: \code{TRUE}.
#' @param compareOffsets Logical value to indicate whether the same heatmap must be drawn twice side by side, once with the normalisation and once on the library sizes alone. Default: \code{FALSE}.
#' @param facetBySet Logical value to indicate whether each region set must get its own heatmap, side by side. Default: \code{FALSE}.
#' @param cluster Logical value to indicate whether the samples must be ordered by hierarchical clustering rather than kept in the order of the object. Default: \code{TRUE}.
#' @param clusteringMethod String with the agglomeration passed to \code{stats::hclust}. Default: \code{"complete"}.
#' @param showDendrogram Logical value to indicate whether the dendrogram of the clustering must be drawn. Default: \code{TRUE}.
#' @param topRegions Numeric value with the number of most variable regions the correlation is computed on. Default: \code{NULL}, all of them.
#' @param excludeDiagonal Logical value to indicate whether the diagonal must be left empty. Default: \code{FALSE}.
#' @param palette Character vector with the colours of the scale. Default: \code{NULL}, \code{viridisLite::mako(100, direction = -1)}.
#' @param limits Numeric vector of length two with the range of the colour scale, either value possibly \code{NA} to take that end from the data. Values outside are drawn at the nearest end rather than dropped, and how many were is reported. Default: \code{NULL}, the range of the values off the diagonal.
#' @param showValues Logical value to indicate whether the correlations must be written in the cells. Default: \code{TRUE}.
#' @param valuesColour String with the colour of the written values. Default: \code{NULL}, black or white on each cell depending on how dark it is.
#' @param digits Numeric value with the number of decimals written. Default: \code{2}.
#' @param title String with the title, written above the first heatmap. Default: \code{NULL}.
#' @param fontSize Numeric value with the font size of the names, the values being written slightly smaller. Default: \code{9}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{Heatmap} built by \code{ComplexHeatmap}, or a \code{HeatmapList} when \code{compareOffsets} or \code{facetBySet} draw several of them, drawn when printed. \code{ComplexHeatmap::draw} opens the layout, for instance \code{draw(x, heatmap_legend_side = "bottom")} or \code{draw(x, ht_gap = grid::unit(5, "mm"))} to push several panels apart. The matrices themselves come out of \code{\link{computeSampleCorrelation}}.
#'
#' @details The scale runs over the values off the diagonal rather than from zero to one, because every sample correlates with itself perfectly and every pair of libraries from the same assay correlates highly. A scale anchored at zero turns the whole matrix one shade and hides the differences that matter. \code{limits} takes that decision back, and either end can be left as \code{NA} to be read from the data: \code{c(NA, 1)} fixes the top at one and lets the bottom follow the values. Values outside \code{limits} are drawn at the nearest end, which hides how far past it they went, so the number of cells concerned is reported.
#'
#' The palette is sequential, since a correlation has a low end and a high end and nothing meaningful in the middle. A diverging scale with white at the centre reads that midpoint as an absence, which on a matrix where everything sits between 0.9 and 1 is exactly wrong.
#'
#' With \code{groupBy}, the mean correlation within a group and between groups is written in the title. Within above between is what a usable experiment looks like; the two being equal says the condition effect is small next to the replicate noise, and that is the answer about whether to block, regardless of what an ordination suggests.
#'
#' The clustering, on one minus the correlation, comes from the first heatmap and is reused by the others, so a comparison across \code{compareOffsets} or \code{facetBySet} shows the values changing rather than the samples moving. The dendrogram on the side is drawn once for the same reason, and so are the sample names, the columns of every panel following the order of the rows. Numeric annotation columns get a grey gradient, and turning one into a factor colours it by level instead, which suits a replicate number.
#'
#' \code{compareOffsets} answers less here than it does on an ordination. A correlation does not see a single factor per sample, so the two panels come out identical unless the normalisation holds one offset per region, as \code{method = "loess"} and offsets supplied from outside do. Whether the scaling factors are driving a grouping is a question for \code{\link{plotRegionPCA}}.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#' plotSampleCorrelation(counts, groupBy = "condition", annotationColumns = c("condition", "sex"))
#'
#' # The same samples before and after the normalisation, on the CpG island promoters only
#' plotSampleCorrelation(counts, set = "promoterCpG", method = "pearson", compareOffsets = TRUE)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{computeSampleCorrelation}}, \code{\link{plotRegionPCA}}, \code{\link{normalizeCounts}}
#'
#' @importFrom stats hclust as.dist as.dendrogram
#' @importFrom grid gpar
#'
#' @export plotSampleCorrelation

plotSampleCorrelation <-
  function(object,
           set = NULL,
           contrast = NULL,
           method = "spearman",
           groupBy = NULL,
           annotationColumns = NULL,
           annotationColours = NULL,
           useOffsets = TRUE,
           compareOffsets = FALSE,
           facetBySet = FALSE,
           cluster = TRUE,
           clusteringMethod = "complete",
           showDendrogram = TRUE,
           topRegions = NULL,
           excludeDiagonal = FALSE,
           palette = NULL,
           limits = NULL,
           showValues = TRUE,
           valuesColour = NULL,
           digits = 2,
           title = NULL,
           fontSize = 9,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE) | !requireNamespace("circlize", quietly = TRUE)) {
      stop("The 'ComplexHeatmap' and 'circlize' packages are needed to draw the correlation heatmap.", call. = FALSE)
    }

    if (is.null(palette)) {
      if (!requireNamespace("viridisLite", quietly = TRUE)) {
        stop("The 'viridisLite' package is needed for the default palette, or pass one through 'palette'.", call. = FALSE)
      }
      palette <- viridisLite::mako(n = 100, direction = -1)
    }

    if (is.null(annotationColumns) & !is.null(groupBy)) {
      annotationColumns <- groupBy
    }

    #-------------------------------#
    # One matrix per panel          #
    #-------------------------------#
    # A matrix computed beforehand is drawn as it is, the panels only make sense on an object with counts
    precomputed <- is.list(object) && all(c("correlation", "samples") %in% names(object))

    matrixList <- list()

    # The transformation is named when the two are being compared, or when the values are not the normalised ones
    nameOffsets <- isTRUE(compareOffsets) | isFALSE(useOffsets)

    if (isTRUE(precomputed)) {
      matrixList[["computed"]] <- object$correlation
      sampleTable <- object$samples
      method <- if (is.null(object$parameters$method)) {method} else {object$parameters$method}
      panelLabels <- if (isFALSE(object$parameters$useOffsets)) {"library size only"} else {""}
    } else {
      counts <- .resolveCounts(object = object, counts = NULL, contrast = contrast)$counts
      .normalisationNotice(counts = counts, useOffsets = useOffsets | isTRUE(compareOffsets), verbose = verbose)

      panelSets <- .panelSets(counts = counts, set = set, facetBySet = facetBySet)
      offsetPanels <- if (isTRUE(compareOffsets)) {c(TRUE, FALSE)} else {useOffsets}
      panelLabels <- character(0)

      for (setName in names(panelSets)) {
        for (offsetFlag in offsetPanels) {
          panelCorrelation <- computeSampleCorrelation(object = counts,
                                                       set = if (setName == "all") {set} else {setName},
                                                       method = method,
                                                       useOffsets = offsetFlag,
                                                       topRegions = topRegions,
                                                       verbose = FALSE)

          panelLabel <- paste(c(if (setName == "all") {NULL} else {setName},
                                if (isTRUE(nameOffsets)) {if (isTRUE(offsetFlag)) {"normalised"} else {"library size only"}} else {NULL}),
                              collapse = ", ")

          matrixList[[length(matrixList) + 1]] <- panelCorrelation$correlation
          panelLabels <- c(panelLabels, panelLabel)
        }
      }

      sampleTable <- .sampleAnnotation(counts = counts)
    }

    sampleOrder <- colnames(matrixList[[1]])
    sampleTable <- sampleTable[match(sampleOrder, sampleTable$sample), , drop = FALSE]

    for (columnName in unique(c(groupBy, annotationColumns))) {
      if (!(columnName %in% colnames(sampleTable))) {
        stop("The column '", columnName, "' is absent from the colData.", call. = FALSE)
      }
    }

    groupVector <- if (is.null(groupBy)) {NULL} else {as.character(sampleTable[[groupBy]])}

    #-------------------------------#
    # Order, scale and annotation   #
    #-------------------------------#
    # The first matrix decides the order, the others follow it so that the values move and not the samples
    sampleClustering <- if (isTRUE(cluster) & length(sampleOrder) > 2) {
      stats::as.dendrogram(stats::hclust(d = stats::as.dist(1 - matrixList[[1]]), method = clusteringMethod))
    } else {
      FALSE
    }

    offDiagonal <- unlist(lapply(matrixList, function(correlationMatrix) {correlationMatrix[row(correlationMatrix) != col(correlationMatrix)]}))
    drawnValues <- if (isTRUE(excludeDiagonal)) {offDiagonal} else {unlist(matrixList)}

    scaleLimits <- .resolveScaleLimits(values = offDiagonal, limits = limits, drawnValues = drawnValues, label = "correlation")

    # Two samples leave a single off-diagonal value, which is not a range to build a scale on
    if (diff(scaleLimits) == 0) {scaleLimits <- scaleLimits + c(-0.01, 0.01)}

    colourFunction <- circlize::colorRamp2(breaks = seq(scaleLimits[1], scaleLimits[2], length.out = length(palette)), colors = palette)

    annotationColourList <- .annotationColourList(sampleTable = sampleTable,
                                                  annotationColumns = annotationColumns,
                                                  annotationColours = annotationColours)

    legendTitle <- paste0(toupper(substring(method, 1, 1)), substring(method, 2), "\ncorrelation")

    #-------------------------------#
    # One heatmap per panel         #
    #-------------------------------#
    heatmapList <- list()

    for (i in seq_along(matrixList)) {
      panelMatrix <- matrixList[[i]]
      if (isTRUE(excludeDiagonal)) {diag(panelMatrix) <- NA}

      panelTitle <- paste0(panelLabels[i], .correlationSummary(correlationMatrix = matrixList[[i]], groupVector = groupVector, digits = digits))
      panelTitle <- sub("^\n", "", panelTitle)
      if (i == 1 & !is.null(title)) {
        panelTitle <- if (panelTitle == "") {title} else {paste(title, panelTitle, sep = "\n")}
      }

      # Legends and side annotations once, on the first heatmap, the others would only repeat them
      topAnnotation <- NULL
      leftAnnotation <- NULL
      if (length(annotationColumns) > 0) {
        topAnnotation <- ComplexHeatmap::HeatmapAnnotation(df = sampleTable[, annotationColumns, drop = FALSE],
                                                           col = annotationColourList,
                                                           show_legend = i == 1,
                                                           show_annotation_name = i == length(matrixList),
                                                           annotation_name_gp = grid::gpar(fontsize = fontSize),
                                                           simple_anno_size = grid::unit(3, "mm"))
        if (i == 1) {
          leftAnnotation <- ComplexHeatmap::rowAnnotation(df = sampleTable[, annotationColumns, drop = FALSE],
                                                          col = annotationColourList,
                                                          show_legend = FALSE,
                                                          show_annotation_name = FALSE,
                                                          simple_anno_size = grid::unit(3, "mm"))
        }
      }

      heatmapList[[i]] <- ComplexHeatmap::Heatmap(matrix = panelMatrix,
                                                  name = paste0("correlation_", i),
                                                  col = colourFunction,
                                                  na_col = "grey90",
                                                  cluster_rows = sampleClustering,
                                                  cluster_columns = sampleClustering,
                                                  show_row_dend = isTRUE(showDendrogram) & i == 1,
                                                  show_column_dend = isTRUE(showDendrogram),
                                                  show_row_names = i == length(matrixList),
                                                  show_column_names = length(matrixList) == 1,
                                                  row_names_gp = grid::gpar(fontsize = fontSize),
                                                  column_names_gp = grid::gpar(fontsize = fontSize),
                                                  top_annotation = topAnnotation,
                                                  left_annotation = leftAnnotation,
                                                  column_title = if (panelTitle == "") {NULL} else {panelTitle},
                                                  column_title_gp = grid::gpar(fontsize = fontSize + 1),
                                                  show_heatmap_legend = i == 1,
                                                  heatmap_legend_param = list(title = legendTitle,
                                                                              title_gp = grid::gpar(fontsize = fontSize),
                                                                              labels_gp = grid::gpar(fontsize = fontSize - 1)),
                                                  rect_gp = grid::gpar(col = "white", lwd = 0.5),
                                                  cell_fun = if (isTRUE(showValues)) {.correlationCellFunction(panelMatrix = panelMatrix,
                                                                                                               digits = digits,
                                                                                                               fontSize = fontSize,
                                                                                                               valuesColour = valuesColour)} else {NULL},
                                                  border = TRUE)
    }

    if (length(heatmapList) == 1) {
      return(heatmapList[[1]])
    }

    return(Reduce(f = `+`, x = heatmapList))
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
      message(outsideCount, " ", label, "s fall outside the scale and are drawn at its ends, ",
              "so how far past they went is not visible.")
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
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
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
#' @description Builds the table of the points of an ordination, with the groups driving the colour and the shape and the labels.
#'
#' @param sampleAnnotation Data.frame with a \code{sample} column and the annotation of the samples.
#' @param colourBy String with a column of the annotation, or \code{NULL}.
#' @param shapeBy String with a column of the annotation, or \code{NULL}.
#' @param labelBy String with a column of the annotation, \code{"sample"}, or \code{NULL}.
#'
#' @return A data.frame with the \code{sample}, \code{colour.group}, \code{shape.group} and \code{point.label} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.sampleTable <-
  function(sampleAnnotation,
           colourBy = NULL,
           shapeBy = NULL,
           labelBy = "sample") {

    colTable <- as.data.frame(sampleAnnotation, stringsAsFactors = FALSE)

    for (columnName in c(colourBy, shapeBy, labelBy)) {
      if (!is.null(columnName) && columnName != "sample" && !(columnName %in% colnames(colTable))) {
        stop("The column '", columnName, "' is absent from the colData.", call. = FALSE)
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




#' @title .normalisationNotice
#'
#' @description Says when the normalised values were asked for but the object carries no normalisation, in which case the library sizes are used.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param useOffsets Logical value indicating whether the normalised values were asked for.
#' @param verbose Logical value to indicate whether the message must be printed.
#'
#' @return Nothing, it prints a message at most.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.normalisationNotice <-
  function(counts,
           useOffsets = TRUE,
           verbose = TRUE) {

    if (isTRUE(useOffsets) & isTRUE(verbose) && is.null(.fitOffsets(counts = counts, useOffsets = TRUE, verbose = FALSE))) {
      message("No normalisation is stored in the object, the samples are scaled by their library sizes alone. Run normalizeCounts() first.")
    }

    return(invisible(TRUE))
  } # END function




#' @title .sampleAnnotation
#'
#' @description Returns the \code{colData} of a counts object as a data.frame, with the sample names in a first \code{sample} column.
#'
#' @param counts \code{RegionSetDE.counts} object.
#'
#' @return A data.frame with one row per sample, in the order of the object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr relocate
#'
#' @keywords internal

.sampleAnnotation <-
  function(counts) {

    sampleTable <- as.data.frame(SummarizedExperiment::colData(counts), optional = TRUE)
    sampleTable$sample <- colnames(counts)
    rownames(sampleTable) <- NULL

    return(dplyr::relocate(sampleTable, "sample"))
  } # END function




#' @title .annotationColourList
#'
#' @description Gives the colours of the annotation columns of a heatmap, from the user when given, a grey gradient for the numeric columns, and a palette per categorical column otherwise, the first one being the palette the other plots of the package give to the groups.
#'
#' @param sampleTable Data.frame with the annotation of the samples.
#' @param annotationColumns Character vector with the columns drawn.
#' @param annotationColours Named list with the colours given by the user, or \code{NULL}.
#'
#' @return A named list with the colours of every annotation column, as \code{ComplexHeatmap} takes them.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom grDevices hcl.colors
#'
#' @keywords internal

.annotationColourList <-
  function(sampleTable,
           annotationColumns,
           annotationColours = NULL) {

    # Blue then red first, as in the other plots, then palettes that do not repeat those hues
    paletteList <- list(c("#2166AC", "#B2182B", "#1B7837", "#E08214", "#762A83", "#01665E"),
                        grDevices::hcl.colors(8, palette = "Dark 3"),
                        grDevices::hcl.colors(8, palette = "Set 2"),
                        grDevices::hcl.colors(8, palette = "Pastel 1"))

    colourList <- list()

    for (k in seq_along(annotationColumns)) {
      columnName <- annotationColumns[k]

      if (!is.null(annotationColours[[columnName]])) {
        colourList[[columnName]] <- annotationColours[[columnName]]
        next
      }

      columnValues <- sampleTable[[columnName]]

      if (is.numeric(columnValues)) {
        # A constant column still needs two distinct breaks to build a gradient
        valueRange <- range(columnValues, na.rm = TRUE)
        if (valueRange[1] == valueRange[2]) {valueRange <- valueRange + c(-0.5, 0.5)}
        colourList[[columnName]] <- circlize::colorRamp2(breaks = valueRange, colors = c("grey92", "grey25"))
        next
      }

      columnLevels <- if (is.factor(columnValues)) {levels(droplevels(columnValues))} else {sort(unique(as.character(columnValues[!is.na(columnValues)])))}
      basePalette <- paletteList[[((k - 1) %% length(paletteList)) + 1]]

      levelColours <- if (length(columnLevels) <= length(basePalette)) {
        basePalette[seq_along(columnLevels)]
      } else {
        grDevices::hcl.colors(length(columnLevels), palette = "Dark 3")
      }

      names(levelColours) <- columnLevels
      colourList[[columnName]] <- levelColours
    }

    return(colourList)
  } # END function




#' @title .correlationCellFunction
#'
#' @description Builds the function writing the correlations in the cells of a heatmap, in black or white depending on how dark the cell is.
#'
#' @param panelMatrix Numeric matrix drawn by the heatmap.
#' @param digits Numeric value with the number of decimals written.
#' @param fontSize Numeric value with the font size of the names of the heatmap.
#' @param valuesColour String with a colour for every value, or \code{NULL}.
#'
#' @return A function with the signature \code{ComplexHeatmap} expects for \code{cell_fun}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom grid grid.text gpar
#'
#' @keywords internal

.correlationCellFunction <-
  function(panelMatrix,
           digits = 2,
           fontSize = 9,
           valuesColour = NULL) {

    # The matrix is captured here, a loop over the panels would otherwise leave every heatmap with the last one
    force(panelMatrix)

    return(function(j, i, x, y, width, height, fill) {
      cellValue <- panelMatrix[i, j]

      if (!is.na(cellValue)) {
        grid::grid.text(label = format(round(cellValue, digits), nsmall = digits),
                        x = x,
                        y = y,
                        gp = grid::gpar(fontsize = fontSize * 0.8,
                                        col = if (is.null(valuesColour)) {.contrastColour(colour = fill)} else {valuesColour}))
      }
    })
  } # END function
