# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title plotPeakUpset
#'
#' @description Draws an UpSet plot of the regions of a \code{RegionSetDE} object built from peaks, showing which groups have a consensus peak on each region, or which samples have a peak on it. It replaces the Venn diagram, which stops being readable beyond three or four sets.
#'
#' @param regionSet \code{RegionSetDE} object returned by \code{\link{loadConsensusPeaks}}.
#' @param by String indicating the sets of the plot, either \code{"group"}, one set per group and its consensus, or \code{"sample"}, one set per sample and its peaks. Default: \code{"group"}.
#' @param set Character vector with the region sets to describe. Default: \code{NULL}, all of them.
#' @param topIntersections Numeric value with the largest number of intersections drawn, the biggest ones. Default: \code{30}.
#' @param groupOrder Character vector with the order of the groups. Default: \code{NULL}, alphabetical.
#' @param groupColours Character vector with one colour per group, named after the groups or in the order of \code{groupOrder}. With \code{by = "sample"} the samples of a group take shades of its colour. Default: \code{NULL}, the palette of the other plots of the package.
#' @param title String with the title. Default: \code{NULL}, the number of regions described.
#'
#' @return An UpSet plot, as a \code{Heatmap} object of \code{ComplexHeatmap} with one row per set and one column per intersection, drawn when printed and open to \code{ComplexHeatmap::draw} for the layout options.
#'
#' @details The rows are the regions of the object, the ones that will be counted, and a region counts in a set when it overlaps the consensus of the group or a peak of the sample. After \code{regionMode = "replace"} the rows are the regions of the user, and those overlapping no peak at all are reported in the title rather than as an empty intersection. The intersections are exclusive, as in any UpSet plot: a region is counted once, in the combination of sets holding it and in no other.
#'
#' @examples
#' if (requireNamespace("consensusRegions", quietly = TRUE) & requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#'   peakFiles <- list.files(system.file("extdata", package = "consensusRegions"),
#'                           pattern = "rep[0-9]\\.narrowPeak$", full.names = TRUE)
#'
#'   sampleSheet <- loadSampleSheet(data.frame(sample = c("A_1", "A_2", "A_3", "B_1", "B_2"),
#'                                             bam = c("A_1.bam", "A_2.bam", "A_3.bam", "B_1.bam", "B_2.bam"),
#'                                             peaks = peakFiles[c(1, 2, 3, 2, 3)],
#'                                             condition = c("A", "A", "A", "B", "B")),
#'                                  checkFiles = FALSE, verbose = FALSE)
#'
#'   regions <- loadConsensusPeaks(sampleSheet, groupBy = "condition", verbose = FALSE)
#'
#'   plotPeakUpset(regions, by = "group")
#'   plotPeakUpset(regions, by = "sample")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadConsensusPeaks}}, \code{\link{consensusData}}
#'
#' @importFrom GenomicRanges GRangesList granges
#' @importFrom IRanges overlapsAny
#' @importFrom BiocGenerics unique
#' @importFrom methods is
#' @importFrom grid gpar unit grid.rect grid.points grid.lines
#'
#' @export plotPeakUpset

plotPeakUpset <-
  function(regionSet,
           by = "group",
           set = NULL,
           topIntersections = 30,
           groupOrder = NULL,
           groupColours = NULL,
           title = NULL) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      stop("The 'ComplexHeatmap' package is needed to draw the UpSet plot.", call. = FALSE)
    }

    if (!methods::is(regionSet, "RegionSetDE")) {
      stop("The 'regionSet' parameter must be a RegionSetDE object built by loadConsensusPeaks().", call. = FALSE)
    }

    consensusList <- consensusData(regionSet)

    by <- as.character(by[1])
    if (!(by %in% c("group", "sample"))) {
      stop("The 'by' parameter must be either 'group' or 'sample'.", call. = FALSE)
    }

    if (!is.numeric(topIntersections) | length(topIntersections) != 1 || is.na(topIntersections) || topIntersections < 1) {
      stop("The 'topIntersections' parameter must be a positive number.", call. = FALSE)
    }

    #-------------------------------#
    # Regions to describe           #
    #-------------------------------#
    regionList <- as.list(regionSet@regions)

    if (!is.null(set)) {
      absentSets <- setdiff(set, names(regionList))
      if (length(absentSets) > 0) {
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
      }
      regionList <- regionList[set]
    }

    # A region shared by two sets of the user is one stretch of genome, counted once
    regionRanges <- BiocGenerics::unique(GenomicRanges::granges(unlist(GenomicRanges::GRangesList(unname(regionList)), use.names = FALSE)))

    #-------------------------------#
    # Sets and their colours        #
    #-------------------------------#
    groupLevels <- .resolveGroupOrder(groupOrder = groupOrder, groupLevels = names(consensusList$groups))
    groupPalette <- .upsetGroupColours(groupLevels = groupLevels, groupColours = groupColours)

    if (by == "group") {
      referenceList <- consensusList$groups[groupLevels]
      setColours <- groupPalette
    } else {
      sampleTable <- consensusList$samples
      sampleTable$group <- factor(sampleTable$group, levels = groupLevels)
      sampleTable <- sampleTable[order(sampleTable$group), , drop = FALSE]

      referenceList <- as.list(consensusList$peaks)[sampleTable$sample]
      setColours <- .groupShades(colTable = sampleTable, baseColours = unname(groupPalette))
    }

    #-------------------------------#
    # Membership of every region    #
    #-------------------------------#
    # matrix() keeps the shape when there is a single region, which vapply would flatten
    membershipMatrix <- matrix(vapply(referenceList,
                                      function(referenceRanges) {IRanges::overlapsAny(regionRanges, referenceRanges, ignore.strand = TRUE)},
                                      logical(length(regionRanges))),
                               nrow = length(regionRanges),
                               dimnames = list(NULL, names(referenceList)))

    regionsWithoutPeak <- sum(rowSums(membershipMatrix) == 0)

    if (regionsWithoutPeak == nrow(membershipMatrix)) {
      stop("None of the regions overlaps a peak, there is nothing to draw.", call. = FALSE)
    }

    #-------------------------------#
    # Combinations and their sizes  #
    #-------------------------------#
    # Every region falls in exactly one combination, written as a 0/1 code over the sets
    peakRegions <- membershipMatrix[rowSums(membershipMatrix) > 0, , drop = FALSE]
    regionCodes <- apply(peakRegions, 1, function(x) {paste(as.integer(x), collapse = "")})
    combinationSizes <- table(regionCodes)

    combinationCodes <- names(combinationSizes)
    combinationDegree <- vapply(strsplit(combinationCodes, ""), function(x) {sum(as.integer(x))}, numeric(1))

    # Largest first, and among equal sizes the combinations shared by more sets
    combinationOrder <- order(as.numeric(combinationSizes), combinationDegree, decreasing = TRUE)

    # Past a few dozen intersections the plot is unreadable, the largest ones carry the message
    combinationOrder <- combinationOrder[seq_len(min(length(combinationOrder), topIntersections))]
    combinationCodes <- combinationCodes[combinationOrder]
    combinationSizes <- as.numeric(combinationSizes)[combinationOrder]

    # One row per set, one column per combination
    combinationMatrix <- matrix(as.integer(unlist(strsplit(combinationCodes, ""))),
                                nrow = ncol(membershipMatrix),
                                dimnames = list(colnames(membershipMatrix), combinationCodes))

    setSizes <- colSums(membershipMatrix)

    #-------------------------------#
    # Draw                          #
    #-------------------------------#
    if (is.null(title)) {
      title <- paste(format(length(regionRanges), big.mark = ",", trim = TRUE), "regions, sets by", by)
      if (regionsWithoutPeak > 0) {
        title <- paste0(title, " (", format(regionsWithoutPeak, big.mark = ",", trim = TRUE), " without any peak)")
      }
    }

    # Rows in alternating shades, one dot per set and combination, a line joining the sets of a combination
    rowShades <- rep(c("grey96", "white"), length.out = nrow(combinationMatrix))

    upsetLayer <- function(j, i, x, y, w, h, fill) {
      nRows <- round(1 / as.numeric(h[1]))
      nColumns <- round(1 / as.numeric(w[1]))
      subMatrix <- matrix(ComplexHeatmap::pindex(combinationMatrix, i, j), nrow = nRows, byrow = FALSE)

      for (k in seq_len(nRows)) {
        grid::grid.rect(y = 1 - (k - 1) / nRows, height = 1 / nRows, just = "top",
                        gp = grid::gpar(fill = rowShades[k], col = NA))
      }

      grid::grid.points(x, y, size = grid::unit(3, "mm"), pch = 16,
                        gp = grid::gpar(col = ifelse(ComplexHeatmap::pindex(combinationMatrix, i, j) == 1, "grey20", "#CCCCCC")))

      for (k in seq_len(nColumns)) {
        if (sum(subMatrix[, k]) >= 2) {
          memberRows <- range(which(subMatrix[, k] > 0))
          grid::grid.lines(x = rep(k - 0.5, 2) / nColumns,
                           y = (nRows - memberRows + 0.5) / nRows,
                           gp = grid::gpar(col = "grey20", lwd = 2))
        }
      }
    }

    # The numbers stay horizontal above the bars, ComplexHeatmap tilts them by default on a column annotation
    topAnnotation <- ComplexHeatmap::HeatmapAnnotation(intersection_size = ComplexHeatmap::anno_barplot(combinationSizes,
                                                                                                       border = FALSE,
                                                                                                       add_numbers = TRUE,
                                                                                                       numbers_rot = 0,
                                                                                                       numbers_gp = grid::gpar(fontsize = 7),
                                                                                                       gp = grid::gpar(fill = "grey30", col = NA),
                                                                                                       height = grid::unit(3, "cm")),
                                                       annotation_label = "Intersection\nsize",
                                                       annotation_name_side = "left",
                                                       annotation_name_rot = 0,
                                                       annotation_name_gp = grid::gpar(fontsize = 9))

    rightAnnotation <- ComplexHeatmap::rowAnnotation(set_size = ComplexHeatmap::anno_barplot(setSizes,
                                                                                             border = FALSE,
                                                                                             add_numbers = TRUE,
                                                                                             numbers_gp = grid::gpar(fontsize = 7),
                                                                                             gp = grid::gpar(fill = setColours[names(referenceList)], col = NA),
                                                                                             width = grid::unit(2, "cm")),
                                                     annotation_label = "Set size",
                                                     annotation_name_side = "bottom",
                                                     annotation_name_gp = grid::gpar(fontsize = 9))

    # A plain matrix, not a comb_mat: on R-devel aperm() keeps the class, and apply() over a comb_mat
    # calls the comb_mat subsetting method again and again until the stack overflows
    upsetPlot <- ComplexHeatmap::Heatmap(combinationMatrix,
                                         col = c("0" = "#CCCCCC", "1" = "grey20"),
                                         rect_gp = grid::gpar(type = "none"),
                                         layer_fun = upsetLayer,
                                         cluster_rows = FALSE,
                                         cluster_columns = FALSE,
                                         show_heatmap_legend = FALSE,
                                         show_column_names = FALSE,
                                         row_names_side = "left",
                                         row_names_gp = grid::gpar(fontsize = 9),
                                         top_annotation = topAnnotation,
                                         right_annotation = rightAnnotation,
                                         column_title = title,
                                         column_title_gp = grid::gpar(fontsize = 10))

    return(upsetPlot)
  } # END function




#' @title .upsetGroupColours
#'
#' @description Gives one colour per group, from the user when given and from the palette of the package otherwise.
#'
#' @param groupLevels Character vector with the groups, ordered.
#' @param groupColours Character vector with the colours of the user, named or in the order of the groups, or \code{NULL}.
#'
#' @return A named character vector with one colour per group.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.upsetGroupColours <-
  function(groupLevels,
           groupColours = NULL) {

    # Blue then red first, the same families the other plots of the package give to the groups
    if (is.null(groupColours)) {
      groupColours <- c("#2166AC", "#B2182B", "#1B7837", "#E08214", "#762A83", "#01665E")
    }

    if (!is.null(names(groupColours))) {
      absentGroups <- setdiff(groupLevels, names(groupColours))
      if (length(absentGroups) > 0) {
        stop("The following groups have no colour in 'groupColours': ", paste(absentGroups, collapse = ", "), ".", call. = FALSE)
      }
      groupColours <- groupColours[groupLevels]
    } else {
      groupColours <- rep(groupColours, length.out = length(groupLevels))
    }

    names(groupColours) <- groupLevels
    return(groupColours)
  } # END function
