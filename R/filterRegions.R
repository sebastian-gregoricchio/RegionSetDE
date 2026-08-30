#' @title filterRegions
#'
#' @description Removes the rows of a \code{RegionSetDE.counts} object that carry too little signal to say anything about a contrast. The decision is taken on the average abundance alone, never on the variance or on a fold change, so that the rows kept are independent of the comparison that will be run on them afterwards.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param method String with the criterion, one of \code{"background"}, \code{"abundance"}, \code{"proportion"}, \code{"byExpr"} and \code{"manual"}. Default: \code{"background"}.
#' @param foldChange Numeric value with the enrichment over the background a region must reach, on the linear scale. Only for \code{method = "background"}. Default: \code{2}.
#' @param minCount Numeric value with the number of reads a region of reference width must carry, on average across the samples. Only for \code{method = "abundance"}. Default: \code{10}.
#' @param proportion Numeric value between 0 and 1 with the fraction of rows to keep. Only for \code{method = "proportion"}. Default: \code{0.5}.
#' @param design Design matrix or formula passed to \code{edgeR::filterByExpr}. Only for \code{method = "byExpr"}. Default: \code{NULL}.
#' @param group String with the name of a \code{colData} column holding the experimental groups, passed to \code{edgeR::filterByExpr}. Only for \code{method = "byExpr"}. Default: \code{NULL}.
#' @param keep Logical vector with one value per row, or a vector of row positions. Only for \code{method = "manual"}. Default: \code{NULL}.
#' @param byWidth Logical value to indicate whether the abundance must be brought to a common width before being compared to the threshold. Default: \code{TRUE}.
#' @param referenceWidth Numeric value with the width the abundances are scaled to. Default: \code{NULL}, the median width of the rows.
#' @param bySet Logical value to indicate whether the threshold of \code{method = "proportion"} must be computed inside each region set rather than over the whole object. Default: \code{TRUE}.
#' @param widthStrata Numeric value with the number of width strata used by \code{method = "proportion"} when \code{byWidth = TRUE}. Default: \code{5}.
#' @param wholeRegion Logical value to indicate whether a tiled region must be kept in full as soon as one of its tiles passes. Ignored when the counts were not tiled. Default: \code{FALSE}.
#' @param minTiles Numeric value with the number of tiles a region must retain to survive. Ignored when the counts were not tiled. Default: \code{1}.
#' @param assay String with the name of the assay holding the values used to compute the abundance. Default: \code{"counts"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.counts} object with fewer rows.
#'
#' @details The filter runs before \code{\link{fitRegions}} because the dispersion trend is fitted across the rows of the object: leave in a few thousand regions that see two reads each and the trend describes them rather than the regions being tested. It also feeds the multiple testing correction, since every row that survives costs power at the adjustment step.
#'
#' Width is the reason a plain count threshold does not work on arbitrary region sets. A 40 kb domain accumulates more reads than a 400 bp promoter window at the same signal density, so a threshold in reads keeps every broad region and drops every narrow one, whatever the biology. With \code{byWidth = TRUE} the abundance of each row is divided by its width and brought back to \code{referenceWidth}, and the threshold then reads as "this many reads in a region of that width". This matters here more than in a window based analysis, where every window has the same size by construction.
#'
#' The methods differ in what they compare against. \code{"background"} needs the bins from \code{\link{countBackground}} and keeps the regions whose abundance exceeds, by \code{foldChange}, what a stretch of genome of the same width would carry: an absolute statement, in the sense that it does not depend on which other regions were loaded. \code{"abundance"} is the simpler version of the same idea, with the threshold given directly in reads. \code{"proportion"} keeps the strongest fraction of the rows and is relative by construction, which is why \code{bySet} defaults to \code{TRUE}: filtering the sets together lets a uniformly weak set be removed entirely, and a set that no longer has rows cannot be tested at the set level later. \code{"byExpr"} hands the decision to \code{edgeR::filterByExpr}, which reads the group sizes from the design.
#'
#' On a tiled object the filter applies to the tiles. A region that loses every tile disappears, and the count of regions lost this way is reported. \code{wholeRegion = TRUE} keeps all the tiles of a region as soon as one of them passes, which preserves the span of the region at the cost of carrying the empty tiles through the fit; the Simes combination in \code{\link{testRegions}} pays for those tiles in multiplicity, so the default leaves them out.
#'
#' @examples
#' \dontrun{
#' counts <- filterRegions(counts, method = "background", foldChange = 3)
#'
#' # Ten reads in a region of 1 kb, whatever the actual widths are
#' counts <- filterRegions(counts, method = "abundance", minCount = 10, referenceWidth = 1000)
#'
#' # Keep the strongest half of each set, within width strata
#' counts <- filterRegions(counts, method = "proportion", proportion = 0.5, bySet = TRUE)
#'
#' counts <- filterRegions(counts, method = "byExpr", group = "condition")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countBackground}}, \code{\link{normalizeCounts}}, \code{\link{fitRegions}}
#'
#' @importFrom SummarizedExperiment assay assayNames colData rowData rowRanges
#' @importFrom S4Vectors metadata
#' @importFrom BiocGenerics width
#' @importFrom edgeR aveLogCPM filterByExpr
#' @importFrom dplyr mutate group_by ungroup filter bind_rows count n
#' @importFrom rlang .data
#' @importFrom stats median quantile model.matrix
#' @importFrom methods is
#'
#' @export filterRegions

filterRegions <-
  function(counts,
           method = "background",
           foldChange = 2,
           minCount = 10,
           proportion = 0.5,
           design = NULL,
           group = NULL,
           keep = NULL,
           byWidth = TRUE,
           referenceWidth = NULL,
           bySet = TRUE,
           widthStrata = 5,
           wholeRegion = FALSE,
           minTiles = 1,
           assay = "counts",
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    if (!(method %in% c("background", "abundance", "proportion", "byExpr", "manual"))) {
      stop("The 'method' parameter must be one of 'background', 'abundance', 'proportion', 'byExpr', 'manual'.", call. = FALSE)
    }

    if (!(assay %in% SummarizedExperiment::assayNames(counts))) {
      stop(paste0("The assay '", assay, "' is absent from the object."), call. = FALSE)
    }

    if (proportion <= 0 | proportion > 1) {
      stop("The 'proportion' parameter must lie between 0 and 1.", call. = FALSE)
    }

    if (nrow(counts) == 0) {
      stop("The 'counts' object holds no row.", call. = FALSE)
    }

    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, assay))

    librarySizes <- SummarizedExperiment::colData(counts)$library.size
    if (is.null(librarySizes) | any(is.na(librarySizes))) {
      librarySizes <- colSums(countMatrix)
    }

    rowTable <- as.data.frame(SummarizedExperiment::rowData(counts))
    rowWidths <- BiocGenerics::width(SummarizedExperiment::rowRanges(counts))

    if (is.null(referenceWidth)) {
      referenceWidth <- stats::median(rowWidths)
    }

    #-------------------------------#
    # Abundance of every row        #
    #-------------------------------#
    rowAbundance <- .regionAbundance(countMatrix = countMatrix,
                                     librarySizes = librarySizes,
                                     rowWidths = if (isTRUE(byWidth)) {rowWidths} else {NULL},
                                     referenceWidth = referenceWidth)

    filterTable <- dplyr::mutate(rowTable,
                                 row.index = seq_len(nrow(counts)),
                                 row.width = rowWidths,
                                 abundance = rowAbundance)

    #-------------------------------#
    # Apply the criterion           #
    #-------------------------------#
    if (method == "background") {
      backgroundBins <- S4Vectors::metadata(counts)$background

      if (is.null(backgroundBins)) {
        stop("No background bin is stored in the object, run 'countBackground' before filtering with this method.", call. = FALSE)
      }

      binMatrix <- as.matrix(SummarizedExperiment::assay(backgroundBins, 1))
      binWidths <- BiocGenerics::width(SummarizedExperiment::rowRanges(backgroundBins))

      # The bins are brought to the same width as the regions, so the comparison is one of signal density
      binAbundance <- .regionAbundance(countMatrix = binMatrix,
                                       librarySizes = librarySizes,
                                       rowWidths = binWidths,
                                       referenceWidth = referenceWidth)

      backgroundLevel <- stats::median(binAbundance, na.rm = TRUE)

      # Without byWidth the regions sit on their own width and the background has to follow them there
      if (isFALSE(byWidth)) {
        backgroundLevel <- backgroundLevel + log2(rowWidths / referenceWidth)
      }

      filterTable <- dplyr::mutate(filterTable, enrichment = .data$abundance - backgroundLevel)
      filterTable <- dplyr::mutate(filterTable, keep.row = .data$enrichment >= log2(foldChange))

      thresholdMessage <- paste0("enrichment over background >= ", round(log2(foldChange), 2),
                                 " log2 (background at ", round(stats::median(backgroundLevel), 2), " log2 CPM)")

    } else if (method == "abundance") {
      # The threshold is written in reads and turned into the same scale as the abundances
      abundanceThreshold <- edgeR::aveLogCPM(y = minCount, lib.size = mean(librarySizes))

      filterTable <- dplyr::mutate(filterTable, keep.row = .data$abundance >= abundanceThreshold)

      thresholdMessage <- paste0("abundance >= ", round(abundanceThreshold, 2), " log2 CPM (", minCount, " reads",
                                 if (isTRUE(byWidth)) {paste0(" over ", referenceWidth, " bp")} else {""}, ")")

    } else if (method == "proportion") {
      #-------------------------------#
      # Strata for the relative cut   #
      #-------------------------------#
      filterTable <- dplyr::mutate(filterTable,
                                   set.stratum = if (isTRUE(bySet)) {as.character(.data$region.set)} else {"all"},
                                   width.stratum = if (isTRUE(byWidth)) {.widthStratum(rowWidths, widthStrata)} else {"all"},
                                   stratum = paste(.data$set.stratum, .data$width.stratum, sep = "|"))

      # slice_max would break the ties arbitrarily, a quantile keeps every row sitting exactly on the cut
      filterTable <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(filterTable, .data$stratum),
                                                  stratum.threshold = stats::quantile(.data$abundance, probs = 1 - proportion, na.rm = TRUE)))

      filterTable <- dplyr::mutate(filterTable, keep.row = .data$abundance >= .data$stratum.threshold)

      thresholdMessage <- paste0("top ", round(proportion * 100), "% of ",
                                 length(unique(filterTable$stratum)), " strata")

    } else if (method == "byExpr") {
      designObject <- design
      if (inherits(design, "formula")) {
        designObject <- stats::model.matrix(object = design, data = as.data.frame(SummarizedExperiment::colData(counts)))
      }

      groupVector <- NULL
      if (!is.null(group)) {
        if (!(group %in% colnames(SummarizedExperiment::colData(counts)))) {
          stop(paste0("The column '", group, "' is absent from the colData."), call. = FALSE)
        }
        groupVector <- SummarizedExperiment::colData(counts)[[group]]
      }

      if (is.null(designObject) & is.null(groupVector)) {
        stop("The 'byExpr' method needs either 'design' or 'group'.", call. = FALSE)
      }

      keptRows <- edgeR::filterByExpr(y = countMatrix,
                                      design = designObject,
                                      group = groupVector,
                                      lib.size = librarySizes)

      filterTable <- dplyr::mutate(filterTable, keep.row = as.logical(keptRows))
      thresholdMessage <- "edgeR::filterByExpr"

    } else {
      if (is.null(keep)) {
        stop("The 'manual' method needs the 'keep' parameter.", call. = FALSE)
      }

      if (is.logical(keep)) {
        if (length(keep) != nrow(counts)) {
          stop("The 'keep' parameter must have one value per row.", call. = FALSE)
        }
        keptRows <- keep
      } else {
        keep <- as.integer(keep)
        if (any(is.na(keep)) | any(keep < 1) | any(keep > nrow(counts))) {
          stop("The 'keep' parameter contains positions outside the range of the rows.", call. = FALSE)
        }
        keptRows <- seq_len(nrow(counts)) %in% keep
      }

      filterTable <- dplyr::mutate(filterTable, keep.row = keptRows)
      thresholdMessage <- "manual selection"
    }

    filterTable$keep.row[is.na(filterTable$keep.row)] <- FALSE

    #-------------------------------#
    # Tiles back to their region    #
    #-------------------------------#
    isTiled <- counts@counting.level == "tile"
    lostRegions <- 0

    if (isTRUE(isTiled)) {
      filterTable <- dplyr::mutate(filterTable, region.key = paste(.data$region.set, .data$region.id, sep = "|"))

      filterTable <- dplyr::ungroup(dplyr::mutate(dplyr::group_by(filterTable, .data$region.key),
                                                  tiles.kept = sum(.data$keep.row),
                                                  tiles.total = dplyr::n()))

      if (isTRUE(wholeRegion)) {
        filterTable <- dplyr::mutate(filterTable, keep.row = .data$tiles.kept >= minTiles)
      } else {
        filterTable <- dplyr::mutate(filterTable, keep.row = .data$keep.row & .data$tiles.kept >= minTiles)
      }

      regionSummary <- dplyr::ungroup(dplyr::filter(dplyr::group_by(filterTable, .data$region.key), dplyr::row_number() == 1))
      lostRegions <- sum(regionSummary$tiles.kept < minTiles)
    }

    keptIndex <- filterTable$row.index[filterTable$keep.row]

    if (length(keptIndex) == 0) {
      stop("No row passes the filter, loosen the threshold.", call. = FALSE)
    }

    #-------------------------------#
    # Subset and log                #
    #-------------------------------#
    filteredCounts <- counts[keptIndex, ]

    inputCounts <- dplyr::count(rowTable, .data$region.set, name = "n.input")
    outputCounts <- dplyr::count(as.data.frame(SummarizedExperiment::rowData(filteredCounts)), .data$region.set, name = "n.output")

    stepLog <- dplyr::mutate(dplyr::left_join(inputCounts, outputCounts, by = "region.set"), step = "filterRegions")
    stepLog$n.output[is.na(stepLog$n.output)] <- 0

    filteredCounts@filtering.log <- dplyr::bind_rows(counts@filtering.log, stepLog)

    filteredCounts@parameters <- c(filteredCounts@parameters,
                                   list(filterRegions = list(method = method,
                                                             foldChange = foldChange,
                                                             minCount = minCount,
                                                             proportion = proportion,
                                                             byWidth = byWidth,
                                                             referenceWidth = referenceWidth,
                                                             bySet = bySet,
                                                             wholeRegion = wholeRegion,
                                                             minTiles = minTiles,
                                                             threshold = thresholdMessage)))

    #-------------------------------#
    # Report                        #
    #-------------------------------#
    if (isTRUE(verbose)) {
      message(paste0("Filter: ", thresholdMessage, "."))
      message(paste0("Kept ", length(keptIndex), " out of ", nrow(counts), " ", counts@counting.level, "s (",
                     round(100 * length(keptIndex) / nrow(counts), 1), "%)."))

      if (lostRegions > 0) {
        message(paste0(lostRegions, " regions lost every tile and have been removed."))
      }

      for (i in seq_len(nrow(stepLog))) {
        message(paste0("  ", stepLog$region.set[i], ": ", stepLog$n.output[i], " / ", stepLog$n.input[i]))
      }

      # A set reduced to a handful of rows cannot carry a set level test afterwards
      emptiedSets <- stepLog$region.set[stepLog$n.output < 10]
      if (length(emptiedSets) > 0) {
        warning(paste0("The following region sets keep fewer than 10 rows: ", paste(emptiedSets, collapse = ", "), "."), call. = FALSE)
      }
    }

    return(filteredCounts)
  } # END function




#' @title .regionAbundance
#'
#' @description Computes the average abundance of every row, optionally brought to a common width so that regions of different sizes can be compared to the same threshold.
#'
#' @param countMatrix Numeric matrix of counts.
#' @param librarySizes Numeric vector with the library sizes.
#' @param rowWidths Numeric vector with the width of every row, or \code{NULL} to leave the abundances on their own width.
#' @param referenceWidth Numeric value with the width the abundances are scaled to.
#' @param priorCount Numeric value with the prior count added before taking the logarithm. Default: \code{2}.
#'
#' @return A numeric vector of log2 CPM values.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom edgeR aveLogCPM
#'
#' @keywords internal

.regionAbundance <-
  function(countMatrix,
           librarySizes,
           rowWidths = NULL,
           referenceWidth = NULL,
           priorCount = 2) {

    rowAbundance <- edgeR::aveLogCPM(y = countMatrix, lib.size = librarySizes, prior.count = priorCount)

    if (is.null(rowWidths)) {
      return(rowAbundance)
    }

    if (any(rowWidths <= 0)) {
      stop("Some rows have a zero or negative width.", call. = FALSE)
    }

    # Reads accumulate with the length of the interval, dividing by the width turns the abundance into a density
    return(rowAbundance - log2(rowWidths / referenceWidth))
  } # END function




#' @title .widthStratum
#'
#' @description Splits a vector of widths into strata of comparable size, used to keep the relative filters from favouring the wide regions.
#'
#' @param rowWidths Numeric vector with the widths.
#' @param strataNumber Numeric value with the number of strata.
#'
#' @return A character vector with the stratum of every row.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats quantile
#'
#' @keywords internal

.widthStratum <-
  function(rowWidths,
           strataNumber = 5) {

    strataNumber <- as.integer(strataNumber[1])
    if (is.na(strataNumber) | strataNumber < 1) {
      stop("The 'widthStrata' parameter must be a positive integer.", call. = FALSE)
    }

    # A set of regions all of the same width gives one breakpoint, the strata then collapse to a single level
    breakPoints <- unique(stats::quantile(rowWidths, probs = seq(0, 1, length.out = strataNumber + 1), na.rm = TRUE))

    if (length(breakPoints) < 2) {
      return(rep("all", length(rowWidths)))
    }

    return(as.character(cut(x = rowWidths, breaks = breakPoints, include.lowest = TRUE, labels = FALSE)))
  } # END function
