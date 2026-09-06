#' @title estimateNullDispersion
#'
#' @description Estimates the biological variation between samples from a collection of rows assumed not to respond to the contrast, so that a design with no replicates has a dispersion to be tested against. The rows are usually the background bins, which cover the genome and should carry no treatment effect, but any region set believed to be invariant works the same way.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param source String with where the null rows come from, one of \code{"background"} (the bins stored by \code{\link{countBackground}}), \code{"regionSet"} (a set of the object) and \code{"supplied"}. Default: \code{"background"}.
#' @param regionSets Character vector with the names of the sets used as null rows. Only for \code{source = "regionSet"}. Default: \code{NULL}.
#' @param index Integer vector with the positions of the null rows. Only for \code{source = "supplied"}. Default: \code{NULL}.
#' @param samples Character vector with the samples the estimate is computed on. Default: \code{NULL}, all of them.
#' @param minCount Numeric value with the average count a null row must carry to be used. Default: \code{10}.
#' @param maxRows Numeric value with the number of null rows kept, drawn deterministically when there are more. Default: \code{50000}.
#' @param holdout Numeric value between 0 and 1 with the fraction of the null rows left out of the estimate, so that \code{\link{checkNullCalibration}} has rows the dispersion has not already been fitted to. Ignored when \code{\link{normalizeCounts}} has already split the background bins, whose split is reused instead. Default: \code{0.5}.
#' @param subset Integer vector restricting the null rows to a subset of the ones \code{source} selects. Default: \code{NULL}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with the \code{dispersion}, its square root as \code{bcv}, the \code{source}, the number of rows it was computed on, the samples used, \code{holdout.index}, the rows kept aside for the calibration check, and \code{holdout.type}, saying which steps those rows were kept out of.
#'
#' @details Without replicates nothing in the data measures how much two libraries differ for reasons that have nothing to do with the treatment, and every engine either refuses to fit or invents an answer. The way out is to assume that some rows do not respond, and to read the variation across those rows as the variation that would be seen between replicates. The estimate is a common dispersion fitted under an intercept-only model, which is exactly what \code{edgeR} recommends for an experiment with no replication, with the housekeeping genes replaced here by the background bins.
#'
#' That assumption is the whole estimate, so it is worth being deliberate about it. Background bins are the safest choice, since a treatment that changed the genome-wide average would have broken the normalisation long before it reached this point. A region set chosen because it looked flat in the data is the unsafe one: picking rows for their small fold change and then measuring the spread of those fold changes gives a dispersion biased towards zero, and p-values that follow it down.
#'
#' A plausible number is not a replicate, and the distinction is not a technicality. What this estimate measures is how much the two observed libraries differ from each other over rows assumed not to respond. What a replicated experiment measures is how much two animals, two donors or two cultures differ. The second contains the first and is usually several times larger, and no amount of null rows recovers it, because the variation that was never sampled is not in the data at all. Everything downstream is conditional on the assumption that the variability of the null rows stands in for the variability that was not observed, which is why \code{\link{checkNullCalibration}} exists and why it should be run before any of the output is believed. Treat the result as exploratory and say in the methods that it rests on this estimate.
#'
#' Checking the estimate against the rows it was fitted to would say nothing, so half the null rows are held out by default and travel back in \code{holdout.index} for that check to use. How much that holdout is worth depends on what came before it. Held out here alone, the rows are independent of the dispersion but not of the normalisation, since the background bins have already been through \code{csaw::normFactors} in \code{\link{normalizeCounts}}: they are dispersion holdouts and \code{holdout.type} says so. Splitting the bins earlier, with \code{normalizeCounts(backgroundHoldout = )}, keeps the same rows out of both steps, and that split is picked up here rather than a second one being drawn on top of it.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#' # The background bins cover the genome and carry no strain effect
#' nullDispersion <- estimateNullDispersion(counts, source = "background",
#'                                          verbose = FALSE)
#' nullDispersion
#'
#' # Any set believed to be invariant works the same way, as long as it has signal
#' fromRegions <- estimateNullDispersion(counts,
#'                                       source = "regionSet",
#'                                       regionSets = "geneBody",
#'                                       verbose = FALSE)
#' fromRegions
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{fitRegions}}, \code{\link{checkNullCalibration}}, \code{\link{countBackground}}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom edgeR DGEList scaleOffset estimateGLMCommonDisp aveLogCPM
#' @importFrom methods is
#'
#' @export estimateNullDispersion

estimateNullDispersion <-
  function(counts,
           source = "background",
           regionSets = NULL,
           index = NULL,
           samples = NULL,
           minCount = 10,
           maxRows = 50000,
           holdout = 0.5,
           subset = NULL,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    if (!(source %in% c("background", "regionSet", "supplied"))) {
      stop("The 'source' parameter must be one of 'background', 'regionSet', 'supplied'.", call. = FALSE)
    }

    if (!is.null(samples)) {
      counts <- selectSamples(counts = counts, samples = samples, dropNormalization = FALSE, verbose = FALSE)
    }

    if (ncol(counts) < 2) {
      stop("At least two samples are needed to measure the variation between them.", call. = FALSE)
    }

    if (holdout < 0 | holdout >= 1) {
      stop("The 'holdout' parameter must lie between 0 and 1, and leave something to estimate from.", call. = FALSE)
    }

    nullObject <- .nullMatrix(counts = counts, source = source, regionSets = regionSets,
                              index = index, subset = subset, minCount = minCount, maxRows = maxRows)

    if (length(nullObject$kept.rows) < 100) {
      stop("Only ", length(nullObject$kept.rows), " null rows reach ", minCount,
           " counts, which is too few to estimate a dispersion from.", call. = FALSE)
    }

    #-------------------------------#
    # Split off the check rows      #
    #-------------------------------#
    # Checking a dispersion against the rows it was fitted to would be calibrated by construction
    estimationRows <- seq_along(nullObject$kept.rows)
    holdoutIndex <- integer(0)
    holdoutType <- "none"

    storedHoldout <- if (source == "background") {S4Vectors::metadata(counts)$background.holdout} else {NULL}

    if (!is.null(storedHoldout) & length(storedHoldout) > 0) {
      # normalizeCounts already kept these bins out of the factors, splitting again would undo that
      holdoutIndex <- intersect(nullObject$kept.rows, as.integer(storedHoldout))
      estimationRows <- which(!(nullObject$kept.rows %in% holdoutIndex))
      holdoutType <- "normalization and dispersion"

      if (length(estimationRows) < 100) {
        stop("Only ", length(estimationRows), " null rows are left outside the normalisation holdout, ",
             "lower 'backgroundHoldout' in normalizeCounts.", call. = FALSE)
      }

    } else if (holdout > 0) {
      estimationRows <- .thinIndex(n = length(nullObject$kept.rows),
                                   maxPoints = ceiling(length(nullObject$kept.rows) * (1 - holdout)))
      holdoutIndex <- nullObject$kept.rows[setdiff(seq_along(nullObject$kept.rows), estimationRows)]
      holdoutType <- "dispersion"
    }

    countMatrix <- nullObject$counts[estimationRows, , drop = FALSE]
    offsetMatrix <- nullObject$offset[estimationRows, , drop = FALSE]

    #-------------------------------#
    # Variation under the null      #
    #-------------------------------#
    dgeList <- edgeR::DGEList(counts = countMatrix, lib.size = nullObject$library.size)
    dgeList <- edgeR::scaleOffset(y = dgeList, offset = offsetMatrix)

    # An intercept-only design leaves the whole between-sample variation in the residual, which is the point
    nullDesign <- matrix(data = 1, nrow = ncol(countMatrix), ncol = 1, dimnames = list(colnames(countMatrix), "(Intercept)"))

    commonDispersion <- edgeR::estimateGLMCommonDisp(y = dgeList, design = nullDesign)$common.dispersion

    if (!is.finite(commonDispersion) | commonDispersion <= 0) {
      stop("The dispersion came out zero or undefined, which usually means the null rows carry no variation at all.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message("Dispersion estimated on ", nrow(countMatrix), " ", source, " rows over ", ncol(countMatrix),
              " samples: ", signif(commonDispersion, 3), " (BCV ", signif(sqrt(commonDispersion), 3), ")",
              if (length(holdoutIndex) > 0) {paste0(", ", length(holdoutIndex), " rows held out for the check")} else {""}, ".")
      message("This is an assumption, not a replicate: it measures how these two libraries differ, not how two ",
              "biological samples would. Run checkNullCalibration before reading the p-values.")

      if (holdoutType == "dispersion") {
        message("The held-out rows are holdouts of the dispersion alone: they were already used to estimate the ",
                "normalisation factors. Split them earlier with normalizeCounts(backgroundHoldout = ) for a check ",
                "that sits outside the whole chain.")
      }
    }

    return(list(dispersion = commonDispersion,
                bcv = sqrt(commonDispersion),
                source = source,
                n.rows = nrow(countMatrix),
                holdout.index = holdoutIndex,
                holdout.type = holdoutType,
                samples = colnames(counts)))
  } # END function




#' @title .nullMatrix
#'
#' @description Collects the counts and the offsets of the rows a null estimate is computed on.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param source String with where the null rows come from.
#' @param regionSets Character vector with the sets used as null rows, or \code{NULL}.
#' @param index Integer vector with the positions of the null rows, or \code{NULL}.
#' @param subset Integer vector restricting the selected rows further, or \code{NULL}.
#' @param minCount Numeric value with the average count a row must carry.
#' @param maxRows Numeric value with the number of rows kept.
#'
#' @return A list with the filtered \code{counts} and \code{offset} matrices, the \code{library.size} vector, the \code{abundance} of the kept rows and \code{kept.rows}, their positions in the unfiltered selection.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom S4Vectors metadata
#' @importFrom edgeR aveLogCPM
#'
#' @keywords internal

.nullMatrix <-
  function(counts,
           source = "background",
           regionSets = NULL,
           index = NULL,
           subset = NULL,
           minCount = 10,
           maxRows = 50000) {

    librarySizes <- SummarizedExperiment::colData(counts)$library.size
    if (is.null(librarySizes) | any(is.na(librarySizes))) {
      librarySizes <- colSums(as.matrix(SummarizedExperiment::assay(counts, 1)))
    }

    #-------------------------------#
    # The rows                      #
    #-------------------------------#
    if (source == "background") {
      backgroundBins <- S4Vectors::metadata(counts)$background

      if (is.null(backgroundBins)) {
        stop("No background bin is stored in the object, run 'countBackground' first.", call. = FALSE)
      }
      countMatrix <- as.matrix(SummarizedExperiment::assay(backgroundBins, 1))

    } else if (source == "regionSet") {
      if (is.null(regionSets)) {
        stop("The 'regionSet' source needs the names of the sets in 'regionSets'.", call. = FALSE)
      }

      rowSets <- as.character(SummarizedExperiment::rowData(counts)$region.set)
      absentSets <- setdiff(regionSets, unique(rowSets))
      if (length(absentSets) > 0) {
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
      }
      countMatrix <- as.matrix(SummarizedExperiment::assay(counts, 1))[rowSets %in% regionSets, , drop = FALSE]

    } else {
      if (is.null(index)) {
        stop("The 'supplied' source needs the row positions in 'index'.", call. = FALSE)
      }
      index <- as.integer(index)
      if (any(is.na(index)) | any(index < 1) | any(index > nrow(counts))) {
        stop("The 'index' parameter contains positions outside the range of the rows.", call. = FALSE)
      }
      countMatrix <- as.matrix(SummarizedExperiment::assay(counts, 1))[index, , drop = FALSE]
    }

    #-------------------------------#
    # The offsets                   #
    #-------------------------------#
    # The null rows have to be normalised the way the regions are, or the estimate absorbs the depth differences
    scalingFactors <- SummarizedExperiment::colData(counts)$scaling.factor

    offsetMatrix <- if (is.null(scalingFactors) | all(is.na(scalingFactors))) {
      matrix(data = rep(log(librarySizes), each = nrow(countMatrix)), nrow = nrow(countMatrix), ncol = ncol(countMatrix))
    } else {
      matrix(data = rep(log(scalingFactors), each = nrow(countMatrix)), nrow = nrow(countMatrix), ncol = ncol(countMatrix))
    }

    dimnames(countMatrix) <- list(rownames(countMatrix), colnames(counts))

    #-------------------------------#
    # Restrict and filter           #
    #-------------------------------#
    if (!is.null(subset)) {
      subset <- as.integer(subset)
      if (any(is.na(subset)) | any(subset < 1) | any(subset > nrow(countMatrix))) {
        stop("The 'subset' parameter contains positions outside the selected null rows.", call. = FALSE)
      }
      countMatrix <- countMatrix[subset, , drop = FALSE]
      offsetMatrix <- offsetMatrix[subset, , drop = FALSE]
    }

    # A row seen a handful of times carries Poisson noise and nothing else, and it would drag the estimate up
    rowAbundance <- edgeR::aveLogCPM(y = countMatrix, lib.size = librarySizes)
    abundanceThreshold <- edgeR::aveLogCPM(y = minCount, lib.size = mean(librarySizes))

    keptRows <- which(rowAbundance >= abundanceThreshold)
    keptRows <- keptRows[.thinIndex(n = length(keptRows), maxPoints = maxRows)]

    return(list(counts = countMatrix[keptRows, , drop = FALSE],
                offset = offsetMatrix[keptRows, , drop = FALSE],
                abundance = rowAbundance[keptRows],
                kept.rows = keptRows,
                library.size = librarySizes))
  } # END function
