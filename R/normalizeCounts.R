#' @title normalizeCounts
#'
#' @description Estimates the scaling factors of a \code{RegionSetDE.counts} object and stores them together with a normalised assay. The factors can be computed from the counts themselves, from the background bins collected by \code{\link{countBackground}}, from a spike-in, or supplied by the user.
#'
#' @param counts \code{RegionSetDE.counts} object returned by \code{\link{countReads}}, \code{\link{countBigwig}} or \code{\link{loadCounts}}.
#' @param method String with the method used to estimate the factors, one among \code{"TMM"}, \code{"TMMwsp"}, \code{"RLE"}, \code{"upperQuartile"}, \code{"librarySize"}, \code{"background"}, \code{"loess"}, \code{"spikeIn"}, \code{"manual"} or \code{"none"}. Default: \code{"TMM"}.
#' @param scalingFactors Numeric vector with one factor per sample, required by the \code{"manual"} method. Named vectors are matched to the sample names and may cover samples absent from the object, unnamed ones must follow the column order. Default: \code{NULL}.
#' @param factorType String declaring how the values of \code{scalingFactors} must be applied, either \code{"division"} when the counts have to be divided by them, as for the size factors of DESeq2, or \code{"multiplication"} when they have to be multiplied, as for the scale factors of deeptools and of most spike-in protocols. Default: \code{"division"}.
#' @param spikeInCounts Numeric vector with the number of reads assigned to the exogenous genome in each sample, required by the \code{"spikeIn"} method. Default: \code{NULL}.
#' @param useRegionSets Character vector with the names of the region sets used to estimate the factors. Default: \code{NULL}, all of them.
#' @param minCount Numeric value with the minimum total count required for a region to take part in the estimation. Default: \code{1}.
#' @param referenceSample String or numeric position of the sample used as reference by the \code{"TMM"} and \code{"TMMwsp"} methods. Default: \code{NULL}, chosen by edgeR.
#' @param normalizedAssay String with the name given to the normalised assay. Default: \code{"norm.counts"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return The input \code{RegionSetDE.counts} object with three additions: the \code{norm.factor} and \code{scaling.factor} columns in the \code{colData}, the normalised values in the assay named after \code{normalizedAssay}, and, for the \code{"loess"} method, the log offsets in the \code{offset} assay.
#'
#' @details The raw counts are never overwritten, the normalisation only adds an assay and the factors beside it, so that the testing functions can keep working on the counts and their offsets.
#'
#' Whatever the method, \code{scaling.factor} always holds a divisor: the normalised values are the raw counts divided by it, and a sample sequenced more deeply than the others therefore gets a factor above one. This is the convention of the size factors of DESeq2 and the opposite of the scale factors written by \code{bamCoverage} or by the spike-in pipelines, which are meant to multiply the signal. Factors coming from those tools must be declared with \code{factorType = "multiplication"} and are inverted on the way in. The factors are centred so that their mean is one, which keeps the normalised values on the scale of the raw counts instead of collapsing them to fractions.
#'
#' The choice of the method matters more than usual on region sets. \code{"TMM"} and \code{"RLE"} assume that most of the regions do not change, which is reasonable for a catalogue of thousands of peaks but not for a handful of hand-picked ones, and not for a mark that is globally redistributed by the treatment. \code{"background"} sidesteps that assumption by estimating the factors on the genome wide bins, where the signal of the experiment is diluted, and is the safest option when a global shift is expected. \code{"spikeIn"} relies on the exogenous genome alone and ignores the regions altogether. \code{"loess"} corrects a bias that changes with the abundance, which no single factor per sample can describe, so it returns a matrix of offsets rather than a vector and leaves \code{scaling.factor} empty.
#'
#' @examples
#' \dontrun{
#' counts <- normalizeCounts(counts, method = "TMM")
#'
#' counts <- normalizeCounts(counts, method = "background")
#'
#' counts <- normalizeCounts(counts,
#'                           method = "spikeIn",
#'                           spikeInCounts = c(dmso_rep1 = 412553, dmso_rep2 = 388120,
#'                                             combo_rep1 = 501233, combo_rep2 = 470981))
#'
#' # Factors written by deeptools multiply the signal, so they are declared as such
#' counts <- normalizeCounts(counts,
#'                           method = "manual",
#'                           scalingFactors = c(0.81, 0.94, 1.12, 1.05),
#'                           factorType = "multiplication")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countReads}}, \code{\link{countBackground}}
#'
#' @importFrom edgeR calcNormFactors
#' @importFrom csaw normFactors normOffsets
#' @importFrom SummarizedExperiment assay assay<- assayNames colData colData<- rowData
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export normalizeCounts

normalizeCounts <-
  function(counts,
           method = "TMM",
           scalingFactors = NULL,
           factorType = "division",
           spikeInCounts = NULL,
           useRegionSets = NULL,
           minCount = 1,
           referenceSample = NULL,
           normalizedAssay = "norm.counts",
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    method <- method[1]
    if (!(method %in% c("TMM", "TMMwsp", "RLE", "upperQuartile", "librarySize", "background", "loess", "spikeIn", "manual", "none"))) {
      stop("The 'method' parameter must be one among 'TMM', 'TMMwsp', 'RLE', 'upperQuartile', 'librarySize', 'background', 'loess', 'spikeIn', 'manual' or 'none'.", call. = FALSE)
    }

    factorType <- tolower(factorType[1])
    if (!(factorType %in% c("division", "multiplication"))) {
      stop("The 'factorType' parameter must be either 'division' or 'multiplication'.", call. = FALSE)
    }

    # Silently ignoring the factors would leave the user convinced that they have been applied
    if (!is.null(scalingFactors) & method != "manual") {
      stop("The 'scalingFactors' parameter is used only by the 'manual' method, set 'method' accordingly.", call. = FALSE)
    }

    if (is.null(scalingFactors) & method == "manual") {
      stop("The 'manual' method requires the 'scalingFactors' parameter.", call. = FALSE)
    }

    if (is.null(spikeInCounts) & method == "spikeIn") {
      stop("The 'spikeIn' method requires the 'spikeInCounts' parameter.", call. = FALSE)
    }

    countMatrix <- SummarizedExperiment::assay(counts, "counts")
    sampleNumber <- ncol(counts)

    #-------------------------------------#
    # Sequencing depth of each sample     #
    #-------------------------------------#
    librarySizes <- SummarizedExperiment::colData(counts)$library.size

    # bigWig and external matrices come without a depth, the signal collected in the regions is the closest stand-in
    if (is.null(librarySizes) | any(is.na(librarySizes))) {
      librarySizes <- as.numeric(colSums(countMatrix))
      if (isTRUE(verbose) & method %in% c("TMM", "TMMwsp", "RLE", "upperQuartile", "librarySize")) {
        message("The library sizes are unknown, the column sums of the counts are used in their place.")
      }
    }

    #-----------------------------------#
    # Regions used for the estimation   #
    #-----------------------------------#
    rowTable <- data.frame(row.index = seq_len(nrow(counts)),
                           region.set = as.character(SummarizedExperiment::rowData(counts)$region.set),
                           total.count = as.numeric(rowSums(countMatrix)),
                           stringsAsFactors = FALSE)

    estimationRows <- dplyr::filter(rowTable, .data$total.count >= minCount)

    if (!is.null(useRegionSets)) {
      absentSets <- setdiff(useRegionSets, unique(rowTable$region.set))
      if (length(absentSets) > 0) {
        stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }
      estimationRows <- dplyr::filter(estimationRows, .data$region.set %in% useRegionSets)
    }

    if (nrow(estimationRows) == 0 & method %in% c("TMM", "TMMwsp", "RLE", "upperQuartile", "loess")) {
      stop("No region passes 'minCount', the factors cannot be estimated.", call. = FALSE)
    }

    # A handful of regions gives a factor that reflects those regions rather than the experiment
    if (nrow(estimationRows) < 100 & method %in% c("TMM", "TMMwsp", "RLE", "upperQuartile", "loess")) {
      warning(paste0("Only ", nrow(estimationRows), " regions are used to estimate the factors, consider 'background' or 'spikeIn' normalisation."), call. = FALSE)
    }

    #-------------------------------#
    # Estimate the scaling factors  #
    #-------------------------------#
    normFactorVector <- rep(NA_real_, sampleNumber)
    offsetMatrix <- NULL

    if (method == "none") {
      normFactorVector <- rep(1, sampleNumber)
      scalingFactorVector <- rep(1, sampleNumber)

    } else if (method == "manual") {
      scalingFactors <- as.numeric(.orderSampleValues(values = scalingFactors, sampleNames = colnames(counts), parameterName = "scalingFactors"))

      # Everything downstream divides, the multiplicative factors are inverted here and only here
      scalingFactorVector <- if (factorType == "division") {scalingFactors} else {1 / scalingFactors}

    } else if (method == "spikeIn") {
      spikeInCounts <- as.numeric(.orderSampleValues(values = spikeInCounts, sampleNames = colnames(counts), parameterName = "spikeInCounts"))

      # More exogenous reads means more material sequenced, the endogenous signal is scaled down accordingly
      scalingFactorVector <- spikeInCounts / mean(spikeInCounts)

    } else if (method == "librarySize") {
      scalingFactorVector <- librarySizes / mean(librarySizes)

    } else if (method %in% c("TMM", "TMMwsp", "RLE", "upperQuartile")) {
      edgeRMethod <- c("TMM" = "TMM", "TMMwsp" = "TMMwsp", "RLE" = "RLE", "upperQuartile" = "upperquartile")[method]

      normFactorVector <- edgeR::calcNormFactors(object = countMatrix[estimationRows$row.index, , drop = FALSE],
                                                 lib.size = librarySizes,
                                                 method = as.character(edgeRMethod),
                                                 refColumn = .resolveSampleIndex(sample = referenceSample, sampleNames = colnames(counts)))

      scalingFactorVector <- (librarySizes * normFactorVector) / mean(librarySizes * normFactorVector)

    } else if (method == "background") {
      backgroundBins <- S4Vectors::metadata(counts)$background

      if (is.null(backgroundBins)) {
        stop("No background bin is stored in the object, run 'countBackground' before normalising with this method.", call. = FALSE)
      }

      # The bins carry mostly noise, the factors they give describe the depth rather than the biology of the regions
      normFactorVector <- as.numeric(csaw::normFactors(object = backgroundBins, se.out = FALSE))
      scalingFactorVector <- (librarySizes * normFactorVector) / mean(librarySizes * normFactorVector)

    } else {
      # A bias that changes along the abundance needs one value per region and per sample, not one per sample
      offsetMatrix <- csaw::normOffsets(object = countMatrix, se.out = FALSE)
      scalingFactorVector <- rep(NA_real_, sampleNumber)
    }

    if (method != "loess") {
      if (any(!is.finite(scalingFactorVector)) | any(scalingFactorVector <= 0)) {
        stop("The scaling factors must be finite and strictly positive.", call. = FALSE)
      }
    }

    #-------------------------------#
    # Store factors and values      #
    #-------------------------------#
    if (is.null(offsetMatrix)) {
      normalizedMatrix <- sweep(countMatrix, MARGIN = 2, STATS = scalingFactorVector, FUN = "/")
    } else {
      # Centring the offsets on the row keeps the normalised values around the magnitude of the raw ones
      normalizedMatrix <- countMatrix / exp(offsetMatrix - rowMeans(offsetMatrix))
      SummarizedExperiment::assay(counts, "offset") <- offsetMatrix
    }

    SummarizedExperiment::colData(counts)$norm.factor <- normFactorVector
    SummarizedExperiment::colData(counts)$scaling.factor <- scalingFactorVector
    SummarizedExperiment::assay(counts, normalizedAssay) <- normalizedMatrix

    S4Vectors::metadata(counts)$normalization <- list(method = method,
                                                      factor.type = "division",
                                                      normalized.assay = normalizedAssay,
                                                      offsets = !is.null(offsetMatrix))

    counts@parameters <- c(counts@parameters,
                           list(normalizeCounts = list(method = method,
                                                       factorType = factorType,
                                                       useRegionSets = useRegionSets,
                                                       minCount = minCount,
                                                       referenceSample = referenceSample,
                                                       n.regions.used = nrow(estimationRows))))

    if (isTRUE(verbose)) {
      if (method == "loess") {
        message("Done. A matrix of log offsets has been stored in the 'offset' assay, no single factor per sample applies.")
      } else {
        message(paste0("Done. The counts have been divided by the scaling factors (",
                       paste(round(range(scalingFactorVector), 3), collapse = " - "), ")."))
      }
    }

    return(counts)
  } # END function




#' @title .orderSampleValues
#'
#' @description Checks a vector of user supplied values against the samples of a counts object and puts it in the order of the columns.
#'
#' @param values Numeric vector, named or not.
#' @param sampleNames Character vector with the sample names, in the order of the columns.
#' @param parameterName String with the name of the parameter, used in the error messages.
#'
#' @return The numeric vector, ordered as the columns of the object.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.orderSampleValues <-
  function(values,
           sampleNames,
           parameterName) {

    if (!is.numeric(values)) {
      stop(paste0("The '", parameterName, "' parameter must be numeric."), call. = FALSE)
    }

    # Names are the only way to be sure of the pairing, an unnamed vector is trusted to follow the columns
    if (!is.null(names(values))) {
      if (any(duplicated(names(values)))) {
        stop(paste0("The '", parameterName, "' parameter contains duplicated sample names."), call. = FALSE)
      }

      absentSamples <- setdiff(sampleNames, names(values))
      if (length(absentSamples) > 0) {
        stop(paste0("The '", parameterName, "' parameter is missing the following samples: ", paste(absentSamples, collapse = ", "), "."), call. = FALSE)
      }

      # Values for samples absent from the object are dropped, so one table can serve several objects
      values <- values[sampleNames]
    } else {
      if (length(values) != length(sampleNames)) {
        stop(paste0("The '", parameterName, "' parameter must have one value per sample, or be named after the samples."), call. = FALSE)
      }
    }

    if (any(is.na(values))) {
      stop(paste0("The '", parameterName, "' parameter contains missing values."), call. = FALSE)
    }

    return(values)
  } # END function




#' @title .resolveSampleIndex
#'
#' @description Turns a sample name or position into a column index.
#'
#' @param sample String or numeric position identifying a sample. Default: \code{NULL}.
#' @param sampleNames Character vector with the sample names, in the order of the columns.
#'
#' @return The column index, or \code{NULL} when no sample is given.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resolveSampleIndex <-
  function(sample = NULL,
           sampleNames) {

    if (is.null(sample)) {
      return(NULL)
    }

    if (is.character(sample)) {
      if (!(sample[1] %in% sampleNames)) {
        stop(paste0("The sample '", sample[1], "' is absent from the object."), call. = FALSE)
      }
      return(match(sample[1], sampleNames))
    }

    sampleIndex <- as.integer(sample[1])
    if (is.na(sampleIndex) | sampleIndex < 1 | sampleIndex > length(sampleNames)) {
      stop("The reference sample is outside the range of the samples.", call. = FALSE)
    }

    return(sampleIndex)
  } # END function
