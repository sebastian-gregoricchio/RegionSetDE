# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title makeGreylist
#'
#' @description Builds a greylist from input libraries: the stretches of genome where an input carries far more fragments than the rest of its genome leads to expect, such as copy number gains of the cell line, collapsed repeats or regions that stick to any immunoprecipitation. Each input is judged against its own coverage, and the regions flagged are pooled over the inputs, ready for \code{\link{applyGreylist}}.
#'
#' @param inputFiles Character vector with the paths of the input BAM files, or the data.frame returned by \code{\link{loadSampleSheet}}, in which case every distinct file of its \code{input} column is used once and named after \code{input.id}.
#' @param inputNames Character vector with the names of the inputs. Default: \code{NULL}, the \code{input.id} of the sample sheet or the file names.
#' @param binSize Numeric value with the width of the windows, in base pairs. Consecutive windows overlap by half of it. Default: \code{1024}.
#' @param quantile Numeric value with the quantile of the fitted negative binomial above which a window is flagged. Default: \code{0.99}.
#' @param maxGap Numeric value, in base pairs: two flagged windows separated by a shorter gap are merged into one region. Default: \code{16384}.
#' @param minInputs Numeric value with the number of inputs that must flag a position for it to enter the greylist. Default: \code{1}, a position flagged by any input.
#' @param excludeChromosomes Character vector with the chromosomes left out, written in either naming style, as in \code{\link{countReads}}. Default: \code{NULL}, none.
#' @param pairedEnd Logical value, one logical value per BAM file, or the string \code{"auto"} to read the layout from the files themselves. Default: \code{"auto"}.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended. Default: \code{150}.
#' @param maxFragmentLength Numeric value with the maximum length accepted for a paired-end fragment. Default: \code{1000}.
#' @param minMapq Numeric value with the minimum mapping quality of a read. Default: \code{20}.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded. Default: \code{TRUE}.
#' @param nThreads Number of threads. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{GRanges} with one element per greylisted region, carrying \code{n.inputs}, the number of inputs flagging it, and \code{inputs}, their names. Its \code{metadata} holds \code{thresholds}, a data.frame with, for every input, the fragments counted, the mean and size of the fitted negative binomial, the threshold and how much of the genome was flagged, and \code{parameters}, the arguments of the call.
#'
#' @details The procedure follows GreyListChIP, which is also what DiffBind uses for its greylists. The genome is cut into windows of \code{binSize} base pairs overlapping by half, the fragments of each input are counted once in each window holding their centre, and a negative binomial is fitted to the counts of every input. Windows above the \code{quantile} of that distribution are flagged, and flagged windows separated by less than \code{maxGap} are merged into one region, the stretch between them included.
#'
#' Two things differ from GreyListChIP. The distribution is fitted by maximum likelihood on all the windows, rather than on 100 bootstrap samples of 30,000 windows, so the threshold does not depend on the random seed. And the fragments are read with the same filters as in \code{\link{countReads}}: mapping quality, duplicates and proper pairs.
#'
#' Empty windows enter the fit, as they do in GreyListChIP. A chromosome without any read, such as chrY in a female sample or an unused contig, adds a block of zeros that widens the fitted distribution and raises the threshold, so it is worth leaving out through \code{excludeChromosomes}.
#'
#' A greylist describes the input, not the chromatin. Marks like H3K9me3 sit largely on satellites and other repeats, which is precisely where inputs pile up, so a greylist can take genuine signal away with the artefacts. Check how much of each set \code{\link{applyGreylist}} removes, and consider \code{trimRegions = TRUE}, which cuts the greylisted stretch out of a broad domain rather than dropping the whole domain.
#'
#' @references Brown G. GreyListChIP: Grey Lists -- Mask Artefact Regions Based on ChIP Inputs. Bioconductor package.
#'
#' @examples
#' # The alignment shipped with Rsamtools stands in for an input library
#' inputFile <- system.file("extdata", "ex1.bam", package = "Rsamtools")
#'
#' greylist <- makeGreylist(inputFile, binSize = 200, maxGap = 200)
#' greylist
#'
#' # The threshold and the fit behind it
#' S4Vectors::metadata(greylist)$thresholds
#'
#' # From a sample sheet, each distinct input once
#' sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#' greylist <- makeGreylist(sampleSheet, binSize = 10000, verbose = FALSE)
#'
#' peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)
#' peakRegions <- applyGreylist(peakRegions, greylist = greylist, verbose = FALSE)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{applyGreylist}}, \code{\link{loadSampleSheet}}, \code{\link{applyBlacklist}}
#'
#' @importFrom Rsamtools scanBamHeader testPairedEndBam
#' @importFrom IRanges reduce
#' @importFrom S4Vectors metadata<-
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom BiocGenerics width
#' @importFrom stats qnbinom qpois
#' @importFrom dplyr distinct filter
#' @importFrom rlang .data
#'
#' @export makeGreylist

makeGreylist <-
  function(inputFiles,
           inputNames = NULL,
           binSize = 1024,
           quantile = 0.99,
           maxGap = 16384,
           minInputs = 1,
           excludeChromosomes = NULL,
           pairedEnd = "auto",
           fragmentLength = 150,
           maxFragmentLength = 1000,
           minMapq = 20,
           removeDuplicates = TRUE,
           nThreads = 1,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    # A sample sheet brings its inputs, each distinct file once whatever the number of samples pointing to it
    if (is.data.frame(inputFiles)) {
      if (!("input" %in% colnames(inputFiles))) {
        stop("The sample sheet has no 'input' column.", call. = FALSE)
      }

      inputTable <- dplyr::distinct(dplyr::filter(inputFiles, !is.na(.data$input)), .data$input, .keep_all = TRUE)

      if (nrow(inputTable) == 0) {
        stop("No sample of the sheet has an input.", call. = FALSE)
      }

      if (is.null(inputNames) & "input.id" %in% colnames(inputTable)) {inputNames <- inputTable$input.id}
      inputFiles <- inputTable$input
    }

    if (!is.character(inputFiles) | length(inputFiles) == 0) {
      stop("The 'inputFiles' parameter must be a character vector with at least one BAM file, or a sample sheet.", call. = FALSE)
    }

    if (anyDuplicated(inputFiles) > 0) {
      stop("The following input files are listed more than once: ",
           paste(unique(basename(inputFiles[duplicated(inputFiles)])), collapse = ", "), ".", call. = FALSE)
    }

    if (is.null(inputNames)) {
      inputNames <- make.unique(sub("\\.bam$", "", basename(inputFiles), ignore.case = TRUE), sep = "_")
    }

    if (length(inputNames) != length(inputFiles) | anyDuplicated(inputNames) > 0 | anyNA(inputNames)) {
      stop("The 'inputNames' parameter must hold one distinct name per input file.", call. = FALSE)
    }

    binSize <- as.integer(binSize[1])
    if (is.na(binSize) | binSize < 2) {
      stop("The 'binSize' parameter must be an integer of at least 2.", call. = FALSE)
    }

    if (!is.numeric(quantile) | length(quantile) != 1 | anyNA(quantile) | isTRUE(quantile <= 0) | isTRUE(quantile >= 1)) {
      stop("The 'quantile' parameter must be a single value between 0 and 1.", call. = FALSE)
    }

    maxGap <- as.integer(maxGap[1])
    if (is.na(maxGap) | maxGap < 0) {
      stop("The 'maxGap' parameter must be a non-negative integer.", call. = FALSE)
    }

    minInputs <- as.integer(minInputs[1])
    if (is.na(minInputs) | minInputs < 1 | minInputs > length(inputFiles)) {
      stop("The 'minInputs' parameter must lie between 1 and the number of inputs (", length(inputFiles), ").", call. = FALSE)
    }

    if (!is.null(excludeChromosomes) & !is.character(excludeChromosomes)) {
      stop("The 'excludeChromosomes' parameter must be a character vector with chromosome names.", call. = FALSE)
    }

    missingFiles <- inputFiles[!file.exists(inputFiles)]
    if (length(missingFiles) > 0) {
      stop("The following input files do not exist: ", paste(missingFiles, collapse = ", "), ".", call. = FALSE)
    }

    if (any(!.hasBamIndex(inputFiles))) {
      stop("The following input files are not indexed: ", paste(basename(inputFiles[!.hasBamIndex(inputFiles)]), collapse = ", "), ".", call. = FALSE)
    }

    if (identical(pairedEnd, "auto")) {
      pairedEnd <- vapply(inputFiles, Rsamtools::testPairedEndBam, logical(1), USE.NAMES = FALSE)
    }

    if (length(pairedEnd) == 1) {pairedEnd <- rep(pairedEnd, length(inputFiles))}

    if (!is.logical(pairedEnd) | length(pairedEnd) != length(inputFiles) | anyNA(pairedEnd)) {
      stop("The 'pairedEnd' parameter must be 'auto', a single logical value, or one logical value per input file.", call. = FALSE)
    }

    #-------------------------------#
    # Windows over the genome       #
    #-------------------------------#
    chromosomeLengths <- Rsamtools::scanBamHeader(inputFiles[1])[[1]]$targets

    # A name in the other naming style would leave the chromosome in and go unnoticed
    excludeChromosomes <- .matchChromosomeNames(chromosomeNames = excludeChromosomes, targetSeqlevels = names(chromosomeLengths))

    chromosomeLengths <- chromosomeLengths[!(names(chromosomeLengths) %in% excludeChromosomes)]

    if (length(chromosomeLengths) == 0) {
      stop("Every chromosome of the BAM files is listed in 'excludeChromosomes', no window is left to count.", call. = FALSE)
    }

    windowRanges <- .greylistWindows(chromosomeLengths = chromosomeLengths, binSize = binSize)

    if (isTRUE(verbose)) {
      message("Counting ", length(inputFiles), " input libraries over ", format(length(windowRanges), big.mark = ",", trim = TRUE),
              " windows of ", binSize, " bp...")
    }

    # Each fragment counts once in every window holding its centre, two windows since they overlap by half
    countList <- .countBamFragments(bamFiles = inputFiles,
                                    ranges = windowRanges,
                                    pairedEnd = pairedEnd,
                                    fragmentLength = fragmentLength,
                                    maxFragmentLength = maxFragmentLength,
                                    minMapq = minMapq,
                                    removeDuplicates = removeDuplicates,
                                    excludeChromosomes = excludeChromosomes,
                                    fullLibrarySize = FALSE,
                                    countMode = "bin",
                                    nThreads = nThreads)

    #-------------------------------#
    # Threshold of every input      #
    #-------------------------------#
    thresholdList <- list()
    flaggedList <- list()

    for (i in seq_along(inputFiles)) {
      windowCounts <- countList$counts[, i]

      if (sum(windowCounts) == 0) {
        stop("The input '", inputNames[i], "' has no fragment in the windows, check 'excludeChromosomes' and the read filters.", call. = FALSE)
      }

      fittedDistribution <- .fitNegativeBinomial(windowCounts = windowCounts)

      # Without overdispersion the negative binomial collapses on the Poisson
      countThreshold <- if (is.finite(fittedDistribution$size)) {
        stats::qnbinom(quantile, size = fittedDistribution$size, mu = fittedDistribution$mu)
      } else {
        stats::qpois(quantile, lambda = fittedDistribution$mu)
      }

      flaggedWindows <- windowRanges[windowCounts > countThreshold]
      flaggedRegions <- IRanges::reduce(flaggedWindows, min.gapwidth = max(maxGap, 1L))

      flaggedList[[inputNames[i]]] <- flaggedRegions
      thresholdList[[i]] <- data.frame(input = inputNames[i],
                                       file = inputFiles[i],
                                       fragments = countList$library.size[i],
                                       mean = fittedDistribution$mu,
                                       size = fittedDistribution$size,
                                       threshold = countThreshold,
                                       windows = length(windowRanges),
                                       flagged.windows = length(flaggedWindows),
                                       regions = length(flaggedRegions),
                                       greylisted.bp = sum(as.numeric(BiocGenerics::width(flaggedRegions))),
                                       stringsAsFactors = FALSE)

      if (isTRUE(verbose)) {
        message("  ", inputNames[i], ": more than ", countThreshold, " fragments per window (mean ", round(fittedDistribution$mu, 1), "), ",
                format(length(flaggedWindows), big.mark = ",", trim = TRUE), " windows flagged, ",
                format(length(flaggedRegions), big.mark = ",", trim = TRUE), " regions, ",
                round(thresholdList[[i]]$greylisted.bp / 1e6, 2), " Mb.")
      }
    }

    thresholdTable <- do.call(what = rbind, args = thresholdList)

    #-------------------------------#
    # Pool the inputs               #
    #-------------------------------#
    greylist <- .poolGreylists(flaggedList = flaggedList, minInputs = minInputs, genomeInfo = GenomeInfoDb::seqinfo(windowRanges))

    S4Vectors::metadata(greylist) <- list(thresholds = thresholdTable,
                                          parameters = list(inputFiles = inputFiles,
                                                            inputNames = inputNames,
                                                            binSize = binSize,
                                                            quantile = quantile,
                                                            maxGap = maxGap,
                                                            minInputs = minInputs,
                                                            excludeChromosomes = excludeChromosomes,
                                                            pairedEnd = pairedEnd,
                                                            fragmentLength = fragmentLength,
                                                            maxFragmentLength = maxFragmentLength,
                                                            minMapq = minMapq,
                                                            removeDuplicates = removeDuplicates))

    if (isTRUE(verbose)) {
      genomeLength <- sum(as.numeric(chromosomeLengths))
      greylistLength <- sum(as.numeric(BiocGenerics::width(greylist)))

      message("Greylist: ", format(length(greylist), big.mark = ",", trim = TRUE), " regions, ", round(greylistLength / 1e6, 2), " Mb (",
              signif(100 * greylistLength / genomeLength, 2), "% of the genome), flagged by at least ", minInputs, " of ", length(inputFiles), " inputs.")
    }

    return(greylist)
  } # END function




#' @title .greylistWindows
#'
#' @description Cuts the chromosomes into windows of fixed width overlapping by half, the last ones trimmed at the end of each chromosome.
#'
#' @param chromosomeLengths Named numeric vector with the length of every chromosome.
#' @param binSize Integer with the width of the windows.
#'
#' @return A \code{GRanges} with the windows and the chromosome lengths in its \code{seqinfo}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @keywords internal

.greylistWindows <-
  function(chromosomeLengths,
           binSize) {

    stepSize <- binSize %/% 2L

    windowsPerChromosome <- as.integer(ceiling(chromosomeLengths / stepSize))
    windowStarts <- as.integer(unlist(lapply(windowsPerChromosome, function(n) {(seq_len(n) - 1L) * stepSize + 1L}), use.names = FALSE))
    windowEnds <- pmin(windowStarts + binSize - 1L, rep(as.integer(chromosomeLengths), times = windowsPerChromosome))

    return(GenomicRanges::GRanges(seqnames = rep(names(chromosomeLengths), times = windowsPerChromosome),
                                  ranges = IRanges::IRanges(start = windowStarts, end = windowEnds),
                                  seqlengths = chromosomeLengths))
  } # END function




#' @title .fitNegativeBinomial
#'
#' @description Fits a negative binomial to window counts by maximum likelihood. The mean is estimated by the sample mean, which is its maximum likelihood estimate whatever the size, and the size by maximising the likelihood over its logarithm, computed once per distinct count value.
#'
#' @param windowCounts Integer vector with the count of every window.
#'
#' @return A list with \code{mu} and \code{size}, the latter infinite when the counts are not overdispersed.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats var dnbinom optimize
#'
#' @keywords internal

.fitNegativeBinomial <-
  function(windowCounts) {

    meanCount <- mean(windowCounts)
    varianceCount <- stats::var(windowCounts)

    if (!is.finite(varianceCount) | varianceCount <= meanCount) {
      return(list(mu = meanCount, size = Inf))
    }

    # Millions of windows share a few thousand distinct values, so the likelihood is summed over those
    countFrequencies <- tabulate(as.integer(windowCounts) + 1L)
    observedValues <- which(countFrequencies > 0) - 1L
    countFrequencies <- countFrequencies[countFrequencies > 0]

    logLikelihood <- function(logSize) {
      sum(countFrequencies * stats::dnbinom(observedValues, size = exp(logSize), mu = meanCount, log = TRUE))
    }

    # The moments give the scale of the size, the likelihood is searched well beyond it on both sides
    momentSize <- meanCount^2 / (varianceCount - meanCount)
    fittedSize <- stats::optimize(f = logLikelihood, interval = log(momentSize) + c(-10, 10), maximum = TRUE)

    return(list(mu = meanCount, size = exp(fittedSize$maximum)))
  } # END function




#' @title .poolGreylists
#'
#' @description Pools the regions flagged in every input, keeping the positions flagged by at least \code{minInputs} of them, and records which inputs flagged each region.
#'
#' @param flaggedList Named list of \code{GRanges}, one per input, each without overlaps within itself.
#' @param minInputs Integer with the number of inputs that must flag a position.
#' @param genomeInfo \code{Seqinfo} of the windows.
#'
#' @return A \code{GRanges} with the \code{n.inputs} and \code{inputs} metadata columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRanges GRangesList findOverlaps
#' @importFrom IRanges coverage
#' @importFrom S4Vectors mcols<- DataFrame queryHits subjectHits
#' @importFrom dplyr group_by summarise n_distinct
#' @importFrom rlang .data
#'
#' @keywords internal

.poolGreylists <-
  function(flaggedList,
           minInputs,
           genomeInfo) {

    stackedRanges <- unlist(GenomicRanges::GRangesList(flaggedList), use.names = FALSE)
    stackedRanges$input <- rep(names(flaggedList), times = lengths(flaggedList))

    # Each input is free of overlaps within itself, so the coverage counts how many inputs flag a position
    inputCoverage <- IRanges::coverage(stackedRanges)
    pooledList <- IRanges::slice(inputCoverage, lower = minInputs, rangesOnly = TRUE)

    greylist <- GenomicRanges::GRanges(seqnames = rep(names(pooledList), times = lengths(pooledList)),
                                       ranges = unlist(pooledList, use.names = FALSE),
                                       seqinfo = genomeInfo)

    if (length(greylist) == 0) {
      S4Vectors::mcols(greylist) <- S4Vectors::DataFrame(n.inputs = integer(0), inputs = character(0))
      return(greylist)
    }

    #-------------------------------#
    # Inputs behind every region    #
    #-------------------------------#
    overlapHits <- GenomicRanges::findOverlaps(greylist, stackedRanges)
    hitTable <- data.frame(region = S4Vectors::queryHits(overlapHits),
                           input = stackedRanges$input[S4Vectors::subjectHits(overlapHits)],
                           stringsAsFactors = FALSE)

    inputTable <- dplyr::summarise(dplyr::group_by(hitTable, .data$region),
                                   n.inputs = dplyr::n_distinct(.data$input),
                                   inputs = paste(sort(unique(.data$input)), collapse = ","),
                                   .groups = "drop")

    greylist$n.inputs <- inputTable$n.inputs[match(seq_along(greylist), inputTable$region)]
    greylist$inputs <- inputTable$inputs[match(seq_along(greylist), inputTable$region)]

    return(greylist)
  } # END function
