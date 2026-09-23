# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title countGreenlist
#'
#' @description Counts the reads falling in the regions of a CUT&RUN or CUT&Tag greenlist, the places where the background of the protocol is reproducible enough to measure how much material was sequenced. The counts are stored in the metadata of the object, where \code{\link{normalizeCounts}} picks them up with \code{method = "greenlist"}.
#'
#' @param counts \code{RegionSetDE.counts} object returned by \code{\link{countReads}}.
#' @param greenlist \code{GRanges} with the greenlist regions, typically from \code{\link{loadGreenlist}}, or the path to a BED file holding them.
#' @param bamFiles Character vector with the paths of the BAM files, in the same order as the samples of \code{counts}. Default: \code{NULL}, the files recorded by \code{\link{countReads}} are reused.
#' @param minCount Numeric value with the minimum total count required to keep a greenlist region. Default: \code{1}.
#' @param pairedEnd Logical value, or one logical value per BAM file, indicating whether the reads must be counted as proper pairs. Default: \code{NULL}, the layouts resolved at the counting step.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended. Default: \code{NULL}, the value used at the counting step.
#' @param maxFragmentLength Numeric value with the maximum insert size accepted for a pair. Default: \code{NULL}, the value used at the counting step.
#' @param minMapq Numeric value with the minimum mapping quality of a read. Default: \code{NULL}, the value used at the counting step.
#' @param removeDuplicates Logical value indicating whether the duplicated reads must be discarded. Default: \code{NULL}, the value used at the counting step.
#' @param nThreads Number of threads. The files are cut into pieces of at most 50 Mb, shared among the threads. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return The input \code{RegionSetDE.counts} object with the greenlist counts stored as a \code{RangedSummarizedExperiment} in \code{metadata(counts)$greenlist}.
#'
#' @details The read filters are taken from \code{\link{countReads}} unless they are given here, for the same reason as in \code{\link{countBackground}}: a reference counted with another mapping quality or duplicate policy describes a library that is not the one under study.
#'
#' Each fragment is counted once, in the region holding its centre, as the background bins do, so the totals stay a share of the library and two neighbouring regions never claim the same fragment. The list is merged beforehand for the same reason. Greenlist regions lying on chromosomes absent from the BAM files are dropped, and how many were is reported, which is what catches a list built for another assembly before it quietly halves the counts.
#'
#' @examples
#' sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#' peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)
#' counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)
#'
#' # A published list, loadGreenlist("hg38", assay = "cuttag"), holds too few reads in these small
#' # slices of chromosome 19, so a few stretches of the window stand in for it here
#' greenlist <- GenomicRanges::GRanges("19", IRanges::IRanges(start = seq(46.5e6, 57.5e6, by = 1e6), width = 2e5))
#'
#' counts <- countGreenlist(counts, greenlist = greenlist, verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "greenlist", verbose = FALSE)
#'
#' SummarizedExperiment::colData(counts)$scaling.factor
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadGreenlist}}, \code{\link{normalizeCounts}}, \code{\link{countBackground}}
#'
#' @importFrom Rsamtools scanBamHeader
#' @importFrom SummarizedExperiment SummarizedExperiment assay
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom BiocGenerics width
#' @importFrom IRanges reduce
#' @importFrom S4Vectors metadata metadata<- DataFrame
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export countGreenlist

countGreenlist <-
  function(counts,
           greenlist,
           bamFiles = NULL,
           minCount = 1,
           pairedEnd = NULL,
           fragmentLength = NULL,
           maxFragmentLength = NULL,
           minMapq = NULL,
           removeDuplicates = NULL,
           nThreads = 1,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    if (is.character(greenlist) & length(greenlist) == 1) {
      greenlist <- .readRegionFile(filePath = greenlist, header = FALSE, asGRanges = TRUE)
    }

    if (!methods::is(greenlist, "GRanges")) {
      stop("The 'greenlist' parameter must be a GRanges or the path to a BED file, see loadGreenlist().", call. = FALSE)
    }

    #-----------------------------------#
    # Recover the counting parameters   #
    #-----------------------------------#
    # A reference counted under a different read filter does not describe the libraries being normalised
    countingParameters <- counts@parameters$countReads

    if (is.null(countingParameters) & is.null(bamFiles)) {
      stop("No BAM file is recorded in the object, provide them through the 'bamFiles' parameter.", call. = FALSE)
    }

    if (is.null(bamFiles)) {bamFiles <- countingParameters$bamFiles}
    if (is.null(pairedEnd)) {pairedEnd <- if (is.null(countingParameters$pairedEnd)) {FALSE} else {countingParameters$pairedEnd}}
    if (is.null(fragmentLength)) {fragmentLength <- if (is.null(countingParameters$fragmentLength)) {150} else {countingParameters$fragmentLength}}
    if (is.null(maxFragmentLength)) {maxFragmentLength <- if (is.null(countingParameters$maxFragmentLength)) {1000} else {countingParameters$maxFragmentLength}}
    if (is.null(minMapq)) {minMapq <- if (is.null(countingParameters$minMapq)) {20} else {countingParameters$minMapq}}
    if (is.null(removeDuplicates)) {removeDuplicates <- if (is.null(countingParameters$removeDuplicates)) {TRUE} else {countingParameters$removeDuplicates}}

    if (length(bamFiles) != ncol(counts)) {
      stop("The number of BAM files does not match the number of samples of the counts object.", call. = FALSE)
    }

    if (length(pairedEnd) == 1) {pairedEnd <- rep(pairedEnd, length(bamFiles))}

    if (!is.logical(pairedEnd) | length(pairedEnd) != length(bamFiles) | any(is.na(pairedEnd))) {
      stop("The 'pairedEnd' parameter must be a single logical value, or one logical value per BAM file.", call. = FALSE)
    }

    missingFiles <- bamFiles[!file.exists(bamFiles)]
    if (length(missingFiles) > 0) {
      stop("The following BAM files do not exist: ", paste(missingFiles, collapse = ", "), ".", call. = FALSE)
    }

    #-----------------------------------#
    # Bring the list onto the BAM files #
    #-----------------------------------#
    # Overlapping entries would otherwise count the same fragment twice
    greenlistRanges <- IRanges::reduce(greenlist, ignore.strand = TRUE)

    bamSeqlevels <- names(Rsamtools::scanBamHeader(bamFiles[1])[[1]]$targets)
    greenlistRanges <- .matchSeqlevels(x = greenlistRanges, targetSeqlevels = bamSeqlevels, fileName = bamFiles[1], verbose = verbose)

    absentRegions <- !(as.character(GenomeInfoDb::seqnames(greenlistRanges)) %in% bamSeqlevels)

    if (all(absentRegions)) {
      stop("No greenlist region lies on a chromosome of the BAM files. The list and the alignments are on different assemblies.", call. = FALSE)
    }

    if (any(absentRegions) & isTRUE(verbose)) {
      message(sum(absentRegions), " greenlist regions out of ", length(greenlistRanges),
              " lie on chromosomes absent from the BAM files and are dropped.")
    }

    greenlistRanges <- greenlistRanges[!absentRegions]

    #-------------------------------#
    # Count over the greenlist      #
    #-------------------------------#
    if (isTRUE(verbose)) {
      message("Counting reads over ", format(length(greenlistRanges), big.mark = ","), " greenlist regions across ",
              length(bamFiles), " samples...")
    }

    greenlistCounts <- .countBamFragments(bamFiles = bamFiles,
                                          ranges = greenlistRanges,
                                          pairedEnd = pairedEnd,
                                          fragmentLength = fragmentLength[1],
                                          maxFragmentLength = maxFragmentLength[1],
                                          minMapq = minMapq,
                                          removeDuplicates = removeDuplicates,
                                          excludeChromosomes = NULL,
                                          fullLibrarySize = FALSE,
                                          countMode = "bin",
                                          nThreads = nThreads)

    greenlistExperiment <-
      SummarizedExperiment::SummarizedExperiment(assays = list(counts = greenlistCounts$counts),
                                                 rowRanges = greenlistRanges,
                                                 colData = S4Vectors::DataFrame(bam.files = bamFiles,
                                                                                totals = as.numeric(colSums(greenlistCounts$counts)),
                                                                                row.names = colnames(counts)))

    colnames(greenlistExperiment) <- colnames(counts)

    #-------------------------------#
    # Filter and store              #
    #-------------------------------#
    regionTable <- data.frame(region.index = seq_len(nrow(greenlistExperiment)),
                              total.count = as.numeric(rowSums(SummarizedExperiment::assay(greenlistExperiment, "counts"))),
                              stringsAsFactors = FALSE)

    keptRegions <- dplyr::filter(regionTable, .data$total.count >= minCount)

    if (nrow(keptRegions) == 0) {
      stop("No greenlist region carries a read, the factors cannot be estimated from them.", call. = FALSE)
    }

    greenlistExperiment <- greenlistExperiment[keptRegions$region.index, ]

    S4Vectors::metadata(counts)$greenlist <- greenlistExperiment

    # Assigned rather than appended, so that counting twice leaves one record and not two
    counts@parameters$countGreenlist <- list(bamFiles = bamFiles,
                                             n.regions = nrow(greenlistExperiment),
                                             covered.bp = sum(as.numeric(BiocGenerics::width(greenlistRanges))),
                                             minCount = minCount,
                                             pairedEnd = pairedEnd,
                                             fragmentLength = fragmentLength,
                                             maxFragmentLength = maxFragmentLength,
                                             minMapq = minMapq,
                                             removeDuplicates = removeDuplicates,
                                             greenlist = S4Vectors::metadata(greenlist))

    if (isTRUE(verbose)) {
      message("Done. ", format(nrow(greenlistExperiment), big.mark = ","), " greenlist regions kept.")

      # The share of the library sitting on the greenlist is the number worth watching between samples
      librarySizes <- SummarizedExperiment::colData(counts)$library.size

      if (!is.null(librarySizes) && !any(is.na(librarySizes))) {
        greenlistFraction <- 100 * colSums(SummarizedExperiment::assay(greenlistExperiment, "counts")) / librarySizes

        message("They hold between ", format(round(min(greenlistFraction), 2), nsmall = 2), " and ",
                format(round(max(greenlistFraction), 2), nsmall = 2), "% of the libraries.")
      }
    }

    return(counts)
  } # END function




#' @title .greenlistMatrix
#'
#' @description Returns the matrix of greenlist counts the factors are estimated from, either the one stored by \code{countGreenlist} or the one handed to \code{normalizeCounts}.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param greenlistCounts Numeric matrix with one row per region and one column per sample, a numeric vector with one total per sample, or \code{NULL}.
#'
#' @return A numeric matrix with one column per sample.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#'
#' @keywords internal

.greenlistMatrix <-
  function(counts,
           greenlistCounts = NULL) {

    #-------------------------------#
    # Counted elsewhere             #
    #-------------------------------#
    if (!is.null(greenlistCounts)) {
      if (is.matrix(greenlistCounts) | is.data.frame(greenlistCounts)) {
        greenlistMatrix <- as.matrix(greenlistCounts)
      } else {
        # A single total per sample is a matrix of one row, and the estimators below work on it
        greenlistMatrix <- matrix(as.numeric(.orderSampleValues(values = greenlistCounts,
                                                                sampleNames = colnames(counts),
                                                                parameterName = "greenlistCounts")),
                                  nrow = 1)
        colnames(greenlistMatrix) <- colnames(counts)
      }

      if (ncol(greenlistMatrix) != ncol(counts)) {
        stop("The 'greenlistCounts' parameter must hold one column per sample of the object.", call. = FALSE)
      }

      # Named columns are matched to the samples, unnamed ones must already follow their order
      if (!is.null(colnames(greenlistMatrix)) && all(colnames(counts) %in% colnames(greenlistMatrix))) {
        greenlistMatrix <- greenlistMatrix[, colnames(counts), drop = FALSE]
      }

      return(greenlistMatrix)
    }

    #-------------------------------#
    # Counted by countGreenlist     #
    #-------------------------------#
    greenlistExperiment <- S4Vectors::metadata(counts)$greenlist

    if (is.null(greenlistExperiment)) {
      stop("No greenlist count is stored in the object, run 'countGreenlist' before normalising with this method, or pass the counts through 'greenlistCounts'.", call. = FALSE)
    }

    return(as.matrix(SummarizedExperiment::assay(greenlistExperiment, "counts")))
  } # END function




#' @title .greenlistFactors
#'
#' @description Turns a matrix of greenlist counts into one scaling factor per sample.
#'
#' @param greenlistMatrix Numeric matrix with one row per greenlist region and one column per sample.
#' @param estimator String with the estimator, one among \code{"medianRatio"}, \code{"TMM"} and \code{"sum"}.
#'
#' @return A numeric vector with one factor per sample, centred on one.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom edgeR normLibSizes
#' @importFrom stats median
#'
#' @keywords internal

.greenlistFactors <-
  function(greenlistMatrix,
           estimator = "medianRatio") {

    sampleTotals <- as.numeric(colSums(greenlistMatrix))

    if (any(sampleTotals == 0)) {
      stop("Some samples carry no read over the greenlist, their factors cannot be estimated from it.", call. = FALSE)
    }

    #-------------------------------#
    # The three estimators          #
    #-------------------------------#
    # One total per sample leaves no region to take ratios over, the total is the only factor it can give
    if (estimator == "sum" | nrow(greenlistMatrix) == 1) {
      return(sampleTotals / mean(sampleTotals))
    }

    if (estimator == "TMM") {
      normFactorVector <- edgeR::normLibSizes(object = greenlistMatrix, lib.size = sampleTotals, method = "TMM")
      return((sampleTotals * normFactorVector) / mean(sampleTotals * normFactorVector))
    }

    # Median of ratios to the geometric mean of the regions, the size factor of DESeq2
    logCounts <- log(greenlistMatrix)
    logCounts[!is.finite(logCounts)] <- NA

    referenceRow <- rowMeans(logCounts)
    usableRows <- is.finite(referenceRow)

    if (sum(usableRows) == 0) {
      stop("No greenlist region carries a read in every sample, the median of ratios has nothing to work on. Use estimator = 'sum'.", call. = FALSE)
    }

    # A median taken over a handful of regions is not a median of anything stable
    if (sum(usableRows) < 10) {
      warning("Only ", sum(usableRows), " greenlist regions carry a read in every sample, the factors rest on those alone.", call. = FALSE)
    }

    logRatios <- logCounts[usableRows, , drop = FALSE] - referenceRow[usableRows]
    sizeFactorVector <- exp(apply(logRatios, MARGIN = 2, FUN = stats::median, na.rm = TRUE))

    return(as.numeric(sizeFactorVector / mean(sizeFactorVector)))
  } # END function
