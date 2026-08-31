#' @title countBackground
#'
#' @description Counts the reads falling in large genome wide bins, which provide the background distribution used to estimate composition-aware normalisation factors. The bins overlapping the counted regions are removed by default, so that the normalisation is not driven by the signal under study. The result is stored in the metadata of the counts object, where the normalisation step retrieves it.
#'
#' @param counts \code{RegionSetDE.counts} object returned by \code{\link{countReads}}.
#' @param bamFiles Character vector with the paths of the BAM files, in the same order as the samples of \code{counts}. Default: \code{NULL}, the files recorded by \code{\link{countReads}} are reused.
#' @param binSize Numeric value with the width of the bins, in base pairs. Default: \code{10000}.
#' @param excludeRegions Logical value indicating whether the bins overlapping the regions of \code{counts} must be discarded. Default: \code{TRUE}.
#' @param minCount Numeric value with the minimum total count required to keep a bin. Default: \code{1}.
#' @param restrictChromosomes Character vector with the chromosomes to read from the BAM files. Default: \code{NULL}, the value used at the counting step.
#' @param pairedEnd Logical value, or one logical value per BAM file, indicating whether the reads must be counted as proper pairs. Default: \code{NULL}, the layouts resolved at the counting step.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended. Default: \code{NULL}, the value used at the counting step.
#' @param maxFragmentLength Numeric value with the maximum insert size accepted for a pair. Default: \code{NULL}, the value used at the counting step.
#' @param minMapq Numeric value with the minimum mapping quality of a read. Default: \code{NULL}, the value used at the counting step.
#' @param removeDuplicates Logical value indicating whether the duplicated reads must be discarded. Default: \code{NULL}, the value used at the counting step.
#' @param nThreads Number of threads used to process the files in parallel. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return The input \code{RegionSetDE.counts} object with the bin counts stored as a \code{RangedSummarizedExperiment} in \code{metadata(counts)$background}.
#'
#' @details Bins of ten kilobases or more are wide enough that most of them carry background reads only, and their counts therefore track the amount of sequencing spent outside the regions of interest. Reusing the read parameters of \code{\link{countReads}} matters here: bins counted with a different mapping quality or duplicate policy would return factors that do not apply to the region counts. The parameters are taken from the object unless they are given explicitly.
#'
#' @examples
#' # The example counts already carry their background bins
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' backgroundBins <- S4Vectors::metadata(counts)$background
#' backgroundBins
#'
#' \dontrun{
#' # Recomputing them needs the BAM files the object was counted from
#' counts <- countBackground(counts, binSize = 10000, nThreads = 4)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countReads}}
#'
#' @importFrom csaw readParam windowCounts
#' @importFrom SummarizedExperiment assay rowRanges colData colData<- SummarizedExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqinfo seqlevels
#' @importFrom IRanges overlapsAny
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export countBackground

countBackground <-
  function(counts,
           bamFiles = NULL,
           binSize = 10000,
           excludeRegions = TRUE,
           minCount = 1,
           restrictChromosomes = NULL,
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

    binSize <- as.integer(binSize[1])
    if (is.na(binSize) | binSize < 1) {
      stop("The 'binSize' parameter must be a positive integer.", call. = FALSE)
    }

    #-----------------------------------#
    # Recover the counting parameters   #
    #-----------------------------------#
    # Bins counted under a different read filter would give factors that do not apply to the regions
    countingParameters <- counts@parameters$countReads

    if (is.null(countingParameters) & is.null(bamFiles)) {
      stop("No BAM file is recorded in the object, provide them through the 'bamFiles' parameter.", call. = FALSE)
    }

    # ifelse would keep the first element only, and the layout is one value per file
    if (is.null(bamFiles)) {bamFiles <- countingParameters$bamFiles}
    if (is.null(pairedEnd)) {pairedEnd <- if (is.null(countingParameters$pairedEnd)) {FALSE} else {countingParameters$pairedEnd}}
    if (is.null(fragmentLength)) {fragmentLength <- if (is.null(countingParameters$fragmentLength)) {150} else {countingParameters$fragmentLength}}
    if (is.null(maxFragmentLength)) {maxFragmentLength <- if (is.null(countingParameters$maxFragmentLength)) {1000} else {countingParameters$maxFragmentLength}}
    if (is.null(minMapq)) {minMapq <- if (is.null(countingParameters$minMapq)) {20} else {countingParameters$minMapq}}
    if (is.null(removeDuplicates)) {removeDuplicates <- if (is.null(countingParameters$removeDuplicates)) {TRUE} else {countingParameters$removeDuplicates}}
    if (is.null(restrictChromosomes)) {restrictChromosomes <- countingParameters$restrictChromosomes}

    if (length(bamFiles) != ncol(counts)) {
      stop("The number of BAM files does not match the number of samples of the counts object.", call. = FALSE)
    }

    if (length(pairedEnd) == 1) {pairedEnd <- rep(pairedEnd, length(bamFiles))}

    if (!is.logical(pairedEnd) | length(pairedEnd) != length(bamFiles) | any(is.na(pairedEnd))) {
      stop("The 'pairedEnd' parameter must be a single logical value, or one logical value per BAM file.", call. = FALSE)
    }

    missingFiles <- bamFiles[!file.exists(bamFiles)]
    if (length(missingFiles) > 0) {
      stop(paste0("The following BAM files do not exist: ", paste(missingFiles, collapse = ", "), "."), call. = FALSE)
    }

    parallelParam <- .makeParallelParam(nThreads = nThreads)

    #-------------------#
    # Count in the bins #
    #-------------------#
    if (isTRUE(verbose)) {
      message(paste0("Counting reads in ", format(binSize, big.mark = ","), " bp bins across ", length(bamFiles), " samples..."))
    }

    backgroundList <-
      lapply(unique(pairedEnd),
             function(layout) {
               readParameters <- csaw::readParam(pe = ifelse(layout, "both", "none"),
                                                 max.frag = maxFragmentLength,
                                                 dedup = removeDuplicates,
                                                 minq = minMapq,
                                                 restrict = restrictChromosomes)

               csaw::windowCounts(bam.files = bamFiles[pairedEnd == layout],
                                  bin = TRUE,
                                  width = binSize,
                                  filter = 0,
                                  ext = ifelse(layout, NA, as.integer(fragmentLength[1])),
                                  param = readParameters,
                                  BPPARAM = parallelParam)
             })

    # The bins come from the BAM headers, files aligned against different chromosome sizes cannot share them
    if (length(backgroundList) > 1) {
      if (!identical(GenomeInfoDb::seqinfo(backgroundList[[1]]), GenomeInfoDb::seqinfo(backgroundList[[2]])) |
          nrow(backgroundList[[1]]) != nrow(backgroundList[[2]])) {
        stop("The paired-end and single-end BAM files do not share the same chromosome sizes, their background bins cannot be merged.", call. = FALSE)
      }
    }

    # The passes are stacked by layout, this index puts the samples back in the order of the counts object
    columnOrder <- order(unlist(lapply(unique(pairedEnd), function(layout) {which(pairedEnd == layout)})))

    backgroundCounts <-
      SummarizedExperiment::SummarizedExperiment(assays = list(counts = do.call(cbind, lapply(backgroundList, function(x) {SummarizedExperiment::assay(x, "counts")}))[, columnOrder, drop = FALSE]),
                                                 rowRanges = SummarizedExperiment::rowRanges(backgroundList[[1]]),
                                                 colData = do.call(rbind, lapply(backgroundList, SummarizedExperiment::colData))[columnOrder, , drop = FALSE],
                                                 metadata = S4Vectors::metadata(backgroundList[[1]]))

    colnames(backgroundCounts) <- colnames(counts)

    # The bins inherit the names of the BAM headers, without this the overlap with the regions would find nothing
    backgroundCounts <- .matchSeqlevels(x = backgroundCounts,
                                        targetSeqlevels = GenomeInfoDb::seqlevels(counts),
                                        fileName = bamFiles[1],
                                        verbose = verbose)

    #-------------------#
    # Filter the bins   #
    #-------------------#
    binTable <- data.frame(bin.index = seq_len(nrow(backgroundCounts)),
                           total.count = as.numeric(rowSums(SummarizedExperiment::assay(backgroundCounts, "counts"))),
                           overlaps.regions = IRanges::overlapsAny(SummarizedExperiment::rowRanges(backgroundCounts),
                                                                   SummarizedExperiment::rowRanges(counts)),
                           stringsAsFactors = FALSE)

    keptBins <- dplyr::filter(binTable, .data$total.count >= minCount)

    # A bin holding one of the regions carries its signal, which is exactly what the background must not contain
    if (isTRUE(excludeRegions)) {
      keptBins <- dplyr::filter(keptBins, !.data$overlaps.regions)
    }

    if (nrow(keptBins) == 0) {
      stop("No background bin survived the filters, lower 'minCount' or reduce 'binSize'.", call. = FALSE)
    }

    backgroundCounts <- backgroundCounts[keptBins$bin.index, ]

    #---------------------------#
    # Store it in the object    #
    #---------------------------#
    S4Vectors::metadata(counts)$background <- backgroundCounts

    counts@parameters <- c(counts@parameters,
                           list(countBackground = list(bamFiles = bamFiles,
                                                       binSize = binSize,
                                                       excludeRegions = excludeRegions,
                                                       minCount = minCount,
                                                       pairedEnd = pairedEnd,
                                                       fragmentLength = fragmentLength,
                                                       maxFragmentLength = maxFragmentLength,
                                                       minMapq = minMapq,
                                                       removeDuplicates = removeDuplicates,
                                                       restrictChromosomes = restrictChromosomes)))

    if (isTRUE(verbose)) {
      message(paste0("Done. ", format(nrow(backgroundCounts), big.mark = ","), " background bins retained out of ",
                     format(nrow(binTable), big.mark = ","), "."))
    }

    return(counts)
  } # END function
