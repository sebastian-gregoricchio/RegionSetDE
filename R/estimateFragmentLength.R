# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title estimateFragmentLength
#'
#' @description Estimates the fragment length of single-end libraries from the strand cross-correlation of their reads, the same signal phantompeakqualtools and MACS2 read the fragment size from. Paired-end libraries need no estimate, and for them the median insert size is reported instead.
#'
#' @param bamFiles Character vector with the paths of the BAM files. Default: \code{NULL}, taken from \code{sampleSheet}.
#' @param regions Regions the reads are collected from: a \code{RegionSetDE} object, a \code{GRangesList}, a named list of \code{GRanges} or a single \code{GRanges}. Default: \code{NULL}, the three longest chromosomes of the BAM files.
#' @param sampleSheet Data.frame returned by \code{\link{loadSampleSheet}}, or the path to a sample sheet, providing the BAM files and the sample names. Default: \code{NULL}.
#' @param sampleNames Character vector with the sample names. Default: \code{NULL}, the BAM file names are used.
#' @param pairedEnd Logical value, one logical value per BAM file, or the string \code{"auto"} to read the layout from the files themselves. Default: \code{"auto"}.
#' @param maxDistance Numeric value with the longest fragment considered, in base pairs. Default: \code{600}.
#' @param minMapq Numeric value with the minimum mapping quality of a read. Default: \code{20}.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded. Default: \code{TRUE}.
#' @param discardRegions \code{GRanges} with regions whose reads must be ignored, for instance a blacklist. Default: \code{NULL}.
#' @param maxReads Numeric value with the maximum number of forward reads per sample entering the cross-correlation. Beyond it the reads are thinned evenly along the genome. Default: \code{1e5}.
#' @param nThreads Number of threads, one file per thread. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with three elements. \code{table} is a data.frame with one row per sample: \code{sample}, \code{paired.end}, \code{read.length}, \code{fragment.length}, \code{method} (\code{"cross-correlation"} or \code{"insert size"}) and \code{n.reads}, the number of reads or pairs the estimate rests on. \code{profile} holds the cross-correlation of the single-end samples, as the number of read pairs found at each distance and the same number relative to its mean. \code{plot} draws that profile, with the estimate marked on each sample.
#'
#' @details Each fragment leaves a read on the forward strand at its left end and a read on the reverse strand at its right end, so the distance between the 5' ends of forward and reverse reads piles up at the fragment length. The function counts, for every forward read, the reverse reads found at each distance up to \code{maxDistance}, and takes the distance with the most pairs. A count of pairs at a given distance is the numerator of the cross-correlation at that shift, the rest of the Pearson formula being the same for every shift, so the two peak at the same place.
#'
#' A second peak sits at the read length. It comes from mappability rather than from the fragments, since the reads of a region that cannot be mapped go missing on both strands at once, and on libraries with little enrichment it can be the higher of the two. The search therefore starts 15 bp past the read length. A library whose fragments are shorter than that gets an estimate at the start of the search window, and for such a library extending the reads changes little anyway.
#'
#' Only forward reads whose 5' end falls within the regions are used, while reverse reads are collected up to \code{maxDistance} beyond them, so every forward read has the same range of distances to look across and the counts carry no edge effect. The regions of an analysis, a catalogue of peaks for instance, hold far more fragments than the genome at large and give a sharper peak. A value of the estimate sitting at \code{maxDistance} means the peak was not reached, and a warning says so.
#'
#' A fragment length computed elsewhere, by phantompeakqualtools or by the peak caller, can be given to \code{\link{countReads}} as a column of the sample sheet instead.
#'
#' @examples
#' # A single-end library of 200 bp fragments piled on a few sites
#' set.seed(1)
#' fragmentStarts <- sort(unlist(lapply(seq(20000, 180000, by = 20000), function(site) {site + round(rnorm(400, 0, 60))})))
#' readStrand <- sample(c("+", "-"), length(fragmentStarts), replace = TRUE)
#' readStart <- ifelse(readStrand == "+", fragmentStarts, fragmentStarts + 200 - 50)
#' readOrder <- order(readStart)
#'
#' samRecords <- paste(paste0("read", seq_along(readOrder)), ifelse(readStrand[readOrder] == "+", 0, 16), "chrT",
#'                     as.integer(readStart[readOrder]), 60, "50M", "*", 0, 0, "*", "*", sep = "\t")
#' samFile <- file.path(tempdir(), "fragments.sam")
#' writeLines(c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chrT\tLN:200000", samRecords), samFile)
#' bamFile <- Rsamtools::asBam(samFile, destination = file.path(tempdir(), "fragments"), overwrite = TRUE, indexDestination = TRUE)
#'
#' fragmentEstimate <- estimateFragmentLength(bamFile, verbose = FALSE)
#' fragmentEstimate$table
#' fragmentEstimate$plot
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countReads}}
#'
#' @importFrom Rsamtools scanBamHeader testPairedEndBam
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges IRanges reduce
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr bind_rows filter
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line geom_vline facet_wrap labs
#' @importFrom methods is
#'
#' @export estimateFragmentLength

estimateFragmentLength <-
  function(bamFiles = NULL,
           regions = NULL,
           sampleSheet = NULL,
           sampleNames = NULL,
           pairedEnd = "auto",
           maxDistance = 600,
           minMapq = 20,
           removeDuplicates = TRUE,
           discardRegions = NULL,
           maxReads = 1e5,
           nThreads = 1,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!is.null(sampleSheet)) {
      if (!is.null(bamFiles)) {
        stop("Give the BAM files either directly or through 'sampleSheet', not both.", call. = FALSE)
      }
      if (is.character(sampleSheet) & length(sampleSheet) == 1) {
        sampleSheet <- loadSampleSheet(sampleSheet, verbose = verbose)
      }
      if (!is.data.frame(sampleSheet) || !all(c("sample", "bam") %in% colnames(sampleSheet))) {
        stop("The sample sheet needs the 'sample' and 'bam' columns.", call. = FALSE)
      }
      bamFiles <- as.character(sampleSheet$bam)
      sampleNames <- as.character(sampleSheet$sample)
    }

    if (!is.character(bamFiles) | length(bamFiles) == 0) {
      stop("The 'bamFiles' parameter must be a character vector with at least one BAM file.", call. = FALSE)
    }

    missingFiles <- bamFiles[!file.exists(bamFiles)]
    if (length(missingFiles) > 0) {
      stop("The following BAM files do not exist: ", paste(missingFiles, collapse = ", "), ".", call. = FALSE)
    }

    hasIndex <- .hasBamIndex(bamFiles)
    if (any(!hasIndex)) {
      stop("The following BAM files are not indexed: ", paste(basename(bamFiles[!hasIndex]), collapse = ", "), ".", call. = FALSE)
    }

    if (is.null(sampleNames)) {
      sampleNames <- sub("\\.bam$", "", basename(bamFiles), ignore.case = TRUE)
    }

    if (length(sampleNames) != length(bamFiles) | any(duplicated(sampleNames))) {
      stop("The 'sampleNames' parameter must hold one unique name per BAM file.", call. = FALSE)
    }

    if (!is.numeric(maxDistance) | length(maxDistance) != 1 || maxDistance < 50) {
      stop("The 'maxDistance' parameter must be a number of at least 50.", call. = FALSE)
    }
    maxDistance <- as.integer(maxDistance)

    if (!is.numeric(maxReads) | length(maxReads) != 1 || maxReads < 100) {
      stop("The 'maxReads' parameter must be a number of at least 100.", call. = FALSE)
    }

    if (!is.null(discardRegions) & !methods::is(discardRegions, "GRanges")) {
      stop("The 'discardRegions' parameter must be a GRanges object.", call. = FALSE)
    }

    if (identical(pairedEnd, "auto")) {
      pairedEnd <- vapply(bamFiles, Rsamtools::testPairedEndBam, logical(1), USE.NAMES = FALSE)
    }
    if (length(pairedEnd) == 1) {pairedEnd <- rep(pairedEnd, length(bamFiles))}

    if (!is.logical(pairedEnd) | length(pairedEnd) != length(bamFiles) | any(is.na(pairedEnd))) {
      stop("The 'pairedEnd' parameter must be 'auto', a single logical value, or one logical value per BAM file.", call. = FALSE)
    }

    #------------------------#
    # Regions to read        #
    #------------------------#
    chromosomeLengths <- Rsamtools::scanBamHeader(bamFiles[1])[[1]]$targets
    bamSeqlevels <- names(chromosomeLengths)

    if (is.null(regions)) {
      # Without regions the longest chromosomes stand for the genome
      longestChromosomes <- names(sort(chromosomeLengths, decreasing = TRUE))[seq_len(min(3, length(chromosomeLengths)))]
      readRegions <- GenomicRanges::GRanges(seqnames = longestChromosomes,
                                            ranges = IRanges::IRanges(start = 1L, end = as.integer(chromosomeLengths[longestChromosomes])))
    } else {
      readRegions <- .regionsAsGRanges(regions)
      readRegions <- .matchSeqlevels(x = readRegions, targetSeqlevels = bamSeqlevels, fileName = bamFiles[1], verbose = FALSE)
    }

    readRegions <- readRegions[as.character(GenomeInfoDb::seqnames(readRegions)) %in% bamSeqlevels]

    if (length(readRegions) == 0) {
      stop("No region lies on a chromosome of the BAM files.", call. = FALSE)
    }

    readRegions <- IRanges::reduce(readRegions, ignore.strand = TRUE)

    if (!is.null(discardRegions)) {
      discardRegions <- .matchSeqlevels(x = discardRegions, targetSeqlevels = bamSeqlevels, fileName = bamFiles[1], verbose = FALSE)
    }

    if (isTRUE(verbose) & any(!pairedEnd)) {
      message("Estimating the fragment length of ", sum(!pairedEnd), " single-end sample(s) from their strand cross-correlation...")
    }

    if (isTRUE(verbose) & any(pairedEnd)) {
      message("Reading the insert size of ", sum(pairedEnd), " paired-end sample(s)...")
    }

    #------------------------#
    # One file per thread    #
    #------------------------#
    estimateList <-
      BiocParallel::bplapply(seq_along(bamFiles),
                             function(fileIndex) {
                               if (pairedEnd[fileIndex]) {
                                 .insertSizeEstimate(bamFile = bamFiles[fileIndex],
                                                     regions = readRegions,
                                                     chromosomeLengths = chromosomeLengths,
                                                     minMapq = minMapq,
                                                     removeDuplicates = removeDuplicates,
                                                     discardRegions = discardRegions)
                               } else {
                                 .crossCorrelationEstimate(bamFile = bamFiles[fileIndex],
                                                           regions = readRegions,
                                                           chromosomeLengths = chromosomeLengths,
                                                           maxDistance = maxDistance,
                                                           minMapq = minMapq,
                                                           removeDuplicates = removeDuplicates,
                                                           discardRegions = discardRegions,
                                                           maxReads = maxReads)
                               }
                             },
                             BPPARAM = .makeParallelParam(nThreads = nThreads, tasks = length(bamFiles)))

    #------------------------#
    # Assemble the output    #
    #------------------------#
    estimateTable <- data.frame(sample = sampleNames,
                                paired.end = pairedEnd,
                                read.length = vapply(estimateList, function(x) {as.numeric(x$read.length)}, numeric(1)),
                                fragment.length = vapply(estimateList, function(x) {as.numeric(x$fragment.length)}, numeric(1)),
                                method = ifelse(pairedEnd, "insert size", "cross-correlation"),
                                n.reads = vapply(estimateList, function(x) {as.numeric(x$n.reads)}, numeric(1)),
                                stringsAsFactors = FALSE)

    profileTable <- dplyr::bind_rows(lapply(which(!pairedEnd),
                                            function(fileIndex) {
                                              profile <- estimateList[[fileIndex]]$profile
                                              if (is.null(profile)) {return(NULL)}
                                              profile$sample <- sampleNames[fileIndex]
                                              return(profile[, c("sample", "distance", "pairs", "relative.pairs")])
                                            }))

    # Problems are gathered and reported once, a warning per sample would bury the table
    noReads <- estimateTable$sample[is.na(estimateTable$fragment.length)]
    if (length(noReads) > 0) {
      warning("No usable read was found for: ", paste(noReads, collapse = ", "), ". Their fragment length is NA.", call. = FALSE)
    }

    atBoundary <- estimateTable$sample[!estimateTable$paired.end & !is.na(estimateTable$fragment.length) & estimateTable$fragment.length >= maxDistance - 5]
    if (length(atBoundary) > 0) {
      warning("The cross-correlation of ", paste(atBoundary, collapse = ", "),
              " peaks at the end of the search window: raise 'maxDistance', or give the fragment length directly.", call. = FALSE)
    }

    fewReads <- estimateTable$sample[!estimateTable$paired.end & !is.na(estimateTable$fragment.length) & estimateTable$n.reads < 1000]
    if (length(fewReads) > 0) {
      warning("Fewer than 1,000 forward reads entered the estimate of: ", paste(fewReads, collapse = ", "),
              ". The value is unreliable, widen 'regions' or give the fragment length directly.", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message("Fragment lengths: ", paste(estimateTable$sample, round(estimateTable$fragment.length), sep = " ", collapse = ", "), ".")
    }

    #------------------------#
    # Profile plot           #
    #------------------------#
    profilePlot <- NULL
    if (nrow(profileTable) > 0) {
      # The panels follow the order of the samples rather than the alphabet
      profileTable$sample <- factor(profileTable$sample, levels = sampleNames)

      markTable <- dplyr::filter(estimateTable, !.data$paired.end, !is.na(.data$fragment.length))

      profilePlot <-
        ggplot2::ggplot(data = profileTable, mapping = ggplot2::aes(x = .data$distance, y = .data$relative.pairs)) +
        ggplot2::geom_line(colour = "grey30", linewidth = 0.4) +
        ggplot2::geom_vline(data = markTable, mapping = ggplot2::aes(xintercept = .data$fragment.length), colour = "#D1495B", linetype = 2) +
        ggplot2::geom_vline(data = markTable, mapping = ggplot2::aes(xintercept = .data$read.length), colour = "grey60", linetype = 3) +
        ggplot2::facet_wrap(~ sample, scales = "free_y") +
        ggplot2::labs(x = "Distance between the ends of forward and reverse reads (bp)",
                      y = "Read pairs, relative to the mean",
                      title = "Strand cross-correlation",
                      subtitle = "Dashed: fragment length estimate; dotted: read length") +
        .regionSetTheme()
    }

    return(list(table = estimateTable,
                profile = profileTable,
                plot = profilePlot))
  } # END function




#' @title .crossCorrelationEstimate
#'
#' @description Estimates the fragment length of one single-end BAM file from the distances between the 5' ends of its forward and reverse reads.
#'
#' @param bamFile String with the path of the BAM file.
#' @param regions \code{GRanges} with the regions, named after the chromosomes of the BAM file and without overlaps.
#' @param chromosomeLengths Named numeric vector with the length of every chromosome of the BAM file.
#' @param maxDistance Integer value with the longest distance considered.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param discardRegions \code{GRanges} with the regions whose reads must be ignored, or \code{NULL}.
#' @param maxReads Numeric value with the maximum number of forward reads used.
#'
#' @return A list with \code{read.length}, \code{fragment.length}, \code{n.reads}, the number of forward reads used, and \code{profile}, a data.frame with the number of pairs at every distance.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges IRanges findOverlaps overlapsAny
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom BiocGenerics width strand
#' @importFrom GenomeInfoDb seqnames
#' @importFrom stats median
#' @importFrom utils modifyList
#'
#' @keywords internal

.crossCorrelationEstimate <-
  function(bamFile,
           regions,
           chromosomeLengths,
           maxDistance,
           minMapq,
           removeDuplicates,
           discardRegions,
           maxReads) {

    emptyEstimate <- list(read.length = NA_real_, fragment.length = NA_real_, n.reads = 0, profile = NULL)

    # Reverse reads are needed up to maxDistance past the regions, the forward ones only inside them
    reads <- .windowFragments(bamFile = bamFile,
                              windows = regions,
                              isPairedEnd = FALSE,
                              fragmentLength = 1L,
                              maxFragmentLength = maxDistance,
                              minMapq = minMapq,
                              removeDuplicates = removeDuplicates,
                              chromosomeLengths = chromosomeLengths,
                              discardRegions = discardRegions,
                              padding = maxDistance,
                              readsOnly = TRUE)

    if (length(reads) == 0) {return(emptyEstimate)}

    readLength <- stats::median(BiocGenerics::width(reads))
    readStrand <- as.character(BiocGenerics::strand(reads))
    fivePrime <- S4Vectors::mcols(reads)$five.prime

    forwardReads <- reads[readStrand == "+" & IRanges::overlapsAny(GenomicRanges::GRanges(GenomeInfoDb::seqnames(reads), IRanges::IRanges(fivePrime, width = 1)),
                                                                     regions, ignore.strand = TRUE)]
    reverseReads <- reads[readStrand == "-"]

    if (length(forwardReads) == 0 | length(reverseReads) == 0) {return(emptyEstimate)}

    #--------------------------#
    # Thin the forward reads   #
    #--------------------------#
    # Every k-th read along the genome, which keeps the estimate reproducible without touching the random seed
    forwardReads <- forwardReads[order(as.character(GenomeInfoDb::seqnames(forwardReads)), S4Vectors::mcols(forwardReads)$five.prime)]
    if (length(forwardReads) > maxReads) {
      forwardReads <- forwardReads[seq(1, length(forwardReads), by = ceiling(length(forwardReads) / maxReads))]
    }

    #--------------------------#
    # Pairs at every distance  #
    #--------------------------#
    # A fragment of length d leaves its reverse 5' end at d - 1 bases from the forward one
    pairCounts <- numeric(maxDistance)
    forwardChromosomes <- as.character(GenomeInfoDb::seqnames(forwardReads))
    reverseChromosomes <- as.character(GenomeInfoDb::seqnames(reverseReads))

    for (chromosome in unique(forwardChromosomes)) {
      forwardEnds <- S4Vectors::mcols(forwardReads)$five.prime[forwardChromosomes == chromosome]
      reverseEnds <- IRanges::IRanges(start = S4Vectors::mcols(reverseReads)$five.prime[reverseChromosomes == chromosome], width = 1)
      if (length(reverseEnds) == 0) {next}

      # In chunks, so that dense peaks cannot blow up the number of hits held at once
      for (chunk in split(seq_along(forwardEnds), ceiling(seq_along(forwardEnds) / 20000))) {
        pairHits <- IRanges::findOverlaps(IRanges::IRanges(start = forwardEnds[chunk], width = maxDistance), reverseEnds)
        pairDistance <- BiocGenerics::start(reverseEnds)[S4Vectors::subjectHits(pairHits)] - forwardEnds[chunk][S4Vectors::queryHits(pairHits)] + 1L
        pairCounts <- pairCounts + tabulate(pairDistance, nbins = maxDistance)
      }
    }

    if (sum(pairCounts) == 0) {return(utils::modifyList(emptyEstimate, list(read.length = readLength)))}

    #--------------------------#
    # Peak past the read length#
    #--------------------------#
    # A running mean over 11 bp keeps a single noisy distance from winning, and the phantom peak is left behind
    smoothedCounts <- as.numeric(stats::filter(pairCounts, rep(1 / 11, 11), sides = 2))
    searchStart <- min(as.integer(readLength + 16), maxDistance)
    searchRange <- seq(searchStart, maxDistance)
    searchValues <- smoothedCounts[searchRange]

    fragmentLength <- if (all(is.na(searchValues))) {searchStart} else {searchRange[which.max(searchValues)]}

    return(list(read.length = readLength,
                fragment.length = fragmentLength,
                n.reads = length(forwardReads),
                profile = data.frame(distance = seq_len(maxDistance),
                                     pairs = pairCounts,
                                     relative.pairs = pairCounts / mean(pairCounts))))
  } # END function




#' @title .insertSizeEstimate
#'
#' @description Reports the median insert size of the proper pairs of one paired-end BAM file over a set of regions.
#'
#' @param bamFile String with the path of the BAM file.
#' @param regions \code{GRanges} with the regions, named after the chromosomes of the BAM file.
#' @param chromosomeLengths Named numeric vector with the length of every chromosome of the BAM file.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param discardRegions \code{GRanges} with the regions whose reads must be ignored, or \code{NULL}.
#'
#' @return A list with \code{read.length} (\code{NA}), \code{fragment.length} and \code{n.reads}, the number of pairs.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocGenerics width
#' @importFrom stats median
#'
#' @keywords internal

.insertSizeEstimate <-
  function(bamFile,
           regions,
           chromosomeLengths,
           minMapq,
           removeDuplicates,
           discardRegions) {

    fragments <- .windowFragments(bamFile = bamFile,
                                  windows = regions,
                                  isPairedEnd = TRUE,
                                  fragmentLength = NA_integer_,
                                  maxFragmentLength = 2000L,
                                  minMapq = minMapq,
                                  removeDuplicates = removeDuplicates,
                                  chromosomeLengths = chromosomeLengths,
                                  discardRegions = discardRegions)

    if (length(fragments) == 0) {
      return(list(read.length = NA_real_, fragment.length = NA_real_, n.reads = 0))
    }

    return(list(read.length = NA_real_,
                fragment.length = stats::median(BiocGenerics::width(fragments)),
                n.reads = length(fragments)))
  } # END function




#' @title .regionsAsGRanges
#'
#' @description Pools the regions of a \code{RegionSetDE} object, a \code{GRangesList} or a list of \code{GRanges} into a single \code{GRanges}, without metadata.
#'
#' @param regions A \code{RegionSetDE} object, a \code{GRangesList}, a list of \code{GRanges} or a \code{GRanges}.
#'
#' @return A \code{GRanges}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRangesList granges
#' @importFrom methods is
#'
#' @keywords internal

.regionsAsGRanges <-
  function(regions) {

    if (methods::is(regions, "RegionSetDE")) {
      regions <- regions@regions
    }

    if (methods::is(regions, "GRanges")) {
      return(GenomicRanges::granges(regions))
    }

    if (is.list(regions) & !methods::is(regions, "GRangesList")) {
      if (!all(vapply(regions, function(x) {methods::is(x, "GRanges")}, logical(1)))) {
        stop("The regions must be a RegionSetDE object, a GRangesList, a list of GRanges or a GRanges.", call. = FALSE)
      }
      regions <- GenomicRanges::GRangesList(lapply(regions, GenomicRanges::granges))
    }

    if (!methods::is(regions, "GRangesList")) {
      stop("The regions must be a RegionSetDE object, a GRangesList, a list of GRanges or a GRanges.", call. = FALSE)
    }

    return(GenomicRanges::granges(unlist(regions, use.names = FALSE)))
  } # END function
