# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title .bamIndexPath
#'
#' @description Looks for the index of a BAM file under the four names one can go under: \code{file.bam.bai}, \code{file.bai}, \code{file.bam.csi} and \code{file.csi}.
#'
#' @param bamFile String with the path of the BAM file.
#'
#' @return String with the path of the first index found, \code{NA} when the file has none.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.bamIndexPath <-
  function(bamFile) {

    indexPaths <- c(paste0(bamFile, ".bai"),
                    sub("\\.bam$", ".bai", bamFile, ignore.case = TRUE),
                    paste0(bamFile, ".csi"),
                    sub("\\.bam$", ".csi", bamFile, ignore.case = TRUE))

    indexFound <- indexPaths[file.exists(indexPaths)]

    if (length(indexFound) == 0) {
      return(NA_character_)
    }

    indexFound[1]
  }



#' @title .bamWithIndex
#'
#' @description Opens a BAM file with its index, handing the path over explicitly. Rsamtools finds a BAI on its own and leaves a CSI alone, so a file indexed with \code{samtools index -c} cannot be read by region unless it is told where the index is. CSI is not an exotic case: BAI cannot address a contig longer than 512 Mb at all, which rules it out for several plant and amphibian assemblies.
#'
#' @param bamFile String with the path of the BAM file.
#'
#' @return A \code{BamFile} carrying the index that was found, or the path unchanged when there is none.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools BamFile
#'
#' @keywords internal

.bamWithIndex <-
  function(bamFile) {

    indexPath <- .bamIndexPath(bamFile)

    if (is.na(indexPath)) {
      return(bamFile)
    }

    Rsamtools::BamFile(file = bamFile, index = indexPath)
  }



#' @title .countBamFragments
#'
#' @description Counts the fragments of a group of BAM files over a set of ranges. The chromosomes are cut into pieces of at most 50 Mb, and the pieces of all the files are shared among the threads, so that even a single file keeps every thread busy. Paired-end fragments are rebuilt from the first mate of each proper pair, whose position, mate position and template length (TLEN) give the start and the width of the fragment: the two reads never have to be matched. Single-end reads are extended to the fragment length from their 5' end.
#'
#' @param bamFiles Character vector with the paths of the BAM files, all sharing the same header.
#' @param ranges \code{GRanges} with the ranges to count, named after the chromosomes of the BAM files. The strand is ignored.
#' @param pairedEnd Logical vector with one value per BAM file.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended. Default: \code{150}.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment. Default: \code{1000}.
#' @param minMapq Numeric value with the minimum mapping quality of a read. For paired-end data the mate is checked through its \code{MQ} tag, when the file carries it. Default: \code{20}.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded. Default: \code{TRUE}.
#' @param excludeChromosomes Character vector with the chromosomes left out of the library sizes, named as in the BAM files. The ranges lying on them are still counted. Default: \code{NULL}, none.
#' @param discardRegions \code{GRanges} with the regions whose reads must be ignored, named after the chromosomes of the BAM files. A fragment is dropped when one of its reads starts inside them. Default: \code{NULL}.
#' @param fullLibrarySize Logical value: \code{TRUE} reads every chromosome that is not excluded to compute the library sizes, \code{FALSE} only the chromosomes carrying ranges, which is faster but leaves the library sizes partial. Default: \code{TRUE}.
#' @param countMode String with the way a fragment is assigned to the ranges: \code{"overlap"} counts it in every range it overlaps, \code{"bin"} counts it once, at its centre for paired-end data and at the 5' end of the read for single-end data, as csaw does for genome wide bins. Default: \code{"overlap"}.
#' @param pieceLength Numeric value with the maximum length of the stretch of genome read by a single job, in base pairs. Default: \code{5e7}.
#' @param nThreads Number of threads. Default: \code{1}.
#'
#' @return A list with three elements: \code{counts}, an integer matrix with one row per range and one column per file; \code{library.size}, the number of fragments that went through the filters on the chromosomes read and not excluded; \code{mate.mapq.found}, telling for each paired-end file whether the \code{MQ} tag was found (\code{NA} for single-end files).
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools scanBamHeader
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges ranges
#' @importFrom BiocParallel bplapply
#' @importFrom dplyr filter mutate arrange desc
#' @importFrom rlang .data
#'
#' @keywords internal

.countBamFragments <-
  function(bamFiles,
           ranges,
           pairedEnd,
           fragmentLength = 150,
           maxFragmentLength = 1000,
           minMapq = 20,
           removeDuplicates = TRUE,
           excludeChromosomes = NULL,
           discardRegions = NULL,
           fullLibrarySize = TRUE,
           countMode = "overlap",
           pieceLength = 5e7,
           nThreads = 1) {

    #--------------------------#
    # Chromosomes to be read   #
    #--------------------------#
    # The pieces are cut once for all the files, which only works if they share the same chromosomes
    targetList <- lapply(Rsamtools::scanBamHeader(bamFiles), function(header) {header$targets})
    sameHeader <- vapply(targetList, identical, logical(1), targetList[[1]])

    if (any(!sameHeader)) {
      stop("The BAM files do not share the same chromosomes and lengths: ", paste(basename(bamFiles[!sameHeader]), collapse = ", "),
           " differ from ", basename(bamFiles[1]), ".", call. = FALSE)
    }

    # A chromosome is read when it carries ranges to count, or, for the full library sizes, whenever it is not excluded
    rangeChromosomes <- as.character(GenomeInfoDb::seqnames(ranges))

    chromosomeTable <- data.frame(chromosome = names(targetList[[1]]),
                                  length = as.numeric(targetList[[1]]),
                                  stringsAsFactors = FALSE)

    chromosomeTable <- dplyr::mutate(chromosomeTable,
                                     has.ranges = .data$chromosome %in% rangeChromosomes,
                                     in.library = !(.data$chromosome %in% excludeChromosomes))

    chromosomeTable <- dplyr::filter(chromosomeTable, .data$has.ranges | (isTRUE(fullLibrarySize) & .data$in.library))

    countMatrix <- matrix(0L, nrow = length(ranges), ncol = length(bamFiles))
    librarySizes <- numeric(length(bamFiles))
    mateMapqFound <- ifelse(pairedEnd, FALSE, NA)

    if (nrow(chromosomeTable) == 0) {
      return(list(counts = countMatrix, library.size = librarySizes, mate.mapq.found = mateMapqFound))
    }

    #--------------------------#
    # Cut the genome in pieces #
    #--------------------------#
    # Long chromosomes are cut, so that a single file still spreads over all the threads.
    # Integer coordinates keep the piece names identical to the ones scanBam gives back (5e+07 would not match).
    pieceLength <- as.integer(pieceLength)
    piecesPerChromosome <- as.integer(ceiling(chromosomeTable$length / pieceLength))

    pieceTable <- data.frame(chromosome = rep(chromosomeTable$chromosome, times = piecesPerChromosome),
                             chromosome.length = as.integer(rep(chromosomeTable$length, times = piecesPerChromosome)),
                             in.library = rep(chromosomeTable$in.library, times = piecesPerChromosome),
                             start = as.integer(unlist(lapply(piecesPerChromosome, function(n) {(seq_len(n) - 1L) * pieceLength + 1L}))),
                             stringsAsFactors = FALSE)

    pieceTable <- dplyr::mutate(pieceTable,
                                end = pmin(.data$start + pieceLength - 1L, .data$chromosome.length),
                                width = .data$end - .data$start + 1L)
    pieceTable <- dplyr::arrange(pieceTable, dplyr::desc(.data$width))

    # Short pieces are pooled until they add up to a full one, otherwise thousands of contigs would become thousands of jobs
    pieceJob <- integer(nrow(pieceTable))
    currentJob <- 1L
    coveredLength <- 0

    for (i in seq_len(nrow(pieceTable))) {
      if (coveredLength >= pieceLength) {
        currentJob <- currentJob + 1L
        coveredLength <- 0
      }
      pieceJob[i] <- currentJob
      coveredLength <- coveredLength + pieceTable$width[i]
    }

    pieceTable$job <- pieceJob

    #--------------------------#
    # One job per file and     #
    # group of pieces          #
    #--------------------------#
    rangeIndexList <- split(seq_along(ranges), factor(rangeChromosomes, levels = unique(rangeChromosomes)))

    discardList <- NULL
    if (!is.null(discardRegions)) {
      if (length(discardRegions) > 0) {
        discardList <- split(IRanges::ranges(discardRegions), as.character(GenomeInfoDb::seqnames(discardRegions)))
      }
    }

    # The biggest pieces go first and the files alternate, so that no thread is left alone with a long job at the end
    jobTable <- expand.grid(file.index = seq_along(bamFiles), job = unique(pieceJob), stringsAsFactors = FALSE)

    jobList <-
      lapply(seq_len(nrow(jobTable)),
             function(j) {
               jobPieces <- dplyr::filter(pieceTable, .data$job == jobTable$job[j])
               jobChromosomes <- unique(jobPieces$chromosome)
               jobRangeIndex <- unlist(rangeIndexList[intersect(jobChromosomes, names(rangeIndexList))], use.names = FALSE)
               if (is.null(jobRangeIndex)) {jobRangeIndex <- integer(0)}

               list(file.index = jobTable$file.index[j],
                    pieces = jobPieces,
                    range.index = jobRangeIndex,
                    ranges = IRanges::ranges(ranges[jobRangeIndex]),
                    range.chromosome = rangeChromosomes[jobRangeIndex],
                    discard = if (is.null(discardList)) {NULL} else {discardList[intersect(jobChromosomes, names(discardList))]})
             })

    #--------------------------#
    # Read and collect         #
    #--------------------------#
    jobResults <- .bplapplySameLibraries(jobList,
                                         .readBamFragments,
                                         bamFiles = bamFiles,
                                         pairedEnd = pairedEnd,
                                         fragmentLength = fragmentLength,
                                         maxFragmentLength = maxFragmentLength,
                                         minMapq = minMapq,
                                         removeDuplicates = removeDuplicates,
                                         countMode = countMode,
                                         BPPARAM = .makeParallelParam(nThreads = nThreads, tasks = length(jobList)))

    # A chromosome cut in several pieces sends its ranges to several jobs, whose counts add up
    for (jobResult in jobResults) {
      fileIndex <- jobResult$file.index
      countMatrix[jobResult$range.index, fileIndex] <- countMatrix[jobResult$range.index, fileIndex] + jobResult$counts
      librarySizes[fileIndex] <- librarySizes[fileIndex] + jobResult$total.fragments
      if (isTRUE(jobResult$mate.mapq.found)) {mateMapqFound[fileIndex] <- TRUE}
    }

    return(list(counts = countMatrix, library.size = librarySizes, mate.mapq.found = mateMapqFound))
  } # END function




#' @title .readBamFragments
#'
#' @description Reads the fragments of one BAM file over a group of pieces of genome and counts them over the ranges of the same chromosomes. It is the job run by every thread of \code{.countBamFragments}.
#'
#' @param job List describing the job: \code{file.index}, the \code{pieces} table, the \code{ranges} to count with their \code{range.chromosome} and \code{range.index}, and the \code{discard} regions of its chromosomes.
#' @param bamFiles Character vector with the paths of the BAM files.
#' @param pairedEnd Logical vector with one value per BAM file.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param countMode String, either \code{"overlap"} or \code{"bin"}, see \code{.countBamFragments}.
#'
#' @return A list with the \code{file.index} and \code{range.index} of the job, the \code{counts} of its ranges, the \code{total.fragments} that went through the filters and \code{mate.mapq.found}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges countOverlaps overlapsAny
#'
#' @keywords internal

.readBamFragments <-
  function(job,
           bamFiles,
           pairedEnd,
           fragmentLength,
           maxFragmentLength,
           minMapq,
           removeDuplicates,
           countMode) {

    isPairedEnd <- pairedEnd[job$file.index]
    pieces <- job$pieces
    checkMapq <- isTRUE(!is.na(minMapq[1]) & minMapq[1] > 0)

    jobCounts <- integer(length(job$range.index))
    totalFragments <- 0
    mateMapqFound <- FALSE

    #--------------------------#
    # Read the pieces          #
    #--------------------------#
    # Of a proper pair only the first mate is read, it already knows where the fragment starts and how long it is
    readParameters <- Rsamtools::ScanBamParam(what = if (isPairedEnd) {c("pos", "mpos", "isize")} else {c("pos", "strand", "cigar")},
                                              tag = if (isPairedEnd & checkMapq) {"MQ"} else {character(0)},
                                              which = GenomicRanges::GRanges(seqnames = pieces$chromosome,
                                                                             ranges = IRanges::IRanges(start = pieces$start, end = pieces$end)),
                                              mapqFilter = if (checkMapq) {as.integer(minMapq)} else {NA_integer_},
                                              flag = Rsamtools::scanBamFlag(isPaired = if (isPairedEnd) {TRUE} else {NA},
                                                                            isProperPair = if (isPairedEnd) {TRUE} else {NA},
                                                                            isFirstMateRead = if (isPairedEnd) {TRUE} else {NA},
                                                                            isUnmappedQuery = FALSE,
                                                                            isSecondaryAlignment = FALSE,
                                                                            isSupplementaryAlignment = FALSE,
                                                                            isDuplicate = if (isTRUE(removeDuplicates)) {FALSE} else {NA}))

    # All the pieces go in a single call: querying an open BamFile a second time returns no reads
    readList <- Rsamtools::scanBam(.bamWithIndex(bamFiles[job$file.index]), param = readParameters)
    readList <- readList[paste0(pieces$chromosome, ":", pieces$start, "-", pieces$end)]

    for (i in seq_len(nrow(pieces))) {
      reads <- readList[[i]]

      # A read belongs to the piece holding its start, the neighbouring piece would count it again otherwise
      keep <- reads$pos >= pieces$start[i] & reads$pos <= pieces$end[i]

      #--------------------------#
      # Fragment boundaries      #
      #--------------------------#
      if (isPairedEnd) {
        keep <- keep & reads$isize != 0 & abs(reads$isize) <= maxFragmentLength

        # The quality of the mate is only known through the MQ tag, without it the mate is let through
        mateMapq <- reads$tag$MQ
        if (checkMapq & !is.null(mateMapq)) {
          mateMapqFound <- TRUE
          keep <- keep & (is.na(mateMapq) | mateMapq >= minMapq)
        }

        readStarts <- list(reads$pos, reads$mpos)
      } else {
        readEnd <- reads$pos + .cigarReferenceWidth(reads$cigar) - 1L
        onMinus <- as.character(reads$strand) == "-"
        readStarts <- list(reads$pos)
      }

      # A read starting in a discarded region takes its whole fragment away
      discardHere <- job$discard[[pieces$chromosome[i]]]
      if (!is.null(discardHere)) {
        for (readStart in readStarts) {
          keep <- keep & !IRanges::overlapsAny(IRanges::IRanges(start = readStart, width = 1), discardHere)
        }
      }

      if (isPairedEnd) {
        fragmentStart <- pmin(reads$pos, reads$mpos)[keep]
        fragmentEnd <- fragmentStart + abs(reads$isize[keep]) - 1L
      } else {
        fragmentStart <- ifelse(onMinus, readEnd - fragmentLength + 1L, reads$pos)[keep]
        fragmentEnd <- ifelse(onMinus, readEnd, reads$pos + fragmentLength - 1L)[keep]
      }

      fragmentStart <- pmax(fragmentStart, 1L)
      fragmentEnd <- pmin(fragmentEnd, pieces$chromosome.length[i])

      # An excluded chromosome is read for its ranges only, its fragments stay out of the library size
      if (pieces$in.library[i]) {
        totalFragments <- totalFragments + length(fragmentStart)
      }

      #--------------------------#
      # Count over the ranges    #
      #--------------------------#
      rangeSlot <- which(job$range.chromosome == pieces$chromosome[i])

      if (length(rangeSlot) > 0 & length(fragmentStart) > 0) {
        if (countMode == "bin") {
          # One count per fragment, at its centre, or at the 5' end of a single read
          countPoint <- if (isPairedEnd) {fragmentStart + (fragmentEnd - fragmentStart + 1L) %/% 2L} else {ifelse(onMinus, readEnd, reads$pos)[keep]}
          fragmentRanges <- IRanges::IRanges(start = countPoint, width = 1)
        } else {
          fragmentRanges <- IRanges::IRanges(start = fragmentStart, end = fragmentEnd)
        }

        jobCounts[rangeSlot] <- jobCounts[rangeSlot] + IRanges::countOverlaps(job$ranges[rangeSlot], fragmentRanges)
      }
    }

    return(list(file.index = job$file.index,
                range.index = job$range.index,
                counts = jobCounts,
                total.fragments = totalFragments,
                mate.mapq.found = mateMapqFound))
  } # END function




#' @title .cigarReferenceWidth
#'
#' @description Returns the number of reference bases covered by each alignment, from its CIGAR string.
#'
#' @param cigar Character vector with the CIGAR strings.
#'
#' @return An integer vector with one width per alignment.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.cigarReferenceWidth <-
  function(cigar) {
    # Most reads are a single match, only the others need their operations added up
    referenceWidth <- suppressWarnings(as.integer(sub("M$", "", cigar)))
    complexCigar <- which(is.na(referenceWidth))

    if (length(complexCigar) > 0) {
      operationList <- regmatches(cigar[complexCigar], gregexpr("[0-9]+[MDN=X]", cigar[complexCigar]))
      referenceWidth[complexCigar] <- vapply(operationList, function(operations) {sum(as.integer(sub("[MDN=X]$", "", operations)))}, integer(1))
    }

    return(referenceWidth)
  } # END function
