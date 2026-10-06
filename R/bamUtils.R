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



#' @title .bamIsPairedEnd
#'
#' @description Tells whether BAM files hold paired-end reads, from the flags of their first records. In a paired-end library every read carries the paired flag, so the head of the file answers for all of it. \code{Rsamtools::testPairedEndBam} stops at the first paired read as well, but on a single-end file it finds none and reads to the last record, a million at a time, printing the running total.
#'
#' @param bamFiles Character vector with the paths of the BAM files.
#' @param nRecords Numeric value with the number of records read from the head of each file. Default: \code{1e5}.
#'
#' @return A logical vector with one value per file, \code{FALSE} for a file without records.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools BamFile scanBam ScanBamParam bamFlagTest
#'
#' @keywords internal

.bamIsPairedEnd <-
  function(bamFiles,
           nRecords = 1e5) {

    return(vapply(bamFiles,
                  function(bamFile) {
                    # No index is needed to read a file from its first record
                    headFlags <- Rsamtools::scanBam(Rsamtools::BamFile(file = bamFile, yieldSize = as.integer(nRecords)),
                                                    param = Rsamtools::ScanBamParam(what = "flag"))[[1]]$flag

                    return(any(Rsamtools::bamFlagTest(headFlags, "isPaired")))
                  },
                  logical(1),
                  USE.NAMES = FALSE))
  } # END function




#' @title .bamChromosomeMap
#'
#' @description Brings the chromosomes of a group of BAM files under one set of names, those of the first file or of a reference given by the caller. Files aligned to the same assembly do not always name its chromosomes alike, \code{chr1} in one header and \code{1} in the next, and the ranges being counted can only be written one way. Every file is then read under its own names and reported under the common ones.
#'
#' @param bamFiles Character vector with the paths of the BAM files.
#' @param referenceTargets Named vector with the chromosome lengths the names are taken from. Default: \code{NULL}, the header of the first file.
#' @param referenceLabel String naming the reference in the error messages. Default: \code{NULL}, the name of the first file.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{FALSE}.
#'
#' @return A list with \code{lengths}, a named integer vector with the length of every chromosome found in the reference or in a file, under the common names; \code{names}, a list with one named character vector per file, giving the name each of its chromosomes carries in the file (the names of the vector are the common ones); \code{shared}, the common names of the chromosomes every file has; and \code{renamed}, a logical vector telling which files are read under names of their own.
#'
#' @details Two files that give different lengths to the same chromosome are not on the same assembly, and the counting stops there. A contig that some files lack, under any name, is no reason to stop: scaffolds and decoys differ between builds of the same assembly, and the files without it simply hold no read there.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools scanBamHeader
#' @importFrom stats setNames
#' @importFrom utils head
#'
#' @keywords internal

.bamChromosomeMap <-
  function(bamFiles,
           referenceTargets = NULL,
           referenceLabel = NULL,
           verbose = FALSE) {

    targetList <- lapply(bamFiles, function(bamFile) {Rsamtools::scanBamHeader(bamFile)[[1]]$targets})

    if (is.null(referenceTargets)) {
      referenceTargets <- targetList[[1]]
      if (is.null(referenceLabel)) {referenceLabel <- paste0("'", basename(bamFiles[1]), "'")}
    }

    if (is.null(referenceLabel)) {referenceLabel <- "the reference"}

    chromosomeLengths <- stats::setNames(as.integer(referenceTargets), names(referenceTargets))
    fileNames <- vector(mode = "list", length = length(bamFiles))
    renamedFiles <- logical(length(bamFiles))

    #-------------------------------#
    # One file at a time            #
    #-------------------------------#
    for (fileIndex in seq_along(bamFiles)) {
      fileSeqlevels <- names(targetList[[fileIndex]])
      fileLengths <- as.integer(targetList[[fileIndex]])

      translatedNames <- .translateChromosomeNames(chromosomeNames = fileSeqlevels, targetSeqlevels = names(chromosomeLengths))

      # A header listing both chrM and MT would send two chromosomes to one name, the second keeps its own
      exactNames <- !is.na(translatedNames) & translatedNames == fileSeqlevels
      collidingNames <- !is.na(translatedNames) & !exactNames &
        (translatedNames %in% translatedNames[exactNames] | duplicated(ifelse(exactNames, NA_character_, translatedNames)))
      translatedNames[collidingNames] <- NA_character_

      commonNames <- ifelse(is.na(translatedNames), fileSeqlevels, translatedNames)
      knownNames <- commonNames %in% names(chromosomeLengths)

      if (!any(knownNames)) {
        stop("The chromosome names of '", basename(bamFiles[fileIndex]), "' cannot be reconciled with the ones of ", referenceLabel,
             ": ", paste(utils::head(fileSeqlevels, 3), collapse = ", "), " against ",
             paste(utils::head(names(chromosomeLengths), 3), collapse = ", "), ".", call. = FALSE)
      }

      # The same chromosome with two lengths means two assemblies, and the regions would land elsewhere
      differentLengths <- which(knownNames & chromosomeLengths[commonNames] != fileLengths)

      if (length(differentLengths) > 0) {
        stop("The BAM files are not aligned to the same assembly: '", basename(bamFiles[fileIndex]), "' and ", referenceLabel,
             " give different lengths to ", paste(utils::head(commonNames[differentLengths], 3), collapse = ", "), ".", call. = FALSE)
      }

      chromosomeLengths <- c(chromosomeLengths, stats::setNames(fileLengths[!knownNames], commonNames[!knownNames]))
      fileNames[[fileIndex]] <- stats::setNames(fileSeqlevels, commonNames)
      renamedFiles[fileIndex] <- any(commonNames != fileSeqlevels)
    }

    sharedNames <- Reduce(f = intersect, x = lapply(fileNames, names))

    #-------------------------------#
    # What the user should know     #
    #-------------------------------#
    if (isTRUE(verbose) & any(renamedFiles)) {
      message("The BAM files do not name their chromosomes alike. The names of ", referenceLabel, " (",
              paste(utils::head(names(referenceTargets), 2), collapse = ", "), ") stand for all of them; files read under their own names (",
              paste(utils::head(fileNames[[which(renamedFiles)[1]]], 2), collapse = ", "), "): ", sum(renamedFiles), " of ", length(bamFiles), ".")
    }

    partialNames <- setdiff(unique(unlist(lapply(fileNames, names), use.names = FALSE)), sharedNames)

    if (isTRUE(verbose) & length(partialNames) > 0) {
      message("Chromosomes or contigs missing from at least one BAM file under any name: ", length(partialNames), " (",
              paste(utils::head(partialNames, 3), collapse = ", "), if (length(partialNames) > 3) {", ..."} else {""},
              "). Their reads enter the library sizes of the files that have them, list them in 'excludeChromosomes' to leave them out.")
    }

    return(list(lengths = chromosomeLengths,
                names = fileNames,
                shared = sharedNames,
                renamed = renamedFiles))
  } # END function




#' @title .countBamFragments
#'
#' @description Counts the fragments of a group of BAM files over a set of ranges. The chromosomes are cut into pieces of at most 50 Mb, and the pieces of all the files are shared among the threads, so that even a single file keeps every thread busy. Paired-end fragments are rebuilt from the first mate of each proper pair, whose position, mate position and template length (TLEN) give the start and the width of the fragment: the two reads never have to be matched. Single-end reads are extended to the fragment length from their 5' end.
#'
#' @param bamFiles Character vector with the paths of the BAM files, all aligned to the same assembly. Their chromosomes may be named in different styles, see \code{.bamChromosomeMap}.
#' @param ranges \code{GRanges} with the ranges to count, named after the chromosomes of the first BAM file, or of \code{referenceTargets}. The strand is ignored.
#' @param pairedEnd Logical vector with one value per BAM file.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended, or one value per BAM file. Default: \code{150}.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment. Default: \code{1000}.
#' @param minMapq Numeric value with the minimum mapping quality of a read. For paired-end data the mate is checked through its \code{MQ} tag, when the file carries it. Default: \code{20}.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded. Default: \code{TRUE}.
#' @param excludeChromosomes Character vector with the chromosomes left out of the library sizes, named as in the BAM files. The ranges lying on them are still counted. Default: \code{NULL}, none.
#' @param discardRegions \code{GRanges} with the regions whose reads must be ignored, named after the chromosomes of the BAM files. A fragment is dropped when one of its reads starts inside them. Default: \code{NULL}.
#' @param fullLibrarySize Logical value: \code{TRUE} reads every chromosome that is not excluded to compute the library sizes, \code{FALSE} only the chromosomes carrying ranges, which is faster but leaves the library sizes partial. Default: \code{TRUE}.
#' @param countMode String with the way a fragment is assigned to the ranges: \code{"overlap"} counts it in every range it overlaps, \code{"bin"} counts it once, at its centre for paired-end data and at the 5' end of the read for single-end data, as csaw does for genome wide bins. Default: \code{"overlap"}.
#' @param pieceLength Numeric value with the maximum length of the stretch of genome read by a single job, in base pairs. Default: \code{5e7}.
#' @param referenceTargets Named vector with the chromosome lengths whose names the ranges, the excluded chromosomes and the discarded regions are written in. Default: \code{NULL}, the header of the first BAM file.
#' @param progressBar Logical value to indicate whether a progress bar must be drawn, one step for every job that comes back. Default: \code{FALSE}.
#' @param nThreads Number of threads. Default: \code{1}.
#'
#' @return A list with three elements: \code{counts}, an integer matrix with one row per range and one column per file; \code{library.size}, the number of fragments that went through the filters on the chromosomes read and not excluded; \code{mate.mapq.found}, telling for each paired-end file whether the \code{MQ} tag was found (\code{NA} for single-end files).
#'
#' @author Sebastian Gregoricchio
#'
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
           referenceTargets = NULL,
           progressBar = FALSE,
           nThreads = 1) {

    #--------------------------#
    # Chromosomes to be read   #
    #--------------------------#
    # The pieces are cut once for all the files, under one set of names. Each file then reads them under
    # its own, which is what lets alignments named chr1 and alignments named 1 be counted together
    bamChromosomes <- .bamChromosomeMap(bamFiles = bamFiles, referenceTargets = referenceTargets)

    # A chromosome is read when it carries ranges to count, or, for the full library sizes, whenever it is not excluded
    rangeChromosomes <- as.character(GenomeInfoDb::seqnames(ranges))

    chromosomeTable <- data.frame(chromosome = names(bamChromosomes$lengths),
                                  length = as.numeric(bamChromosomes$lengths),
                                  stringsAsFactors = FALSE)

    # A chromosome of the reference that no file has would only make empty jobs
    chromosomeTable <- dplyr::filter(chromosomeTable, .data$chromosome %in% unlist(lapply(bamChromosomes$names, names), use.names = FALSE))

    chromosomeTable <- dplyr::mutate(chromosomeTable,
                                     has.ranges = .data$chromosome %in% rangeChromosomes,
                                     in.library = !(.data$chromosome %in% excludeChromosomes))

    chromosomeTable <- dplyr::filter(chromosomeTable, .data$has.ranges | (isTRUE(fullLibrarySize) & .data$in.library))

    countMatrix <- matrix(0L, nrow = length(ranges), ncol = length(bamFiles))
    librarySizes <- numeric(length(bamFiles))
    mateMapqFound <- ifelse(pairedEnd, FALSE, NA)

    # One length per file, so that libraries of different fragment sizes are each extended to their own
    fragmentLength <- rep_len(as.integer(round(fragmentLength)), length(bamFiles))

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

               # The name each chromosome carries in this file, NA when the file does not have it
               jobPieces$file.chromosome <- unname(bamChromosomes$names[[jobTable$file.index[j]]][jobPieces$chromosome])
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
    jobResults <- BiocParallel::bplapply(jobList,
                                         .readBamFragments,
                                         bamFiles = bamFiles,
                                         pairedEnd = pairedEnd,
                                         fragmentLength = fragmentLength,
                                         maxFragmentLength = maxFragmentLength,
                                         minMapq = minMapq,
                                         removeDuplicates = removeDuplicates,
                                         countMode = countMode,
                                         BPPARAM = .makeParallelParam(nThreads = nThreads, tasks = length(jobList), progressBar = progressBar))

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
#' @param job List describing the job: \code{file.index}, the \code{pieces} table, whose \code{file.chromosome} column holds the name each chromosome carries in the file, the \code{ranges} to count with their \code{range.chromosome} and \code{range.index}, and the \code{discard} regions of its chromosomes.
#' @param bamFiles Character vector with the paths of the BAM files.
#' @param pairedEnd Logical vector with one value per BAM file.
#' @param fragmentLength Integer vector with the length to which single-end reads are extended, one value per BAM file.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param countMode String, either \code{"overlap"} or \code{"bin"}, see \code{.countBamFragments}.
#'
#' @return A list with the \code{file.index} and \code{range.index} of the job, the \code{counts} of its ranges, the \code{total.fragments} that went through the filters and \code{mate.mapq.found}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools scanBam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges countOverlaps
#' @importFrom dplyr filter
#' @importFrom rlang .data
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

    # A chromosome the file does not have holds no read, and asking for it would be an error
    pieces <- dplyr::filter(job$pieces, !is.na(.data$file.chromosome))

    jobCounts <- integer(length(job$range.index))
    totalFragments <- 0
    mateMapqFound <- FALSE

    if (nrow(pieces) == 0) {
      return(list(file.index = job$file.index,
                  range.index = job$range.index,
                  counts = jobCounts,
                  total.fragments = totalFragments,
                  mate.mapq.found = mateMapqFound))
    }

    #--------------------------#
    # Read the pieces          #
    #--------------------------#
    # The file is asked under its own chromosome names, everything else speaks the common ones
    readParameters <- .bamReadParameters(which = GenomicRanges::GRanges(seqnames = pieces$file.chromosome,
                                                                        ranges = IRanges::IRanges(start = pieces$start, end = pieces$end)),
                                         isPairedEnd = isPairedEnd,
                                         minMapq = minMapq,
                                         removeDuplicates = removeDuplicates)

    # All the pieces go in a single call: querying an open BamFile a second time returns no reads
    readList <- Rsamtools::scanBam(.bamWithIndex(bamFiles[job$file.index]), param = readParameters)
    readList <- readList[paste0(pieces$file.chromosome, ":", pieces$start, "-", pieces$end)]

    for (i in seq_len(nrow(pieces))) {
      reads <- readList[[i]]

      # A read belongs to the piece holding its start, the neighbouring piece would count it again otherwise
      pieceFragments <- .fragmentsFromReads(reads = reads,
                                            keep = reads$pos >= pieces$start[i] & reads$pos <= pieces$end[i],
                                            isPairedEnd = isPairedEnd,
                                            fragmentLength = fragmentLength[job$file.index],
                                            maxFragmentLength = maxFragmentLength,
                                            minMapq = minMapq,
                                            chromosomeLength = pieces$chromosome.length[i],
                                            discard = job$discard[[pieces$chromosome[i]]])

      if (isTRUE(pieceFragments$mate.mapq.found)) {mateMapqFound <- TRUE}

      # An excluded chromosome is read for its ranges only, its fragments stay out of the library size
      if (pieces$in.library[i]) {
        totalFragments <- totalFragments + length(pieceFragments$start)
      }

      #--------------------------#
      # Count over the ranges    #
      #--------------------------#
      rangeSlot <- which(job$range.chromosome == pieces$chromosome[i])

      if (length(rangeSlot) > 0 & length(pieceFragments$start) > 0) {
        # In bin mode a fragment counts once, at its centre, or at the 5' end of a single read
        fragmentRanges <- if (countMode == "bin") {
          IRanges::IRanges(start = pieceFragments$point, width = 1)
        } else {
          IRanges::IRanges(start = pieceFragments$start, end = pieceFragments$end)
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




#' @title .bamReadParameters
#'
#' @description Builds the \code{ScanBamParam} shared by every function reading fragments from a BAM file, so that counting, summits and profiles apply the same filters. Of a proper pair only the first mate is read, since its position, the position of its mate and the template length already describe the fragment.
#'
#' @param which \code{GRanges} with the stretches of genome to read.
#' @param isPairedEnd Logical value, \code{TRUE} for a paired-end file.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#'
#' @return A \code{ScanBamParam} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#'
#' @keywords internal

.bamReadParameters <-
  function(which,
           isPairedEnd,
           minMapq,
           removeDuplicates) {

    checkMapq <- isTRUE(!is.na(minMapq[1]) & minMapq[1] > 0)

    return(Rsamtools::ScanBamParam(what = if (isPairedEnd) {c("pos", "mpos", "isize")} else {c("pos", "strand", "cigar")},
                                   tag = if (isPairedEnd & checkMapq) {"MQ"} else {character(0)},
                                   which = which,
                                   mapqFilter = if (checkMapq) {as.integer(minMapq[1])} else {NA_integer_},
                                   flag = Rsamtools::scanBamFlag(isPaired = if (isPairedEnd) {TRUE} else {NA},
                                                                 isProperPair = if (isPairedEnd) {TRUE} else {NA},
                                                                 isFirstMateRead = if (isPairedEnd) {TRUE} else {NA},
                                                                 isUnmappedQuery = FALSE,
                                                                 isSecondaryAlignment = FALSE,
                                                                 isSupplementaryAlignment = FALSE,
                                                                 isDuplicate = if (isTRUE(removeDuplicates)) {FALSE} else {NA})))
  } # END function




#' @title .fragmentsFromReads
#'
#' @description Turns the records returned by \code{scanBam} into fragments. Paired-end fragments are rebuilt from the first mate of each proper pair, single-end reads are extended from their 5' end to the fragment length. The records of a pair longer than the maximum, of a mate below the mapping quality, or starting in a discarded region are dropped.
#'
#' @param reads List returned by \code{scanBam} for one stretch of genome, or several of them pasted together.
#' @param keep Logical vector with one value per record, telling which records enter at all.
#' @param isPairedEnd Logical value, \code{TRUE} for a paired-end file.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment.
#' @param minMapq Numeric value with the minimum mapping quality of a read, applied to the mate through its \code{MQ} tag.
#' @param chromosomeLength Numeric value, or one value per record, with the length of the chromosome the fragments are clipped to.
#' @param discard \code{IRanges} with the discarded regions of the chromosome, or \code{NULL}. Only for records of a single chromosome.
#'
#' @return A list with the \code{start}, the \code{end} and the \code{point} of every fragment, the point being its centre for paired-end data and the 5' end of the read for single-end data, the \code{index} of the record each fragment comes from, and \code{mate.mapq.found}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges IRanges overlapsAny
#'
#' @keywords internal

.fragmentsFromReads <-
  function(reads,
           keep,
           isPairedEnd,
           fragmentLength,
           maxFragmentLength,
           minMapq,
           chromosomeLength,
           discard = NULL) {

    checkMapq <- isTRUE(!is.na(minMapq[1]) & minMapq[1] > 0)
    mateMapqFound <- FALSE
    keep[is.na(keep)] <- FALSE

    #--------------------------#
    # Filters of the records   #
    #--------------------------#
    if (isPairedEnd) {
      keep <- keep & reads$isize != 0 & abs(reads$isize) <= maxFragmentLength

      # The quality of the mate is only known through the MQ tag, without it the mate is let through
      mateMapq <- reads$tag$MQ
      if (checkMapq & !is.null(mateMapq)) {
        mateMapqFound <- TRUE
        keep <- keep & (is.na(mateMapq) | mateMapq >= minMapq[1])
      }

      readStarts <- list(reads$pos, reads$mpos)
    } else {
      readEnd <- reads$pos + .cigarReferenceWidth(reads$cigar) - 1L
      onMinus <- as.character(reads$strand) == "-"
      readStarts <- list(reads$pos)
    }

    # A read starting in a discarded region takes its whole fragment away
    if (!is.null(discard)) {
      for (readStart in readStarts) {
        keep <- keep & !IRanges::overlapsAny(IRanges::IRanges(start = readStart, width = 1), discard)
      }
    }

    #--------------------------#
    # Fragment boundaries      #
    #--------------------------#
    if (isPairedEnd) {
      fragmentStart <- pmin(reads$pos, reads$mpos)[keep]
      fragmentEnd <- fragmentStart + abs(reads$isize[keep]) - 1L
    } else {
      fragmentLength <- as.integer(fragmentLength[1])
      fragmentStart <- ifelse(onMinus, readEnd - fragmentLength + 1L, reads$pos)[keep]
      fragmentEnd <- ifelse(onMinus, readEnd, reads$pos + fragmentLength - 1L)[keep]
    }

    if (length(chromosomeLength) > 1) {chromosomeLength <- chromosomeLength[keep]}

    fragmentStart <- as.integer(pmax(fragmentStart, 1L))
    fragmentEnd <- as.integer(pmin(fragmentEnd, chromosomeLength))

    # The centre of a pair, once clipped to the chromosome, or the 5' end of a single read
    fragmentPoint <- if (isPairedEnd) {fragmentStart + (fragmentEnd - fragmentStart + 1L) %/% 2L} else {ifelse(onMinus, readEnd, reads$pos)[keep]}

    return(list(start = fragmentStart,
                end = fragmentEnd,
                point = as.integer(fragmentPoint),
                index = which(keep),
                mate.mapq.found = mateMapqFound))
  } # END function




#' @title .windowFragments
#'
#' @description Reads the fragments of one BAM file that overlap a set of windows, with the same filters as the counting. The windows are widened by the longest fragment accepted, so that a fragment reaching a window from outside is not lost, and merged, so that every record is read once.
#'
#' @param bamFile String with the path of the BAM file.
#' @param windows \code{GRanges} with the windows, named as in \code{chromosomeLengths}.
#' @param isPairedEnd Logical value, \code{TRUE} for a paired-end file.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param chromosomeLengths Named numeric vector with the length of every chromosome, as returned by \code{.bamChromosomeMap}. The file may name the same chromosomes in another style, it is read under its own names and the fragments come back under these.
#' @param discardRegions \code{GRanges} with the regions whose reads must be ignored, named as in \code{chromosomeLengths}. Default: \code{NULL}.
#' @param padding Numeric value with the number of base pairs added on both sides of the windows. Default: \code{NULL}, the longest fragment accepted.
#' @param readsOnly Logical value: \code{TRUE} returns the reads as they are, with their strand and 5' end, instead of the fragments. Only for single-end files. Default: \code{FALSE}.
#'
#' @return A \code{GRanges} with one element per fragment, or per read when \code{readsOnly = TRUE}, carrying the 5' end of the read in the \code{five.prime} column.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom Rsamtools scanBam
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom IRanges IRanges ranges reduce
#'
#' @keywords internal

.windowFragments <-
  function(bamFile,
           windows,
           isPairedEnd,
           fragmentLength,
           maxFragmentLength,
           minMapq,
           removeDuplicates,
           chromosomeLengths,
           discardRegions = NULL,
           padding = NULL,
           readsOnly = FALSE) {

    emptyResult <- GenomicRanges::GRanges(seqnames = character(0), ranges = IRanges::IRanges(), strand = character(0), five.prime = integer(0))

    # The names the chromosomes carry in this file, which may differ from the ones the windows are written in
    fileChromosomes <- .bamChromosomeMap(bamFiles = bamFile, referenceTargets = chromosomeLengths)$names[[1]]

    windows <- windows[as.character(GenomeInfoDb::seqnames(windows)) %in% names(fileChromosomes)]
    if (length(windows) == 0) {return(emptyResult)}

    #--------------------------#
    # Windows to read          #
    #--------------------------#
    # Merged after the widening, so that the windows never share a record and each read is taken once
    if (is.null(padding)) {padding <- max(c(maxFragmentLength, fragmentLength), na.rm = TRUE)}
    windowChromosomes <- as.character(GenomeInfoDb::seqnames(windows))

    readWindows <- GenomicRanges::GRanges(seqnames = windowChromosomes,
                                          ranges = IRanges::IRanges(start = as.integer(pmax(BiocGenerics::start(windows) - padding, 1)),
                                                                    end = as.integer(pmin(BiocGenerics::end(windows) + padding, chromosomeLengths[windowChromosomes]))))
    readWindows <- IRanges::reduce(readWindows, ignore.strand = TRUE)

    fileWindowChromosomes <- unname(fileChromosomes[as.character(GenomeInfoDb::seqnames(readWindows))])

    readParameters <- .bamReadParameters(which = GenomicRanges::GRanges(seqnames = fileWindowChromosomes,
                                                                        ranges = IRanges::ranges(readWindows)),
                                         isPairedEnd = isPairedEnd,
                                         minMapq = minMapq,
                                         removeDuplicates = removeDuplicates)

    readList <- Rsamtools::scanBam(.bamWithIndex(bamFile), param = readParameters)
    readList <- readList[paste0(fileWindowChromosomes, ":", BiocGenerics::start(readWindows), "-", BiocGenerics::end(readWindows))]

    #--------------------------#
    # One chromosome at a time #
    #--------------------------#
    # The records of a chromosome are pooled, a read belonging to the window holding its start
    readWindowChromosomes <- as.character(GenomeInfoDb::seqnames(readWindows))
    discardList <- if (is.null(discardRegions) || length(discardRegions) == 0) {NULL} else {
      split(IRanges::ranges(discardRegions), as.character(GenomeInfoDb::seqnames(discardRegions)))
    }

    fragmentList <- list()

    for (chromosome in unique(readWindowChromosomes)) {
      windowIndex <- which(readWindowChromosomes == chromosome)
      chromosomeReads <- readList[windowIndex]

      recordNumber <- vapply(chromosomeReads, function(reads) {length(reads$pos)}, integer(1))
      if (sum(recordNumber) == 0) {next}

      pooledReads <- list(pos = unlist(lapply(chromosomeReads, `[[`, "pos"), use.names = FALSE))
      if (isPairedEnd) {
        pooledReads$mpos <- unlist(lapply(chromosomeReads, `[[`, "mpos"), use.names = FALSE)
        pooledReads$isize <- unlist(lapply(chromosomeReads, `[[`, "isize"), use.names = FALSE)
        mateMapq <- lapply(chromosomeReads, function(reads) {reads$tag$MQ})
        if (!all(vapply(mateMapq, is.null, logical(1)))) {
          pooledReads$tag <- list(MQ = unlist(lapply(seq_along(mateMapq), function(j) {
            if (is.null(mateMapq[[j]])) {rep(NA_integer_, recordNumber[j])} else {mateMapq[[j]]}
          }), use.names = FALSE))
        }
      } else {
        pooledReads$strand <- unlist(lapply(chromosomeReads, function(reads) {as.character(reads$strand)}), use.names = FALSE)
        pooledReads$cigar <- unlist(lapply(chromosomeReads, `[[`, "cigar"), use.names = FALSE)
      }

      windowStart <- rep(BiocGenerics::start(readWindows)[windowIndex], times = recordNumber)
      windowEnd <- rep(BiocGenerics::end(readWindows)[windowIndex], times = recordNumber)
      inWindow <- pooledReads$pos >= windowStart & pooledReads$pos <= windowEnd

      chromosomeFragments <- .fragmentsFromReads(reads = pooledReads,
                                                 keep = inWindow,
                                                 isPairedEnd = isPairedEnd,
                                                 fragmentLength = if (isTRUE(readsOnly)) {1L} else {fragmentLength},
                                                 maxFragmentLength = maxFragmentLength,
                                                 minMapq = minMapq,
                                                 chromosomeLength = chromosomeLengths[[chromosome]],
                                                 discard = discardList[[chromosome]])

      if (length(chromosomeFragments$start) == 0) {next}

      if (isTRUE(readsOnly)) {
        # The read itself, with its strand, is what a cross-correlation looks at
        readIndex <- chromosomeFragments$index
        readStart <- pooledReads$pos[readIndex]
        readEnd <- readStart + .cigarReferenceWidth(pooledReads$cigar[readIndex]) - 1L
        fragmentList[[chromosome]] <- GenomicRanges::GRanges(seqnames = chromosome,
                                                             ranges = IRanges::IRanges(start = readStart, end = readEnd),
                                                             strand = pooledReads$strand[readIndex],
                                                             five.prime = chromosomeFragments$point)
      } else {
        fragmentList[[chromosome]] <- GenomicRanges::GRanges(seqnames = chromosome,
                                                             ranges = IRanges::IRanges(start = chromosomeFragments$start, end = chromosomeFragments$end),
                                                             strand = "*",
                                                             five.prime = chromosomeFragments$point)
      }
    }

    if (length(fragmentList) == 0) {return(emptyResult)}

    return(do.call(what = c, args = unname(fragmentList)))
  } # END function




#' @title .windowCoverage
#'
#' @description Reads the fragments of one BAM file around a set of windows and returns their coverage, chromosome by chromosome.
#'
#' @param bamFile String with the path of the BAM file.
#' @param windows \code{GRanges} with the windows, named as in \code{chromosomeLengths}.
#' @param chromosomeLengths Named numeric vector with the length of every chromosome, as returned by \code{.bamChromosomeMap}.
#' @param ... Read filters passed to \code{.windowFragments}.
#'
#' @return A list with \code{coverage}, an \code{RleList} with one element per chromosome of \code{chromosomeLengths}, and \code{fragments}, the number of fragments read.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges coverage
#' @importFrom GenomeInfoDb seqlevels<- seqlengths<-
#'
#' @keywords internal

.windowCoverage <-
  function(bamFile,
           windows,
           chromosomeLengths,
           ...) {

    fragments <- .windowFragments(bamFile = bamFile, windows = windows, chromosomeLengths = chromosomeLengths, ...)

    # Every chromosome gets an element, so that a window on a chromosome without fragments reads zero instead of failing
    GenomeInfoDb::seqlevels(fragments) <- names(chromosomeLengths)
    GenomeInfoDb::seqlengths(fragments) <- as.integer(chromosomeLengths)

    return(list(coverage = GenomicRanges::coverage(fragments),
                fragments = length(fragments)))
  } # END function




#' @title .bamSummits
#'
#' @description Finds, in every BAM file, the summit of the fragment pileup within each region: the middle of the highest stretch of coverage, and its height.
#'
#' @param bamFiles Character vector with the paths of the BAM files.
#' @param regions \code{GRanges} with the regions, named after the chromosomes of the BAM files.
#' @param pairedEnd Logical vector with one value per BAM file.
#' @param fragmentLength Numeric vector with the length to which single-end reads are extended, one value per BAM file.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param discardRegions \code{GRanges} with the regions whose reads must be ignored, or \code{NULL}.
#' @param progressBar Logical value to indicate whether a progress bar must be drawn, one step for every file. Default: \code{FALSE}.
#' @param nThreads Number of threads, one file per thread. Default: \code{1}.
#'
#' @return A list with \code{position} and \code{height}, two matrices with one row per region and one column per file, and \code{fragments}, the number of fragments read in each file.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom BiocParallel bplapply
#' @importFrom GenomeInfoDb seqnames
#' @importFrom IRanges Views viewMaxs viewRangeMaxs ranges
#' @importFrom BiocGenerics start end
#'
#' @keywords internal

.bamSummits <-
  function(bamFiles,
           regions,
           pairedEnd,
           fragmentLength,
           maxFragmentLength,
           minMapq,
           removeDuplicates,
           discardRegions = NULL,
           progressBar = FALSE,
           nThreads = 1) {

    # The names of the first file stand for all of them, each file being read under its own
    chromosomeLengths <- .bamChromosomeMap(bamFiles = bamFiles)$lengths
    fragmentLength <- rep_len(fragmentLength, length(bamFiles))

    regionChromosomes <- as.character(GenomeInfoDb::seqnames(regions))
    regionIndexList <- split(seq_along(regions), factor(regionChromosomes, levels = unique(regionChromosomes)))
    regionIndexList <- regionIndexList[names(regionIndexList) %in% names(chromosomeLengths)]

    summitList <-
      BiocParallel::bplapply(seq_along(bamFiles),
                             function(fileIndex) {
                               fileCoverage <- .windowCoverage(bamFile = bamFiles[fileIndex],
                                                               windows = regions,
                                                               chromosomeLengths = chromosomeLengths,
                                                               isPairedEnd = pairedEnd[fileIndex],
                                                               fragmentLength = fragmentLength[fileIndex],
                                                               maxFragmentLength = maxFragmentLength,
                                                               minMapq = minMapq,
                                                               removeDuplicates = removeDuplicates,
                                                               discardRegions = discardRegions)

                               summitPosition <- rep(NA_integer_, length(regions))
                               summitHeight <- rep(0, length(regions))

                               # The middle of the highest plateau, a flat top would otherwise pull the summit to its left edge
                               for (chromosome in names(regionIndexList)) {
                                 regionIndex <- regionIndexList[[chromosome]]
                                 coverageViews <- IRanges::Views(fileCoverage$coverage[[chromosome]], IRanges::ranges(regions[regionIndex]))
                                 plateauRanges <- IRanges::viewRangeMaxs(coverageViews)
                                 summitHeight[regionIndex] <- as.numeric(IRanges::viewMaxs(coverageViews))
                                 summitPosition[regionIndex] <- as.integer((BiocGenerics::start(plateauRanges) + BiocGenerics::end(plateauRanges)) %/% 2L)
                               }

                               summitPosition[summitHeight == 0] <- NA_integer_

                               list(position = summitPosition, height = summitHeight, fragments = fileCoverage$fragments)
                             },
                             BPPARAM = .makeParallelParam(nThreads = nThreads, tasks = length(bamFiles), progressBar = progressBar))

    return(list(position = matrix(vapply(summitList, `[[`, integer(length(regions)), "position"), nrow = length(regions)),
                height = matrix(vapply(summitList, `[[`, numeric(length(regions)), "height"), nrow = length(regions)),
                fragments = vapply(summitList, `[[`, numeric(1), "fragments")))
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
    simpleCigar <- grepl("^[0-9]+M$", cigar)
    referenceWidth <- rep(NA_integer_, length(cigar))
    referenceWidth[simpleCigar] <- as.integer(sub("M$", "", cigar[simpleCigar]))
    complexCigar <- which(!simpleCigar)

    if (length(complexCigar) > 0) {
      operationList <- regmatches(cigar[complexCigar], gregexpr("[0-9]+[MDN=X]", cigar[complexCigar]))
      referenceWidth[complexCigar] <- vapply(operationList, function(operations) {sum(as.integer(sub("[MDN=X]$", "", operations)))}, integer(1))
    }

    return(referenceWidth)
  } # END function
