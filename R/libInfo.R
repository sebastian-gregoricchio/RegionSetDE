# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title libInfo
#'
#' @description Summarises the libraries of a counts object in one table: how many reads each BAM file holds, how many fragments went through the read filters, how many of those fall in the regions, and the fraction of reads in regions (FRiP) that follows from the two.
#'
#' @param counts \code{RegionSetDE.counts} object returned by \code{\link{countReads}}, or a \code{RegionSetDE.fit}, whose counts are used.
#' @param annotationColumns Character vector with the columns of the \code{colData} to add after the sample names, for instance \code{"condition"}. Default: \code{NULL}, none.
#' @param bamFiles Character vector with the BAM files, in the order of the samples. Default: \code{NULL}, the files recorded by \code{\link{countReads}}.
#'
#' @return A data.frame with one row per sample: \code{sample}, the \code{annotationColumns}, \code{paired.end}, \code{fragment.length} (the length single-end reads were extended to, \code{NA} for paired-end samples), \code{bam.reads} (every record of the BAM file), \code{bam.mapped} (the mapped ones), \code{library.size} (the fragments that went through the filters of the counting), \code{reads.in.regions} and \code{FRiP}, the ratio of the last two. When \code{\link{countReads}} discarded the reads of a blacklist, of a greylist or of other regions, \code{discarded.reads} follows \code{library.size} with the fragments each sample lost to them. When the inputs were counted, \code{input.id}, \code{input.library.size}, \code{input.in.regions} and \code{input.FRiP} give the same numbers for the input of each sample.
#'
#' @details The three counts measure different things and are not expected to agree. \code{bam.reads} and \code{bam.mapped} come from the index of each file, so they are read in an instant, and they count alignment records: a paired-end fragment is two of them. \code{library.size} and \code{reads.in.regions} come from the counting, where a paired-end fragment counts once and only after the mapping quality, duplicate and proper pair filters. On paired-end data \code{bam.mapped} is therefore about twice \code{library.size}, less what the filters removed.
#'
#' \code{reads.in.regions} counts a fragment once per region it overlaps, the way the counting does. A region shared by several sets is counted once, but a fragment lying across two neighbouring regions, or two tiles of the same region, is counted in both, so on a tiled object the FRiP comes out slightly high.
#'
#' The FRiP is computed on \code{library.size}, the reads that went through the same filters as the counts. Computed on \code{bam.reads} it would mix fragments with alignment records and mapped with filtered reads.
#'
#' \code{discarded.reads} counts the fragments that passed every other filter and lie on the discarded regions, on the chromosomes entering the library sizes. They are in neither \code{library.size} nor \code{reads.in.regions}, and a sample losing a much larger share of its library than the others is worth a look.
#'
#' The FRiP of the input is the share of the input library falling in the regions, which is what the regions would collect with no enrichment at all. A sample whose FRiP sits close to the one of its input carries little signal in the regions, however many reads it has.
#'
#' @examples
#' sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#' peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)
#' counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)
#'
#' libInfo(counts, annotationColumns = "condition")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countReads}}, \code{\link{countTable}}
#'
#' @importFrom Rsamtools idxstatsBam
#' @importFrom SummarizedExperiment colData rowRanges assay assayNames
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom dplyr mutate bind_cols
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export libInfo

libInfo <-
  function(counts,
           annotationColumns = NULL,
           bamFiles = NULL) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (methods::is(counts, "RegionSetDE.fit")) {
      counts <- fitCounts(counts)
    }

    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts or a RegionSetDE.fit object.", call. = FALSE)
    }

    # Signal read from bigWig files has no reads to put a fraction on
    if (!is.null(counts@parameters$countBigwig)) {
      stop("The object was counted from bigWig files: it holds signal rather than reads, and has no library to summarise.", call. = FALSE)
    }

    sampleTable <- as.data.frame(SummarizedExperiment::colData(counts), optional = TRUE)

    missingColumns <- setdiff(annotationColumns, colnames(sampleTable))
    if (length(missingColumns) > 0) {
      stop("The following 'annotationColumns' are not in the colData: ", paste(missingColumns, collapse = ", "), ".", call. = FALSE)
    }

    if (is.null(bamFiles)) {
      bamFiles <- counts@parameters$countReads$bamFiles
    }

    if (!is.null(bamFiles) && length(bamFiles) != ncol(counts)) {
      stop("The number of BAM files does not match the number of samples of the counts object.", call. = FALSE)
    }

    #------------------------#
    # Reads in the BAM files #
    #------------------------#
    # The index already knows how many records every chromosome holds, the file itself is never read
    if (is.null(bamFiles)) {
      bamReads <- rep(NA_real_, ncol(counts))
      bamMapped <- rep(NA_real_, ncol(counts))
    } else {
      missingFiles <- bamFiles[!file.exists(bamFiles)]
      if (length(missingFiles) > 0) {
        stop("The following BAM files do not exist: ", paste(missingFiles, collapse = ", "), ".", call. = FALSE)
      }

      indexStats <- lapply(bamFiles, function(bamFile) {Rsamtools::idxstatsBam(.bamWithIndex(bamFile))})
      bamMapped <- vapply(indexStats, function(stats) {sum(as.numeric(stats$mapped))}, numeric(1))
      bamReads <- bamMapped + vapply(indexStats, function(stats) {sum(as.numeric(stats$unmapped))}, numeric(1))
    }

    #------------------------#
    # Reads in the regions   #
    #------------------------#
    # A region shared by several sets is one stretch of genome, it enters the sum once
    regionRanges <- SummarizedExperiment::rowRanges(counts)
    regionKey <- paste0(as.character(GenomeInfoDb::seqnames(regionRanges)), ":",
                        BiocGenerics::start(regionRanges), "-", BiocGenerics::end(regionRanges))

    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, "counts"))
    readsInRegions <- colSums(countMatrix[!duplicated(regionKey), , drop = FALSE])

    librarySize <- if ("library.size" %in% colnames(sampleTable)) {as.numeric(sampleTable$library.size)} else {rep(NA_real_, ncol(counts))}
    pairedEnd <- if ("paired.end" %in% colnames(sampleTable)) {as.logical(sampleTable$paired.end)} else {rep(NA, ncol(counts))}
    fragmentLength <- if ("fragment.length" %in% colnames(sampleTable)) {as.numeric(sampleTable$fragment.length)} else {rep(NA_real_, ncol(counts))}

    #------------------------#
    # Assemble the table     #
    #------------------------#
    infoTable <- data.frame(sample = colnames(counts),
                            paired.end = pairedEnd,
                            fragment.length = fragmentLength,
                            bam.reads = bamReads,
                            bam.mapped = bamMapped,
                            library.size = librarySize,
                            reads.in.regions = as.numeric(readsInRegions),
                            stringsAsFactors = FALSE)

    infoTable <- dplyr::mutate(infoTable, FRiP = round(.data$reads.in.regions / .data$library.size, 4))

    # What the blacklist, the greylist and the other discarded regions took out of each library
    if ("discarded.reads" %in% colnames(sampleTable)) {
      infoTable <- dplyr::mutate(infoTable, discarded.reads = as.numeric(sampleTable$discarded.reads), .after = "library.size")
    }

    #------------------------#
    # The inputs             #
    #------------------------#
    # The same fraction on the input tells how much of the FRiP the regions would collect with no enrichment
    if ("input" %in% SummarizedExperiment::assayNames(counts)) {
      inputMatrix <- as.matrix(SummarizedExperiment::assay(counts, "input"))

      infoTable <- dplyr::mutate(infoTable,
                                 input.id = if ("input.id" %in% colnames(sampleTable)) {as.character(sampleTable$input.id)} else {NA_character_},
                                 input.library.size = if ("input.library.size" %in% colnames(sampleTable)) {as.numeric(sampleTable$input.library.size)} else {NA_real_},
                                 input.in.regions = as.numeric(colSums(inputMatrix[!duplicated(regionKey), , drop = FALSE])),
                                 input.FRiP = round(.data$input.in.regions / .data$input.library.size, 4))
    }

    if (length(annotationColumns) > 0) {
      infoTable <- dplyr::bind_cols(infoTable[, "sample", drop = FALSE],
                                    sampleTable[, annotationColumns, drop = FALSE],
                                    infoTable[, setdiff(colnames(infoTable), "sample"), drop = FALSE])
    }

    rownames(infoTable) <- NULL
    return(infoTable)
  } # END function
