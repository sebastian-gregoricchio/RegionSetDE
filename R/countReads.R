# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title countReads
#'
#' @description Counts the reads of a group of BAM files over the regions of a \code{RegionSetDE} object. Paired-end data are counted as fragments, while single-end reads are extended to the fragment length before the overlap is evaluated. The regions can be recentred on the summit of the signal, as DiffBind does with the peaks of a consensus, or cut into tiles of fixed width, in which case each tile becomes a row of the resulting object. The input libraries of the samples, when there are any, are counted over the same rows and stored beside the counts.
#'
#' @param regionSet \code{RegionSetDE} object returned by \code{\link{loadRegions}} or \code{\link{loadConsensusPeaks}}, a named \code{GRangesList}, or a single \code{GRanges}, which is loaded by \code{\link{loadRegions}} as one set named \code{regions}.
#' @param bamFiles Character vector with the paths of the BAM files. Each file must be indexed, and all of them must be aligned to the same assembly, whose chromosomes they may name in different styles (\code{chr1} in some files and \code{1} in others). Default: \code{NULL}, taken from \code{sampleSheet}, or from the sample sheet a consensus was built from by \code{\link{loadConsensusPeaks}}.
#' @param sampleSheet Data.frame returned by \code{\link{loadSampleSheet}}, or the path to a sample sheet, providing the BAM files, the input files, the sample names and the annotation in one go. Default: \code{NULL}.
#' @param sampleNames Character vector with the sample names. Default: \code{NULL}, the BAM file names are used.
#' @param sampleMetadata Data.frame with the sample annotation, stored in the \code{colData}. When it contains a \code{sample} column the rows are matched by name, otherwise they must follow the order of \code{bamFiles}. Default: \code{NULL}.
#' @param keepMetadata Logical value to indicate whether the metadata columns carried by the regions must be kept in the \code{rowData}, harmonised across the sets. Default: \code{TRUE}.
#' @param regionId String with the name of a metadata column holding the region identifiers, for instance a gene name. It must hold a different value for every region of every set. Default: \code{NULL}, the names of the ranges, and their coordinates when they are unnamed.
#' @param tileWidth Numeric value with the width of the tiles, in base pairs. Cannot be combined with \code{summits}. Default: \code{NULL}, one row per region.
#' @param partialTiles Logical value: \code{TRUE} keeps the trailing tile of each region even when narrower than \code{tileWidth}, \code{FALSE} discards it together with the regions narrower than a single tile. Default: \code{TRUE}.
#' @param pairedEnd Logical value, one logical value per BAM file, or the string \code{"auto"} to read the layout from the files themselves. Default: \code{"auto"}.
#' @param fragmentLength Length to which single-end reads are extended. Either a number applied to every single-end sample; one number per BAM file, matched to the sample names when the vector is named; the name of a column of the sample sheet or of \code{sampleMetadata} holding one value per sample, such as the fragment length computed by phantompeakqualtools; or \code{"auto"}, to estimate it from the reads with \code{\link{estimateFragmentLength}}. Ignored for paired-end samples, whose fragments are rebuilt from the pairs. Default: \code{150}.
#' @param maxFragmentLength Numeric value with the maximum length accepted for a paired-end fragment. Applied to the paired-end samples only. Default: \code{1000}.
#' @param minMapq Numeric value with the minimum mapping quality of a read. Default: \code{20}.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded. Default: \code{TRUE}.
#' @param excludeChromosomes Character vector with the chromosomes left out of the library sizes, written in either naming style, \code{chrM} and \code{MT} both reaching the mitochondrial genome of a file naming it either way, for instance the mitochondrial genome, chrY or the unplaced and alternative contigs. A name matching no chromosome once converted raises a warning, since it would leave that chromosome inside the library sizes without a word. The regions lying on them are still counted: to leave those out as well, filter the regions when loading them. The contigs can be collected from the BAM header, e.g. \code{grep("_|EBV", names(Rsamtools::scanBamHeader(bamFile)[[1]]$targets), value = TRUE)}. Default: \code{NULL}, every chromosome enters the library sizes.
#' @param discardRegions \code{GRanges} with regions whose reads must be ignored, for instance a blacklist. A fragment is dropped when one of its reads starts inside them. Default: \code{NULL}.
#' @param fullLibrarySize Logical value: \code{TRUE} reads every chromosome that is not excluded, even those without any region, so that the library sizes cover the whole library; \code{FALSE} reads only the chromosomes carrying regions, which is much faster for a few regions but leaves library sizes that must not be used for normalisation. Default: \code{TRUE}.
#' @param inputFiles Character vector with the path of the input BAM file of every sample, \code{NA} for a sample without input. One input can serve several samples and is counted once. Default: \code{NULL}, the \code{input} column of the sample sheet or of \code{sampleMetadata}, when there is one.
#' @param countInput Logical value to indicate whether the input files must be counted. Default: \code{TRUE}.
#' @param summits Numeric value with the half width of the regions recentred on their summit, which then span \code{2 * summits + 1} bp, as with the \code{summits} argument of \code{DiffBind::dba.count}. \code{0} locates the summits and stores them without moving the regions. Default: \code{NULL}, the regions are counted as they are.
#' @param summitSource String with where the summits are taken from: \code{"reads"}, the highest point of the fragment pileup of the samples, or \code{"peaks"}, the summits written by the peak caller in the narrowPeak files of a consensus built by \code{\link{loadConsensusPeaks}}. Default: \code{"reads"}.
#' @param nThreads Number of threads. The files are cut into pieces of at most 50 Mb, shared among the threads, so even a single file benefits from several of them. Default: \code{1}.
#' @param progressBar Logical value to indicate whether a progress bar must be drawn while the files are read. It advances with the pieces of the files as the threads hand them back, and it is drawn only when \code{verbose = TRUE}. Default: \code{interactive()}, which keeps it out of scripts and rendered documents.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.counts} object with one row per region, or per tile, and one column per sample. The library sizes are stored in the \code{library.size} column of the \code{colData}, the set membership in the \code{region.set} column of the \code{rowData}, and the length the single-end reads were extended to in \code{fragment.length} (\code{NA} for paired-end samples). With inputs, the \code{input} assay holds for every sample the counts of its input over the same rows (\code{NA} for a sample without one), and the \code{colData} gains \code{input.id} and \code{input.library.size}. With \code{summits}, the \code{rowData} gains \code{summit}, the position of the summit of every region.
#'
#' @details Regions shared by several sets are counted only once and the values are then copied to every set they belong to, which keeps the running time proportional to the number of distinct regions.
#'
#' A paired-end fragment is counted in every region it overlaps, including the regions it spans with both reads outside them. The fragment is rebuilt from the first mate of each proper pair, whose position and template length (TLEN) give its start and width, so the two reads never have to be matched in memory. The pairs therefore have to be flagged as proper by the aligner, and those longer than \code{maxFragmentLength} are dropped. The mapping quality of the second mate is read from the \code{MQ} tag, which \code{samtools fixmate} and Picard write; on files without it only the first mate is checked, and a message says so.
#'
#' Counts and library sizes are kept apart. The counts only need the chromosomes carrying regions, and every region is counted, wherever it lies. The library size of a sample is the number of fragments that went through the same filters as the counts, on every chromosome of the BAM files except those in \code{excludeChromosomes}. Leaving out the mitochondrial genome matters in ATAC-seq, where its share of the reads changes from sample to sample. With \code{fullLibrarySize = FALSE} only the chromosomes carrying regions are read, and the library sizes are partial: the counts do not change, but the library sizes are not usable for normalisation, and \code{\link{normalizeCounts}} warns when a method relies on them.
#'
#' Paired-end and single-end samples can be mixed in the same call, paired-end libraries being counted as fragments and single-end ones as reads extended to their fragment length, so that both end up with one count per sequenced fragment. Forcing a paired-end file through the single-end path counts each mate on its own and nearly doubles its values, while the opposite mistake finds no pair and returns a column of zeros, which is why the layout is read from the files by default. The resolved layout of each sample is stored in the \code{paired.end} column of the \code{colData}.
#'
#' Single-end libraries rarely share the same fragment length, and a read extended too far spills into the neighbouring regions while one extended too little misses the centre of its own. The length can come from the pipeline that produced the files, as a column of the sample sheet, or from \code{fragmentLength = "auto"}, which runs the strand cross-correlation of \code{\link{estimateFragmentLength}} over the regions being counted. \code{\link{countBackground}} and \code{\link{countGreenlist}} reuse the lengths chosen here, sample by sample.
#'
#' The inputs are counted with the same filters as the samples, over the same rows, and a single-end input is extended to the mean fragment length of the samples it serves. Nothing downstream subtracts them: the counts of a sample stay the reads of that sample, as the count models need. The input assay is there to check the enrichment of the regions, and to see whether a change between conditions also shows in the inputs, which points to copy number rather than to binding. \code{\link{libInfo}} reports the inputs beside the libraries.
#'
#' A consensus of peaks has regions of every width, and a wide region collects more background than a narrow one around the same summit. With \code{summits}, every region is replaced by a window of fixed width centred on its summit, which is what DiffBind does by default. With \code{summitSource = "reads"} the summit is the middle of the highest stretch of fragment pileup in each sample, averaged over the samples with weights proportional to the height of their pileup, scaled by their depth, so that the samples carrying signal decide where it sits. With \code{summitSource = "peaks"} it is the average of the summits the peak caller wrote for the peaks overlapping the region, weighted by their significance, which needs no BAM file but works only for narrowPeak files. A region with neither reads nor peaks keeps its midpoint. The regions keep their identifiers, so each window can be traced back to the region it came from, and windows of neighbouring regions may overlap, as in DiffBind.
#'
#' Regions and BAM files do not need to share the same chromosome naming style, and neither do the BAM files among themselves. The names of the first file are the reference. The regions, \code{discardRegions} and \code{excludeChromosomes} are brought to them chromosome by chromosome, for the counting only, and every other file is read under its own names. UCSC regions can then be counted on Ensembl alignments, or on a mix of the two, and the object still comes back with the names of the input sets. The same holds for the inputs. What the files cannot differ in is the assembly: two files giving different lengths to the same chromosome are refused. Contigs that some files lack under any name, as scaffolds and decoys often do between two builds of one assembly, hold no read in those files and enter the library sizes of the others, unless they are listed in \code{excludeChromosomes}.
#'
#' A single \code{GRanges} is taken as one set of regions. It goes through \code{\link{loadRegions}} with its order and its chromosome names kept, the regions with identical coordinates collapsed into one, and the set is called \code{regions}. To give the set another name, to sort it or to split it into several sets, call \code{\link{loadRegions}} or \code{\link{splitLoadRegions}} first.
#'
#' @examples
#' # The peaks of one sample of the AR example stand in for a region set
#' sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#' peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)
#'
#' # The sheet brings the BAM files, the inputs, the sample names and the annotation
#' counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)
#' counts
#' head(SummarizedExperiment::assay(counts, "input"), 3)
#'
#' # The same files given one by one, with the annotation as a table and no input
#' counts <- countReads(peakRegions,
#'                      bamFiles = sampleSheet$bam,
#'                      sampleNames = sampleSheet$sample,
#'                      sampleMetadata = sampleSheet[, c("sample", "condition")],
#'                      verbose = FALSE)
#'
#' # Windows of 401 bp centred on the summit of the reads, as DiffBind counts a consensus
#' summitCounts <- countReads(peakRegions, sampleSheet = sampleSheet, summits = 200, verbose = FALSE)
#' head(SummarizedExperiment::rowRanges(summitCounts), 3)
#'
#' # Tiles of 100 bp, one row each
#' tiledCounts <- countReads(peakRegions, sampleSheet = sampleSheet, tileWidth = 100, verbose = FALSE)
#' head(SummarizedExperiment::rowData(tiledCounts), 3)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{estimateFragmentLength}}, \code{\link{countBigwig}}, \code{\link{loadCounts}}, \code{\link{countBackground}}, \code{\link{libInfo}}
#'
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom IRanges ranges ranges<-
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom SummarizedExperiment assay<-
#' @importFrom dplyr mutate n_distinct
#' @importFrom methods is
#' @importFrom utils head
#'
#' @export countReads

countReads <-
  function(regionSet,
           bamFiles = NULL,
           sampleSheet = NULL,
           sampleNames = NULL,
           sampleMetadata = NULL,
           tileWidth = NULL,
           keepMetadata = TRUE,
           regionId = NULL,
           partialTiles = TRUE,
           pairedEnd = "auto",
           fragmentLength = 150,
           maxFragmentLength = 1000,
           minMapq = 20,
           removeDuplicates = TRUE,
           excludeChromosomes = NULL,
           discardRegions = NULL,
           fullLibrarySize = TRUE,
           inputFiles = NULL,
           countInput = TRUE,
           summits = NULL,
           summitSource = "reads",
           nThreads = 1,
           progressBar = interactive(),
           verbose = TRUE) {

    startTime <- Sys.time()

    # The bar belongs to the messages: a call asked to be silent draws nothing
    showProgress <- isTRUE(verbose) & isTRUE(progressBar)

    #------------------------#
    # Check of the arguments #
    #------------------------#
    # A single GRanges becomes a RegionSetDE with one set, the object everything below is written for
    regionSet <- .asRegionSet(regionSet = regionSet, verbose = verbose)

    # A sample sheet brings the files, the names and the annotation in one go
    sheetInput <- .sheetCountingInput(regionSet = regionSet,
                                      sampleSheet = sampleSheet,
                                      files = bamFiles,
                                      sampleNames = sampleNames,
                                      sampleMetadata = sampleMetadata,
                                      fileField = "bam",
                                      verbose = verbose)

    bamFiles <- sheetInput$files
    sampleNames <- sheetInput$sampleNames
    sampleMetadata <- sheetInput$sampleMetadata

    if (!is.character(bamFiles) | length(bamFiles) == 0) {
      stop("The 'bamFiles' parameter must be a character vector with at least one BAM file.", call. = FALSE)
    }

    missingFiles <- bamFiles[!file.exists(bamFiles)]
    if (length(missingFiles) > 0) {
      stop("The following BAM files do not exist: ", paste(missingFiles, collapse = ", "), ".", call. = FALSE)
    }

    # Every chromosome is reached through the index, a file without one cannot be read by pieces
    hasIndex <- .hasBamIndex(bamFiles)

    if (any(!hasIndex)) {
      stop("The following BAM files are not indexed: ", paste(basename(bamFiles[!hasIndex]), collapse = ", "), ".", call. = FALSE)
    }

    if (!(is.numeric(fragmentLength) | (is.character(fragmentLength) & length(fragmentLength) == 1))) {
      stop("The 'fragmentLength' parameter must be a number, one number per BAM file, the name of a column of the sample annotation, or 'auto'.", call. = FALSE)
    }

    if (!is.numeric(maxFragmentLength) | maxFragmentLength[1] < 1) {
      stop("The 'maxFragmentLength' parameter must be a positive number.", call. = FALSE)
    }

    if (!is.logical(fullLibrarySize) | length(fullLibrarySize) != 1 | any(is.na(fullLibrarySize))) {
      stop("The 'fullLibrarySize' parameter must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.logical(countInput) | length(countInput) != 1 | any(is.na(countInput))) {
      stop("The 'countInput' parameter must be TRUE or FALSE.", call. = FALSE)
    }

    if (!is.null(summits)) {
      if (!is.numeric(summits) | length(summits) != 1 || is.na(summits) || summits < 0) {
        stop("The 'summits' parameter must be NULL, 0, or a positive number of base pairs.", call. = FALSE)
      }

      if (!is.null(tileWidth)) {
        stop("The 'summits' and 'tileWidth' parameters cannot be used together: recentre the regions or tile them.", call. = FALSE)
      }

      summitSource <- tolower(as.character(summitSource[1]))
      if (!(summitSource %in% c("reads", "peaks"))) {
        stop("The 'summitSource' parameter must be either 'reads' or 'peaks'.", call. = FALSE)
      }

      if (summitSource == "peaks" && !(methods::is(regionSet, "RegionSetDE") && length(regionSet@consensus) > 0)) {
        stop("summitSource = 'peaks' needs the regions returned by loadConsensusPeaks(), whose peak calls carry the summits. Use summitSource = 'reads' for any other region set.", call. = FALSE)
      }
    }

    # The flags of the first records are enough to tell a paired library from a single-end one
    if (identical(pairedEnd, "auto")) {
      pairedEnd <- .bamIsPairedEnd(bamFiles = bamFiles)
    }

    if (length(pairedEnd) == 1) {pairedEnd <- rep(pairedEnd, length(bamFiles))}

    if (!is.logical(pairedEnd) | length(pairedEnd) != length(bamFiles) | any(is.na(pairedEnd))) {
      stop("The 'pairedEnd' parameter must be 'auto', a single logical value, or one logical value per BAM file.", call. = FALSE)
    }

    if (!is.null(discardRegions) & !methods::is(discardRegions, "GRanges")) {
      stop("The 'discardRegions' parameter must be a GRanges object.", call. = FALSE)
    }

    if (!is.null(excludeChromosomes) & !is.character(excludeChromosomes)) {
      stop("The 'excludeChromosomes' parameter must be a character vector with chromosome names.", call. = FALSE)
    }

    # The files may name the chromosomes in different styles: those of the first file stand for all of them
    bamTargets <- .bamChromosomeMap(bamFiles = bamFiles, verbose = verbose)$lengths
    bamSeqlevels <- names(bamTargets)

    # 'chrM' against a BAM naming it 'MT' would exclude nothing, and the library sizes would carry
    # the mitochondrial reads into the normalisation without a word
    excludeChromosomes <- .matchChromosomeNames(chromosomeNames = excludeChromosomes, targetSeqlevels = bamSeqlevels)

    #------------------------#
    # Samples and regions    #
    #------------------------#
    sampleTable <- .buildSampleTable(files = bamFiles,
                                     sampleNames = sampleNames,
                                     sampleMetadata = sampleMetadata,
                                     fileColumn = "bam.file",
                                     extensionPattern = "\\.bam$")

    allRegions <- .flattenRegionSets(regionSet = regionSet,
                                     tileWidth = tileWidth,
                                     keepMetadata = keepMetadata,
                                     regionId = regionId,
                                     partialTiles = partialTiles,
                                     verbose = verbose)

    # Overlapping sets share regions, counting them once and copying the values back is much cheaper.
    # The strand is left out of the key because the counting ignores it.
    regionKey <- paste0(as.character(GenomeInfoDb::seqnames(allRegions)), ":",
                        BiocGenerics::start(allRegions), "-",
                        BiocGenerics::end(allRegions))

    uniqueRegions <- allRegions[!duplicated(regionKey)]
    expansionIndex <- match(regionKey, regionKey[!duplicated(regionKey)])

    #---------------------------------#
    # Chromosome names of the samples #
    #---------------------------------#
    # Only the copies used for the counting are renamed, the returned object keeps the style of the region sets
    countingRegions <- .matchSeqlevels(x = uniqueRegions, targetSeqlevels = bamSeqlevels, fileName = bamFiles[1], verbose = verbose)

    if (!is.null(discardRegions)) {
      discardRegions <- .matchSeqlevels(x = discardRegions, targetSeqlevels = bamSeqlevels, fileName = bamFiles[1], verbose = FALSE)
    }

    # A region on a chromosome missing from the files keeps a count of zero, better to say it than to let it pass for an empty region
    regionChromosomes <- as.character(GenomeInfoDb::seqnames(countingRegions))
    absentRegions <- sum(!(regionChromosomes %in% bamSeqlevels))
    excludedRegions <- sum(regionChromosomes %in% excludeChromosomes)

    if (isTRUE(verbose) & absentRegions > 0) {
      message(absentRegions, " regions lie on chromosomes absent from the BAM files and get zero counts.")
    }

    if (isTRUE(verbose) & excludedRegions > 0) {
      message(excludedRegions, " regions lie on chromosomes listed in 'excludeChromosomes': they are counted, but their chromosomes stay out of the library sizes.")
    }

    #------------------------#
    # Fragment lengths       #
    #------------------------#
    # One length per sample, NA for the paired-end ones whose fragments come from the pairs
    fragmentInfo <- .resolveFragmentLength(fragmentLength = fragmentLength,
                                           sampleTable = sampleTable,
                                           pairedEnd = pairedEnd,
                                           bamFiles = bamFiles,
                                           regions = countingRegions,
                                           minMapq = minMapq,
                                           removeDuplicates = removeDuplicates,
                                           discardRegions = discardRegions,
                                           nThreads = nThreads,
                                           verbose = verbose)

    sampleFragmentLength <- fragmentInfo$values

    #------------------------#
    # Summits                #
    #------------------------#
    # The regions are moved before any read is counted, so that samples and inputs see the same windows
    if (!is.null(summits)) {
      if (isTRUE(verbose) & summitSource == "reads") {
        message("Locating the summits of ", length(uniqueRegions), " regions on the pileup of ", length(bamFiles), " samples...")
      }

      summitPosition <- if (summitSource == "reads") {
        .readSummits(bamFiles = bamFiles,
                     regions = countingRegions,
                     pairedEnd = pairedEnd,
                     fragmentLength = sampleFragmentLength,
                     maxFragmentLength = maxFragmentLength[1],
                     minMapq = minMapq,
                     removeDuplicates = removeDuplicates,
                     discardRegions = discardRegions,
                     progressBar = showProgress,
                     nThreads = nThreads)
      } else {
        .peakSummits(regions = uniqueRegions, peakList = regionSet@consensus$peaks)
      }

      withoutSummit <- is.na(summitPosition)
      summitPosition[withoutSummit] <- as.integer((BiocGenerics::start(uniqueRegions)[withoutSummit] + BiocGenerics::end(uniqueRegions)[withoutSummit]) %/% 2L)

      if (summits > 0) {
        recentredRanges <- .recentreRanges(summitPosition = summitPosition,
                                           halfWidth = summits,
                                           chromosomeLength = bamTargets[regionChromosomes])

        # The names live in the ranges, they are put back so that the rows keep their identifiers
        recentredRows <- recentredRanges[expansionIndex]
        names(recentredRows) <- names(allRegions)
        names(recentredRanges) <- names(uniqueRegions)

        IRanges::ranges(countingRegions) <- recentredRanges
        IRanges::ranges(uniqueRegions) <- recentredRanges
        IRanges::ranges(allRegions) <- recentredRows
      }

      S4Vectors::mcols(allRegions)$summit <- summitPosition[expansionIndex]

      if (isTRUE(verbose) & summits > 0) {
        message("Regions recentred on the summit of the ", if (summitSource == "reads") {"reads"} else {"peak calls"}, ", ", 2 * summits + 1, " bp wide.")
      }

      if (isTRUE(verbose) & summits == 0) {
        message("Summits located, the regions keep their coordinates.")
      }

      if (isTRUE(verbose) & any(withoutSummit)) {
        message(sum(withoutSummit), " regions without ", summitSource, " were centred on their midpoint.")
      }
    }

    #----------------#
    # Read the files #
    #----------------#
    if (isTRUE(verbose)) {
      message("Counting reads in ", length(bamFiles), " samples over ", length(uniqueRegions), " unique regions (",
              length(allRegions), " rows, ", dplyr::n_distinct(S4Vectors::mcols(allRegions)$region.set), " sets)...")
    }

    if (isTRUE(verbose) & length(unique(pairedEnd)) > 1) {
      message("Mixed layouts: ", sum(pairedEnd), " paired-end and ", sum(!pairedEnd), " single-end samples.")
    }

    fragmentCounts <- .countBamFragments(bamFiles = bamFiles,
                                         ranges = countingRegions,
                                         pairedEnd = pairedEnd,
                                         fragmentLength = sampleFragmentLength,
                                         maxFragmentLength = maxFragmentLength[1],
                                         minMapq = minMapq,
                                         removeDuplicates = removeDuplicates,
                                         excludeChromosomes = excludeChromosomes,
                                         discardRegions = discardRegions,
                                         fullLibrarySize = fullLibrarySize,
                                         countMode = "overlap",
                                         progressBar = showProgress,
                                         nThreads = nThreads)

    countMatrix <- fragmentCounts$counts
    storage.mode(countMatrix) <- "double"

    # Without the MQ tag the quality of the second mate is unknown, the pair then rests on the first mate alone
    missingMateMapq <- which(pairedEnd & !fragmentCounts$mate.mapq.found)

    if (isTRUE(verbose) & length(missingMateMapq) > 0 & isTRUE(minMapq > 0)) {
      message("No MQ tag in: ", paste(sampleTable$sample[missingMateMapq], collapse = ", "),
              ". For these samples 'minMapq' is applied to the first mate of each pair only.")
    }

    # An empty library breaks every normalisation downstream, and usually points to a wrong layout or to excluding too much
    emptyLibraries <- which(fragmentCounts$library.size == 0)
    if (length(emptyLibraries) > 0) {
      warning("No fragment entered the library size of: ", paste(sampleTable$sample[emptyLibraries], collapse = ", "),
              ". Check 'pairedEnd', 'excludeChromosomes' and the read filters.", call. = FALSE)
    }

    #------------------------#
    # Inputs                 #
    #------------------------#
    inputPaths <- if (isTRUE(countInput)) {.resolveInputFiles(inputFiles = inputFiles, sampleTable = sampleTable)} else {NULL}
    inputCounts <- NULL

    # The names of the sheet belong to its own inputs, files given here are named after themselves
    if (!is.null(inputPaths)) {
      inputCounts <- .countInputFiles(inputPaths = inputPaths,
                                      inputIds = if (is.null(inputFiles) & "input.id" %in% colnames(sampleTable)) {as.character(sampleTable$input.id)} else {NULL},
                                      sampleFragmentLength = sampleFragmentLength,
                                      fragmentLength = fragmentLength,
                                      bamTargets = bamTargets,
                                      regions = countingRegions,
                                      maxFragmentLength = maxFragmentLength[1],
                                      minMapq = minMapq,
                                      removeDuplicates = removeDuplicates,
                                      excludeChromosomes = excludeChromosomes,
                                      discardRegions = discardRegions,
                                      fullLibrarySize = fullLibrarySize,
                                      nThreads = nThreads,
                                      progressBar = showProgress,
                                      verbose = verbose)
    }

    #-------------------------#
    # Assemble the object     #
    #-------------------------#
    # The library sizes are the fragments surviving the filters, which is the denominator the normalisation expects
    sampleTable <- dplyr::mutate(sampleTable,
                                 paired.end = pairedEnd,
                                 fragment.length = sampleFragmentLength,
                                 library.size = fragmentCounts$library.size)

    if (!is.null(inputCounts)) {
      sampleTable <- dplyr::mutate(sampleTable,
                                   input = inputPaths,
                                   input.id = inputCounts$input.id,
                                   input.library.size = inputCounts$library.size)
    }

    countMatrix <- countMatrix[expansionIndex, , drop = FALSE]

    newParameters <- list(countReads = list(bamFiles = bamFiles,
                                            tileWidth = tileWidth,
                                            partialTiles = partialTiles,
                                            pairedEnd = pairedEnd,
                                            fragmentLength = sampleFragmentLength,
                                            fragmentLengthSource = fragmentInfo$source,
                                            maxFragmentLength = maxFragmentLength,
                                            minMapq = minMapq,
                                            removeDuplicates = removeDuplicates,
                                            excludeChromosomes = excludeChromosomes,
                                            fullLibrarySize = fullLibrarySize,
                                            inputFiles = inputPaths,
                                            summits = summits,
                                            summitSource = if (is.null(summits)) {NULL} else {summitSource}))

    # A tiled object has to say so it is tiled, otherwise testRegions treats every tile as a region
    # and the combination step that puts the tiles back together never runs
    counts <- .newCountsObject(countMatrix = countMatrix,
                               regions = allRegions,
                               sampleTable = sampleTable,
                               provenance = .provenanceSlots(regionSet),
                               countingLevel = if (is.null(tileWidth)) {"region"} else {"tile"},
                               newParameters = newParameters,
                               metadataList = list(signal.type = "reads", count.like = TRUE))

    # The input of every sample over the same rows, a sample without input keeps a column of NA
    if (!is.null(inputCounts)) {
      inputMatrix <- inputCounts$counts[expansionIndex, , drop = FALSE]
      dimnames(inputMatrix) <- dimnames(countMatrix)
      SummarizedExperiment::assay(counts, "input", withDimnames = FALSE) <- inputMatrix
    }

    if (isTRUE(verbose)) {
      message("Done in ", .elapsedTime(startTime), ". Library sizes: ",
              paste(round(range(sampleTable$library.size) / 1e6, 1), collapse = " - "), " million fragments.")

      if (isFALSE(fullLibrarySize)) {
        message("The library sizes cover only the chromosomes carrying regions: count again with 'fullLibrarySize = TRUE' before a normalisation that relies on them.")
      }
    }

    return(counts)
  } # END function




#' @title .resolveFragmentLength
#'
#' @description Works out the length every single-end sample is extended to, from a single value, one value per sample, a column of the sample annotation, or an estimate made on the reads.
#'
#' @param fragmentLength Value given to \code{countReads}.
#' @param sampleTable Data.frame with one row per sample, as built by \code{.buildSampleTable}.
#' @param pairedEnd Logical vector with one value per sample.
#' @param bamFiles Character vector with the paths of the BAM files.
#' @param regions \code{GRanges} with the regions, named after the chromosomes of the BAM files, used by the estimate.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param discardRegions \code{GRanges} with the regions whose reads must be ignored, or \code{NULL}.
#' @param nThreads Number of threads.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A list with \code{values}, a numeric vector with one length per sample, \code{NA} for the paired-end ones, and \code{source}, a string saying where the lengths came from.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resolveFragmentLength <-
  function(fragmentLength,
           sampleTable,
           pairedEnd,
           bamFiles,
           regions,
           minMapq,
           removeDuplicates,
           discardRegions,
           nThreads,
           verbose) {

    sampleNumber <- nrow(sampleTable)
    singleEnd <- !pairedEnd

    #-------------------------------#
    # Where the values come from    #
    #-------------------------------#
    if (is.numeric(fragmentLength)) {
      if (length(fragmentLength) == 1) {
        lengthValues <- rep(as.numeric(fragmentLength), sampleNumber)
        lengthSource <- "fixed"
      } else {
        lengthValues <- as.numeric(.orderSampleValues(values = fragmentLength, sampleNames = sampleTable$sample, parameterName = "fragmentLength"))
        lengthSource <- "per sample"
      }
    } else if (identical(fragmentLength, "auto")) {
      lengthValues <- rep(NA_real_, sampleNumber)
      lengthSource <- "cross-correlation"

      if (any(singleEnd)) {
        if (isTRUE(verbose)) {
          message("Estimating the fragment length of ", sum(singleEnd), " single-end samples from the reads over the regions...")
        }

        fragmentEstimate <- estimateFragmentLength(bamFiles = bamFiles[singleEnd],
                                                   regions = regions,
                                                   sampleNames = sampleTable$sample[singleEnd],
                                                   pairedEnd = FALSE,
                                                   minMapq = minMapq,
                                                   removeDuplicates = removeDuplicates,
                                                   discardRegions = discardRegions,
                                                   nThreads = nThreads,
                                                   verbose = FALSE)

        lengthValues[singleEnd] <- fragmentEstimate$table$fragment.length
      }
    } else {
      if (!(fragmentLength %in% colnames(sampleTable))) {
        stop("The 'fragmentLength' parameter names no column of the sample annotation: '", fragmentLength,
             "'. Give a number, a column of the sample sheet, or 'auto'.", call. = FALSE)
      }

      lengthValues <- suppressWarnings(as.numeric(as.character(sampleTable[[fragmentLength]])))
      lengthSource <- paste("column", fragmentLength)
    }

    #-------------------------------#
    # Checks on the single-end ones #
    #-------------------------------#
    lengthValues[pairedEnd] <- NA_real_
    invalidLengths <- singleEnd & (is.na(lengthValues) | lengthValues < 1)

    if (any(invalidLengths)) {
      stop("No valid fragment length for the single-end samples: ", paste(sampleTable$sample[invalidLengths], collapse = ", "),
           ". Give a positive number for each of them.", call. = FALSE)
    }

    lengthValues <- round(lengthValues)

    if (isTRUE(verbose) & any(singleEnd) & lengthSource != "fixed") {
      message("Single-end reads extended to (", lengthSource, "): ",
              paste(sampleTable$sample[singleEnd], lengthValues[singleEnd], sep = " ", collapse = ", "), " bp.")
    }

    return(list(values = lengthValues, source = lengthSource))
  } # END function




#' @title .resolveInputFiles
#'
#' @description Returns the input file of every sample, given directly or taken from the \code{input} column of the sample annotation.
#'
#' @param inputFiles Character vector given to \code{countReads}, or \code{NULL}.
#' @param sampleTable Data.frame with one row per sample.
#'
#' @return A character vector with one path per sample, \code{NA} where there is no input, or \code{NULL} when no sample has one.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resolveInputFiles <-
  function(inputFiles,
           sampleTable) {

    if (is.null(inputFiles)) {
      if (!("input" %in% colnames(sampleTable))) {return(NULL)}
      inputFiles <- as.character(sampleTable$input)
    }

    if (!is.character(inputFiles) & !all(is.na(inputFiles))) {
      stop("The 'inputFiles' parameter must be a character vector with one path per sample, NA for the samples without input.", call. = FALSE)
    }

    inputFiles <- as.character(inputFiles)
    if (length(inputFiles) == 1) {inputFiles <- rep(inputFiles, nrow(sampleTable))}

    if (length(inputFiles) != nrow(sampleTable)) {
      stop("The 'inputFiles' parameter must hold one path per sample, NA for the samples without input.", call. = FALSE)
    }

    inputFiles[!is.na(inputFiles) & inputFiles == ""] <- NA_character_

    if (all(is.na(inputFiles))) {return(NULL)}

    return(inputFiles)
  } # END function




#' @title .countInputFiles
#'
#' @description Counts every distinct input file once over the regions, with the read filters of the samples, and spreads the counts over the samples they serve.
#'
#' @param inputPaths Character vector with the input of every sample, \code{NA} where there is none.
#' @param inputIds Character vector with the name of the input of every sample, or \code{NULL} to name them after the files.
#' @param sampleFragmentLength Numeric vector with the fragment length of every sample, \code{NA} for the paired-end ones.
#' @param fragmentLength Value given to \code{countReads}, used to tell a single length from lengths varying by sample.
#' @param bamTargets Named vector with the chromosome lengths of the sample BAM files, whose names the regions are written in.
#' @param regions \code{GRanges} with the regions, named after the chromosomes of the BAM files.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param excludeChromosomes Character vector with the chromosomes left out of the library sizes.
#' @param discardRegions \code{GRanges} with the regions whose reads must be ignored, or \code{NULL}.
#' @param fullLibrarySize Logical value indicating whether the library sizes cover every chromosome.
#' @param nThreads Number of threads.
#' @param progressBar Logical value to indicate whether a progress bar must be drawn while the inputs are read. Default: \code{FALSE}.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A list with \code{counts}, a matrix with one column per sample, \code{NA} for the samples without input; \code{library.size}, the library size of the input of every sample; and \code{input.id}, its name.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.countInputFiles <-
  function(inputPaths,
           inputIds,
           sampleFragmentLength,
           fragmentLength,
           bamTargets,
           regions,
           maxFragmentLength,
           minMapq,
           removeDuplicates,
           excludeChromosomes,
           discardRegions,
           fullLibrarySize,
           nThreads,
           progressBar = FALSE,
           verbose) {

    distinctInputs <- unique(inputPaths[!is.na(inputPaths)])

    #-------------------------------#
    # Files usable at all           #
    #-------------------------------#
    missingInputs <- distinctInputs[!file.exists(distinctInputs)]
    if (length(missingInputs) > 0) {
      stop("The following input files do not exist: ", paste(missingInputs, collapse = ", "), ". Set countInput = FALSE to skip them.", call. = FALSE)
    }

    unindexedInputs <- distinctInputs[!.hasBamIndex(distinctInputs)]
    if (length(unindexedInputs) > 0) {
      stop("The following input files are not indexed: ", paste(basename(unindexedInputs), collapse = ", "), ".", call. = FALSE)
    }

    # Inputs aligned against another assembly would place the regions somewhere else, and are refused here.
    # Inputs naming the chromosomes in another style than the samples are fine, the counting reads them under their own names
    .bamChromosomeMap(bamFiles = distinctInputs, referenceTargets = bamTargets, referenceLabel = "the sample BAM files")

    #-------------------------------#
    # Layout and fragment length    #
    #-------------------------------#
    inputPairedEnd <- .bamIsPairedEnd(bamFiles = distinctInputs)

    # A single-end input is extended as the samples it serves, which are the fragments it has to be compared with
    inputFragmentLength <- vapply(distinctInputs,
                                  function(inputFile) {
                                    if (is.numeric(fragmentLength) & length(fragmentLength) == 1) {return(as.numeric(fragmentLength))}
                                    servedLengths <- sampleFragmentLength[!is.na(inputPaths) & inputPaths == inputFile]
                                    if (all(is.na(servedLengths))) {return(NA_real_)}
                                    return(round(mean(servedLengths, na.rm = TRUE)))
                                  },
                                  numeric(1), USE.NAMES = FALSE)

    # Only paired-end samples behind a single-end input: the input has to speak for itself
    needEstimate <- !inputPairedEnd & is.na(inputFragmentLength)
    if (any(needEstimate)) {
      inputEstimate <- estimateFragmentLength(bamFiles = distinctInputs[needEstimate],
                                              regions = regions,
                                              pairedEnd = FALSE,
                                              minMapq = minMapq,
                                              removeDuplicates = removeDuplicates,
                                              discardRegions = discardRegions,
                                              verbose = FALSE)
      inputFragmentLength[needEstimate] <- inputEstimate$table$fragment.length
      inputFragmentLength[needEstimate & is.na(inputFragmentLength)] <- 150
    }

    if (isTRUE(verbose)) {
      message("Counting ", length(distinctInputs), " input file(s) serving ", sum(!is.na(inputPaths)), " samples...")
    }

    #-------------------------------#
    # Count and spread              #
    #-------------------------------#
    inputFragments <- .countBamFragments(bamFiles = distinctInputs,
                                         ranges = regions,
                                         pairedEnd = inputPairedEnd,
                                         fragmentLength = inputFragmentLength,
                                         maxFragmentLength = maxFragmentLength,
                                         minMapq = minMapq,
                                         removeDuplicates = removeDuplicates,
                                         excludeChromosomes = excludeChromosomes,
                                         discardRegions = discardRegions,
                                         fullLibrarySize = fullLibrarySize,
                                         countMode = "overlap",
                                         referenceTargets = bamTargets,
                                         progressBar = progressBar,
                                         nThreads = nThreads)

    inputIndex <- match(inputPaths, distinctInputs)

    inputMatrix <- inputFragments$counts[, inputIndex, drop = FALSE]
    storage.mode(inputMatrix) <- "double"
    inputMatrix[, is.na(inputIndex)] <- NA_real_

    if (is.null(inputIds) || all(is.na(inputIds))) {
      inputIds <- .inputIdentifiers(inputPaths = inputPaths)
    }
    inputIds[is.na(inputPaths)] <- NA_character_

    return(list(counts = inputMatrix,
                library.size = inputFragments$library.size[inputIndex],
                input.id = inputIds))
  } # END function




#' @title .readSummits
#'
#' @description Places the summit of every region from the fragment pileup of the samples: the summit of each sample, averaged with weights proportional to the height of its pileup over its depth.
#'
#' @param bamFiles Character vector with the paths of the BAM files.
#' @param regions \code{GRanges} with the regions, named after the chromosomes of the BAM files.
#' @param pairedEnd Logical vector with one value per BAM file.
#' @param fragmentLength Numeric vector with the fragment length of every sample, \code{NA} for the paired-end ones.
#' @param maxFragmentLength Numeric value with the maximum length of a paired-end fragment.
#' @param minMapq Numeric value with the minimum mapping quality of a read.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded.
#' @param discardRegions \code{GRanges} with the regions whose reads must be ignored, or \code{NULL}.
#' @param progressBar Logical value to indicate whether a progress bar must be drawn, one step for every file. Default: \code{FALSE}.
#' @param nThreads Number of threads.
#'
#' @return An integer vector with the summit of every region, \code{NA} where no sample has a fragment.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.readSummits <-
  function(bamFiles,
           regions,
           pairedEnd,
           fragmentLength,
           maxFragmentLength,
           minMapq,
           removeDuplicates,
           discardRegions,
           progressBar = FALSE,
           nThreads) {

    summitData <- .bamSummits(bamFiles = bamFiles,
                              regions = regions,
                              pairedEnd = pairedEnd,
                              fragmentLength = fragmentLength,
                              maxFragmentLength = maxFragmentLength,
                              minMapq = minMapq,
                              removeDuplicates = removeDuplicates,
                              discardRegions = discardRegions,
                              progressBar = progressBar,
                              nThreads = nThreads)

    # Heights over depth, so that a deeper library does not decide on its own where the summit sits
    sampleDepth <- pmax(summitData$fragments, 1)
    summitWeight <- sweep(summitData$height, MARGIN = 2, STATS = sampleDepth, FUN = "/")
    summitWeight[is.na(summitData$position)] <- 0

    weightedPosition <- summitData$position
    weightedPosition[is.na(weightedPosition)] <- 0

    weightTotal <- rowSums(summitWeight)
    summitPosition <- as.integer(round(rowSums(summitWeight * weightedPosition) / weightTotal))
    summitPosition[weightTotal == 0] <- NA_integer_

    return(summitPosition)
  } # END function




#' @title .peakSummits
#'
#' @description Places the summit of every region from the peak calls overlapping it: the summits written in the narrowPeak files, averaged with weights proportional to their significance.
#'
#' @param regions \code{GRanges} with the regions.
#' @param peakList \code{GRangesList} with the peaks of every sample, as kept by \code{\link{loadConsensusPeaks}}.
#'
#' @return An integer vector with the summit of every region, \code{NA} where no peak overlaps it.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom BiocGenerics start end width
#' @importFrom dplyr filter group_by summarise
#' @importFrom rlang .data
#'
#' @keywords internal

.peakSummits <-
  function(regions,
           peakList) {

    allPeaks <- unlist(peakList, use.names = FALSE)
    peakMetadata <- S4Vectors::mcols(allPeaks)

    if (!("peak" %in% colnames(peakMetadata)) || all(is.na(peakMetadata$peak) | peakMetadata$peak < 0)) {
      warning("The peak files carry no summit, as a narrowPeak file would in its tenth column: the peak midpoints are used instead.", call. = FALSE)
    }

    #-------------------------------#
    # Summit of every peak          #
    #-------------------------------#
    # The offset of narrowPeak counts from the 0-based start, which is the 1-based start minus one
    summitOffset <- if ("peak" %in% colnames(peakMetadata)) {as.numeric(peakMetadata$peak)} else {rep(NA_real_, length(allPeaks))}
    peakSummit <- ifelse(!is.na(summitOffset) & summitOffset >= 0,
                         BiocGenerics::start(allPeaks) + summitOffset,
                         BiocGenerics::start(allPeaks) + BiocGenerics::width(allPeaks) %/% 2)

    # A missing significance must not take a peak out of the average
    peakWeight <- if ("negLog10P" %in% colnames(peakMetadata)) {as.numeric(peakMetadata$negLog10P)} else {rep(1, length(allPeaks))}
    peakWeight[is.na(peakWeight) | peakWeight <= 0] <- 1

    #-------------------------------#
    # Average over each region      #
    #-------------------------------#
    peakHits <- IRanges::findOverlaps(regions, allPeaks, ignore.strand = TRUE)

    memberTable <- data.frame(region = S4Vectors::queryHits(peakHits),
                              region.start = BiocGenerics::start(regions)[S4Vectors::queryHits(peakHits)],
                              region.end = BiocGenerics::end(regions)[S4Vectors::queryHits(peakHits)],
                              summit = peakSummit[S4Vectors::subjectHits(peakHits)],
                              weight = peakWeight[S4Vectors::subjectHits(peakHits)])

    # A summit outside the region belongs to the part of the peak hanging over its edge
    memberTable <- dplyr::filter(memberTable, .data$summit >= .data$region.start, .data$summit <= .data$region.end)

    summitTable <- dplyr::summarise(dplyr::group_by(memberTable, .data$region),
                                    summit = round(sum(.data$summit * .data$weight) / sum(.data$weight)),
                                    .groups = "drop")

    summitPosition <- rep(NA_integer_, length(regions))
    summitPosition[summitTable$region] <- as.integer(summitTable$summit)

    return(summitPosition)
  } # END function




#' @title .recentreRanges
#'
#' @description Builds windows of fixed width around the summits, clipped to the chromosome ends.
#'
#' @param summitPosition Integer vector with the summits.
#' @param halfWidth Numeric value with the number of base pairs on each side of the summit.
#' @param chromosomeLength Numeric vector with the length of the chromosome of every summit, \code{NA} when unknown.
#'
#' @return An \code{IRanges} with one window per summit.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges IRanges
#'
#' @keywords internal

.recentreRanges <-
  function(summitPosition,
           halfWidth,
           chromosomeLength = NULL) {

    halfWidth <- as.integer(round(halfWidth))
    windowStart <- pmax(summitPosition - halfWidth, 1L)
    windowEnd <- summitPosition + halfWidth

    if (!is.null(chromosomeLength)) {
      chromosomeLength <- as.numeric(chromosomeLength)
      clipEnd <- !is.na(chromosomeLength) & windowEnd > chromosomeLength
      windowEnd[clipEnd] <- as.integer(chromosomeLength[clipEnd])
    }

    return(IRanges::IRanges(start = as.integer(windowStart), end = as.integer(windowEnd)))
  } # END function
