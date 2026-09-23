# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title countReads
#'
#' @description Counts the reads of a group of BAM files over the regions of a \code{RegionSetDE} object. Paired-end data are counted as fragments, while single-end reads are extended to the expected fragment length before the overlap is evaluated. The regions can be cut into tiles of fixed width, in which case each tile becomes a row of the resulting object.
#'
#' @param regionSet \code{RegionSetDE} object returned by \code{\link{loadRegions}}, or a named \code{GRangesList}.
#' @param bamFiles Character vector with the paths of the BAM files. Each file must be indexed, and all of them must share the same header. Default: \code{NULL}, taken from \code{sampleSheet}, or from the sample sheet a consensus was built from by \code{\link{loadConsensusPeaks}}.
#' @param sampleSheet Data.frame returned by \code{\link{loadSampleSheet}}, or the path to a sample sheet, providing the BAM files, the sample names and the annotation in one go. Default: \code{NULL}.
#' @param sampleNames Character vector with the sample names. Default: \code{NULL}, the BAM file names are used.
#' @param sampleMetadata Data.frame with the sample annotation, stored in the \code{colData}. When it contains a \code{sample} column the rows are matched by name, otherwise they must follow the order of \code{bamFiles}. Default: \code{NULL}.
#' @param keepMetadata Logical value to indicate whether the metadata columns carried by the regions must be kept in the \code{rowData}, harmonised across the sets. Default: \code{TRUE}.
#' @param regionId String with the name of a metadata column holding the region identifiers, for instance a gene name. It must hold a different value for every region of every set. Default: \code{NULL}, the names of the ranges, and their coordinates when they are unnamed.
#' @param tileWidth Numeric value with the width of the tiles, in base pairs. Default: \code{NULL}, one row per region.
#' @param partialTiles Logical value: \code{TRUE} keeps the trailing tile of each region even when narrower than \code{tileWidth}, \code{FALSE} discards it together with the regions narrower than a single tile. Default: \code{TRUE}.
#' @param pairedEnd Logical value, one logical value per BAM file, or the string \code{"auto"} to read the layout from the files themselves. Default: \code{"auto"}.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended. Applied to the single-end samples only. Default: \code{150}.
#' @param maxFragmentLength Numeric value with the maximum length accepted for a paired-end fragment. Applied to the paired-end samples only. Default: \code{1000}.
#' @param minMapq Numeric value with the minimum mapping quality of a read. Default: \code{20}.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded. Default: \code{TRUE}.
#' @param excludeChromosomes Character vector with the chromosomes left out of the library sizes, written in either naming style, \code{chrM} and \code{MT} both reaching the mitochondrial genome of a file naming it either way, for instance the mitochondrial genome, chrY or the unplaced and alternative contigs. A name matching no chromosome once converted raises a warning, since it would leave that chromosome inside the library sizes without a word. The regions lying on them are still counted: to leave those out as well, filter the regions when loading them. The contigs can be collected from the BAM header, e.g. \code{grep("_|EBV", names(Rsamtools::scanBamHeader(bamFile)[[1]]$targets), value = TRUE)}. Default: \code{NULL}, every chromosome enters the library sizes.
#' @param discardRegions \code{GRanges} with regions whose reads must be ignored, for instance a blacklist. A fragment is dropped when one of its reads starts inside them. Default: \code{NULL}.
#' @param fullLibrarySize Logical value: \code{TRUE} reads every chromosome that is not excluded, even those without any region, so that the library sizes cover the whole library; \code{FALSE} reads only the chromosomes carrying regions, which is much faster for a few regions but leaves library sizes that must not be used for normalisation. Default: \code{TRUE}.
#' @param nThreads Number of threads. The files are cut into pieces of at most 50 Mb, shared among the threads, so even a single file benefits from several of them. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.counts} object with one row per region, or per tile, and one column per sample. The library sizes are stored in the \code{library.size} column of the \code{colData}, the set membership in the \code{region.set} column of the \code{rowData}.
#'
#' @details Regions shared by several sets are counted only once and the values are then copied to every set they belong to, which keeps the running time proportional to the number of distinct regions.
#'
#' A paired-end fragment is counted in every region it overlaps, including the regions it spans with both reads outside them. The fragment is rebuilt from the first mate of each proper pair, whose position and template length (TLEN) give its start and width, so the two reads never have to be matched in memory. The pairs therefore have to be flagged as proper by the aligner, and those longer than \code{maxFragmentLength} are dropped. The mapping quality of the second mate is read from the \code{MQ} tag, which \code{samtools fixmate} and Picard write; on files without it only the first mate is checked, and a message says so.
#'
#' Counts and library sizes are kept apart. The counts only need the chromosomes carrying regions, and every region is counted, wherever it lies. The library size of a sample is the number of fragments that went through the same filters as the counts, on every chromosome of the BAM files except those in \code{excludeChromosomes}. Leaving out the mitochondrial genome matters in ATAC-seq, where its share of the reads changes from sample to sample. With \code{fullLibrarySize = FALSE} only the chromosomes carrying regions are read, and the library sizes are partial: the counts do not change, but the library sizes are not usable for normalisation, and \code{\link{normalizeCounts}} warns when a method relies on them.
#'
#' Paired-end and single-end samples can be mixed in the same call, paired-end libraries being counted as fragments and single-end ones as reads extended to \code{fragmentLength}, so that both end up with one count per sequenced fragment. Forcing a paired-end file through the single-end path counts each mate on its own and nearly doubles its values, while the opposite mistake finds no pair and returns a column of zeros, which is why the layout is read from the files by default. The resolved layout of each sample is stored in the \code{paired.end} column of the \code{colData}.
#'
#' Regions and BAM files do not need to share the same chromosome naming style. When no chromosome is shared, the regions are converted to the style of the files for the counting only, so that UCSC regions can be counted on Ensembl alignments and the object still comes back with the names of the input sets.
#'
#' @examples
#' # The peaks of one sample of the AR example stand in for a region set
#' sampleSheet <- loadExampleData("peakSheet", verbose = FALSE)
#' peakRegions <- loadRegions(list(peaks = sampleSheet$peaks[7]), genomeAssembly = "hg38", verbose = FALSE)
#'
#' # The sheet brings the BAM files, the sample names and the annotation
#' counts <- countReads(peakRegions, sampleSheet = sampleSheet, verbose = FALSE)
#' counts
#'
#' # The same files given one by one, with the annotation as a table
#' counts <- countReads(peakRegions,
#'                      bamFiles = sampleSheet$bam,
#'                      sampleNames = sampleSheet$sample,
#'                      sampleMetadata = sampleSheet[, c("sample", "condition")],
#'                      verbose = FALSE)
#'
#' # Tiles of 100 bp, one row each
#' tiledCounts <- countReads(peakRegions, sampleSheet = sampleSheet, tileWidth = 100, verbose = FALSE)
#' head(SummarizedExperiment::rowData(tiledCounts), 3)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countBigwig}}, \code{\link{loadCounts}}, \code{\link{countBackground}}
#'
#' @importFrom Rsamtools scanBamHeader testPairedEndBam
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end
#' @importFrom S4Vectors mcols
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
           nThreads = 1,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
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

    if (!is.numeric(fragmentLength) | fragmentLength[1] < 1) {
      stop("The 'fragmentLength' parameter must be a positive number.", call. = FALSE)
    }

    if (!is.numeric(maxFragmentLength) | maxFragmentLength[1] < 1) {
      stop("The 'maxFragmentLength' parameter must be a positive number.", call. = FALSE)
    }

    if (!is.logical(fullLibrarySize) | length(fullLibrarySize) != 1 | any(is.na(fullLibrarySize))) {
      stop("The 'fullLibrarySize' parameter must be TRUE or FALSE.", call. = FALSE)
    }

    # The flags of the first records are enough to tell a paired library from a single-end one
    if (identical(pairedEnd, "auto")) {
      pairedEnd <- vapply(bamFiles, Rsamtools::testPairedEndBam, logical(1), USE.NAMES = FALSE)
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

    bamSeqlevels <- names(Rsamtools::scanBamHeader(bamFiles[1])[[1]]$targets)

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
                                         fragmentLength = fragmentLength[1],
                                         maxFragmentLength = maxFragmentLength[1],
                                         minMapq = minMapq,
                                         removeDuplicates = removeDuplicates,
                                         excludeChromosomes = excludeChromosomes,
                                         discardRegions = discardRegions,
                                         fullLibrarySize = fullLibrarySize,
                                         countMode = "overlap",
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

    #-------------------------#
    # Assemble the object     #
    #-------------------------#
    # The library sizes are the fragments surviving the filters, which is the denominator the normalisation expects
    sampleTable <- dplyr::mutate(sampleTable, paired.end = pairedEnd, library.size = fragmentCounts$library.size)

    countMatrix <- countMatrix[expansionIndex, , drop = FALSE]

    newParameters <- list(countReads = list(bamFiles = bamFiles,
                                            tileWidth = tileWidth,
                                            partialTiles = partialTiles,
                                            pairedEnd = pairedEnd,
                                            fragmentLength = fragmentLength,
                                            maxFragmentLength = maxFragmentLength,
                                            minMapq = minMapq,
                                            removeDuplicates = removeDuplicates,
                                            excludeChromosomes = excludeChromosomes,
                                            fullLibrarySize = fullLibrarySize))

    # A tiled object has to say so it is tiled, otherwise testRegions treats every tile as a region
    # and the combination step that puts the tiles back together never runs
    counts <- .newCountsObject(countMatrix = countMatrix,
                               regions = allRegions,
                               sampleTable = sampleTable,
                               provenance = .provenanceSlots(regionSet),
                               countingLevel = if (is.null(tileWidth)) {"region"} else {"tile"},
                               newParameters = newParameters,
                               metadataList = list(signal.type = "reads", count.like = TRUE))

    if (isTRUE(verbose)) {
      message("Done. Library sizes: ",
              paste(round(range(sampleTable$library.size) / 1e6, 1), collapse = " - "), " million fragments.")

      if (isFALSE(fullLibrarySize)) {
        message("The library sizes cover only the chromosomes carrying regions: count again with 'fullLibrarySize = TRUE' before a normalisation that relies on them.")
      }
    }

    return(counts)
  } # END function
