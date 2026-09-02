#' @title countReads
#'
#' @description Counts the reads of a group of BAM files over the regions of a \code{RegionSetDE} object. Single-end reads are extended to the expected fragment length before the overlap is evaluated, while paired-end data are counted as fragments. The regions can be cut into tiles of fixed width, in which case each tile becomes a row of the resulting object.
#'
#' @param regionSet \code{RegionSetDE} object returned by \code{\link{loadRegions}}, or a named \code{GRangesList}.
#' @param bamFiles Character vector with the paths of the BAM files. Each file must be indexed.
#' @param sampleNames Character vector with the sample names. Default: \code{NULL}, the BAM file names are used.
#' @param sampleMetadata Data.frame with the sample annotation, stored in the \code{colData}. When it contains a \code{sample} column the rows are matched by name, otherwise they must follow the order of \code{bamFiles}. Default: \code{NULL}.
#' @param keepMetadata Logical value to indicate whether the metadata columns carried by the regions must be kept in the \code{rowData}, harmonised across the sets. Default: \code{TRUE}.
#' @param regionId String with the name of a metadata column holding the region identifiers, for instance a gene name. It must hold a different value for every region of every set. Default: \code{NULL}, the names of the ranges, and their coordinates when they are unnamed.
#' @param tileWidth Numeric value with the width of the tiles, in base pairs. Default: \code{NULL}, one row per region.
#' @param partialTiles Logical value: \code{TRUE} keeps the trailing tile of each region even when narrower than \code{tileWidth}, \code{FALSE} discards it together with the regions narrower than a single tile. Default: \code{TRUE}.
#' @param pairedEnd Logical value, one logical value per BAM file, or the string \code{"auto"} to read the layout from the files themselves. Default: \code{"auto"}.
#' @param fragmentLength Numeric value with the length to which single-end reads are extended. Applied to the single-end samples only. Default: \code{150}.
#' @param maxFragmentLength Numeric value with the maximum insert size accepted for a pair. Applied to the paired-end samples only. Default: \code{1000}.
#' @param minMapq Numeric value with the minimum mapping quality of a read. Default: \code{20}.
#' @param removeDuplicates Logical value indicating whether the reads flagged as duplicates must be discarded. Default: \code{TRUE}.
#' @param restrictChromosomes Character vector with the chromosomes to read from the BAM files. Default: \code{NULL}, all of them.
#' @param discardRegions \code{GRanges} with regions whose reads must be ignored, for instance a blacklist. Default: \code{NULL}.
#' @param nThreads Number of threads used to process the files in parallel. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.counts} object with one row per region, or per tile, and one column per sample. The library sizes are stored in the \code{library.size} column of the \code{colData}, the set membership in the \code{region.set} column of the \code{rowData}.
#'
#' @details Regions shared by several sets are counted only once and the values are then copied to every set they belong to, which keeps the running time proportional to the number of distinct regions.
#'
#' Paired-end and single-end samples can be mixed in the same call. Each layout is counted in a separate pass, paired-end libraries as fragments and single-end ones as reads extended to \code{fragmentLength}, so that both end up with one count per sequenced fragment. Forcing a paired-end file through the single-end path counts each mate on its own and nearly doubles its values, while the opposite mistake keeps only the proper pairs and returns a column of zeros, which is why the layout is read from the files by default. The resolved layout of each sample is stored in the \code{paired.end} column of the \code{colData}.
#'
#' Regions and BAM files do not need to share the same chromosome naming style. When no chromosome is shared, the regions are converted to the style of the files for the counting only, so that UCSC regions can be counted on Ensembl alignments and the object still comes back with the names of the input sets.
#'
#' @examples
#' \dontrun{
#' counts <- countReads(regions,
#'                      bamFiles = list.files("bam", pattern = "\\.bam$", full.names = TRUE),
#'                      sampleMetadata = data.frame(sample = c("ctrl1", "ctrl2", "treat1", "treat2"),
#'                                                  condition = c("ctrl", "ctrl", "treat", "treat")),
#'                      pairedEnd = TRUE,
#'                      nThreads = 4)
#'
#' countsTiled <- countReads(regions, bamFiles = bamPaths, tileWidth = 500)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countBigwig}}, \code{\link{loadCounts}}, \code{\link{countBackground}}
#'
#' @importFrom csaw readParam regionCounts
#' @importFrom Rsamtools scanBamHeader testPairedEndBam
#' @importFrom GenomicRanges GRanges
#' @importFrom SummarizedExperiment assay colData
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start end strand
#' @importFrom S4Vectors mcols
#' @importFrom dplyr mutate n_distinct
#' @importFrom methods is
#'
#' @export countReads

countReads <-
  function(regionSet,
           bamFiles,
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
           restrictChromosomes = NULL,
           discardRegions = NULL,
           nThreads = 1,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!is.character(bamFiles) | length(bamFiles) == 0) {
      stop("The 'bamFiles' parameter must be a character vector with at least one BAM file.", call. = FALSE)
    }

    missingFiles <- bamFiles[!file.exists(bamFiles)]
    if (length(missingFiles) > 0) {
      stop("The following BAM files do not exist: ", paste(missingFiles, collapse = ", "), ".", call. = FALSE)
    }

    # Counting without an index would read each file from the beginning for every region
    hasIndex <-
      vapply(bamFiles,
             function(bamFile) {
               any(file.exists(c(paste0(bamFile, ".bai"),
                                 paste0(bamFile, ".csi"),
                                 sub("\\.bam$", ".bai", bamFile, ignore.case = TRUE))))
             },
             logical(1))

    if (any(!hasIndex)) {
      stop("The following BAM files are not indexed: ", paste(basename(bamFiles[!hasIndex]), collapse = ", "), ".", call. = FALSE)
    }

    if (!is.numeric(fragmentLength) | fragmentLength[1] < 1) {
      stop("The 'fragmentLength' parameter must be a positive number.", call. = FALSE)
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

    parallelParam <- .makeParallelParam(nThreads = nThreads)

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

    # Overlapping sets share regions, counting them once and copying the values back is much cheaper
    regionKey <- paste0(as.character(GenomeInfoDb::seqnames(allRegions)), ":",
                        BiocGenerics::start(allRegions), "-",
                        BiocGenerics::end(allRegions), ":",
                        as.character(BiocGenerics::strand(allRegions)))

    uniqueRegions <- allRegions[!duplicated(regionKey)]
    expansionIndex <- match(regionKey, regionKey[!duplicated(regionKey)])

    #---------------------------------#
    # Chromosome names of the samples #
    #---------------------------------#
    # Only the copy handed to csaw is renamed, the returned object keeps the style of the region sets
    bamSeqlevels <- names(Rsamtools::scanBamHeader(bamFiles[1])[[1]]$targets)
    countingRegions <- .matchSeqlevels(x = uniqueRegions, targetSeqlevels = bamSeqlevels, fileName = bamFiles[1], verbose = verbose)

    #----------------#
    # Read the files #
    #----------------#
    if (isTRUE(verbose)) {
      message("Counting reads in ", length(bamFiles), " samples over ", length(uniqueRegions), " unique regions (",
              length(allRegions), " rows, ", dplyr::n_distinct(S4Vectors::mcols(allRegions)$region.set), " sets)...")
    }

    if (isTRUE(verbose) & length(unique(pairedEnd)) > 1) {
      message("Mixed layouts: ", sum(pairedEnd), " paired-end and ", sum(!pairedEnd),
              " single-end samples, counted in two separate passes.")
    }

    countMatrix <- matrix(0, nrow = length(uniqueRegions), ncol = length(bamFiles))
    librarySizes <- numeric(length(bamFiles))

    # A readParam object describes one layout only, so the paired and single-end files go through their own pass
    for (layout in unique(pairedEnd)) {
      layoutIndex <- which(pairedEnd == layout)

      readParameters <- csaw::readParam(pe = ifelse(layout, "both", "none"),
                                        max.frag = maxFragmentLength,
                                        dedup = removeDuplicates,
                                        minq = minMapq,
                                        restrict = restrictChromosomes,
                                        discard = if (is.null(discardRegions)) {GenomicRanges::GRanges()} else {discardRegions})

      rawCounts <- csaw::regionCounts(bam.files = bamFiles[layoutIndex],
                                      regions = countingRegions,
                                      ext = ifelse(layout, NA, as.integer(fragmentLength[1])),
                                      param = readParameters,
                                      BPPARAM = parallelParam)

      countMatrix[, layoutIndex] <- SummarizedExperiment::assay(rawCounts, "counts")
      librarySizes[layoutIndex] <- as.numeric(SummarizedExperiment::colData(rawCounts)$totals)
    }

    #-------------------------#
    # Assemble the object     #
    #-------------------------#
    # The totals are the reads surviving the filters, which is the denominator the normalisation expects
    sampleTable <- dplyr::mutate(sampleTable, paired.end = pairedEnd, library.size = librarySizes)

    countMatrix <- countMatrix[expansionIndex, , drop = FALSE]

    newParameters <- list(countReads = list(bamFiles = bamFiles,
                                            tileWidth = tileWidth,
                                            partialTiles = partialTiles,
                                            pairedEnd = pairedEnd,
                                            fragmentLength = fragmentLength,
                                            maxFragmentLength = maxFragmentLength,
                                            minMapq = minMapq,
                                            removeDuplicates = removeDuplicates,
                                            restrictChromosomes = restrictChromosomes))

    counts <- .newCountsObject(countMatrix = countMatrix,
                               regions = allRegions,
                               sampleTable = sampleTable,
                               provenance = .provenanceSlots(regionSet),
                               countingLevel = "region",
                               newParameters = newParameters,
                               metadataList = list(signal.type = "reads"))

    if (isTRUE(verbose)) {
      message("Done. Library sizes: ",
              paste0(round(range(sampleTable$library.size) / 1e6, 1), collapse = " - "), " million reads.")
    }

    return(counts)
  } # END function
