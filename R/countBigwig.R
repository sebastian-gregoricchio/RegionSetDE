# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title countBigwig
#'
#' @description Summarises the signal of a group of bigWig files over the regions of a \code{RegionSetDE} object. Useful when the BAM files are not available, or when the coverage has been produced by an external pipeline. The regions can be cut into tiles of fixed width, in which case each tile becomes a row of the resulting object.
#'
#' @param regionSet \code{RegionSetDE} object returned by \code{\link{loadRegions}}, or a named \code{GRangesList}.
#' @param bigwigFiles Character vector with the paths of the bigWig files. Default: \code{NULL}, taken from \code{sampleSheet}, or from the sample sheet a consensus was built from by \code{\link{loadConsensusPeaks}}.
#' @param sampleSheet Data.frame returned by \code{\link{loadSampleSheet}}, or the path to a sample sheet, providing the bigWig files, the sample names and the annotation in one go. Default: \code{NULL}.
#' @param sampleNames Character vector with the sample names. Default: \code{NULL}, the bigWig file names are used.
#' @param sampleMetadata Data.frame with the sample annotation, stored in the \code{colData}. When it contains a \code{sample} column the rows are matched by name, otherwise they must follow the order of \code{bigwigFiles}. Default: \code{NULL}.
#' @param keepMetadata Logical value to indicate whether the metadata columns carried by the regions must be kept in the \code{rowData}, harmonised across the sets. Default: \code{TRUE}.
#' @param regionId String with the name of a metadata column holding the region identifiers, for instance a gene name. It must hold a different value for every region of every set. Default: \code{NULL}, the names of the ranges, and their coordinates when they are unnamed.
#' @param tileWidth Numeric value with the width of the tiles, in base pairs. Default: \code{NULL}, one row per region.
#' @param partialTiles Logical value: \code{TRUE} keeps the trailing tile of each region even when narrower than \code{tileWidth}, \code{FALSE} discards it together with the regions narrower than a single tile. Default: \code{TRUE}.
#' @param summaryFunction String indicating how the per-base values are collapsed into a single value per region, one among \code{"sum"}, \code{"mean"}, \code{"max"} or \code{"min"}. Default: \code{"sum"}.
#' @param missingAsZero Logical value indicating whether the positions not covered by the bigWig must be treated as zeros rather than as missing values. Default: \code{TRUE}.
#' @param countLike Logical value with which you assert that the bigWig holds count-like values, meaning raw, unnormalised coverage that a count model may legitimately be applied to. Default: \code{FALSE}.
#' @param roundValues Logical value indicating whether the summarised values must be rounded to integers. Rounding is a formatting step and does not turn coverage into counts. Default: \code{FALSE}.
#' @param nThreads Number of threads used to process the files in parallel. Default: \code{1}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.counts} object with one row per region, or per tile, and one column per sample. The \code{signal.type} entry of the metadata is set to \code{"bigwig"} and \code{count.like} to whatever was declared through \code{countLike}, which is what \code{\link{fitRegions}} reads to decide which engines it will run.
#'
#' @details A bigWig holds coverage, not reads, so the library sizes cannot be recovered from it: the \code{library.size} column of the \code{colData} is left as \code{NA} and the total signal falling in the regions is reported in \code{total.signal} instead. Normalisation factors must therefore be supplied externally, or estimated from a background bigWig.
#'
#' Coverage is not a count and rounding it does not make it one. A value of 12.72 becomes 13 and looks like a count, but the negative binomial likelihood that \code{edgeR} and \code{DESeq2} are built on describes the number of fragments falling in an interval, and a rounded coverage value is not that number. The gap is widest for files carrying an already normalised signal, CPM, RPKM, RPGC, fold enrichment over input, where the values have been divided by a factor that the count model then has no way of knowing about. It does not close entirely even for raw coverage: coverage summed over an interval weights every fragment by how many of its bases fall inside, so it is over-dispersed relative to the fragment count it stands in for.
#'
#' The object therefore records where its values came from, and \code{\link{fitRegions}} refuses the count engines on it unless \code{countLike} says otherwise. The engine to reach for on bigWig input is \code{"limma"}, which models the log2 signal directly and asks nothing of the values that they cannot supply. Setting \code{countLike = TRUE} is an assertion about the files, not a setting: it says these are raw, unnormalised coverage tracks and the count model is close enough for the purpose, and it should be stated in the methods when it is used.
#'
#' @examples
#' # The small bigWig shipped with rtracklayer, with a region set built on its own intervals
#' bigwigFile <- system.file("tests", "test.bw", package = "rtracklayer")
#'
#' if (file.exists(bigwigFile)) {
#'   bigwigRanges <- GenomicRanges::reduce(rtracklayer::import(bigwigFile))
#'   regions <- loadRegions(list(covered = bigwigRanges), seqlevelsStyle = NULL, verbose = FALSE)
#'
#'   signal <- countBigwig(regions,
#'                         bigwigFiles = bigwigFile,
#'                         sampleNames = "example",
#'                         summaryFunction = "mean",
#'                         verbose = FALSE)
#'   SummarizedExperiment::assay(signal)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countReads}}, \code{\link{loadCounts}}
#'
#' @importFrom rtracklayer BigWigFile import
#' @importFrom GenomeInfoDb seqnames seqlevels
#' @importFrom BiocGenerics start end strand width
#' @importFrom BiocParallel bplapply
#' @importFrom S4Vectors mcols
#' @importFrom dplyr mutate n_distinct
#'
#' @export countBigwig

countBigwig <-
  function(regionSet,
           bigwigFiles = NULL,
           sampleSheet = NULL,
           sampleNames = NULL,
           sampleMetadata = NULL,
           tileWidth = NULL,
           keepMetadata = TRUE,
           regionId = NULL,
           partialTiles = TRUE,
           summaryFunction = "sum",
           missingAsZero = TRUE,
           countLike = FALSE,
           roundValues = FALSE,
           nThreads = 1,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    # A sample sheet brings the files, the names and the annotation in one go
    sheetInput <- .sheetCountingInput(regionSet = regionSet,
                                      sampleSheet = sampleSheet,
                                      files = bigwigFiles,
                                      sampleNames = sampleNames,
                                      sampleMetadata = sampleMetadata,
                                      fileField = "bigwig",
                                      verbose = verbose)

    bigwigFiles <- sheetInput$files
    sampleNames <- sheetInput$sampleNames
    sampleMetadata <- sheetInput$sampleMetadata

    if (!is.character(bigwigFiles) | length(bigwigFiles) == 0) {
      stop("The 'bigwigFiles' parameter must be a character vector with at least one bigWig file.", call. = FALSE)
    }

    missingFiles <- bigwigFiles[!file.exists(bigwigFiles)]
    if (length(missingFiles) > 0) {
      stop("The following bigWig files do not exist: ", paste(missingFiles, collapse = ", "), ".", call. = FALSE)
    }

    summaryFunction <- tolower(summaryFunction[1])
    if (!(summaryFunction %in% c("sum", "mean", "max", "min"))) {
      stop("The 'summaryFunction' parameter must be one among 'sum', 'mean', 'max' or 'min'.", call. = FALSE)
    }

    if (!is.logical(countLike) | length(countLike) != 1 | is.na(countLike)) {
      stop("The 'countLike' parameter must be a single logical value.", call. = FALSE)
    }

    # Rounding is what usually gets coverage into a count model unnoticed, so it is said out loud here
    if (isTRUE(roundValues) & !isTRUE(countLike)) {
      warning("The values are being rounded to integers, which does not make them counts. ",
              "Fit them with engine = 'limma', or declare 'countLike = TRUE' if the files hold raw coverage.", call. = FALSE)
    }

    parallelParam <- .makeParallelParam(nThreads = nThreads)

    #---------------------#
    # Samples and regions #
    #---------------------#
    sampleTable <- .buildSampleTable(files = bigwigFiles,
                                     sampleNames = sampleNames,
                                     sampleMetadata = sampleMetadata,
                                     fileColumn = "bigwig.file",
                                     extensionPattern = "\\.(bw|bigwig|bigWig)$")

    allRegions <- .flattenRegionSets(regionSet = regionSet,
                                     tileWidth = tileWidth,
                                     keepMetadata = keepMetadata,
                                     regionId = regionId,
                                     partialTiles = partialTiles,
                                     verbose = verbose)

    # Reading the same interval once per set would multiply the query time by the number of sets
    regionKey <- paste0(as.character(GenomeInfoDb::seqnames(allRegions)), ":",
                        BiocGenerics::start(allRegions), "-",
                        BiocGenerics::end(allRegions), ":",
                        as.character(BiocGenerics::strand(allRegions)))

    uniqueRegions <- allRegions[!duplicated(regionKey)]
    expansionIndex <- match(regionKey, regionKey[!duplicated(regionKey)])
    uniqueWidths <- BiocGenerics::width(uniqueRegions)

    # The query must speak the language of the file, the returned object keeps the style of the region sets
    uniqueRegions <- .matchSeqlevels(x = uniqueRegions,
                                     targetSeqlevels = GenomeInfoDb::seqlevels(rtracklayer::BigWigFile(bigwigFiles[1])),
                                     fileName = bigwigFiles[1],
                                     verbose = verbose)

    #----------------#
    # Read the files #
    #----------------#
    if (isTRUE(verbose)) {
      message("Extracting the signal of ", length(bigwigFiles), " bigWig files over ", length(uniqueRegions), " unique regions (",
              length(allRegions), " rows, ", dplyr::n_distinct(S4Vectors::mcols(allRegions)$region.set), " sets)...")
    }

    signalList <-
      BiocParallel::bplapply(X = bigwigFiles,
                             BPPARAM = parallelParam,
                             FUN = function(bigwigFile) {
                               baseValues <- rtracklayer::import(con = rtracklayer::BigWigFile(bigwigFile),
                                                                 which = uniqueRegions,
                                                                 as = "NumericList")

                               # Uncovered bases are missing in the file, either they count as zeros or they leave the region shorter
                               regionSignal <-
                                 switch(summaryFunction,
                                        "sum" = sum(baseValues, na.rm = TRUE),
                                        "mean" = if (isTRUE(missingAsZero)) {sum(baseValues, na.rm = TRUE) / uniqueWidths} else {mean(baseValues, na.rm = TRUE)},
                                        "max" = max(baseValues, na.rm = TRUE),
                                        "min" = min(baseValues, na.rm = TRUE))

                               regionSignal <- as.numeric(regionSignal)
                               regionSignal[!is.finite(regionSignal)] <- 0
                               return(regionSignal)
                             })

    signalMatrix <- do.call(what = cbind, args = signalList)

    #---------------------#
    # Assemble the object #
    #---------------------#
    # Only a formatting step: it is here for the engines that insist on integers, and it changes nothing
    # about what the values are
    if (isTRUE(roundValues)) {
      signalMatrix <- round(signalMatrix)
    }

    sampleTable <- dplyr::mutate(sampleTable,
                                 library.size = NA_real_,
                                 total.signal = as.numeric(colSums(signalMatrix)))

    signalMatrix <- signalMatrix[expansionIndex, , drop = FALSE]

    newParameters <- list(countBigwig = list(bigwigFiles = bigwigFiles,
                                             tileWidth = tileWidth,
                                             partialTiles = partialTiles,
                                             summaryFunction = summaryFunction,
                                             missingAsZero = missingAsZero,
                                             countLike = countLike,
                                             roundValues = roundValues))

    counts <- .newCountsObject(countMatrix = signalMatrix,
                               regions = allRegions,
                               sampleTable = sampleTable,
                               provenance = .provenanceSlots(regionSet),
                               countingLevel = if (is.null(tileWidth)) {"region"} else {"tile"},
                               newParameters = newParameters,
                               metadataList = list(signal.type = "bigwig", count.like = countLike))

    if (isTRUE(verbose)) {
      message("Done. The library sizes are unknown for bigWig input, supply the normalisation factors at the normalisation step.")
      if (!isTRUE(countLike)) {
        message("These values are coverage rather than counts: fit them with engine = 'limma'.")
      }
    }

    return(counts)
  } # END function
