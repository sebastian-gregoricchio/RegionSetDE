#' @title asDGEList
#'
#' @description Turns a \code{RegionSetDE.counts} object into a \code{DGEList}, with the normalisation stored in the object carried across as offsets and the region annotation kept in the \code{genes} slot.
#'
#' @param counts \code{RegionSetDE.counts} object, or a \code{RegionSetDE.fit}, in which case the list the model was fitted on is returned as it is.
#' @param assay String with the name of the assay holding the counts. Default: \code{"counts"}.
#' @param useOffsets Logical value to indicate whether the normalisation stored in the object must be carried across. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{DGEList}.
#'
#' @details This exists because the conversion is easy to get backwards. The object stores a \code{scaling.factor} that normalised values are divided by, while a generalised linear model wants an offset on the log scale of the effective library size, and the two differ by a sign and a constant. Writing the factors into \code{lib.size} by hand, or exponentiating them the wrong way, produces a \code{DGEList} that fits without complaint and reports fold changes in the wrong direction. Here the factors go through \code{edgeR::scaleOffset}, which puts them on the scale of the library sizes and leaves the coefficients readable, and there is no direction left to get wrong.
#'
#' \code{as(counts, "DGEList")} does the same with the defaults.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#' dgeList <- asDGEList(counts)
#' dgeList$samples
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{asDESeqDataSet}}, \code{\link{fitRegions}}, \code{\link{normalizeCounts}}
#'
#' @importFrom SummarizedExperiment assay assayNames colData rowData rowRanges
#' @importFrom BiocGenerics width
#' @importFrom GenomeInfoDb seqnames
#' @importFrom edgeR DGEList scaleOffset
#' @importFrom methods is
#'
#' @export asDGEList

asDGEList <-
  function(counts,
           assay = "counts",
           useOffsets = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    # A fit already carries the list the model saw, rebuilding it could only make it differ
    if (methods::is(counts, "RegionSetDE.fit")) {
      if (!is.null(counts@fit$dge)) {
        return(counts@fit$dge)
      }
      counts <- counts@counts
    }

    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts or a RegionSetDE.fit object.", call. = FALSE)
    }

    if (!(assay %in% SummarizedExperiment::assayNames(counts))) {
      stop(paste0("The assay '", assay, "' is absent from the object."), call. = FALSE)
    }

    #-------------------------------#
    # Counts, sizes and annotation  #
    #-------------------------------#
    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, assay))

    librarySizes <- SummarizedExperiment::colData(counts)$library.size
    if (is.null(librarySizes) | any(is.na(librarySizes))) {
      librarySizes <- colSums(countMatrix)
    }

    rowRangesObject <- SummarizedExperiment::rowRanges(counts)
    geneTable <- data.frame(as.data.frame(SummarizedExperiment::rowData(counts)),
                            seqnames = as.character(GenomeInfoDb::seqnames(rowRangesObject)),
                            start = BiocGenerics::start(rowRangesObject),
                            end = BiocGenerics::end(rowRangesObject),
                            width = BiocGenerics::width(rowRangesObject),
                            stringsAsFactors = FALSE)

    dgeList <- edgeR::DGEList(counts = countMatrix,
                              samples = as.data.frame(SummarizedExperiment::colData(counts)),
                              genes = geneTable,
                              lib.size = librarySizes)

    #-------------------------------#
    # Normalisation as offsets      #
    #-------------------------------#
    offsetMatrix <- .fitOffsets(counts = counts, useOffsets = useOffsets, verbose = FALSE)

    if (!is.null(offsetMatrix)) {
      dgeList <- edgeR::scaleOffset(y = dgeList, offset = offsetMatrix)

      if (isTRUE(verbose)) {
        message("The normalisation has been carried across as offsets, so do not set 'norm.factors' or 'lib.size' on top of it.")
      }
    } else if (isTRUE(verbose)) {
      message("No offset was carried across, the list holds the raw counts and the library sizes.")
    }

    return(dgeList)
  } # END function




#' @title asDESeqDataSet
#'
#' @description Turns a \code{RegionSetDE.counts} object into a \code{DESeqDataSet}, with the normalisation stored in the object carried across as normalisation factors.
#'
#' @param counts \code{RegionSetDE.counts} object, or a \code{RegionSetDE.fit}.
#' @param design Formula evaluated on the \code{colData}, the same formula written as a string, or a design matrix. Default: \code{NULL}, the design of the fit when one is given and \code{~ 1} otherwise.
#' @param assay String with the name of the assay holding the counts. Default: \code{"counts"}.
#' @param useOffsets Logical value to indicate whether the normalisation stored in the object must be carried across. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{DESeqDataSet}, unfitted.
#'
#' @details The normalisation goes into \code{normalizationFactors} rather than \code{sizeFactors}, since the object may hold one value per region and per sample rather than one per sample. \code{DESeq2} asks those factors to have a geometric mean of one on each row, which is what the row centring here is for; without it the fitted values come out on a scale that has nothing to do with the counts.
#' The object is returned unfitted, so \code{DESeq2::DESeq} still has to be run on it. Nothing stops a design being passed that differs from the one the package used, which is the point of the conversion, but a result obtained that way is a different analysis and not a check of this one.
#' There is no \code{as(counts, "DESeqDataSet")} to go with it, because a coercion has to be registered when the package is built and the class it points at only exists when \code{DESeq2} is installed.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#'
#' if (requireNamespace("DESeq2", quietly = TRUE)) {
#'   deseqData <- asDESeqDataSet(counts, design = ~ condition, verbose = FALSE)
#'   deseqData
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{asDGEList}}, \code{\link{fitRegions}}
#'
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom methods is
#'
#' @export asDESeqDataSet

asDESeqDataSet <-
  function(counts,
           design = NULL,
           assay = "counts",
           useOffsets = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("The 'DESeq2' package is needed for this conversion.", call. = FALSE)
    }

    if (methods::is(counts, "RegionSetDE.fit")) {
      if (is.null(design)) {
        design <- counts@design
      }
      counts <- counts@counts
    }

    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts or a RegionSetDE.fit object.", call. = FALSE)
    }

    if (!(assay %in% SummarizedExperiment::assayNames(counts))) {
      stop(paste0("The assay '", assay, "' is absent from the object."), call. = FALSE)
    }

    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, assay))

    if (any(abs(countMatrix - round(countMatrix)) > 1e-8)) {
      stop("DESeq2 needs integer counts, use the raw 'counts' assay and let the offsets carry the normalisation.", call. = FALSE)
    }
    storage.mode(countMatrix) <- "integer"

    #-------------------------------#
    # Design                        #
    #-------------------------------#
    colTable <- as.data.frame(SummarizedExperiment::colData(counts))
    rownames(colTable) <- colnames(counts)

    designMatrix <- NULL
    if (!is.null(design)) {
      designMatrix <- .buildDesign(design = design, random = NULL, colTable = colTable, engine = "deseq2")$matrix
    }

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countMatrix,
                                          colData = colTable,
                                          design = ~ 1,
                                          rowRanges = SummarizedExperiment::rowRanges(counts))

    #-------------------------------#
    # Normalisation                 #
    #-------------------------------#
    offsetMatrix <- .fitOffsets(counts = counts, useOffsets = useOffsets, verbose = FALSE)

    if (!is.null(offsetMatrix)) {
      # DESeq2 asks the factors to have geometric mean one on each row, which is where the row centring comes from
      DESeq2::normalizationFactors(dds) <- exp(offsetMatrix - rowMeans(offsetMatrix))
    }

    if (isTRUE(verbose)) {
      message(paste0("The object is unfitted. Run DESeq2::DESeq(dds",
                     if (!is.null(designMatrix)) {", full = <design matrix>, betaPrior = FALSE"} else {""}, ")."))

      if (!is.null(designMatrix)) {
        message(paste0("The design matrix is returned as the 'designMatrix' attribute, with the coefficients: ",
                       paste(colnames(designMatrix), collapse = ", "), "."))
      }
    }

    attr(dds, "designMatrix") <- designMatrix
    return(dds)
  } # END function




#' @title Coerce a counts object to a DGEList
#'
#' @description Registers the \code{as(x, "DGEList")} idiom, which calls \code{\link{asDGEList}} with its defaults.
#'
#' @name coerce-RegionSetDE.counts-DGEList
#' @aliases coerce,RegionSetDE.counts,DGEList-method
#'
#' @return A \code{DGEList} holding the counts, the sample metadata, the region annotation and the normalisation as offsets.
#'
#' @examples
#' \dontrun{
#' dgeList <- as(counts, "DGEList")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{asDGEList}}
#'
#' @importFrom methods setAs setOldClass isClass
#'
#' @exportMethod coerce
NULL


# edgeR has carried DGEList as an old-style class, and registering it is what lets setAs point at it.
# The guard is there in case a future version defines it formally, where re-registering would be an error.
if (!methods::isClass("DGEList")) {
  methods::setOldClass("DGEList")
}


setAs(from = "RegionSetDE.counts",
      to = "DGEList",
      def = function(from) {asDGEList(counts = from, verbose = FALSE)})
