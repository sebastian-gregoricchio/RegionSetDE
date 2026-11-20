# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title Methods for the generics of BiocGenerics
#'
#' @description The objects of the package answer the generics other Bioconductor packages define for the same kind of data, so that code written for \code{DESeq2} or \code{edgeR} objects reads them too. \code{counts} returns the raw count matrix, \code{design} the design matrix of a fit, \code{plotMA} the MA plot and \code{plotPCA} the principal component plot. They call the functions of the package of the same purpose, which are documented in full on their own pages.
#'
#' @param object For \code{counts}, a \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results} or \code{RegionSetDE.setResults} object, or one of the two list classes holding several contrasts. For \code{design}, a \code{RegionSetDE.fit} object. For \code{plotMA}, a \code{RegionSetDE.results} or \code{RegionSetDE.resultsList} object, drawn by \code{\link{plotResultsMA}}, or a normalised \code{RegionSetDE.counts} object, drawn by \code{\link{plotSetMA}}. For \code{plotPCA}, any object accepted by \code{\link{plotRegionPCA}}.
#' @param ... For \code{plotMA} and \code{plotPCA}, the arguments of \code{\link{plotResultsMA}}, \code{\link{plotSetMA}} or \code{\link{plotRegionPCA}}, such as \code{set} or \code{contrast}. Not used by \code{counts} and \code{design}.
#'
#' @return \code{counts} returns a numeric matrix with one row per region, or per tile, and one column per sample, holding the raw counts. For a result they are the counts it carries, which needs the test to have been run with \code{carryCounts = TRUE}; \code{\link{countTable}} returns the normalised values or a table with the coordinates. \code{design} returns the matrix of the fixed effects, one row per sample, as \code{\link{designMatrix}} does. \code{plotMA} and \code{plotPCA} return a \code{ggplot} object.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' dim(counts(fit))
#' head(design(fit))
#'
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#' plotMA(results)
#' plotPCA(fit)
#'
#' @author Sebastian Gregoricchio
#'
#' @rdname biocGenerics-methods
#'
#' @importFrom BiocGenerics counts design plotMA plotPCA
#' @importFrom SummarizedExperiment assay
#' @importFrom methods setMethod
#'
#' @export
setMethod(
    f = "counts",
    signature = "RegionSetDE.counts",
    definition = function(object, ...) {
        return(SummarizedExperiment::assay(object, "counts"))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "counts",
    signature = "RegionSetDE.fit",
    definition = function(object, ...) {
        return(counts(fitCounts(object)))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "counts",
    signature = "RegionSetDE.results",
    definition = function(object, ...) {
        return(counts(.carriedCounts(object)))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "counts",
    signature = "RegionSetDE.setResults",
    definition = function(object, ...) {
        return(counts(.carriedCounts(object)))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "counts",
    signature = "RegionSetDE.resultsList",
    definition = function(object, ...) {
        return(counts(.carriedCounts(object)))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "counts",
    signature = "RegionSetDE.setResultsList",
    definition = function(object, ...) {
        return(counts(.carriedCounts(object)))
    }
)


#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "design",
    signature = "RegionSetDE.fit",
    definition = function(object, ...) {
        return(designMatrix(object))
    }
)


#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "plotMA",
    signature = "RegionSetDE.results",
    definition = function(object, ...) {
        return(plotResultsMA(results = object, ...))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "plotMA",
    signature = "RegionSetDE.resultsList",
    definition = function(object, ...) {
        return(plotResultsMA(results = object, ...))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "plotMA",
    signature = "RegionSetDE.counts",
    definition = function(object, ...) {
        return(plotSetMA(counts = object, ...))
    }
)


#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "plotPCA",
    signature = "RegionSetDE.counts",
    definition = function(object, ...) {
        return(plotRegionPCA(object = object, ...))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "plotPCA",
    signature = "RegionSetDE.fit",
    definition = function(object, ...) {
        return(plotRegionPCA(object = object, ...))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "plotPCA",
    signature = "RegionSetDE.results",
    definition = function(object, ...) {
        return(plotRegionPCA(object = object, ...))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "plotPCA",
    signature = "RegionSetDE.resultsList",
    definition = function(object, ...) {
        return(plotRegionPCA(object = object, ...))
    }
)

#' @rdname biocGenerics-methods
#' @export
setMethod(
    f = "plotPCA",
    signature = "RegionSetDE.setResults",
    definition = function(object, ...) {
        return(plotRegionPCA(object = object, ...))
    }
)
