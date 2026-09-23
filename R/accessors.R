# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title regionRanges
#'
#' @description Returns the region sets of a \code{RegionSetDE} object as a named \code{GRangesList}, one element per set.
#'
#' @param object \code{RegionSetDE} object returned by \code{\link{loadRegions}}, \code{\link{splitLoadRegions}} or \code{\link{loadConsensusPeaks}}.
#'
#' @return A named \code{GRangesList} with the regions of every set, together with their metadata columns.
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#'
#' regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
#'                             splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)
#'
#' lengths(regionRanges(regions))
#' head(regionRanges(regions)$promoterCpG, 3)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{regionSetNames}}, \code{\link{resultRanges}}
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export regionRanges

setGeneric(name = "regionRanges", def = function(object) {standardGeneric("regionRanges")})

#' @rdname regionRanges
#' @export
setMethod(f = "regionRanges",
          signature = "RegionSetDE",
          definition = function(object) {
            return(object@regions)
          })




#' @title filteringLog
#'
#' @description Returns the record of the filters applied to a RegionSetDE object, with the number of regions before and after each step. The record travels with the regions into the counts, the fit and the results.
#'
#' @param object Any object of the package: \code{RegionSetDE}, \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results}, \code{RegionSetDE.setResults} or \code{RegionSetDE.setScores}.
#'
#' @return A data.frame with one row per filtering step, empty when no filter has been applied.
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#' exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)
#'
#' regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
#'                             splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE)
#' regions <- applyBlacklist(regions, blacklist = exclusionRegions, verbose = FALSE)
#'
#' filteringLog(regions)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{applyBlacklist}}, \code{\link{applyWhitelist}}, \code{\link{applyGreylist}}
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export filteringLog

setGeneric(name = "filteringLog", def = function(object) {standardGeneric("filteringLog")})

#' @rdname filteringLog
#' @export
setMethod(f = "filteringLog",
          signature = "RegionSetDE.provenance",
          definition = function(object) {
            return(object@filtering.log)
          })




#' @title genomeAssembly
#'
#' @description Returns the genome assembly declared for the regions of a RegionSetDE object.
#'
#' @param object Any object of the package: \code{RegionSetDE}, \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results}, \code{RegionSetDE.setResults} or \code{RegionSetDE.setScores}.
#'
#' @return A string with the assembly, or \code{NULL} when none was declared.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' genomeAssembly(counts)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadRegions}}
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export genomeAssembly

setGeneric(name = "genomeAssembly", def = function(object) {standardGeneric("genomeAssembly")})

#' @rdname genomeAssembly
#' @export
setMethod(f = "genomeAssembly",
          signature = "RegionSetDE.provenance",
          definition = function(object) {
            return(object@genome.assembly)
          })




#' @title countingLevel
#'
#' @description Tells whether the rows of an object are whole regions or tiles of a region.
#'
#' @param object \code{RegionSetDE.counts}, \code{RegionSetDE.fit} or \code{RegionSetDE.results} object.
#'
#' @return Either \code{"region"} or \code{"tile"}.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' countingLevel(counts)
#'
#' fit <- loadExampleData("fit", verbose = FALSE)
#' countingLevel(fit)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{countReads}}, \code{\link{tileTable}}
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export countingLevel

setGeneric(name = "countingLevel", def = function(object) {standardGeneric("countingLevel")})

#' @rdname countingLevel
#' @export
setMethod(f = "countingLevel",
          signature = "RegionSetDE.counts",
          definition = function(object) {
            return(object@counting.level)
          })

#' @rdname countingLevel
#' @export
setMethod(f = "countingLevel",
          signature = "RegionSetDE.fit",
          definition = function(object) {
            return(object@counting.level)
          })

#' @rdname countingLevel
#' @export
setMethod(f = "countingLevel",
          signature = "RegionSetDE.results",
          definition = function(object) {
            return(object@counting.level)
          })




#' @title dispersionInfo
#'
#' @description Returns the summary of the dispersion estimated, or supplied, when a model was fitted: the common value, whether it was held fixed, whether the design had no replicates, and where the null rows came from.
#'
#' @param fit \code{RegionSetDE.fit} object returned by \code{\link{fitRegions}}.
#'
#' @return A named list. Which elements it holds depends on the engine and on how the dispersion was obtained.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' names(dispersionInfo(fit))
#' dispersionInfo(fit)$common
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{fitRegions}}, \code{\link{estimateNullDispersion}}, \code{\link{checkNullCalibration}}
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export dispersionInfo

setGeneric(name = "dispersionInfo", def = function(fit) {standardGeneric("dispersionInfo")})

#' @rdname dispersionInfo
#' @export
setMethod(f = "dispersionInfo",
          signature = "RegionSetDE.fit",
          definition = function(fit) {
            return(fit@dispersion)
          })
