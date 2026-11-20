# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title regionRanges
#'
#' @description Returns the region sets of a \code{RegionSetDE} object as a named \code{GRangesList}, one element per set, or the tested coordinates of a per-region result.
#'
#' @param object \code{RegionSetDE} object returned by \code{\link{loadRegions}}, \code{\link{splitLoadRegions}} or \code{\link{loadConsensusPeaks}}, or a \code{RegionSetDE.results} object.
#'
#' @return For a \code{RegionSetDE} object, a named \code{GRangesList} with the regions of every set, together with their metadata columns. For a \code{RegionSetDE.results} object, a \code{GRanges} with the coordinates of the rows of the results table, in the same order and without the statistics, which \code{\link{resultRanges}} attaches.
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#'
#' regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
#'     splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE
#' )
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

setGeneric(name = "regionRanges", def = function(object) {
    standardGeneric("regionRanges")
})

#' @rdname regionRanges
#' @export
setMethod(
    f = "regionRanges",
    signature = "RegionSetDE",
    definition = function(object) {
        return(object@regions)
    }
)


#' @title filteringLog
#'
#' @description Returns the record of the filters applied to a RegionSetDE object, with the number of regions before and after each step. The record travels with the regions into the counts, the fit and the results.
#'
#' @param object Any object of the package: \code{RegionSetDE}, \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results}, \code{RegionSetDE.setResults} or \code{RegionSetDE.setScores}.
#'
#' @return A data.frame with one row per filtering step, empty when no filter has been applied. The object stores it as a \code{DataFrame}.
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#' exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)
#'
#' regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
#'     splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE
#' )
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

setGeneric(name = "filteringLog", def = function(object) {
    standardGeneric("filteringLog")
})

#' @rdname filteringLog
#' @export
setMethod(
    f = "filteringLog",
    signature = "RegionSetDE.provenance",
    definition = function(object) {
        return(.asTable(object@filtering.log))
    }
)


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

setGeneric(name = "genomeAssembly", def = function(object) {
    standardGeneric("genomeAssembly")
})

#' @rdname genomeAssembly
#' @export
setMethod(
    f = "genomeAssembly",
    signature = "RegionSetDE.provenance",
    definition = function(object) {
        return(object@genome.assembly)
    }
)


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

setGeneric(name = "countingLevel", def = function(object) {
    standardGeneric("countingLevel")
})

#' @rdname countingLevel
#' @export
setMethod(
    f = "countingLevel",
    signature = "RegionSetDE.counts",
    definition = function(object) {
        return(object@counting.level)
    }
)

#' @rdname countingLevel
#' @export
setMethod(
    f = "countingLevel",
    signature = "RegionSetDE.fit",
    definition = function(object) {
        return(object@counting.level)
    }
)

#' @rdname countingLevel
#' @export
setMethod(
    f = "countingLevel",
    signature = "RegionSetDE.results",
    definition = function(object) {
        return(object@counting.level)
    }
)


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

setGeneric(name = "dispersionInfo", def = function(fit) {
    standardGeneric("dispersionInfo")
})

#' @rdname dispersionInfo
#' @export
setMethod(
    f = "dispersionInfo",
    signature = "RegionSetDE.fit",
    definition = function(fit) {
        return(fit@dispersion)
    }
)


#' @title .asTable
#'
#' @description Turns a table stored as a \code{DataFrame} into the data.frame the accessors return, keeping the column names exactly as they are.
#'
#' @param x \code{DataFrame} or data.frame.
#'
#' @return A data.frame.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.asTable <-
    function(x) {
        return(as.data.frame(x, optional = TRUE))
    } # END function


##########################################################################################
###    PROVENANCE, every object of the package
##########################################################################################


#' @title Provenance accessors
#'
#' @description Return what every object of the package carries about the origin of its regions: the lists of regions applied to them, the chromosome naming style, and the parameters of every step that built the object. The record travels from the regions into the counts, the fit and the results, so the same accessors answer at every step.
#'
#' @param object Any object of the package: \code{RegionSetDE}, \code{RegionSetDE.counts}, \code{RegionSetDE.fit}, \code{RegionSetDE.results}, \code{RegionSetDE.setResults}, \code{RegionSetDE.setScores} or \code{RegionSetDE.geneResults}.
#'
#' @return \code{blacklistRanges}, \code{greylistRanges} and \code{whitelistRanges} return a \code{GRanges} with the list applied by \code{\link{applyBlacklist}}, \code{\link{applyGreylist}} or \code{\link{applyWhitelist}}, or by the matching argument of \code{\link{countReads}} and \code{\link{loadConsensusPeaks}}, and \code{NULL} when none was applied. \code{chromosomeStyle} returns a string with the naming style of the chromosomes, such as \code{"UCSC"}, or \code{NA} when the names were kept as they came. \code{parameterLog} returns a named list with one element per step, named after the function, holding the arguments it was called with.
#'
#' @examples
#' regionTable <- loadExampleData("regions", verbose = FALSE)
#' exclusionRegions <- loadExampleData("exclusionRegions", verbose = FALSE)
#'
#' regions <- splitLoadRegions(GenomicRanges::makeGRangesFromDataFrame(regionTable, keep.extra.columns = TRUE),
#'     splitBy = "setName", genomeAssembly = "rn4", verbose = FALSE
#' )
#' regions <- applyBlacklist(regions, blacklist = exclusionRegions, verbose = FALSE)
#'
#' blacklistRanges(regions)
#' is.null(greylistRanges(regions))
#' chromosomeStyle(regions)
#'
#' # The steps that built the example counts
#' counts <- loadExampleData("counts", verbose = FALSE)
#' names(parameterLog(counts))
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{filteringLog}}, \code{\link{genomeAssembly}}
#'
#' @rdname provenanceAccessors
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export blacklistRanges

setGeneric(name = "blacklistRanges", def = function(object) {
    standardGeneric("blacklistRanges")
})

#' @rdname provenanceAccessors
#' @export
setMethod(
    f = "blacklistRanges",
    signature = "RegionSetDE.provenance",
    definition = function(object) {
        return(object@blacklist)
    }
)


#' @rdname provenanceAccessors
#' @export greylistRanges
setGeneric(name = "greylistRanges", def = function(object) {
    standardGeneric("greylistRanges")
})

#' @rdname provenanceAccessors
#' @export
setMethod(
    f = "greylistRanges",
    signature = "RegionSetDE.provenance",
    definition = function(object) {
        return(object@greylist)
    }
)


#' @rdname provenanceAccessors
#' @export whitelistRanges
setGeneric(name = "whitelistRanges", def = function(object) {
    standardGeneric("whitelistRanges")
})

#' @rdname provenanceAccessors
#' @export
setMethod(
    f = "whitelistRanges",
    signature = "RegionSetDE.provenance",
    definition = function(object) {
        return(object@whitelist)
    }
)


#' @rdname provenanceAccessors
#' @export chromosomeStyle
setGeneric(name = "chromosomeStyle", def = function(object) {
    standardGeneric("chromosomeStyle")
})

#' @rdname provenanceAccessors
#' @export
setMethod(
    f = "chromosomeStyle",
    signature = "RegionSetDE.provenance",
    definition = function(object) {
        return(object@seqlevels.style)
    }
)


#' @rdname provenanceAccessors
#' @export parameterLog
setGeneric(name = "parameterLog", def = function(object) {
    standardGeneric("parameterLog")
})

#' @rdname provenanceAccessors
#' @export
setMethod(
    f = "parameterLog",
    signature = "RegionSetDE.provenance",
    definition = function(object) {
        return(object@parameters)
    }
)


##########################################################################################
###    REGIONS of a result
##########################################################################################


#' @rdname regionRanges
#' @export
setMethod(
    f = "regionRanges",
    signature = "RegionSetDE.results",
    definition = function(object) {
        # The coordinates as they were tested, without the statistics resultRanges() attaches to them
        return(object@regions)
    }
)


##########################################################################################
###    UNIVERSE
##########################################################################################


#' @title Universe accessors
#'
#' @description Return the parts of a \code{RegionSetDE.universe} object, the rows every region set is compared against in a competitive test. \code{\link{setUniverse}} reads the universe out of a fit or a set level result.
#'
#' @param universe \code{RegionSetDE.universe} object, returned by \code{\link{makeSetUniverse}} or by \code{\link{setUniverse}}.
#'
#' @return \code{universeIndex} returns a named list with one element per region set, the positions of the rows forming its universe. \code{universeType} returns a string with the way the universe was built, for instance \code{"otherSets"}, \code{"supplied"}, \code{"allRows"} or \code{"pairedSet"}. \code{universeMatching} returns the covariates the comparison rows were matched on, empty when they were not matched. \code{comparisonSets} returns the names of the sets the comparison rows came from. \code{universeDiagnostics} returns a data.frame with, for every set, the median width and abundance of the set and of the rows it is compared against. \code{universeRows} returns the number of rows of the object the positions refer to.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' universe <- setUniverse(fit)
#'
#' universeType(universe)
#' universeMatching(universe)
#' comparisonSets(universe)
#' lengths(universeIndex(universe))
#' universeDiagnostics(universe)
#' universeRows(universe)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{makeSetUniverse}}, \code{\link{plotUniverseMatching}}
#'
#' @rdname universeAccessors
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export universeIndex

setGeneric(name = "universeIndex", def = function(universe) {
    standardGeneric("universeIndex")
})

#' @rdname universeAccessors
#' @export
setMethod(
    f = "universeIndex",
    signature = "RegionSetDE.universe",
    definition = function(universe) {
        return(universe@index)
    }
)


#' @rdname universeAccessors
#' @export universeType
setGeneric(name = "universeType", def = function(universe) {
    standardGeneric("universeType")
})

#' @rdname universeAccessors
#' @export
setMethod(
    f = "universeType",
    signature = "RegionSetDE.universe",
    definition = function(universe) {
        return(universe@type)
    }
)


#' @rdname universeAccessors
#' @export universeMatching
setGeneric(name = "universeMatching", def = function(universe) {
    standardGeneric("universeMatching")
})

#' @rdname universeAccessors
#' @export
setMethod(
    f = "universeMatching",
    signature = "RegionSetDE.universe",
    definition = function(universe) {
        return(universe@matching)
    }
)


#' @rdname universeAccessors
#' @export comparisonSets
setGeneric(name = "comparisonSets", def = function(universe) {
    standardGeneric("comparisonSets")
})

#' @rdname universeAccessors
#' @export
setMethod(
    f = "comparisonSets",
    signature = "RegionSetDE.universe",
    definition = function(universe) {
        return(universe@comparison.sets)
    }
)


#' @rdname universeAccessors
#' @export universeDiagnostics
setGeneric(name = "universeDiagnostics", def = function(universe) {
    standardGeneric("universeDiagnostics")
})

#' @rdname universeAccessors
#' @export
setMethod(
    f = "universeDiagnostics",
    signature = "RegionSetDE.universe",
    definition = function(universe) {
        return(.asTable(universe@diagnostics))
    }
)


#' @rdname universeAccessors
#' @export universeRows
setGeneric(name = "universeRows", def = function(universe) {
    standardGeneric("universeRows")
})

#' @rdname universeAccessors
#' @export
setMethod(
    f = "universeRows",
    signature = "RegionSetDE.universe",
    definition = function(universe) {
        return(universe@n.rows)
    }
)


##########################################################################################
###    FIT
##########################################################################################


#' @title Fit accessors
#'
#' @description Return the parts of a \code{RegionSetDE.fit} object: what the engine returned, the engine itself, the design and its formula, the blocking, the samples, and the universe built for the set level tests. \code{fitEngine} and \code{setUniverse} answer for the results as well, since both travel from the fit into them.
#'
#' @param object \code{RegionSetDE.fit} object returned by \code{\link{fitRegions}}. \code{fitEngine} also takes a \code{RegionSetDE.results}, \code{RegionSetDE.setResults} or \code{RegionSetDE.geneResults} object, and \code{setUniverse} a \code{RegionSetDE.setResults} object.
#'
#' @return \code{fitComponents} returns a named list with everything the engine produced: \code{object}, the fit itself, which \code{\link{fitObject}} returns alone, and the auxiliary objects needed to test it, such as the voom weights or the sample table. \code{fitEngine} returns a string, one among \code{"edgeR"}, \code{"voom"}, \code{"dream"}, \code{"deseq2"} and \code{"limma"}. \code{designMatrix} returns the matrix of the fixed effects, one row per sample, as \code{design} does. \code{designFormula} returns the formula of the design, or \code{NULL} when the design was given as a matrix. \code{blockingInfo} returns a named list with the blocking variable, the random effect formula and the consensus correlation, empty when there is no blocking. \code{fitSamples} returns the names of the samples the model was fitted on. \code{setUniverse} returns the \code{RegionSetDE.universe} object with the rows every set is compared against, read with the \link[=universeIndex]{universe accessors}.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' fitEngine(fit)
#' designFormula(fit)
#' head(designMatrix(fit))
#' fitSamples(fit)
#' names(fitComponents(fit))
#' blockingInfo(fit)
#' setUniverse(fit)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{fitObject}}, \code{\link{fitCounts}}, \code{\link{dispersionInfo}}
#'
#' @rdname fitAccessors
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export fitComponents

setGeneric(name = "fitComponents", def = function(object) {
    standardGeneric("fitComponents")
})

#' @rdname fitAccessors
#' @export
setMethod(
    f = "fitComponents",
    signature = "RegionSetDE.fit",
    definition = function(object) {
        return(object@fit)
    }
)


#' @rdname fitAccessors
#' @export fitEngine
setGeneric(name = "fitEngine", def = function(object) {
    standardGeneric("fitEngine")
})

#' @rdname fitAccessors
#' @export
setMethod(
    f = "fitEngine",
    signature = "RegionSetDE.fit",
    definition = function(object) {
        return(object@engine)
    }
)

#' @rdname fitAccessors
#' @export
setMethod(
    f = "fitEngine",
    signature = "RegionSetDE.results",
    definition = function(object) {
        return(object@engine)
    }
)

#' @rdname fitAccessors
#' @export
setMethod(
    f = "fitEngine",
    signature = "RegionSetDE.setResults",
    definition = function(object) {
        return(object@engine)
    }
)

#' @rdname fitAccessors
#' @export
setMethod(
    f = "fitEngine",
    signature = "RegionSetDE.geneResults",
    definition = function(object) {
        return(object@engine)
    }
)


#' @rdname fitAccessors
#' @export designMatrix
setGeneric(name = "designMatrix", def = function(object) {
    standardGeneric("designMatrix")
})

#' @rdname fitAccessors
#' @export
setMethod(
    f = "designMatrix",
    signature = "RegionSetDE.fit",
    definition = function(object) {
        return(object@design)
    }
)


#' @rdname fitAccessors
#' @export designFormula
setGeneric(name = "designFormula", def = function(object) {
    standardGeneric("designFormula")
})

#' @rdname fitAccessors
#' @export
setMethod(
    f = "designFormula",
    signature = "RegionSetDE.fit",
    definition = function(object) {
        # A design given as a matrix has no formula behind it
        return(object@design.formula$formula)
    }
)


#' @rdname fitAccessors
#' @export blockingInfo
setGeneric(name = "blockingInfo", def = function(object) {
    standardGeneric("blockingInfo")
})

#' @rdname fitAccessors
#' @export
setMethod(
    f = "blockingInfo",
    signature = "RegionSetDE.fit",
    definition = function(object) {
        return(object@blocking)
    }
)


#' @rdname fitAccessors
#' @export fitSamples
setGeneric(name = "fitSamples", def = function(object) {
    standardGeneric("fitSamples")
})

#' @rdname fitAccessors
#' @export
setMethod(
    f = "fitSamples",
    signature = "RegionSetDE.fit",
    definition = function(object) {
        return(object@samples)
    }
)


#' @rdname fitAccessors
#' @export setUniverse
setGeneric(name = "setUniverse", def = function(object) {
    standardGeneric("setUniverse")
})

#' @rdname fitAccessors
#' @export
setMethod(
    f = "setUniverse",
    signature = "RegionSetDE.fit",
    definition = function(object) {
        return(object@universe)
    }
)

#' @rdname fitAccessors
#' @export
setMethod(
    f = "setUniverse",
    signature = "RegionSetDE.setResults",
    definition = function(object) {
        return(object@universe)
    }
)


##########################################################################################
###    CONTRAST of a result
##########################################################################################


#' @title Contrast accessors
#'
#' @description Return how a result was obtained: the coefficients of the contrast, the two groups it compares, the cut-offs its labels were set at, and, for tiled counts, how the tiles were combined into their region.
#'
#' @param results \code{RegionSetDE.results} object returned by \code{\link{testRegions}}. \code{contrastGroups} also takes a \code{RegionSetDE.setResults} object, and \code{thresholdInfo} a \code{RegionSetDE.setResults} or \code{RegionSetDE.geneResults} object.
#'
#' @return \code{contrastVector} returns the coefficients of the contrast over the columns of the design. \code{contrastGroups} returns a list with the \code{column} of the \code{colData} the contrast separates and the two \code{groups} it compares, the first being the level the fold change is positive for, and for a per-region result \code{n.samples}, the number of samples in each; the list is empty when the contrast is not a difference between two levels of one variable. \code{thresholdInfo} returns a named list with the cut-offs, such as \code{FDR} and \code{log2FC}, and the correction method, \code{adjust.method}. \code{tileCombination} returns a list with the \code{method} that combined the tiles of each region, empty when the counts were not tiled.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' contrastVector(results)
#' contrastGroups(results)
#' thresholdInfo(results)
#' tileCombination(results)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{contrastName}}, \code{\link{contrastInfo}}, \code{\link{updateThresholds}}
#'
#' @rdname contrastAccessors
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export contrastVector

setGeneric(name = "contrastVector", def = function(results) {
    standardGeneric("contrastVector")
})

#' @rdname contrastAccessors
#' @export
setMethod(
    f = "contrastVector",
    signature = "RegionSetDE.results",
    definition = function(results) {
        return(results@contrast.vector)
    }
)


#' @rdname contrastAccessors
#' @export contrastGroups
setGeneric(name = "contrastGroups", def = function(results) {
    standardGeneric("contrastGroups")
})

#' @rdname contrastAccessors
#' @export
setMethod(
    f = "contrastGroups",
    signature = "RegionSetDE.results",
    definition = function(results) {
        return(results@contrast.groups)
    }
)

#' @rdname contrastAccessors
#' @export
setMethod(
    f = "contrastGroups",
    signature = "RegionSetDE.setResults",
    definition = function(results) {
        return(results@contrast.groups)
    }
)


#' @rdname contrastAccessors
#' @export thresholdInfo
setGeneric(name = "thresholdInfo", def = function(results) {
    standardGeneric("thresholdInfo")
})

#' @rdname contrastAccessors
#' @export
setMethod(
    f = "thresholdInfo",
    signature = "RegionSetDE.results",
    definition = function(results) {
        return(results@thresholds)
    }
)

#' @rdname contrastAccessors
#' @export
setMethod(
    f = "thresholdInfo",
    signature = "RegionSetDE.setResults",
    definition = function(results) {
        return(results@thresholds)
    }
)

#' @rdname contrastAccessors
#' @export
setMethod(
    f = "thresholdInfo",
    signature = "RegionSetDE.geneResults",
    definition = function(results) {
        return(results@thresholds)
    }
)


#' @rdname contrastAccessors
#' @export tileCombination
setGeneric(name = "tileCombination", def = function(results) {
    standardGeneric("tileCombination")
})

#' @rdname contrastAccessors
#' @export
setMethod(
    f = "tileCombination",
    signature = "RegionSetDE.results",
    definition = function(results) {
        return(results@combination)
    }
)


##########################################################################################
###    SET LEVEL RESULTS
##########################################################################################


#' @title Set level result accessors
#'
#' @description Return the parts of a \code{RegionSetDE.setResults} object beyond its table: the per-region statistics the set level test was computed on, the kind of test and the methods that were run.
#'
#' @param results \code{RegionSetDE.setResults} object returned by \code{\link{testRegionSets}} or \code{\link{testSetContrast}}.
#'
#' @return \code{regionStats} returns a data.frame with one row per region, holding the statistics the set level test read. \code{setTestType} returns \code{"set"} for \code{\link{testRegionSets}} and \code{"setContrast"} for \code{\link{testSetContrast}}. \code{setTestMethods} returns the names of the tests that were run, for instance \code{"camera"} and \code{"fry"}.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' setRes <- testRegionSets(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' setTestType(setRes)
#' setTestMethods(setRes)
#' head(regionStats(setRes), 3)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{resultsTable}}, \code{\link{setUniverse}}, \code{\link{contrastGroups}}
#'
#' @rdname setResultsAccessors
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export regionStats

setGeneric(name = "regionStats", def = function(results) {
    standardGeneric("regionStats")
})

#' @rdname setResultsAccessors
#' @export
setMethod(
    f = "regionStats",
    signature = "RegionSetDE.setResults",
    definition = function(results) {
        return(.asTable(results@regionStats))
    }
)


#' @rdname setResultsAccessors
#' @export setTestType
setGeneric(name = "setTestType", def = function(results) {
    standardGeneric("setTestType")
})

#' @rdname setResultsAccessors
#' @export
setMethod(
    f = "setTestType",
    signature = "RegionSetDE.setResults",
    definition = function(results) {
        return(results@test)
    }
)


#' @rdname setResultsAccessors
#' @export setTestMethods
setGeneric(name = "setTestMethods", def = function(results) {
    standardGeneric("setTestMethods")
})

#' @rdname setResultsAccessors
#' @export
setMethod(
    f = "setTestMethods",
    signature = "RegionSetDE.setResults",
    definition = function(results) {
        return(results@methods)
    }
)


##########################################################################################
###    SET SCORES
##########################################################################################


#' @title Set score accessors
#'
#' @description Return the parts of a \code{RegionSetDE.setScores} object beyond the per-library scores of \code{\link{scoreTable}}: the comparisons between the sets, the samples behind the libraries, and how the signal of a set was summarised.
#'
#' @param scores \code{RegionSetDE.setScores} object returned by \code{\link{scoreRegionSets}}.
#'
#' @return \code{scoreComparisons} returns a data.frame with one row per pair of sets, as \code{resultsTable} does. \code{scoreSamples} returns a data.frame with the \code{colData} of the counts the scores came from. \code{scoreReference} returns what the signal of a set was divided by, one among \code{"background"}, \code{"regions"} and \code{"none"}. \code{scoreAssay} returns the name of the assay the signal was read from. \code{scoreSummary} returns \code{"mean"} or \code{"median"}, the summary of the regions of a set. \code{scorePerBasepair} returns \code{TRUE} when the signal was divided by the width of the regions before being summarised.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' setScores <- scoreRegionSets(counts, verbose = FALSE)
#'
#' scoreComparisons(setScores)
#' head(scoreSamples(setScores), 3)
#' scoreReference(setScores)
#' scoreAssay(setScores)
#' scoreSummary(setScores)
#' scorePerBasepair(setScores)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{scoreTable}}, \code{\link{scoreRegionSets}}
#'
#' @rdname scoreAccessors
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export scoreComparisons

setGeneric(name = "scoreComparisons", def = function(scores) {
    standardGeneric("scoreComparisons")
})

#' @rdname scoreAccessors
#' @export
setMethod(
    f = "scoreComparisons",
    signature = "RegionSetDE.setScores",
    definition = function(scores) {
        return(.asTable(scores@comparisons))
    }
)


#' @rdname scoreAccessors
#' @export scoreSamples
setGeneric(name = "scoreSamples", def = function(scores) {
    standardGeneric("scoreSamples")
})

#' @rdname scoreAccessors
#' @export
setMethod(
    f = "scoreSamples",
    signature = "RegionSetDE.setScores",
    definition = function(scores) {
        return(.asTable(scores@sample.metadata))
    }
)


#' @rdname scoreAccessors
#' @export scoreReference
setGeneric(name = "scoreReference", def = function(scores) {
    standardGeneric("scoreReference")
})

#' @rdname scoreAccessors
#' @export
setMethod(
    f = "scoreReference",
    signature = "RegionSetDE.setScores",
    definition = function(scores) {
        return(scores@reference)
    }
)


#' @rdname scoreAccessors
#' @export scoreAssay
setGeneric(name = "scoreAssay", def = function(scores) {
    standardGeneric("scoreAssay")
})

#' @rdname scoreAccessors
#' @export
setMethod(
    f = "scoreAssay",
    signature = "RegionSetDE.setScores",
    definition = function(scores) {
        return(scores@assay)
    }
)


#' @rdname scoreAccessors
#' @export scoreSummary
setGeneric(name = "scoreSummary", def = function(scores) {
    standardGeneric("scoreSummary")
})

#' @rdname scoreAccessors
#' @export
setMethod(
    f = "scoreSummary",
    signature = "RegionSetDE.setScores",
    definition = function(scores) {
        return(scores@summary)
    }
)


#' @rdname scoreAccessors
#' @export scorePerBasepair
setGeneric(name = "scorePerBasepair", def = function(scores) {
    standardGeneric("scorePerBasepair")
})

#' @rdname scoreAccessors
#' @export
setMethod(
    f = "scorePerBasepair",
    signature = "RegionSetDE.setScores",
    definition = function(scores) {
        return(scores@per.basepair)
    }
)


##########################################################################################
###    GENE RESULTS
##########################################################################################


#' @title Gene result accessors
#'
#' @description Return how a \code{RegionSetDE.geneResults} object was built: the aggregation of the regions of a gene and what weighted them. The table itself comes from \code{\link{resultsTable}}, the links from \code{\link{regionLinks}}.
#'
#' @param results \code{RegionSetDE.geneResults} object returned by \code{\link{aggregateToGenes}}.
#'
#' @return \code{geneMethod} returns the aggregations that were run, \code{"simes"}, \code{"potential"} or both. \code{geneWeighting} returns what weighted the regions of a gene, one among \code{"link.weight"}, \code{"distance"} and \code{"none"}.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' # Links of your own, four regions per gene
#' regionTable <- as.data.frame(SummarizedExperiment::rowRanges(fitCounts(fit)))[, 1:3]
#' geneNumber <- ceiling(nrow(regionTable) / 4)
#'
#' annotatedFit <- annotateRegions(fit,
#'     linkTable = data.frame(regionTable,
#'         gene.id = rep(paste0("gene", seq_len(geneNumber)),
#'             each = 4
#'         )[seq_len(nrow(regionTable))]
#'     ),
#'     verbose = FALSE
#' )
#'
#' results <- testRegions(annotatedFit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#' geneResults <- aggregateToGenes(results, method = c("simes", "potential"), verbose = FALSE)
#'
#' geneMethod(geneResults)
#' geneWeighting(geneResults)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{aggregateToGenes}}, \code{\link{resultsTable}}, \code{\link{regionLinks}}
#'
#' @rdname geneResultsAccessors
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export geneMethod

setGeneric(name = "geneMethod", def = function(results) {
    standardGeneric("geneMethod")
})

#' @rdname geneResultsAccessors
#' @export
setMethod(
    f = "geneMethod",
    signature = "RegionSetDE.geneResults",
    definition = function(results) {
        return(results@method)
    }
)


#' @rdname geneResultsAccessors
#' @export geneWeighting
setGeneric(name = "geneWeighting", def = function(results) {
    standardGeneric("geneWeighting")
})

#' @rdname geneResultsAccessors
#' @export
setMethod(
    f = "geneWeighting",
    signature = "RegionSetDE.geneResults",
    definition = function(results) {
        return(results@weighting)
    }
)


##########################################################################################
###    LISTS of results
##########################################################################################


#' @rdname resultsList-accessors
#'
#' @importFrom BiocGenerics as.list
#'
#' @export
setMethod(
    f = "as.list",
    signature = "RegionSetDE.resultsList",
    definition = function(x, ...) {
        return(x@results)
    }
)

#' @rdname resultsList-accessors
#' @export
setMethod(
    f = "as.list",
    signature = "RegionSetDE.setResultsList",
    definition = function(x, ...) {
        return(x@results)
    }
)


##########################################################################################
###    INTERNAL ACCESSORS, for the slots that may be empty
##########################################################################################


#' @title .storedCounts
#'
#' @description Internal accessor of the counts carried by a fit or a result. Unlike \code{\link{resultCounts}} it does not stop when the result was built without them, so that the functions of the package can check for them first.
#'
#' @param object \code{RegionSetDE.fit}, \code{RegionSetDE.results} or \code{RegionSetDE.setResults} object.
#'
#' @return A \code{RegionSetDE.counts} object, with no row and no column when no counts were carried.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.storedCounts <-
    function(object) {
        return(object@counts)
    } # END function


#' @title .storedTiles
#'
#' @description Internal accessor of the tile level table of a per-region result. Unlike \code{\link{tileTable}} it does not stop when the counts were not tiled.
#'
#' @param results \code{RegionSetDE.results} object.
#'
#' @return A data.frame, empty when the counts were not tiled.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.storedTiles <-
    function(results) {
        return(.asTable(results@tiles))
    } # END function


#' @title .storedConsensus
#'
#' @description Internal accessor of the consensus data of a \code{RegionSetDE} object. Unlike \code{\link{consensusData}} it does not stop when the regions were not built from peaks, and it answers an empty list for any other object.
#'
#' @param object Any object.
#'
#' @return The list returned by \code{\link{consensusData}}, or an empty list.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods is
#'
#' @keywords internal

.storedConsensus <-
    function(object) {
        if (!methods::is(object, "RegionSetDE")) {
            return(list())
        }
        return(object@consensus)
    } # END function


##########################################################################################
###    INTERNAL SETTERS
##########################################################################################


#' @title Internal setters
#'
#' @description Replacement functions used inside the package to update the slots of an object, so that \code{parameterLog(object)$step <- value} reads as the accessor it mirrors. They are not exported: the provenance of an object is written by the functions of the package only. Every one of them goes through \code{.replaceSlots}, which checks the class of the value and turns a data.frame into the \code{DataFrame} its slot holds.
#'
#' @param object An object of one of the classes of the package.
#' @param value The new content of the slot.
#'
#' @return The object with the slot replaced.
#'
#' @author Sebastian Gregoricchio
#'
#' @noRd

`parameterLog<-` <- function(object, value) {
    return(.replaceSlots(object, parameters = value))
}

`filteringLog<-` <- function(object, value) {
    return(.replaceSlots(object, filtering.log = value))
}

`regionRanges<-` <- function(object, value) {
    return(.replaceSlots(object, regions = value))
}

`blacklistRanges<-` <- function(object, value) {
    return(.replaceSlots(object, blacklist = value))
}

`greylistRanges<-` <- function(object, value) {
    return(.replaceSlots(object, greylist = value))
}

`whitelistRanges<-` <- function(object, value) {
    return(.replaceSlots(object, whitelist = value))
}

`consensusData<-` <- function(object, value) {
    return(.replaceSlots(object, consensus = value))
}

`universeIndex<-` <- function(object, value) {
    return(.replaceSlots(object, index = value))
}

`universeRows<-` <- function(object, value) {
    return(.replaceSlots(object, n.rows = value))
}

`thresholdInfo<-` <- function(object, value) {
    return(.replaceSlots(object, thresholds = value))
}

`resultsTable<-` <- function(object, value) {
    return(.replaceSlots(object, results = value))
}
