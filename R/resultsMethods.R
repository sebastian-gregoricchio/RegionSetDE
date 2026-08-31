#' @title topRegions
#'
#' @description Returns the regions that respond most strongly to a contrast, optionally restricted to one region set or to one direction.
#'
#' @param results \code{RegionSetDE.results} or \code{RegionSetDE.resultsList} object.
#' @param n Numeric value with the number of regions to return. When \code{set} holds more than one name, \code{n} regions are returned for each of them. Default: \code{Inf}, every region passing the thresholds.
#' @param set Character vector with the names of the region sets to consider. Default: \code{NULL}, all of them pooled together.
#' @param contrast String with the name of the contrast to take, or its position, when \code{results} holds several of them. Default: \code{NULL}.
#' @param sortBy String with the column driving the ranking, one of \code{"FDR"}, \code{"p.value"}, \code{"log2FC"} and \code{"stat"}. Default: \code{"FDR"}.
#' @param direction String restricting the output to one direction of change, one of \code{"both"}, \code{"up"} and \code{"down"}. Default: \code{"both"}.
#' @param FDR Numeric value with an adjusted p-value cut-off applied before the ranking. Default: \code{NULL}, the threshold stored in the object.
#' @param log2FC Numeric value with an absolute log2 fold change cut-off applied before the ranking. Default: \code{NULL}, the threshold stored in the object.
#' @param level String indicating whether the regions (\code{"region"}) or the tiles (\code{"tile"}) must be ranked. Default: \code{"region"}.
#'
#' @return A data.frame with the selected rows.
#'
#' @details Ranking by \code{"log2FC"} sorts on the effect size alone and returns whatever passes the cut-offs, which on a small object is often a handful of low-count regions with a large and badly estimated fold change. Keep an FDR cut-off in place when doing so.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' # By default only the regions passing the thresholds come back
#' topRegions(results, n = 10)
#'
#' # FDR = 1 ranks everything instead, which is what you want when power is low
#' topRegions(results, n = 10, FDR = 1)
#'
#' # One set at a time, ranked rather than filtered
#' topRegions(results, n = 5, set = "promoterCpG", FDR = 1)
#'
#' # Sorting by effect size instead of significance
#' topRegions(results, n = 5, FDR = 1, sortBy = "log2FC", direction = "down")
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegions}}, \code{\link{plotVolcano}}
#'
#' @importFrom dplyr filter arrange group_by slice_head ungroup desc
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export topRegions

topRegions <-
  function(results,
           n = Inf,
           set = NULL,
           contrast = NULL,
           sortBy = "FDR",
           direction = "both",
           FDR = NULL,
           log2FC = NULL,
           level = "region") {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    results <- .pickResults(results = results, contrast = contrast)

    if (!(sortBy %in% c("FDR", "p.value", "log2FC", "stat"))) {
      stop("The 'sortBy' parameter must be one of 'FDR', 'p.value', 'log2FC', 'stat'.", call. = FALSE)
    }

    if (!(direction %in% c("both", "up", "down"))) {
      stop("The 'direction' parameter must be one of 'both', 'up', 'down'.", call. = FALSE)
    }

    resultTable <- if (level == "tile") {tileTable(results)} else {resultsTable(results)}

    if (nrow(resultTable) == 0) {
      stop("The object holds no row at the requested level.", call. = FALSE)
    }

    #-------------------------------#
    # Filter                        #
    #-------------------------------#
    FDRthreshold <- if (is.null(FDR)) {results@thresholds$FDR} else {FDR}
    log2FCthreshold <- if (is.null(log2FC)) {results@thresholds$log2FC} else {log2FC}

    if (!is.null(set)) {
      absentSets <- setdiff(set, unique(resultTable$region.set))
      if (length(absentSets) > 0) {
        stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }
      resultTable <- dplyr::filter(resultTable, .data$region.set %in% set)
    }

    resultTable <- dplyr::filter(resultTable,
                                 .data$FDR < FDRthreshold,
                                 abs(.data$log2FC) > log2FCthreshold)

    if (direction == "up") {
      resultTable <- dplyr::filter(resultTable, .data$log2FC > 0)
    } else if (direction == "down") {
      resultTable <- dplyr::filter(resultTable, .data$log2FC < 0)
    }

    if (nrow(resultTable) == 0) {
      warning("No region passes the thresholds.", call. = FALSE)
      return(resultTable)
    }

    #-------------------------------#
    # Rank                          #
    #-------------------------------#
    resultTable <- switch(sortBy,
                          "FDR" = dplyr::arrange(resultTable, .data$FDR, .data$p.value),
                          "p.value" = dplyr::arrange(resultTable, .data$p.value),
                          "log2FC" = dplyr::arrange(resultTable, dplyr::desc(abs(.data$log2FC))),
                          "stat" = dplyr::arrange(resultTable, dplyr::desc(abs(.data$stat))))

    # With several sets asked for, the top of each one is wanted, not the top of the pool
    if (is.finite(n)) {
      if (length(set) > 1) {
        resultTable <- dplyr::ungroup(dplyr::slice_head(dplyr::group_by(resultTable, .data$region.set), n = n))
      } else {
        resultTable <- dplyr::slice_head(resultTable, n = n)
      }
    }

    return(as.data.frame(resultTable))
  } # END function




#' @title resultsTable
#'
#' @description Returns the per-region table of a \code{RegionSetDE.results} object.
#'
#' @param results \code{RegionSetDE.results} object.
#'
#' @return A data.frame with one row per region.
#'
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' resultTable <- resultsTable(results)
#' head(resultTable)
#'
#' table(resultTable$region.set, resultTable$diff.status)
#'
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export
setGeneric(name = "resultsTable", def = function(results) {standardGeneric("resultsTable")})

#' @rdname resultsTable
#' @export
setMethod(f = "resultsTable",
          signature = "RegionSetDE.results",
          definition = function(results) {
            return(results@results)
          })




#' @title tileTable
#'
#' @description Returns the per-tile table of a \code{RegionSetDE.results} object, empty when the counts were not tiled.
#'
#' @param results \code{RegionSetDE.results} object.
#'
#' @return A data.frame with one row per tile.
#'
#' @examples
#' # Tiling splits each region into fixed-width windows, kept alongside the region
#' bamFile <- system.file("extdata", "ex1.bam", package = "Rsamtools")
#'
#' exampleRegions <- GenomicRanges::GRanges(
#'   seqnames = rep(c("seq1", "seq2"), each = 2),
#'   ranges = IRanges::IRanges(start = rep(c(1, 800), 2), width = 400))
#'
#' exampleRegions$setName <- rep(c("firstSet", "secondSet"), each = 2)
#'
#' exampleSets <- splitLoadRegions(exampleRegions, splitBy = "setName",
#'                                 seqlevelsStyle = NULL, verbose = FALSE)
#'
#' tiledCounts <- countReads(exampleSets, bamFiles = bamFile,
#'                           sampleNames = "example", tileWidth = 100,
#'                           verbose = FALSE)
#'
#' head(SummarizedExperiment::rowData(tiledCounts))
#'
#' \dontrun{
#' # After testing, the per-tile statistics sit behind the combined ones
#' tiledResults <- testRegions(tiledFit, contrast = c("condition", "SHR", "BN"))
#' head(tileTable(tiledResults))
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export
setGeneric(name = "tileTable", def = function(results) {standardGeneric("tileTable")})

#' @rdname tileTable
#' @export
setMethod(f = "tileTable",
          signature = "RegionSetDE.results",
          definition = function(results) {
            if (nrow(results@tiles) == 0) {
              stop("The counts were not tiled, no tile level table is available.", call. = FALSE)
            }
            return(results@tiles)
          })




#' @title resultRanges
#'
#' @description Returns the coordinates of a \code{RegionSetDE.results} object with the statistics attached as metadata columns, ready to be exported as a BED-like file.
#'
#' @param results \code{RegionSetDE.results} object.
#'
#' @return A \code{GRanges} with one element per region.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric setMethod
#' @importFrom S4Vectors mcols<- DataFrame
#'
#' @export
setGeneric(name = "resultRanges", def = function(results) {standardGeneric("resultRanges")})

#' @rdname resultRanges
#' @export
setMethod(f = "resultRanges",
          signature = "RegionSetDE.results",
          definition = function(results) {
            regionRanges <- results@regions

            # The coordinates are already in the ranges, repeating them in the metadata only makes the export noisy
            metadataColumns <- setdiff(colnames(results@results), c("seqnames", "start", "end", "width"))
            S4Vectors::mcols(regionRanges) <- S4Vectors::DataFrame(results@results[, metadataColumns, drop = FALSE])

            return(regionRanges)
          })




#' @title regionSetNames
#'
#' @description Returns the names of the region sets carried by an object, at any stage of the analysis.
#'
#' @param object \code{RegionSetDE}, \code{RegionSetDE.counts}, \code{RegionSetDE.fit} or \code{RegionSetDE.results} object.
#'
#' @return Character vector with the names of the region sets, in the order they are stored.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' regionSetNames(counts)
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{loadRegions}}, \code{\link{splitLoadRegions}}
#'
#' @export

setGeneric(name = "regionSetNames",
           def = function(object) {standardGeneric("regionSetNames")})


#' @rdname regionSetNames
#' @export

setMethod(f = "regionSetNames",
          signature = "RegionSetDE",
          definition = function(object) {
            return(names(object@regions))
          })


#' @rdname regionSetNames
#'
#' @importFrom SummarizedExperiment rowData
#'
#' @export

setMethod(f = "regionSetNames",
          signature = "RegionSetDE.counts",
          definition = function(object) {
            # After counting the membership lives in the rowData, one label per row
            return(unique(as.character(SummarizedExperiment::rowData(object)$region.set)))
          })


#' @rdname regionSetNames
#' @export

setMethod(f = "regionSetNames",
          signature = "RegionSetDE.fit",
          definition = function(object) {
            return(regionSetNames(fitCounts(object)))
          })


#' @rdname regionSetNames
#' @export

setMethod(f = "regionSetNames",
          signature = "RegionSetDE.results",
          definition = function(object) {
            return(regionSetNames(resultCounts(object)))
          })




#' @title contrastName
#'
#' @description Returns the contrast that produced a \code{RegionSetDE.results} object.
#'
#' @param results \code{RegionSetDE.results} object.
#'
#' @return A string.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#' results <- testRegions(fit, contrast = c("condition", "SHR", "BN"), verbose = FALSE)
#'
#' contrastName(results)
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export
setGeneric(name = "contrastName", def = function(results) {standardGeneric("contrastName")})

#' @rdname contrastName
#' @export
setMethod(f = "contrastName",
          signature = "RegionSetDE.results",
          definition = function(results) {
            return(results@contrast)
          })




#' @title fitCounts
#'
#' @description Returns the \code{RegionSetDE.counts} object on which a model was fitted.
#'
#' @param fit \code{RegionSetDE.fit} object.
#'
#' @return A \code{RegionSetDE.counts} object.
#'
#' @examples
#' fit <- loadExampleData("fit", verbose = FALSE)
#'
#' # The counts the model was fitted on, normalisation and filtering included
#' fitCounts(fit)
#'
#' # The underlying engine object, here an edgeR fit
#' class(fitObject(fit))
#'
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export
setGeneric(name = "fitCounts", def = function(fit) {standardGeneric("fitCounts")})

#' @rdname fitCounts
#' @export
setMethod(f = "fitCounts",
          signature = "RegionSetDE.fit",
          definition = function(fit) {
            return(fit@counts)
          })




#' @title fitObject
#'
#' @description Returns the object produced by the engine, untouched, so that any function of \code{edgeR}, \code{limma}, \code{variancePartition} or \code{DESeq2} can be applied to it.
#'
#' @param fit \code{RegionSetDE.fit} object.
#'
#' @return The engine specific fit object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export
setGeneric(name = "fitObject", def = function(fit) {standardGeneric("fitObject")})

#' @rdname fitObject
#' @export
setMethod(f = "fitObject",
          signature = "RegionSetDE.fit",
          definition = function(fit) {
            return(fit@fit$object)
          })




#' @title resultCounts
#'
#' @description Returns the \code{RegionSetDE.counts} object carried inside a result, the one the contrast was computed on.
#'
#' @param results \code{RegionSetDE.results} object.
#'
#' @return A \code{RegionSetDE.counts} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods setGeneric setMethod
#'
#' @export
setGeneric(name = "resultCounts", def = function(results) {standardGeneric("resultCounts")})

#' @rdname resultCounts
#' @export
setMethod(f = "resultCounts",
          signature = "RegionSetDE.results",
          definition = function(results) {
            if (ncol(results@counts) == 0) {
              stop("The result carries no counts, run testRegions with carryCounts = TRUE.", call. = FALSE)
            }
            return(results@counts)
          })








## ---------------------------------------------------------------------------
## RegionSetDE accessor methods for the set-level and list classes
##
## Each group goes in the file that already holds the corresponding generic.
## The full roxygen block stays on the setGeneric call; these are @rdname stubs,
## which is what makes roxygen write the exportMethods entries.
##
## ---------------------------------------------------------------------------


## === resultsTable ==========================================================

#' @rdname resultsTable
#' @export

setMethod(f = "resultsTable",
          signature = "RegionSetDE.setResults",
          definition = function(results) {
            # One row per region set rather than per region
            return(results@results)
          })


#' @rdname resultsTable
#'
#' @importFrom dplyr bind_rows mutate relocate
#'
#' @export

setMethod(f = "resultsTable",
          signature = "RegionSetDE.resultsList",
          definition = function(results) {
            # The tables of every contrast, stacked and labelled
            tableList <- lapply(results@results,
                                function(singleResult) {
                                  dplyr::relocate(
                                    dplyr::mutate(resultsTable(singleResult),
                                                  contrast = contrastName(singleResult)),
                                    contrast)
                                })

            return(dplyr::bind_rows(tableList))
          })


#' @rdname resultsTable
#'
#' @importFrom dplyr bind_rows mutate relocate
#'
#' @export

setMethod(f = "resultsTable",
          signature = "RegionSetDE.setResultsList",
          definition = function(results) {
            tableList <- lapply(results@results,
                                function(singleResult) {
                                  dplyr::relocate(
                                    dplyr::mutate(resultsTable(singleResult),
                                                  contrast = contrastName(singleResult)),
                                    contrast)
                                })

            return(dplyr::bind_rows(tableList))
          })


## === resultCounts ==========================================================

#' @rdname resultCounts
#' @export

setMethod(f = "resultCounts",
          signature = "RegionSetDE.setResults",
          definition = function(results) {
            return(results@counts)
          })


#' @rdname resultCounts
#' @export

setMethod(f = "resultCounts",
          signature = "RegionSetDE.resultsList",
          definition = function(results) {
            # Every contrast comes from the same fit, so the counts are returned once
            if (length(results@results) == 0) {
              stop("The results list is empty, there are no counts to return.", call. = FALSE)
            }

            return(resultCounts(results@results[[1]]))
          })


#' @rdname resultCounts
#' @export

setMethod(f = "resultCounts",
          signature = "RegionSetDE.setResultsList",
          definition = function(results) {
            if (length(results@results) == 0) {
              stop("The results list is empty, there are no counts to return.", call. = FALSE)
            }

            return(resultCounts(results@results[[1]]))
          })


## === resultRanges ==========================================================
## No method for the set classes: a region set is not a range, and returning the
## ranges of its members under the same name would mean something else.

#' @rdname resultRanges
#' @export

setMethod(f = "resultRanges",
          signature = "RegionSetDE.resultsList",
          definition = function(results) {
            # One GRanges per contrast, kept separate rather than stacked
            rangesList <- lapply(results@results, resultRanges)
            names(rangesList) <- vapply(results@results, contrastName, character(1))

            return(GenomicRanges::GRangesList(rangesList))
          })


## === tileTable =============================================================
## No method for the set classes: set-level results carry no tiles.

#' @rdname tileTable
#'
#' @importFrom dplyr bind_rows mutate relocate
#'
#' @export

setMethod(f = "tileTable",
          signature = "RegionSetDE.resultsList",
          definition = function(results) {
            tableList <- lapply(results@results,
                                function(singleResult) {
                                  singleTable <- tileTable(singleResult)

                                  if (is.null(singleTable) || nrow(singleTable) == 0) {
                                    return(NULL)
                                  }

                                  dplyr::relocate(
                                    dplyr::mutate(singleTable,
                                                  contrast = contrastName(singleResult)),
                                    contrast)
                                })

            return(dplyr::bind_rows(tableList))
          })


## === contrastName ==========================================================

#' @rdname contrastName
#' @export

setMethod(f = "contrastName",
          signature = "RegionSetDE.setResults",
          definition = function(results) {
            return(results@contrast)
          })


#' @rdname contrastName
#' @export

setMethod(f = "contrastName",
          signature = "RegionSetDE.resultsList",
          definition = function(results) {
            # A vector, one entry per contrast held in the list
            return(results@contrasts)
          })


#' @rdname contrastName
#' @export

setMethod(f = "contrastName",
          signature = "RegionSetDE.setResultsList",
          definition = function(results) {
            return(results@contrasts)
          })


## === regionSetNames ========================================================

#' @rdname regionSetNames
#' @export

setMethod(f = "regionSetNames",
          signature = "RegionSetDE.setResults",
          definition = function(object) {
            # Here the sets are the rows of the results table
            return(unique(as.character(object@results$region.set)))
          })


#' @rdname regionSetNames
#' @export

setMethod(f = "regionSetNames",
          signature = "RegionSetDE.resultsList",
          definition = function(object) {
            if (length(object@results) == 0) {
              stop("The results list is empty, there are no region sets to return.", call. = FALSE)
            }

            return(regionSetNames(object@results[[1]]))
          })


#' @rdname regionSetNames
#' @export

setMethod(f = "regionSetNames",
          signature = "RegionSetDE.setResultsList",
          definition = function(object) {
            if (length(object@results) == 0) {
              stop("The results list is empty, there are no region sets to return.", call. = FALSE)
            }

            return(regionSetNames(object@results[[1]]))
          })






