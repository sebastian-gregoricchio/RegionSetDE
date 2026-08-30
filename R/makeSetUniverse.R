#' @title makeSetUniverse
#'
#' @description Builds, for every region set, the universe of its competitive test: the set itself together with the rows it will be compared against. The comparison rows come from the other sets of the object, or from an index of your own, and they can be matched to each set on width and on baseline abundance so that the comparison is not driven by the sets simply being made of different kinds of intervals.
#'
#' Calling this is optional. \code{\link{fitRegions}} builds a matched universe and keeps it in the fit, and \code{matchOn} and \code{universeRatio} cover the common adjustments there. Come here to reach \code{strata} and \code{seed}, or to supply an index of your own through \code{type = "supplied"}.
#'
#' @param object \code{RegionSetDE.fit} or \code{RegionSetDE.counts} object.
#' @param type String with the source of the comparison rows, either \code{"otherSets"} (the rows of every other set) or \code{"supplied"}. Default: \code{"otherSets"}.
#' @param match Character vector with the covariates the comparison rows are matched on, among \code{"width"} and \code{"abundance"}. An empty vector takes every other row as it is. Default: \code{c("width", "abundance")}.
#' @param ratio Numeric value with the number of comparison rows drawn per region of the set, when matching. Default: \code{5}.
#' @param strata Numeric value with the number of strata used per covariate. Default: \code{5}.
#' @param regionSets Character vector with the names of the sets to build a universe for. Default: \code{NULL}, all of them.
#' @param index List with one vector of row positions per set, holding the comparison rows. Only for \code{type = "supplied"}. Default: \code{NULL}.
#' @param seed Numeric value with the seed of the sampling, so that a universe can be rebuilt identically. Default: \code{42}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.universe} object.
#'
#' @details Whether a set responds "more than the rest" depends entirely on what the rest is, so the choice is worth making explicitly rather than inheriting it from whatever happened to be loaded. The comparison rows are drawn from the other sets in the object, which asks whether a set responds differently from the others it was loaded alongside. That is usually the interesting claim when the sets were chosen to be compared with each other, and it is the only claim available when the object holds nothing else.
#'
#' Matching on width and abundance is what keeps the answer from being about the intervals rather than the biology. A set of 40 kb domains has more reads per region than a set of 400 bp promoter windows, and more reads mean a tighter fold change estimate, so an unmatched competitive test can separate the two sets on precision alone. The rows are binned on the quantiles of each covariate and the comparison is drawn within the bins, at \code{ratio} rows per region of the set. When a stratum runs out of eligible rows the whole stratum is taken and the shortfall shows up in the diagnostics, where the medians of the set and of its comparison should sit close together.
#'
#' The abundance used for the matching is the width-adjusted one, the same quantity \code{\link{filterRegions}} thresholds on, so a universe matched on abundance is also matched on signal density rather than on total signal.
#'
#' A single region set leaves nothing to compare against and this function stops. Either load the sets that make the comparison interesting, or add the genome bins as a set of their own before counting.
#'
#' @examples
#' \dontrun{
#' # Not needed for the usual case, fitRegions does this itself
#' fit <- fitRegions(counts, design = ~ replicate + condition)
#'
#' # Finer control over the matching
#' universeObject <- makeSetUniverse(fit, match = c("width", "abundance"), ratio = 10, strata = 8)
#' setRes <- testRegionSets(fit, contrast = "conditionCOMBO", universe = universeObject)
#'
#' # An index of your own
#' universeObject <- makeSetUniverse(fit, type = "supplied", index = list(enhancers = myRows))
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{fitRegions}}, \code{\link{testRegionSets}}, \code{\link{plotUniverseMatching}}, and \code{\link{countBackground}}, which is a different thing entirely: it counts genome bins for the normalisation and has nothing to do with the comparison built here.
#'
#' @importFrom SummarizedExperiment assay colData rowData rowRanges
#' @importFrom BiocGenerics width
#' @importFrom stats median
#' @importFrom dplyr filter mutate
#' @importFrom rlang .data
#' @importFrom methods is new
#'
#' @export makeSetUniverse

makeSetUniverse <-
  function(object,
           type = "otherSets",
           match = c("width", "abundance"),
           ratio = 5,
           strata = 5,
           regionSets = NULL,
           index = NULL,
           seed = 42,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (methods::is(object, "RegionSetDE.fit")) {
      counts <- object@counts
    } else if (methods::is(object, "RegionSetDE.counts")) {
      counts <- object
    } else {
      stop("The 'object' parameter must be a RegionSetDE.fit or a RegionSetDE.counts object.", call. = FALSE)
    }

    if (!(type %in% c("otherSets", "supplied"))) {
      stop("The 'type' parameter must be either 'otherSets' or 'supplied'.", call. = FALSE)
    }

    match <- match[!is.na(match)]
    if (length(match) > 0 & !all(match %in% c("width", "abundance"))) {
      stop("The 'match' parameter must contain only 'width' and 'abundance'.", call. = FALSE)
    }

    #-------------------------------#
    # Row level covariates          #
    #-------------------------------#
    rowTable <- .universeCovariates(counts = counts)

    setNames <- unique(as.character(rowTable$region.set))
    if (!is.null(regionSets)) {
      absentSets <- setdiff(regionSets, setNames)
      if (length(absentSets) > 0) {
        stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }
      setNames <- regionSets
    }

    #-------------------------------#
    # Supplied comparison rows      #
    #-------------------------------#
    if (type == "supplied") {
      if (is.null(index) | is.null(names(index))) {
        stop("The 'supplied' type needs a named list in the 'index' parameter.", call. = FALSE)
      }

      absentSets <- setdiff(setNames, names(index))
      if (length(absentSets) > 0) {
        stop(paste0("The 'index' list is missing the following sets: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }

      universeIndex <- lapply(setNames,
                              function(setName) {
                                setRows <- rowTable$row.index[rowTable$region.set == setName]
                                return(sort(unique(c(setRows, as.integer(index[[setName]])))))
                              })
      names(universeIndex) <- setNames

      return(new(Class = "RegionSetDE.universe",
                 index = universeIndex,
                 type = type,
                 matching = character(0),
                 diagnostics = .universeDiagnostics(rowTable = rowTable, universeIndex = universeIndex, setNames = setNames),
                 n.rows = nrow(counts)))
    }

    #-------------------------------#
    # Other sets, matched or not    #
    #-------------------------------#
    if (length(unique(rowTable$region.set)) < 2) {
      stop(paste0("The object holds a single region set ('", unique(rowTable$region.set)[1],
                  "'), and a competitive test needs something to compare it to."), call. = FALSE)
    }

    set.seed(seed)

    universeIndex <-
      lapply(setNames,
             function(setName) {
               setRows <- dplyr::filter(rowTable, .data$region.set == setName)
               eligibleRows <- dplyr::filter(rowTable, .data$region.set != setName)

               # Without matching every other row is eligible, which is what plain camera does
               comparisonIndex <- if (length(match) == 0) {
                 eligibleRows$row.index
               } else {
                 .matchStrata(setRows = setRows, eligibleRows = eligibleRows,
                              match = match, ratio = ratio, strata = strata)
               }

               return(sort(unique(c(setRows$row.index, comparisonIndex))))
             })

    names(universeIndex) <- setNames

    diagnosticsTable <- .universeDiagnostics(rowTable = rowTable, universeIndex = universeIndex, setNames = setNames)

    if (isTRUE(verbose)) {
      message(paste0("Universe built for ", length(setNames), " sets",
                     if (length(match) > 0) {paste0(", matched on ", paste(match, collapse = " and "))} else {""}, "."))
      print(diagnosticsTable, row.names = FALSE)
    }

    return(new(Class = "RegionSetDE.universe",
               index = universeIndex,
               type = type,
               matching = match,
               diagnostics = diagnosticsTable,
               n.rows = nrow(counts)))
  } # END function




#' @title .resolveUniverse
#'
#' @description Turns the \code{universe} argument of \code{\link{fitRegions}} and \code{\link{testRegionSets}} into a \code{RegionSetDE.universe} object, building it when a keyword was given rather than an object.
#'
#' @param object \code{RegionSetDE.fit} or \code{RegionSetDE.counts} object.
#' @param universe String with a keyword, a \code{RegionSetDE.universe} object, or \code{NULL}.
#' @param matchOn Character vector with the covariates the comparison rows are matched on.
#' @param universeRatio Numeric value with the number of comparison rows drawn per region of the set.
#' @param regionSets Character vector with the sets being tested, or \code{NULL}.
#' @param soft Logical value to indicate whether a universe that cannot be built must return empty rather than stop. Default: \code{FALSE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.universe} object, empty when none could be built under \code{soft}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods is new
#'
#' @keywords internal

.resolveUniverse <-
  function(object,
           universe = "matched",
           matchOn = c("width", "abundance"),
           universeRatio = 5,
           regionSets = NULL,
           soft = FALSE,
           verbose = TRUE) {

    nRows <- if (methods::is(object, "RegionSetDE.fit")) {nrow(object@counts)} else {nrow(object)}

    #-------------------------------#
    # An object is taken as it is   #
    #-------------------------------#
    if (methods::is(universe, "RegionSetDE.universe")) {
      if (universe@n.rows != nRows) {
        stop("The universe was built on an object of a different size, rebuild it on this one.", call. = FALSE)
      }
      return(universe)
    }

    if (is.null(universe) | isFALSE(universe)) {
      return(new(Class = "RegionSetDE.universe"))
    }

    if (!is.character(universe) | length(universe) != 1) {
      stop("The 'universe' parameter must be 'matched', 'all', NULL, or a RegionSetDE.universe object.", call. = FALSE)
    }

    if (!(universe %in% c("matched", "all"))) {
      stop("The 'universe' parameter must be 'matched', 'all', NULL, or a RegionSetDE.universe object.", call. = FALSE)
    }

    #-------------------------------#
    # Build it                      #
    #-------------------------------#
    builtUniverse <- try(makeSetUniverse(object = object,
                                         type = "otherSets",
                                         match = if (universe == "matched") {matchOn} else {character(0)},
                                         ratio = universeRatio,
                                         regionSets = regionSets,
                                         verbose = verbose),
                         silent = TRUE)

    if (inherits(builtUniverse, "try-error")) {
      # A fit is useful without a universe, so failing to build one there is a remark rather than an error
      if (isTRUE(soft)) {
        if (isTRUE(verbose)) {
          message(paste0("No comparison universe was built: ", sub("^Error[^:]*: ", "", builtUniverse[1])))
        }
        return(new(Class = "RegionSetDE.universe"))
      }
      stop(sub("^Error[^:]*: ", "", builtUniverse[1]), call. = FALSE)
    }

    return(builtUniverse)
  } # END function




#' @title .universeCovariates
#'
#' @description Assembles the width and the width-adjusted abundance of every row, the two covariates the matching works on.
#'
#' @param counts \code{RegionSetDE.counts} object.
#'
#' @return A data.frame with one row per row of the object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment assay colData rowData rowRanges
#' @importFrom BiocGenerics width
#' @importFrom dplyr mutate
#' @importFrom stats median
#'
#' @keywords internal

.universeCovariates <-
  function(counts) {

    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, 1))
    rowWidths <- BiocGenerics::width(SummarizedExperiment::rowRanges(counts))

    librarySizes <- SummarizedExperiment::colData(counts)$library.size
    if (is.null(librarySizes) | any(is.na(librarySizes))) {
      librarySizes <- colSums(countMatrix)
    }

    rowTable <- as.data.frame(SummarizedExperiment::rowData(counts))

    return(dplyr::mutate(rowTable,
                         row.index = seq_len(nrow(counts)),
                         row.width = rowWidths,
                         abundance = .regionAbundance(countMatrix = countMatrix,
                                                      librarySizes = librarySizes,
                                                      rowWidths = rowWidths,
                                                      referenceWidth = stats::median(rowWidths))))
  } # END function




#' @title .matchStrata
#'
#' @description Draws comparison rows from the strata occupied by a region set, so that the two share a distribution of width and abundance.
#'
#' @param setRows Data.frame with the rows of the set.
#' @param eligibleRows Data.frame with the rows the comparison can be drawn from.
#' @param match Character vector with the covariates to match on.
#' @param ratio Numeric value with the number of comparison rows drawn per region of the set.
#' @param strata Numeric value with the number of strata per covariate.
#'
#' @return An integer vector with the positions of the comparison rows.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats quantile
#'
#' @keywords internal

.matchStrata <-
  function(setRows,
           eligibleRows,
           match = c("width", "abundance"),
           ratio = 5,
           strata = 5) {

    #-------------------------------#
    # Cut on the strata of the set  #
    #-------------------------------#
    # The breaks come from the set, so its own rows spread evenly over the strata and none of them is left unmatched
    strataAssignment <-
      lapply(match,
             function(covariate) {
               columnName <- if (covariate == "width") {"row.width"} else {"abundance"}
               setValues <- setRows[[columnName]]

               breakPoints <- unique(stats::quantile(setValues, probs = seq(0, 1, length.out = strata + 1), na.rm = TRUE))
               if (length(breakPoints) < 2) {
                 return(list(set = rep("1", nrow(setRows)), eligible = rep("1", nrow(eligibleRows))))
               }

               # The tails are opened so that comparison rows outside the range of the set fall in the closest stratum
               breakPoints[1] <- -Inf
               breakPoints[length(breakPoints)] <- Inf

               return(list(set = as.character(cut(setValues, breaks = breakPoints, labels = FALSE)),
                           eligible = as.character(cut(eligibleRows[[columnName]], breaks = breakPoints, labels = FALSE))))
             })

    setStratum <- do.call(what = paste, args = c(lapply(strataAssignment, function(x) {x$set}), sep = "|"))
    eligibleStratum <- do.call(what = paste, args = c(lapply(strataAssignment, function(x) {x$eligible}), sep = "|"))

    #-------------------------------#
    # Draw within each stratum      #
    #-------------------------------#
    setStratumSize <- table(setStratum)

    drawnIndex <-
      unlist(lapply(names(setStratumSize),
                    function(stratumName) {
                      candidateIndex <- eligibleRows$row.index[eligibleStratum == stratumName]
                      requestedNumber <- ceiling(setStratumSize[[stratumName]] * ratio)

                      if (length(candidateIndex) == 0) {
                        return(integer(0))
                      }

                      # A stratum with too few candidates is taken whole, the shortfall shows up in the diagnostics
                      if (length(candidateIndex) <= requestedNumber) {
                        return(candidateIndex)
                      }

                      return(sample(x = candidateIndex, size = requestedNumber, replace = FALSE))
                    }),
             use.names = FALSE)

    drawnIndex <- unique(drawnIndex)

    if (length(drawnIndex) < 10) {
      stop("Fewer than 10 comparison rows could be matched, loosen 'strata' or drop the matching.", call. = FALSE)
    }

    return(sort(drawnIndex))
  } # END function




#' @title .universeDiagnostics
#'
#' @description Summarises how close the comparison rows of a universe sit to the set they were drawn for.
#'
#' @param rowTable Data.frame with the row level covariates.
#' @param universeIndex List with the universe positions of every set.
#' @param setNames Character vector with the names of the sets.
#'
#' @return A data.frame with one row per set.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats median
#'
#' @keywords internal

.universeDiagnostics <-
  function(rowTable,
           universeIndex,
           setNames) {

    diagnosticsList <-
      lapply(setNames,
             function(setName) {
               setIndex <- rowTable$row.index[rowTable$region.set == setName]

               # The universe holds the set as well, the diagnostics describe what is left once it is taken out
               comparisonIndex <- setdiff(universeIndex[[setName]], setIndex)

               setRows <- rowTable[rowTable$row.index %in% setIndex, , drop = FALSE]
               comparisonRows <- rowTable[rowTable$row.index %in% comparisonIndex, , drop = FALSE]

               return(data.frame(region.set = setName,
                                 n.regions = nrow(setRows),
                                 n.comparison = nrow(comparisonRows),
                                 median.width = stats::median(setRows$row.width),
                                 median.width.comparison = stats::median(comparisonRows$row.width),
                                 median.abundance = round(stats::median(setRows$abundance), 2),
                                 median.abundance.comparison = round(stats::median(comparisonRows$abundance), 2),
                                 stringsAsFactors = FALSE))
             })

    return(do.call(what = rbind, args = diagnosticsList))
  } # END function
