#' @title testRegionSets
#'
#' @description Asks whether a region set responds to a contrast as a whole. Two questions are answered side by side: whether the regions of the set move away from zero, which is a self-contained claim, and whether they move more than the regions they are compared against, which is a competitive one. Both are computed from the per-region statistics of the same fit, so they never disagree with \code{\link{testRegions}} on the design, the offsets or the dispersion.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrast Contrast to test, in the syntax accepted by \code{\link{testRegions}}, or a named list of contrasts to run on the same fit.
#' @param method Character vector with the tests to run, among \code{"camera"} (competitive) and \code{"fry"} (self-contained). Default: \code{c("camera", "fry")}.
#' @param universe What each set is compared against in the competitive test. Default: \code{NULL}, the universe carried by the fit. A \code{RegionSetDE.universe} object, or the strings \code{"matched"} and \code{"all"}, override it and are built here.
#' @param matchOn Character vector with the covariates the comparison rows are matched on, when one has to be built here. Default: \code{c("width", "abundance")}.
#' @param universeRatio Numeric value with the number of comparison rows drawn per region of the set, when one has to be built here. Default: \code{5}.
#' @param interRegionCor Numeric value with the correlation between regions, used to inflate the variance of both the set and the rows it is compared against. Default: \code{NULL}, estimated separately for each of the two from the residuals of the fit, or held at 0.01 when the design leaves no residual to estimate it from.
#' @param useRanks Logical value to indicate whether \code{camera} must work on the ranks rather than on the statistics, which is more robust and less powerful. Default: \code{FALSE}.
#' @param FDR Numeric value with the adjusted p-value cut-off reported in the output. Default: \code{0.05}.
#' @param adjustMethod String with the multiple testing correction across the sets. Default: \code{"BH"}.
#' @param regionSets Character vector with the names of the sets to test. Default: \code{NULL}, all of them.
#' @param carryCounts Logical value to indicate whether the counts must travel inside the result, so that \code{\link{plotSetSignal}} can draw the signal without being handed the counts object again. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.setResults} object, or a \code{RegionSetDE.setResultsList} when \code{contrast} is a named list.
#'
#' @details The effect size, not the p-value, is the primary output here. A set of 30,000 promoters tested as if its regions were independent returns a p-value below anything a computer will print for a mean shift of 0.05 log2, which says nothing about whether the shift matters. The regions of a set are not independent either, since neighbouring elements inside the same domain move together, so the variance of the mean log2 fold change is inflated by the factor \code{1 + (n - 1) * rho}, with \code{rho} estimated from the residuals of the fit through \code{limma::interGeneCorrelation}. The confidence interval in the output carries that inflation; read it before reading the p-value.
#' The two tests answer different questions and the pattern between them is informative. \code{camera} is competitive: it asks whether the regions of the set moved more than the regions they are compared against, and it is invariant to a scaling error affecting every region equally. \code{fry} is self-contained: it asks whether they moved away from zero at all, which a global shift in the mark, or a residual normalisation error, will satisfy for every set at once. When camera separates the sets and fry does not, the sets redistributed the signal between them; when fry is significant everywhere and camera nowhere, everything moved together and the normalisation deserves a second look before the biology does.
#' The comparison universe comes from the fit, which built it once, and travels on into the result, so \code{\link{plotUniverseMatching}} can check the matching afterwards without anything being kept on the side. Passing a \code{RegionSetDE.universe} object, or one of the two keywords, overrides it for this test alone.
#' The interval on \code{delta.log2FC} carries the correlation of both sides. The regions of the comparison are no less correlated than the regions of the set, so treating their mean as if it were known would leave the interval narrower than the data supports, by around a factor of the square root of two when the two sides are of similar size.
#'
#' A fit with no replicates loses the self-contained test. \code{fry} builds a linear model inside each set and needs a residual to measure it against, which a design with one sample per level does not have, so it is dropped with a message and only the competitive test runs. The correlation between regions goes the same way: it is estimated from the residuals of the fit, and without them it falls back to 0.01, the value \code{limma} uses when nothing better is available. That number sets how much the confidence interval is widened, so on such a fit the interval is as assumed as the dispersion is, and \code{interRegionCor} is worth setting by hand from a replicated experiment on the same assay when one exists.
#' The competitive test runs through \code{limma::cameraPR} on the per-region statistics, which is what makes it work identically for the four engines. The self-contained test needs the values themselves and is computed on the log-CPM matrix of the fit; for \code{edgeR} and \code{DESeq2} that matrix is a transformation of the counts rather than the quantity the model was fitted on, so the two are close but not identical, and the competitive test is the one to lead with.
#'
#' @examples
#' \dontrun{
#' fit <- fitRegions(counts, design = ~ replicate + condition, engine = "edgeR")
#'
#' # The universe comes from the fit and travels into the result
#' setRes <- testRegionSets(fit, contrast = "conditionCOMBO")
#'
#' plotUniverseMatching(setRes)
#' plotSetEffect(setRes)
#'
#' # Overriding it for one test
#' setRes <- testRegionSets(fit, contrast = "conditionCOMBO", universe = "all")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testSetContrast}}, \code{\link{makeSetUniverse}}, \code{\link{testRegions}}, \code{\link{plotSetEffect}}
#'
#' @importFrom SummarizedExperiment colData rowData rowRanges
#' @importFrom BiocGenerics width
#' @importFrom limma cameraPR fry
#' @importFrom stats p.adjust median
#' @importFrom dplyr mutate filter arrange
#' @importFrom rlang .data
#' @importFrom methods is new
#'
#' @export testRegionSets

testRegionSets <-
  function(fit,
           contrast,
           method = c("camera", "fry"),
           universe = NULL,
           matchOn = c("width", "abundance"),
           universeRatio = 5,
           interRegionCor = NULL,
           useRanks = FALSE,
           FDR = 0.05,
           adjustMethod = "BH",
           regionSets = NULL,
           carryCounts = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(fit, "RegionSetDE.fit")) {
      stop("The 'fit' parameter must be a RegionSetDE.fit object.", call. = FALSE)
    }

    #-------------------------------#
    # Several contrasts at once     #
    #-------------------------------#
    # One fit, several contrasts: the recursion keeps a single code path for the test itself
    if (is.list(contrast) & !is.data.frame(contrast)) {
      if (is.null(names(contrast)) | any(names(contrast) == "")) {
        names(contrast) <- paste0("contrast", seq_along(contrast))
      }

      # Resolved once here, so an override is not rebuilt for every contrast
      if ("camera" %in% method) {
        universe <- .setUniverseOf(fit = fit, universe = universe, matchOn = matchOn,
                                   universeRatio = universeRatio, regionSets = regionSets, verbose = verbose)
      }

      resultsList <-
        lapply(names(contrast),
               function(contrastName) {
                 if (isTRUE(verbose)) {
                   message("--- ", contrastName, " ---")
                 }
                 return(testRegionSets(fit = fit, contrast = contrast[[contrastName]], method = method,
                                       universe = universe, matchOn = matchOn, universeRatio = universeRatio,
                                       interRegionCor = interRegionCor,
                                       useRanks = useRanks, FDR = FDR, adjustMethod = adjustMethod,
                                       regionSets = regionSets, carryCounts = carryCounts, verbose = verbose))
               })

      names(resultsList) <- names(contrast)

      return(new(Class = "RegionSetDE.setResultsList",
                 results = resultsList,
                 contrasts = names(contrast)))
    }

    method <- unique(method)
    if (!all(method %in% c("camera", "fry"))) {
      stop("The 'method' parameter must contain only 'camera' and 'fry'.", call. = FALSE)
    }

    #-------------------------------#
    # What the design can support   #
    #-------------------------------#
    residualDegrees <- nrow(fit@design) - ncol(fit@design)

    # fry fits a model inside each set, which needs something left over to measure it against
    if ("fry" %in% method & residualDegrees < 1) {
      if (identical(method, "fry")) {
        stop("The design uses ", ncol(fit@design), " coefficients for ", nrow(fit@design),
             " samples, leaving no residual for the self-contained test. Use method = 'camera'.", call. = FALSE)
      }

      method <- setdiff(method, "fry")
      if (isTRUE(verbose)) {
        message("No residual degree of freedom: the self-contained test has been dropped and only camera is run.")
      }
    }

    # Not gated on 'verbose': 0.01 against a measured value of 0.4 is a fortyfold change in every
    # variance, and a design that gains a coefficient can cross this line with nothing else looking different
    if (is.null(interRegionCor) & residualDegrees < 2) {
      interRegionCor <- 0.01
      warning("The correlation between regions cannot be estimated with ", residualDegrees,
              " residual degrees of freedom, and is held at 0.01. Every confidence interval and p-value ",
              "below rests on that number. Set 'interRegionCor' from a replicated experiment on the same ",
              "assay, and state the value in the methods.", call. = FALSE)
    }

    #-------------------------------#
    # Universe of the sets          #
    #-------------------------------#
    if ("camera" %in% method) {
      universe <- .setUniverseOf(fit = fit, universe = universe, matchOn = matchOn,
                                 universeRatio = universeRatio, regionSets = regionSets, verbose = verbose)
    } else {
      universe <- new(Class = "RegionSetDE.universe")
    }

    contrastObject <- .resolveContrast(contrast = contrast,
                                       design = fit@design,
                                       colData = SummarizedExperiment::colData(fit@counts))

    #-------------------------------#
    # Per-region statistics         #
    #-------------------------------#
    regionStats <- .setStatistics(fit = fit, contrastObject = contrastObject)

    setNames <- unique(as.character(regionStats$region.set))
    if (!is.null(regionSets)) {
      absentSets <- setdiff(regionSets, setNames)
      if (length(absentSets) > 0) {
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
      }
      setNames <- regionSets
    }

    if (isTRUE(verbose)) {
      message("Testing ", length(setNames), " region sets for '", contrastObject$label, "'.")
    }

    # Only fry and the correlation read the values themselves, the competitive test reads the statistics
    expressionMatrix <- if ("fry" %in% method | is.null(interRegionCor)) {.expressionMatrix(fit = fit)} else {NULL}

    #-------------------------------#
    # One set at a time             #
    #-------------------------------#
    resultList <-
      lapply(setNames,
             function(setName) {
               setIndex <- which(regionStats$region.set == setName)

               if (length(setIndex) < 2) {
                 stop("The set '", setName, "' holds fewer than 2 regions, a set level test needs more.", call. = FALSE)
               }

               # The universe holds the set as well, the comparison is what is left once it is taken out
               backgroundIndex <- if (length(universe@index) == 0) {
                 setdiff(seq_len(nrow(regionStats)), setIndex)
               } else {
                 setdiff(universe@index[[setName]], setIndex)
               }

               if (length(backgroundIndex) < 10) {
                 stop("The set '", setName, "' has fewer than 10 rows in its universe to be compared against.", call. = FALSE)
               }

               #-------------------------------#
               # Correlation between regions   #
               #-------------------------------#
               setCorrelation <- if (is.null(interRegionCor)) {
                 .interRegionCor(expressionMatrix = expressionMatrix, design = fit@design,
                                 index = setIndex, label = setName)
               } else {
                 interRegionCor
               }

               # The comparison rows are as correlated as the set, and treating their mean as known
               # would make the interval on the difference narrower than the data supports
               universeCorrelation <- if (is.null(interRegionCor)) {
                 .interRegionCor(expressionMatrix = expressionMatrix, design = fit@design,
                                 index = backgroundIndex,
                                 label = paste0("the universe of ", setName))
               } else {
                 interRegionCor
               }

               #-------------------------------#
               # Effect size, with inflation   #
               #-------------------------------#
               effectSize <- .setEffectSize(logFC = regionStats$log2FC,
                                            setIndex = setIndex,
                                            backgroundIndex = backgroundIndex,
                                            correlation = setCorrelation,
                                            backgroundCorrelation = universeCorrelation)

               setRow <- data.frame(region.set = setName,
                                    n.regions = length(setIndex),
                                    n.comparison = length(backgroundIndex),
                                    mean.log2FC = effectSize$mean.set,
                                    median.log2FC = stats::median(regionStats$log2FC[setIndex]),
                                    mean.log2FC.comparison = effectSize$mean.background,
                                    delta.log2FC = effectSize$delta,
                                    CI.lower = effectSize$ci.lower,
                                    CI.upper = effectSize$ci.upper,
                                    inter.region.cor = setCorrelation,
                                    inter.region.cor.universe = universeCorrelation,
                                    median.width = stats::median(regionStats$width[setIndex]),
                                    stringsAsFactors = FALSE)

               #-------------------------------#
               # Competitive test              #
               #-------------------------------#
               if ("camera" %in% method) {
                 # The test is restricted to the universe of the set, which is what makes the matching count
                 universeIndex <- c(setIndex, backgroundIndex)
                 cameraTable <- limma::cameraPR(statistic = regionStats$stat[universeIndex],
                                                index = list(set = seq_along(setIndex)),
                                                use.ranks = useRanks,
                                                inter.gene.cor = setCorrelation,
                                                sort = FALSE)

                 setRow$camera.direction <- as.character(cameraTable$Direction[1])
                 setRow$camera.p <- cameraTable$PValue[1]
               }

               #-------------------------------#
               # Self-contained test           #
               #-------------------------------#
               if ("fry" %in% method) {
                 fryTable <- limma::fry(y = expressionMatrix,
                                        index = list(set = setIndex),
                                        design = fit@design,
                                        contrast = contrastObject$vector,
                                        sort = FALSE)

                 setRow$fry.direction <- as.character(fryTable$Direction[1])
                 setRow$fry.p <- fryTable$PValue[1]
               }

               return(setRow)
             })

    resultTable <- do.call(what = rbind, args = resultList)

    #-------------------------------#
    # Correct across the sets       #
    #-------------------------------#
    if ("camera" %in% method) {
      resultTable$camera.FDR <- stats::p.adjust(resultTable$camera.p, method = adjustMethod)
    }
    if ("fry" %in% method) {
      resultTable$fry.FDR <- stats::p.adjust(resultTable$fry.p, method = adjustMethod)
    }

    resultTable <- dplyr::arrange(resultTable, dplyr::desc(abs(.data$delta.log2FC)))
    rownames(resultTable) <- NULL

    #-------------------------------#
    # Assemble the object           #
    #-------------------------------#
    setResults <- new(Class = "RegionSetDE.setResults",
                      counts = if (isTRUE(carryCounts)) {fit@counts} else {new(Class = "RegionSetDE.counts")},
                      results = resultTable,
                      regionStats = regionStats,
                      contrast = contrastObject$label,
                      contrast.groups = contrastObject[intersect(c("column", "groups"), names(contrastObject))],
                      test = "set",
                      methods = method,
                      universe = universe,
                      engine = fit@engine,
                      thresholds = list(FDR = FDR, adjust.method = adjustMethod),
                      blacklist = fit@blacklist,
                      whitelist = fit@whitelist,
                      genome.assembly = fit@genome.assembly,
                      seqlevels.style = fit@seqlevels.style,
                      filtering.log = fit@filtering.log,
                      parameters = c(fit@parameters,
                                     list(testRegionSets = list(contrast = contrastObject$label,
                                                                method = method,
                                                                universe = universe@type,
                                                                interRegionCor = interRegionCor,
                                                                useRanks = useRanks,
                                                                adjustMethod = adjustMethod,
                                                                carryCounts = carryCounts))))

    if (isTRUE(verbose)) {
      message("Done. The confidence interval on 'delta.log2FC' carries the inflation for the correlation between regions.")
    }

    return(setResults)
  } # END function




#' @title testSetContrast
#'
#' @description Asks whether a contrast affects one region set differently from another. This is the comparison behind questions of the kind "does the treatment reduce the mark more at Polycomb promoters than at active enhancers", and it is the one claim that a global normalisation error cannot manufacture, since a scaling factor that is wrong for one set is wrong for the other in the same way.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrast Contrast to test, in the syntax accepted by \code{\link{testRegions}}, or a named list of contrasts to run on the same fit.
#' @param set1 Character vector with the name, or names, of the first region set. Default: \code{NULL}, every pair of sets is tested.
#' @param set2 Character vector with the name, or names, of the second region set. Default: \code{NULL}.
#' @param interRegionCor Numeric value with the correlation between the regions of a set. Default: \code{NULL}, estimated from the residuals.
#' @param useRanks Logical value to indicate whether the test must work on the ranks rather than on the statistics. Default: \code{FALSE}.
#' @param sharedRegions String with what to do with the regions belonging to both sets, either \code{"drop"} or \code{"stop"}. Default: \code{"drop"}.
#' @param FDR Numeric value with the adjusted p-value cut-off reported in the output. Default: \code{0.05}.
#' @param adjustMethod String with the multiple testing correction across the pairs. Default: \code{"BH"}.
#' @param carryCounts Logical value to indicate whether the counts must travel inside the result. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.setResults} object with one row per pair of sets, or a \code{RegionSetDE.setResultsList} when \code{contrast} is a named list.
#'
#' @details The test restricts the universe to the two sets and runs the competitive test of \code{\link{testRegionSets}} on the first of them, which is exactly a comparison of the first set against the second. The effect size is the difference between the two mean log2 fold changes, with a confidence interval carrying the variance inflation of both sets.
#'
#' A region that belongs to both sets carries the same reads into both sides of the comparison and pulls the difference towards zero. Those regions are removed by default and the number removed is reported; \code{sharedRegions = "stop"} refuses to run instead, which is the safer setting when the overlap is unexpected.
#'
#' @examples
#' \dontrun{
#' setContrast <- testSetContrast(fit, contrast = "conditionCOMBO",
#'                                set1 = "enhancers", set2 = "tss")
#'
#' # Every pair at once
#' allPairs <- testSetContrast(fit, contrast = "conditionCOMBO")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegionSets}}, \code{\link{plotSetEffect}}
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom limma cameraPR
#' @importFrom stats p.adjust median
#' @importFrom dplyr arrange desc
#' @importFrom rlang .data
#' @importFrom methods is new
#'
#' @export testSetContrast

testSetContrast <-
  function(fit,
           contrast,
           set1 = NULL,
           set2 = NULL,
           interRegionCor = NULL,
           useRanks = FALSE,
           sharedRegions = "drop",
           FDR = 0.05,
           adjustMethod = "BH",
           carryCounts = TRUE,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(fit, "RegionSetDE.fit")) {
      stop("The 'fit' parameter must be a RegionSetDE.fit object.", call. = FALSE)
    }

    #-------------------------------#
    # Several contrasts at once     #
    #-------------------------------#
    if (is.list(contrast) & !is.data.frame(contrast)) {
      if (is.null(names(contrast)) | any(names(contrast) == "")) {
        names(contrast) <- paste0("contrast", seq_along(contrast))
      }

      resultsList <-
        lapply(names(contrast),
               function(contrastName) {
                 if (isTRUE(verbose)) {
                   message("--- ", contrastName, " ---")
                 }
                 return(testSetContrast(fit = fit, contrast = contrast[[contrastName]], set1 = set1, set2 = set2,
                                        interRegionCor = interRegionCor, useRanks = useRanks,
                                        sharedRegions = sharedRegions, FDR = FDR, adjustMethod = adjustMethod,
                                        carryCounts = carryCounts, verbose = verbose))
               })

      names(resultsList) <- names(contrast)

      return(new(Class = "RegionSetDE.setResultsList",
                 results = resultsList,
                 contrasts = names(contrast)))
    }

    if (!(sharedRegions %in% c("drop", "stop"))) {
      stop("The 'sharedRegions' parameter must be either 'drop' or 'stop'.", call. = FALSE)
    }

    if (is.null(interRegionCor) & (nrow(fit@design) - ncol(fit@design)) < 2) {
      interRegionCor <- 0.01
      warning("The correlation between regions cannot be estimated from this design, and is held at 0.01. ",
              "Every confidence interval and p-value below rests on that number.", call. = FALSE)
    }

    contrastObject <- .resolveContrast(contrast = contrast,
                                       design = fit@design,
                                       colData = SummarizedExperiment::colData(fit@counts))
    regionStats <- .setStatistics(fit = fit, contrastObject = contrastObject)
    setNames <- unique(as.character(regionStats$region.set))

    #-------------------------------#
    # Pairs to test                 #
    #-------------------------------#
    if (is.null(set1) & is.null(set2)) {
      if (length(setNames) < 2) {
        stop("At least two region sets are needed to contrast them.", call. = FALSE)
      }
      pairTable <- as.data.frame(t(utils::combn(setNames, 2)), stringsAsFactors = FALSE)
      colnames(pairTable) <- c("set.1", "set.2")

    } else {
      if (is.null(set1) | is.null(set2)) {
        stop("Both 'set1' and 'set2' must be given, or neither of them.", call. = FALSE)
      }
      absentSets <- setdiff(c(set1, set2), setNames)
      if (length(absentSets) > 0) {
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
      }
      pairTable <- data.frame(set.1 = paste(set1, collapse = "+"),
                              set.2 = paste(set2, collapse = "+"),
                              stringsAsFactors = FALSE)
    }

    expressionMatrix <- if (is.null(interRegionCor)) {.expressionMatrix(fit = fit)} else {NULL}

    #-------------------------------#
    # One pair at a time            #
    #-------------------------------#
    resultList <-
      lapply(seq_len(nrow(pairTable)),
             function(i) {
               firstNames <- if (is.null(set1)) {pairTable$set.1[i]} else {set1}
               secondNames <- if (is.null(set2)) {pairTable$set.2[i]} else {set2}

               firstIndex <- which(regionStats$region.set %in% firstNames)
               secondIndex <- which(regionStats$region.set %in% secondNames)

               #-------------------------------#
               # Regions sitting in both sets  #
               #-------------------------------#
               # The same reads on both sides of the comparison drag the difference towards zero
               sharedKeys <- intersect(regionStats$region.key.plain[firstIndex], regionStats$region.key.plain[secondIndex])

               if (length(sharedKeys) > 0) {
                 if (sharedRegions == "stop") {
                   stop(length(sharedKeys), " regions belong to both '", pairTable$set.1[i], "' and '",
                        pairTable$set.2[i], "'.", call. = FALSE)
                 }
                 firstIndex <- firstIndex[!(regionStats$region.key.plain[firstIndex] %in% sharedKeys)]
                 secondIndex <- secondIndex[!(regionStats$region.key.plain[secondIndex] %in% sharedKeys)]
               }

               if (length(firstIndex) < 2 | length(secondIndex) < 2) {
                 stop("Fewer than 2 regions are left in one side of the pair '", pairTable$set.1[i],
                      "' versus '", pairTable$set.2[i], "'.", call. = FALSE)
               }

               #-------------------------------#
               # Correlation and effect size   #
               #-------------------------------#
               firstCorrelation <- if (is.null(interRegionCor)) {
                 .interRegionCor(expressionMatrix = expressionMatrix, design = fit@design,
                                 index = firstIndex, label = pairTable$set.1[i])
               } else {interRegionCor}

               secondCorrelation <- if (is.null(interRegionCor)) {
                 .interRegionCor(expressionMatrix = expressionMatrix, design = fit@design,
                                 index = secondIndex, label = pairTable$set.2[i])
               } else {interRegionCor}

               effectSize <- .setEffectSize(logFC = regionStats$log2FC,
                                            setIndex = firstIndex,
                                            backgroundIndex = secondIndex,
                                            correlation = firstCorrelation,
                                            backgroundCorrelation = secondCorrelation)

               #-------------------------------#
               # Competitive test on the pair  #
               #-------------------------------#
               universeIndex <- c(firstIndex, secondIndex)
               cameraTable <- limma::cameraPR(statistic = regionStats$stat[universeIndex],
                                              index = list(set = seq_along(firstIndex)),
                                              use.ranks = useRanks,
                                              inter.gene.cor = firstCorrelation,
                                              sort = FALSE)

               return(data.frame(set.1 = pairTable$set.1[i],
                                 set.2 = pairTable$set.2[i],
                                 n.regions.1 = length(firstIndex),
                                 n.regions.2 = length(secondIndex),
                                 n.shared.dropped = length(sharedKeys),
                                 mean.log2FC.1 = effectSize$mean.set,
                                 mean.log2FC.2 = effectSize$mean.background,
                                 delta.log2FC = effectSize$delta,
                                 CI.lower = effectSize$ci.lower,
                                 CI.upper = effectSize$ci.upper,
                                 inter.region.cor.1 = firstCorrelation,
                                 inter.region.cor.2 = secondCorrelation,
                                 camera.direction = as.character(cameraTable$Direction[1]),
                                 camera.p = cameraTable$PValue[1],
                                 stringsAsFactors = FALSE))
             })

    resultTable <- do.call(what = rbind, args = resultList)
    resultTable$camera.FDR <- stats::p.adjust(resultTable$camera.p, method = adjustMethod)
    resultTable <- dplyr::arrange(resultTable, dplyr::desc(abs(.data$delta.log2FC)))
    rownames(resultTable) <- NULL

    if (isTRUE(verbose)) {
      message("Done. ", nrow(resultTable), " pairs tested for '", contrastObject$label, "'.")
    }

    return(new(Class = "RegionSetDE.setResults",
               counts = if (isTRUE(carryCounts)) {fit@counts} else {new(Class = "RegionSetDE.counts")},
               results = resultTable,
               regionStats = regionStats,
               contrast = contrastObject$label,
               contrast.groups = contrastObject[intersect(c("column", "groups"), names(contrastObject))],
               test = "setContrast",
               methods = "camera",
               universe = new(Class = "RegionSetDE.universe", type = "pairedSet", n.rows = nrow(fit@counts)),
               engine = fit@engine,
               thresholds = list(FDR = FDR, adjust.method = adjustMethod),
               blacklist = fit@blacklist,
               whitelist = fit@whitelist,
               genome.assembly = fit@genome.assembly,
               seqlevels.style = fit@seqlevels.style,
               filtering.log = fit@filtering.log,
               parameters = c(fit@parameters,
                              list(testSetContrast = list(contrast = contrastObject$label,
                                                          set1 = set1,
                                                          set2 = set2,
                                                          sharedRegions = sharedRegions,
                                                          useRanks = useRanks,
                                                          adjustMethod = adjustMethod,
                                                          carryCounts = carryCounts)))))
  } # END function




#' @title .setStatistics
#'
#' @description Runs the per-region contrast and returns the statistics the set level tests are built on, one row per row of the fit.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrastObject List returned by \code{.resolveContrast}.
#'
#' @return A data.frame with the region annotation and the per-region \code{log2FC} and \code{stat} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment rowData rowRanges
#' @importFrom BiocGenerics width
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#'
#' @keywords internal

.setStatistics <-
  function(fit,
           contrastObject) {

    rawTable <- switch(fit@engine,
                       "edgeR" = .testEdgeR(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = 0),
                       "voom" = .testVoom(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = 0),
                       "dream" = .testDream(fit = fit, contrastObject = contrastObject, lfcThreshold = 0, verbose = FALSE),
                       "deseq2" = .testDESeq2(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = 0))

    # The quasi-likelihood F carries no sign, so a signed statistic has to be rebuilt for the competitive test
    if (fit@engine == "edgeR") {
      rawTable$stat <- sign(rawTable$log2FC) * sqrt(pmax(rawTable$stat, 0))
    }

    rowTable <- as.data.frame(SummarizedExperiment::rowData(fit@counts))

    return(dplyr::mutate(rawTable,
                         region.set = as.character(rowTable$region.set),
                         region.id = as.character(rowTable$region.id),
                         region.key.plain = as.character(rowTable$region.id),
                         width = BiocGenerics::width(SummarizedExperiment::rowRanges(fit@counts))))
  } # END function




#' @title .expressionMatrix
#'
#' @description Returns the matrix of log2 values the fit was built on, or the closest transformation of the counts when the engine works on the count scale.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param priorCount Numeric value with the prior count added before taking the logarithm. Default: \code{2}.
#'
#' @return A numeric matrix with one row per region and one column per sample.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom edgeR getOffset
#'
#' @keywords internal

.expressionMatrix <-
  function(fit,
           priorCount = 2) {

    # voom and dream already hold the values the model saw, nothing has to be rebuilt
    if (fit@engine %in% c("voom", "dream")) {
      return(as.matrix(fit@fit$voom$E))
    }

    countMatrix <- as.matrix(SummarizedExperiment::assay(fit@counts, 1))

    offsetMatrix <- .fitOffsets(counts = fit@counts, useOffsets = TRUE, verbose = FALSE)

    if (is.null(offsetMatrix)) {
      librarySizes <- SummarizedExperiment::colData(fit@counts)$library.size
      if (is.null(librarySizes) | any(is.na(librarySizes))) {
        librarySizes <- colSums(countMatrix)
      }
      offsetMatrix <- matrix(data = rep(log(librarySizes), each = nrow(countMatrix)),
                             nrow = nrow(countMatrix), ncol = ncol(countMatrix))
    } else {
      # The offsets are centred on the library sizes so that the values come out on a CPM scale
      offsetMatrix <- offsetMatrix - rowMeans(offsetMatrix) + mean(log(colSums(countMatrix)))
    }

    logMatrix <- log2(countMatrix + priorCount) - offsetMatrix / log(2) + log2(1e6)
    dimnames(logMatrix) <- dimnames(countMatrix)

    return(logMatrix)
  } # END function




#' @title .interRegionCor
#'
#' @description Estimates the correlation between the regions of a set, from the residuals of the design.
#'
#' @param expressionMatrix Numeric matrix of log2 values.
#' @param design Design matrix.
#' @param index Integer vector with the rows of the set.
#' @param label String naming what is being estimated, used in the warning. Default: \code{"the set"}.
#'
#' @return A numeric value.
#'
#' @details When the design leaves fewer than two residual degrees of freedom there is nothing to estimate a correlation from, and the value falls back to the 0.01 that \code{limma} uses in the same situation. That fallback is loud rather than silent, because the difference between 0.01 and a measured 0.4 is a fortyfold change in every variance, and a design that gains a coefficient can cross that line without anything else about the analysis appearing to change.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom limma interGeneCorrelation
#' @importFrom stats var
#'
#' @keywords internal

.interRegionCor <-
  function(expressionMatrix,
           design,
           index,
           label = "the set") {

    residualDegrees <- nrow(design) - ncol(design)

    # Below two residual degrees of freedom limma cannot separate the shared variation from the noise
    if (residualDegrees < 2) {
      warning("The correlation between regions cannot be estimated with ", residualDegrees,
              " residual degrees of freedom, and is held at 0.01.", call. = FALSE)
      return(0.01)
    }

    if (length(index) < 3) {
      warning("Fewer than 3 rows in ", label, ", the correlation between regions is held at 0.01, ",
              "which leaves that set effectively uninflated.", call. = FALSE)
      return(0.01)
    }

    setMatrix <- expressionMatrix[index, , drop = FALSE]

    # A row with no variation gives no residual to correlate and turns the estimate into NaN
    rowVariance <- apply(setMatrix, MARGIN = 1, FUN = stats::var)
    setMatrix <- setMatrix[is.finite(rowVariance) & rowVariance > 0, , drop = FALSE]

    if (nrow(setMatrix) < 3) {
      warning("Fewer than 3 rows of ", label, " vary at all, the correlation between regions is held at 0.01.",
              call. = FALSE)
      return(0.01)
    }

    correlationValue <- try(limma::interGeneCorrelation(y = setMatrix, design = design)$correlation, silent = TRUE)

    if (inherits(correlationValue, "try-error") | !is.finite(correlationValue)) {
      warning("The correlation between the regions of ", label, " could not be computed, and is held at 0.01.",
              call. = FALSE)
      return(0.01)
    }

    # A negative estimate would deflate the variance instead of inflating it, which is not a claim worth making
    return(max(correlationValue, 0))
  } # END function




#' @title .setEffectSize
#'
#' @description Computes the mean log2 fold change of a set, the difference with its background, and a confidence interval inflated for the correlation between the regions.
#'
#' @param logFC Numeric vector with the per-region log2 fold changes.
#' @param setIndex Integer vector with the rows of the set.
#' @param backgroundIndex Integer vector with the rows of the background.
#' @param correlation Numeric value with the correlation between the regions of the set.
#' @param backgroundCorrelation Numeric value with the correlation between the regions of the background. Default: \code{NULL}, the same as the set, since a comparison drawn from the same object is correlated in the same way.
#' @param level Numeric value with the confidence level. Default: \code{0.95}.
#'
#' @return A list with the means, the difference and the bounds of the interval.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats sd qnorm
#'
#' @keywords internal

.setEffectSize <-
  function(logFC,
           setIndex,
           backgroundIndex,
           correlation,
           backgroundCorrelation = NULL,
           level = 0.95) {

    # Zero here would say the comparison mean is known exactly, which it never is
    if (is.null(backgroundCorrelation)) {
      backgroundCorrelation <- correlation
    }

    setValues <- logFC[setIndex]
    backgroundValues <- logFC[backgroundIndex]

    setNumber <- length(setValues)
    backgroundNumber <- length(backgroundValues)

    # Neighbouring regions move together, so the mean of n of them carries the information of far fewer than n
    setVIF <- 1 + (setNumber - 1) * correlation
    backgroundVIF <- 1 + (backgroundNumber - 1) * backgroundCorrelation

    setVariance <- (stats::sd(setValues)^2 / setNumber) * setVIF
    backgroundVariance <- (stats::sd(backgroundValues)^2 / backgroundNumber) * backgroundVIF

    deltaValue <- mean(setValues) - mean(backgroundValues)
    deltaError <- sqrt(setVariance + backgroundVariance)
    criticalValue <- stats::qnorm(1 - (1 - level) / 2)

    return(list(mean.set = mean(setValues),
                mean.background = mean(backgroundValues),
                delta = deltaValue,
                ci.lower = deltaValue - criticalValue * deltaError,
                ci.upper = deltaValue + criticalValue * deltaError))
  } # END function




#' @title .setUniverseOf
#'
#' @description Picks the comparison universe a set level test must use: the one carried by the fit, or one built here when the caller asked for something else.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param universe String with a keyword, a \code{RegionSetDE.universe} object, or \code{NULL} to take the one in the fit.
#' @param matchOn Character vector with the covariates the comparison rows are matched on.
#' @param universeRatio Numeric value with the number of comparison rows drawn per region of the set.
#' @param regionSets Character vector with the sets being tested, or \code{NULL}.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A \code{RegionSetDE.universe} object.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom methods is
#'
#' @keywords internal

.setUniverseOf <-
  function(fit,
           universe = NULL,
           matchOn = c("width", "abundance"),
           universeRatio = 5,
           regionSets = NULL,
           verbose = TRUE) {

    if (methods::is(universe, "RegionSetDE.universe")) {
      return(.resolveUniverse(object = fit, universe = universe, verbose = verbose))
    }

    #-------------------------------#
    # Nothing asked for             #
    #-------------------------------#
    if (is.null(universe)) {
      # fitRegions already built one on these very rows, rebuilding it would give the same answer twice
      if (length(fit@universe@index) > 0) {
        return(fit@universe)
      }
      universe <- "matched"
    }

    return(.resolveUniverse(object = fit,
                            universe = universe,
                            matchOn = matchOn,
                            universeRatio = universeRatio,
                            regionSets = regionSets,
                            soft = FALSE,
                            verbose = verbose))
  } # END function
