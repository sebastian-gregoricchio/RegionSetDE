# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

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
#' @param universeSets Character vector with the names of the sets the comparison rows are drawn from, when one has to be built here. Default: \code{NULL}, every set other than the one being tested.
#' @param effectMethod String with what the confidence interval on the effect is computed from, either \code{"sample"}, which treats the biological samples as the replication, or \code{"region"}, which treats the regions as it. Default: \code{"sample"}.
#' @param interRegionCor Numeric value with the correlation between regions, used to inflate the variance of the region-heterogeneity interval of both the set and the rows it is compared against. Default: \code{NULL}, estimated separately for each of the two from the residuals of the fit, or held at 0.01 when the design leaves no residual to estimate it from.
#' @param tileHandling String with what to do when the fit was built on tiles, either \code{"collapse"}, which averages the tiles of a region back into one row before the set is assembled, or \code{"keep"}, which lets every tile count on its own. Default: \code{"collapse"}.
#' @param overlapPolicy String with what to do about comparison rows overlapping the set in the genome, one among \code{"allow"}, \code{"drop"} and \code{"stop"}. Default: \code{"drop"}.
#' @param useRanks Logical value to indicate whether \code{camera} must work on the ranks rather than on the statistics, which is more robust and less powerful. Default: \code{FALSE}.
#' @param FDR Numeric value with the adjusted p-value cut-off reported in the output. Default: \code{0.05}.
#' @param adjustMethod String with the multiple testing correction across the sets. Default: \code{"BH"}.
#' @param regionSets Character vector with the names of the sets to test. Default: \code{NULL}, all of them.
#' @param carryCounts Logical value to indicate whether the counts must travel inside the result, so that \code{\link{plotSetSignal}} can draw the signal without being handed the counts object again. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.setResults} object, or a \code{RegionSetDE.setResultsList} when \code{contrast} is a named list. The table carries \code{CI.lower} and \code{CI.upper} for the interval selected by \code{effectMethod}, \code{CI.type} naming which one that is, and \code{heterogeneity.CI.lower} and \code{heterogeneity.CI.upper} for the region-level one, always.
#'
#' @details The effect size, not the p-value, is the primary output here. A set of 30,000 promoters tested as if its regions were independent returns a p-value below anything a computer will print for a mean shift of 0.05 log2, which says nothing about whether the shift matters.
#'
#' Two different intervals can be put around that effect and they answer different questions, so both are reported and \code{effectMethod} decides which one is called the confidence interval. The \code{"sample"} interval is the default and is the one to quote as a biological result. One number is computed per library, the mean signal over the set minus the mean signal over its comparison, and those numbers are then run through the design of the experiment. The replication is the biological samples, which is where it comes from in the experiment, and adding regions to a set makes that interval more stable without ever making it narrower than four libraries can support.
#'
#' The \code{"region"} interval is the mean of the per-region log2 fold changes with its variance inflated by \code{1 + (n - 1) * rho}, with \code{rho} estimated from the residuals of the fit through \code{limma::interGeneCorrelation}. It describes how much the effect varies from locus to locus within the set, conditional on these libraries, and it is reported under \code{heterogeneity.CI.lower} and \code{heterogeneity.CI.upper} whatever \code{effectMethod} is set to. It is a useful quantity and it is not a confidence interval on a condition effect: the sampling units behind it are genomic loci, and no number of loci substitutes for the biological replication that was or was not done. Reading it as the second thing rather than the first is the safe habit. Note also that its width does not fall away as the set grows, since \code{sd^2 / n * (1 + (n - 1) * rho)} tends to \code{sd^2 * rho}; it flattens rather than collapsing.
#'
#' The two tests answer different questions and neither of them, alone or in combination, establishes that a set did not change. \code{camera} is competitive: it asks whether the regions of the set moved more than the regions they are compared against, and it is invariant to a scaling error affecting every region equally. \code{fry} is self-contained: it asks whether they moved away from zero at all. Read the four outcomes as evidence and not as mechanism:
#' \itemize{
#'   \item camera significant, fry significant: evidence both that the set moved and that it moved more than its comparison.
#'   \item camera significant, fry not: evidence of a difference relative to the comparison, with the absolute claim left open. Failing to reject the self-contained null is not evidence that the absolute change is zero, and the two tests do not have the same power.
#'   \item camera not significant, fry significant: evidence that the set moved, none that it moved differently from its comparison.
#'   \item neither significant: neither test found evidence, which is not the same as evidence of no effect.
#' }
#' The word for a set gaining what another set lost is redistribution, and it is a claim about two sets rather than about one set and its universe, so it belongs to \code{\link{testSetContrast}}. To argue that a mark did not change globally, an equivalence test against a bounded near-zero interval is what the claim needs; a non-significant \code{fry} is not that. Bear in mind too that any centring normalisation, TMM and background included, removes a genuinely global shift from the data before \code{fry} ever sees it, so the self-contained test is not the place to look for one.
#'
#' The comparison universe comes from the fit, which built it once, and travels on into the result, so \code{\link{plotUniverseMatching}} can check the matching afterwards without anything being kept on the side. Passing a \code{RegionSetDE.universe} object, or one of the two keywords, overrides it for this test alone. Whatever it is, it is made of the other sets loaded into the object, and the competitive p-value is a statement about the set relative to those and not relative to the genome. Load two sets and the test compares them to each other; load four that behave alike and every one of them can come out unremarkable against the other three. \code{\link{makeSetUniverse}} takes \code{universeSets} for choosing that comparison pool explicitly, which is worth doing when the sets were not all picked for the same reason.
#'
#' Both intervals carry the uncertainty of both sides. The regions of the comparison are no less correlated than the regions of the set, so treating their mean as if it were known would leave the interval narrower than the data supports.
#'
#' A set that overlaps its own comparison in the genome shares reads with it and drags the difference towards zero. Overlap is measured on the coordinates rather than on the identifiers, so two sets holding chr1:1000-2000 and chr1:1500-2500 are seen as overlapping even though no region identifier is shared, and \code{overlapPolicy} decides what happens next. The number of comparison rows removed, or left in place, is reported in \code{n.comparison.overlapping}.
#'
#' On a tiled fit the row is a tile, and a set assembled from tiles weights each region by how many tiles it was cut into: a 40 kb domain would count forty times a 2 kb one. That changes the question from the average response of the regions in the set to the average response of the base pairs in it. \code{tileHandling = "collapse"} averages the tiles of a region back together first, which keeps the region as the unit and matches what \code{\link{testRegions}} does at its own level. \code{"keep"} is the base-pair version, and is a deliberate choice rather than a default.
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
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges overlapsAny
#' @importFrom limma cameraPR fry
#' @importFrom stats p.adjust median
#' @importFrom dplyr mutate filter arrange desc
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
           universeSets = NULL,
           effectMethod = "sample",
           interRegionCor = NULL,
           tileHandling = "collapse",
           overlapPolicy = "drop",
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

    effectMethod <- tolower(as.character(effectMethod[1]))
    if (!(effectMethod %in% c("sample", "region"))) {
      stop("The 'effectMethod' parameter must be either 'sample' or 'region'.", call. = FALSE)
    }

    tileHandling <- tolower(as.character(tileHandling[1]))
    if (!(tileHandling %in% c("collapse", "keep"))) {
      stop("The 'tileHandling' parameter must be either 'collapse' or 'keep'.", call. = FALSE)
    }

    overlapPolicy <- tolower(as.character(overlapPolicy[1]))
    if (!(overlapPolicy %in% c("allow", "drop", "stop"))) {
      stop("The 'overlapPolicy' parameter must be one among 'allow', 'drop' or 'stop'.", call. = FALSE)
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
                                   universeRatio = universeRatio, regionSets = regionSets,
                                   universeSets = universeSets, verbose = verbose)
      }

      resultsList <-
        lapply(names(contrast),
               function(contrastName) {
                 if (isTRUE(verbose)) {
                   message("--- ", contrastName, " ---")
                 }
                 return(testRegionSets(fit = fit, contrast = contrast[[contrastName]], method = method,
                                       universe = universe, matchOn = matchOn, universeRatio = universeRatio,
                                       universeSets = universeSets, effectMethod = effectMethod, interRegionCor = interRegionCor,
                                       tileHandling = tileHandling, overlapPolicy = overlapPolicy,
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

    # The sample level interval is a model on one number per library, and one library per coefficient
    # leaves nothing to measure it against
    if (effectMethod == "sample" & residualDegrees < 1) {
      effectMethod <- "region"
      warning("The design leaves no residual degree of freedom, so no interval can be built on the samples. ",
              "The reported interval is the region-heterogeneity one, which describes variation between loci ",
              "and not between biological samples.", call. = FALSE)
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
                                 universeRatio = universeRatio, regionSets = regionSets,
                                 universeSets = universeSets, verbose = verbose)
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

    #-------------------------------#
    # Tiles back into regions       #
    #-------------------------------#
    # A row is the unit the set is averaged over, so on a tiled fit the wide regions would otherwise
    # be counted once per tile and carry the set on their own
    tileMap <- NULL

    if (identical(fit@counting.level, "tile") & tileHandling == "collapse") {
      collapsedList <- .collapseTileStats(regionStats = regionStats)
      tileMap <- collapsedList$map
      regionStats <- collapsedList$stats

      # The universe was indexed on the tiles and has to follow them
      if (length(universe@index) > 0) {
        universe@index <- lapply(universe@index, function(rowIndex) {sort(unique(tileMap[rowIndex]))})
        universe@n.rows <- nrow(regionStats)
      }

      if (isTRUE(verbose)) {
        message(length(tileMap), " tiles collapsed into ", nrow(regionStats),
                " regions, so that the region rather than the base pair is the unit of the set.")
      }
    }

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

    # fry, the correlation and the sample level effect read the values themselves; the competitive test
    # reads the statistics only
    expressionMatrix <- if ("fry" %in% method | is.null(interRegionCor) | effectMethod == "sample") {
      .expressionMatrix(fit = fit)
    } else {
      NULL
    }

    # The values have to sit at the level the statistics do, or the two describe different rows
    if (!is.null(tileMap) & !is.null(expressionMatrix)) {
      expressionMatrix <- .collapseTileMatrix(expressionMatrix = expressionMatrix, tileMap = tileMap)
    }

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

               #-------------------------------#
               # Overlap with the comparison   #
               #-------------------------------#
               # Identifiers say nothing about shared reads: two sets can hold different regions covering
               # the same chromatin, and the difference between them is then partly a comparison with itself
               overlappingRows <- .overlappingRows(regionStats = regionStats,
                                                   setIndex = setIndex,
                                                   comparisonIndex = backgroundIndex)

               if (length(overlappingRows) > 0) {
                 if (overlapPolicy == "stop") {
                   stop(length(overlappingRows), " comparison rows of the set '", setName,
                        "' overlap the set itself in the genome. Set 'overlapPolicy' to 'drop' to remove them.", call. = FALSE)
                 }
                 if (overlapPolicy == "drop") {
                   backgroundIndex <- setdiff(backgroundIndex, overlappingRows)
                 }
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
               # Effect size, both intervals   #
               #-------------------------------#
               # Between loci, conditional on these libraries
               heterogeneity <- .setEffectSize(logFC = regionStats$log2FC,
                                               setIndex = setIndex,
                                               backgroundIndex = backgroundIndex,
                                               correlation = setCorrelation,
                                               backgroundCorrelation = universeCorrelation)

               # Between biological samples, which is where the replication of the experiment actually is
               sampleEffect <- if (effectMethod == "sample") {
                 .sampleSetEffect(expressionMatrix = expressionMatrix,
                                  setIndex = setIndex,
                                  backgroundIndex = backgroundIndex,
                                  design = fit@design,
                                  contrastVector = contrastObject$vector)
               } else {
                 NULL
               }

               setRow <- data.frame(region.set = setName,
                                    n.regions = length(setIndex),
                                    n.comparison = length(backgroundIndex),
                                    n.comparison.overlapping = length(overlappingRows),
                                    mean.log2FC = heterogeneity$mean.set,
                                    median.log2FC = stats::median(regionStats$log2FC[setIndex]),
                                    mean.log2FC.comparison = heterogeneity$mean.background,
                                    delta.log2FC = heterogeneity$delta,
                                    CI.lower = if (is.null(sampleEffect)) {heterogeneity$ci.lower} else {sampleEffect$ci.lower},
                                    CI.upper = if (is.null(sampleEffect)) {heterogeneity$ci.upper} else {sampleEffect$ci.upper},
                                    CI.type = effectMethod,
                                    heterogeneity.CI.lower = heterogeneity$ci.lower,
                                    heterogeneity.CI.upper = heterogeneity$ci.upper,
                                    inter.region.cor = setCorrelation,
                                    inter.region.cor.universe = universeCorrelation,
                                    median.width = stats::median(regionStats$width[setIndex]),
                                    stringsAsFactors = FALSE)

               if (!is.null(sampleEffect)) {
                 setRow$sample.delta.log2FC <- sampleEffect$delta
                 setRow$sample.delta.SE <- sampleEffect$standard.error
                 setRow$sample.delta.df <- sampleEffect$df
                 setRow$sample.delta.p <- sampleEffect$p.value
               }

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
    if ("sample.delta.p" %in% colnames(resultTable)) {
      resultTable$sample.delta.FDR <- stats::p.adjust(resultTable$sample.delta.p, method = adjustMethod)
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
                                                                effectMethod = effectMethod,
                                                                interRegionCor = interRegionCor,
                                                                tileHandling = tileHandling,
                                                                overlapPolicy = overlapPolicy,
                                                                useRanks = useRanks,
                                                                adjustMethod = adjustMethod,
                                                                carryCounts = carryCounts))))

    if (isTRUE(verbose)) {
      if (effectMethod == "sample") {
        message("Done. 'CI.lower' and 'CI.upper' come from the samples through the design; the region-to-region ",
                "spread is in the 'heterogeneity' columns beside them.")
      } else {
        message("Done. 'CI.lower' and 'CI.upper' describe how the effect varies between loci, not between ",
                "biological samples. Use effectMethod = 'sample' for the biological interval.")
      }
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
#' @param effectMethod String with what the confidence interval on the difference is computed from, either \code{"sample"} or \code{"region"}. Default: \code{"sample"}.
#' @param sharedRegions String with what to do with the regions the two sets share in the genome, either \code{"drop"} or \code{"stop"}. Default: \code{"drop"}.
#' @param FDR Numeric value with the adjusted p-value cut-off reported in the output. Default: \code{0.05}.
#' @param adjustMethod String with the multiple testing correction across the pairs. Default: \code{"BH"}.
#' @param carryCounts Logical value to indicate whether the counts must travel inside the result. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.setResults} object with one row per pair of sets, or a \code{RegionSetDE.setResultsList} when \code{contrast} is a named list.
#'
#' @details The test restricts the universe to the two sets and runs the competitive test of \code{\link{testRegionSets}} on the first of them, which is exactly a comparison of the first set against the second. The effect size is the difference between the two mean log2 fold changes, with the interval selected by \code{effectMethod} beside it and the region-heterogeneity one always reported next to it.
#'
#' This is the function for the redistribution question. A set gaining what another set lost is a claim about two sets, and it is the one claim a global normalisation error cannot manufacture, since a scaling factor that is wrong for one set is wrong for the other in the same way. Asking it through the pattern of a competitive and a self-contained test on a single set does not work, because failing to reject a self-contained null is not evidence that the absolute change was zero.
#'
#' A region shared by the two sets carries the same reads into both sides of the comparison and pulls the difference towards zero. Sharing is measured on the genome and not on the identifiers: chr1:1000-2000 in one set and chr1:1500-2500 in the other are half the same chromatin even though neither region identifier appears twice. Overlapping regions are removed from both sides by default and the number removed is reported in \code{n.shared.dropped}; \code{sharedRegions = "stop"} refuses to run instead, which is the safer setting when the overlap is unexpected.
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
#' @importFrom IRanges overlapsAny
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
           effectMethod = "sample",
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
                                        effectMethod = effectMethod, interRegionCor = interRegionCor, useRanks = useRanks,
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

    effectMethod <- tolower(as.character(effectMethod[1]))
    if (!(effectMethod %in% c("sample", "region"))) {
      stop("The 'effectMethod' parameter must be either 'sample' or 'region'.", call. = FALSE)
    }

    if (effectMethod == "sample" & (nrow(fit@design) - ncol(fit@design)) < 1) {
      effectMethod <- "region"
      warning("The design leaves no residual degree of freedom, so no interval can be built on the samples. ",
              "The reported interval is the region-heterogeneity one.", call. = FALSE)
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

    expressionMatrix <- if (is.null(interRegionCor) | effectMethod == "sample") {.expressionMatrix(fit = fit)} else {NULL}

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
               # The same reads on both sides of the comparison drag the difference towards zero, and reads
               # are shared through the coordinates whatever the identifiers happen to say
               sharedFirst <- .overlappingRows(regionStats = regionStats,
                                               setIndex = secondIndex,
                                               comparisonIndex = firstIndex)
               sharedSecond <- .overlappingRows(regionStats = regionStats,
                                                setIndex = firstIndex,
                                                comparisonIndex = secondIndex)

               sharedNumber <- length(sharedFirst) + length(sharedSecond)

               if (sharedNumber > 0) {
                 if (sharedRegions == "stop") {
                   stop(sharedNumber, " regions of '", pairTable$set.1[i], "' and '", pairTable$set.2[i],
                        "' overlap each other in the genome.", call. = FALSE)
                 }
                 firstIndex <- setdiff(firstIndex, sharedFirst)
                 secondIndex <- setdiff(secondIndex, sharedSecond)
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

               heterogeneity <- .setEffectSize(logFC = regionStats$log2FC,
                                               setIndex = firstIndex,
                                               backgroundIndex = secondIndex,
                                               correlation = firstCorrelation,
                                               backgroundCorrelation = secondCorrelation)

               sampleEffect <- if (effectMethod == "sample") {
                 .sampleSetEffect(expressionMatrix = expressionMatrix,
                                  setIndex = firstIndex,
                                  backgroundIndex = secondIndex,
                                  design = fit@design,
                                  contrastVector = contrastObject$vector)
               } else {
                 NULL
               }

               #-------------------------------#
               # Competitive test on the pair  #
               #-------------------------------#
               universeIndex <- c(firstIndex, secondIndex)
               cameraTable <- limma::cameraPR(statistic = regionStats$stat[universeIndex],
                                              index = list(set = seq_along(firstIndex)),
                                              use.ranks = useRanks,
                                              inter.gene.cor = firstCorrelation,
                                              sort = FALSE)

               pairRow <- data.frame(set.1 = pairTable$set.1[i],
                                     set.2 = pairTable$set.2[i],
                                     n.regions.1 = length(firstIndex),
                                     n.regions.2 = length(secondIndex),
                                     n.shared.dropped = sharedNumber,
                                     mean.log2FC.1 = heterogeneity$mean.set,
                                     mean.log2FC.2 = heterogeneity$mean.background,
                                     delta.log2FC = heterogeneity$delta,
                                     CI.lower = if (is.null(sampleEffect)) {heterogeneity$ci.lower} else {sampleEffect$ci.lower},
                                     CI.upper = if (is.null(sampleEffect)) {heterogeneity$ci.upper} else {sampleEffect$ci.upper},
                                     CI.type = effectMethod,
                                     heterogeneity.CI.lower = heterogeneity$ci.lower,
                                     heterogeneity.CI.upper = heterogeneity$ci.upper,
                                     inter.region.cor.1 = firstCorrelation,
                                     inter.region.cor.2 = secondCorrelation,
                                     camera.direction = as.character(cameraTable$Direction[1]),
                                     camera.p = cameraTable$PValue[1],
                                     stringsAsFactors = FALSE)

               if (!is.null(sampleEffect)) {
                 pairRow$sample.delta.log2FC <- sampleEffect$delta
                 pairRow$sample.delta.SE <- sampleEffect$standard.error
                 pairRow$sample.delta.df <- sampleEffect$df
                 pairRow$sample.delta.p <- sampleEffect$p.value
               }

               return(pairRow)
             })

    resultTable <- do.call(what = rbind, args = resultList)
    resultTable$camera.FDR <- stats::p.adjust(resultTable$camera.p, method = adjustMethod)
    if ("sample.delta.p" %in% colnames(resultTable)) {
      resultTable$sample.delta.FDR <- stats::p.adjust(resultTable$sample.delta.p, method = adjustMethod)
    }
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
                                                          effectMethod = effectMethod,
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
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics width start end
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
                       "limma" = .testVoom(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = 0, trend = TRUE),
                       "dream" = .testDream(fit = fit, contrastObject = contrastObject, lfcThreshold = 0, verbose = FALSE),
                       "deseq2" = .testDESeq2(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = 0))

    # The quasi-likelihood F carries no sign, so a signed statistic has to be rebuilt for the competitive test
    if (fit@engine == "edgeR") {
      rawTable$stat <- sign(rawTable$log2FC) * sqrt(pmax(rawTable$stat, 0))
    }

    rowTable <- as.data.frame(SummarizedExperiment::rowData(fit@counts))

    rowRangesObject <- SummarizedExperiment::rowRanges(fit@counts)

    # The coordinates travel along so that the overlap between a set and its comparison can be measured
    # on the genome rather than on the identifiers
    return(dplyr::mutate(rawTable,
                         region.set = as.character(rowTable$region.set),
                         region.id = as.character(rowTable$region.id),
                         region.key.plain = as.character(rowTable$region.id),
                         seqnames = as.character(GenomeInfoDb::seqnames(rowRangesObject)),
                         start = BiocGenerics::start(rowRangesObject),
                         end = BiocGenerics::end(rowRangesObject),
                         width = BiocGenerics::width(rowRangesObject)))
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

    # voom, dream and limma already hold the values the model saw, nothing has to be rebuilt
    if (fit@engine %in% c("voom", "dream", "limma")) {
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
#' @description Computes the mean log2 fold change of a set, the difference with its background, and a region-heterogeneity interval inflated for the correlation between the regions.
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
#' @details The sampling units here are the genomic loci, so the interval describes how much the effect varies from one region of the set to another. It is not a confidence interval on a condition effect, whatever the number of regions: the biological replication of the experiment lives in the samples, and \code{.sampleSetEffect} is where that interval comes from. The two are reported side by side and \code{effectMethod} decides which one \code{\link{testRegionSets}} labels as the confidence interval.
#'
#' Growing the set does not shrink this interval away. With \code{rho} held constant, \code{sd^2 / n * (1 + (n - 1) * rho)} tends to \code{sd^2 * rho} as the set grows, so the width flattens out at a floor set by the correlation rather than falling towards zero. What the inflation does not fix is that \code{sd} is the spread of estimated fold changes, which carries the estimation error of each region along with the genuine variation between them.
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
#' @param universeSets Character vector with the sets the comparison rows are drawn from, or \code{NULL}.
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
           universeSets = NULL,
           verbose = TRUE) {

    if (methods::is(universe, "RegionSetDE.universe")) {
      return(.resolveUniverse(object = fit, universe = universe, verbose = verbose))
    }

    #-------------------------------#
    # Nothing asked for             #
    #-------------------------------#
    if (is.null(universe)) {
      # fitRegions already built one on these very rows, rebuilding it would give the same answer twice,
      # unless the comparison pool is being narrowed here and the stored one was built on all of them
      if (length(fit@universe@index) > 0 & is.null(universeSets)) {
        return(fit@universe)
      }
      universe <- "matched"
    }

    return(.resolveUniverse(object = fit,
                            universe = universe,
                            matchOn = matchOn,
                            universeRatio = universeRatio,
                            regionSets = regionSets,
                            universeSets = universeSets,
                            soft = FALSE,
                            verbose = verbose))
  } # END function




#' @title .statsRanges
#'
#' @description Rebuilds a \code{GRanges} from the coordinate columns of a statistics table.
#'
#' @param regionStats Data.frame returned by \code{.setStatistics}.
#' @param index Integer vector with the rows to take. Default: \code{NULL}, all of them.
#'
#' @return A \code{GRanges}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#'
#' @keywords internal

.statsRanges <-
  function(regionStats,
           index = NULL) {

    if (is.null(index)) {
      index <- seq_len(nrow(regionStats))
    }

    return(GenomicRanges::GRanges(seqnames = regionStats$seqnames[index],
                                  ranges = IRanges::IRanges(start = regionStats$start[index],
                                                            end = regionStats$end[index])))
  } # END function




#' @title .overlappingRows
#'
#' @description Finds the comparison rows that overlap a region set in the genome.
#'
#' @param regionStats Data.frame returned by \code{.setStatistics}.
#' @param setIndex Integer vector with the rows of the set.
#' @param comparisonIndex Integer vector with the rows the set is compared against.
#'
#' @return An integer vector with the positions of the overlapping comparison rows.
#'
#' @details Two region sets can describe the same chromatin without sharing a single identifier, so the check is on the coordinates. A comparison row that overlaps the set carries some of the same reads as the set does, which pulls the difference between the two towards zero however the identifiers were assigned.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom IRanges overlapsAny
#'
#' @keywords internal

.overlappingRows <-
  function(regionStats,
           setIndex,
           comparisonIndex) {

    if (length(comparisonIndex) == 0 | is.null(regionStats$seqnames)) {
      return(integer(0))
    }

    overlapFlag <- IRanges::overlapsAny(query = .statsRanges(regionStats = regionStats, index = comparisonIndex),
                                        subject = .statsRanges(regionStats = regionStats, index = setIndex),
                                        ignore.strand = TRUE)

    return(comparisonIndex[overlapFlag])
  } # END function




#' @title .collapseTileStats
#'
#' @description Averages the tiles of a region back into a single row, so that the region rather than the tile is the unit a set is built from.
#'
#' @param regionStats Data.frame returned by \code{.setStatistics}, one row per tile.
#'
#' @return A list with the collapsed \code{stats} table and \code{map}, an integer vector giving the collapsed row of every tile.
#'
#' @details A set assembled from tiles gives every region a weight equal to the number of tiles it was cut into, which turns the mean over the set into a mean over base pairs. Averaging first restores the region as the unit. The averaged statistic is deliberately the plain mean of the per-tile statistics rather than a combined one: it keeps the scale of the values comparable between a region cut into two tiles and one cut into forty, which is what the competitive test ranks them on.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom dplyr group_by summarise n first
#' @importFrom stats median
#' @importFrom rlang .data
#'
#' @keywords internal

.collapseTileStats <-
  function(regionStats) {

    collapsedTable <-
      dplyr::summarise(dplyr::group_by(regionStats, .data$region.set, .data$region.id),
                       log2FC = mean(.data$log2FC, na.rm = TRUE),
                       average.signal = mean(.data$average.signal, na.rm = TRUE),
                       stat = mean(.data$stat, na.rm = TRUE),
                       p.value = stats::median(.data$p.value, na.rm = TRUE),
                       region.key.plain = dplyr::first(.data$region.key.plain),
                       seqnames = dplyr::first(.data$seqnames),
                       start = min(.data$start, na.rm = TRUE),
                       end = max(.data$end, na.rm = TRUE),
                       n.tiles = dplyr::n(),
                       .groups = "drop")

    collapsedTable <- as.data.frame(collapsedTable, stringsAsFactors = FALSE)

    # The region spans its tiles, so the width follows the collapsed coordinates rather than being summed
    collapsedTable$width <- collapsedTable$end - collapsedTable$start + 1

    # dplyr orders the output on the grouping columns, so the map is built by matching rather than assumed
    tileKey <- paste(regionStats$region.set, regionStats$region.id, sep = "\r")
    collapsedKey <- paste(collapsedTable$region.set, collapsedTable$region.id, sep = "\r")

    return(list(stats = collapsedTable,
                map = match(tileKey, collapsedKey)))
  } # END function




#' @title .collapseTileMatrix
#'
#' @description Averages the rows of a value matrix over the tiles of each region.
#'
#' @param expressionMatrix Numeric matrix with one row per tile.
#' @param tileMap Integer vector giving the collapsed row of every tile, as returned by \code{.collapseTileStats}.
#'
#' @return A numeric matrix with one row per region.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.collapseTileMatrix <-
  function(expressionMatrix,
           tileMap) {

    summedMatrix <- rowsum(x = expressionMatrix, group = tileMap, reorder = TRUE)
    tileNumber <- as.numeric(table(tileMap))

    return(summedMatrix / tileNumber)
  } # END function




#' @title .sampleSetEffect
#'
#' @description Computes one set score per biological sample and runs it through the design, so that the interval on the effect rests on the samples rather than on the genomic loci.
#'
#' @param expressionMatrix Numeric matrix of log2 values, one row per region and one column per sample.
#' @param setIndex Integer vector with the rows of the set.
#' @param backgroundIndex Integer vector with the rows the set is compared against.
#' @param design Design matrix of the fit.
#' @param contrastVector Numeric vector with the contrast, in the columns of the design.
#' @param level Numeric value with the confidence level. Default: \code{0.95}.
#'
#' @return A list with the per-sample scores, the estimate of the contrast on them, its standard error, degrees of freedom, p-value and the bounds of the interval.
#'
#' @details The score of a sample is the mean signal over the set minus the mean signal over its comparison, inside that library. Taking the difference within the sample removes anything that scales the whole library, the sequencing depth and the normalisation factor included, which is what makes the score comparable across samples in the first place. Those scores are then a single response fitted on the design of the experiment, and the interval that comes out has as many degrees of freedom as the design leaves, whether the set holds twenty regions or thirty thousand.
#'
#' This is a different quantity from the interval built on the per-region fold changes, and it is the one that answers the question a reader takes a set-level confidence interval to be answering. The regions still contribute, through the precision of each sample's score, but they are not counted as replicates.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom limma lmFit contrasts.fit
#' @importFrom stats qt pt
#'
#' @keywords internal

.sampleSetEffect <-
  function(expressionMatrix,
           setIndex,
           backgroundIndex,
           design,
           contrastVector,
           level = 0.95) {

    # One number per library: the set against its comparison, both read inside the same sample
    sampleDelta <- colMeans(expressionMatrix[setIndex, , drop = FALSE]) -
      colMeans(expressionMatrix[backgroundIndex, , drop = FALSE])

    emptyResult <- list(sample.delta = sampleDelta,
                        delta = mean(sampleDelta),
                        standard.error = NA_real_,
                        df = NA_real_,
                        p.value = NA_real_,
                        ci.lower = NA_real_,
                        ci.upper = NA_real_)

    if ((nrow(design) - ncol(design)) < 1) {
      return(emptyResult)
    }

    # lmFit on a one-row matrix takes care of the pivoting and of contrasts on a rank-deficient design
    deltaMatrix <- matrix(data = sampleDelta, nrow = 1,
                          dimnames = list("setScore", names(sampleDelta)))

    linearFit <- try(limma::contrasts.fit(fit = limma::lmFit(object = deltaMatrix, design = design),
                                          contrasts = contrastVector),
                     silent = TRUE)

    if (inherits(linearFit, "try-error")) {
      warning("The set score could not be fitted on the design, no sample level interval is reported.", call. = FALSE)
      return(emptyResult)
    }

    estimateValue <- as.numeric(linearFit$coefficients[1])
    standardError <- as.numeric(linearFit$stdev.unscaled[1]) * as.numeric(linearFit$sigma[1])
    degreesFreedom <- as.numeric(linearFit$df.residual[1])

    # A design that fits the scores exactly leaves no spread to build an interval from
    if (!is.finite(standardError) | standardError <= 0) {
      emptyResult$delta <- estimateValue
      return(emptyResult)
    }

    criticalValue <- stats::qt(1 - (1 - level) / 2, df = degreesFreedom)

    return(list(sample.delta = sampleDelta,
                delta = estimateValue,
                standard.error = standardError,
                df = degreesFreedom,
                p.value = 2 * stats::pt(-abs(estimateValue / standardError), df = degreesFreedom),
                ci.lower = estimateValue - criticalValue * standardError,
                ci.upper = estimateValue + criticalValue * standardError))
  } # END function
