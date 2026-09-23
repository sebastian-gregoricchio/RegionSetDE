# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title testRegions
#'
#' @description Tests a contrast on a \code{RegionSetDE.fit} object and returns one row per region. When the counts were tiled, every tile is tested on its own and the p-values are then combined back to the region, so that the region stays the unit of inference even though the signal was measured at a finer scale.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrast Contrast to test, given in one of four ways. A character vector of length three, \code{c("column", "groupA", "groupB")}, naming a column of the \code{colData} and two of its levels, which is the form to reach for when the design uses a reference level. A string with the name of a design column, e.g. \code{"conditionCOMBO"}. A string written as an expression over the design columns, e.g. \code{"conditionCOMBO - conditionEPZ"}. Or a numeric vector with one coefficient per column of the design. A named list of any of these runs every contrast on the same fit and returns a \code{RegionSetDE.resultsList}.
#' @param combine Logical value to indicate whether the tile level p-values must be combined into one value per region. Ignored when the counts were not tiled. Default: \code{TRUE}.
#' @param combineMethod String with the method used to combine the tiles into their region. \code{"simes"}, through \code{csaw::combineTests}, asks whether any part of the region changes, and a single strong tile is enough. \code{"holm-min"}, through \code{csaw::minimalTests}, asks for several tiles to change together, three of them or 40\% of the region when that is more, and all of them in a region shorter than three tiles, which suits broad domains where one tile moving on its own is more likely noise than biology. Default: \code{"simes"}.
#' @param lfcThreshold Numeric value with the log2 fold change against which the null hypothesis is tested. A value above zero moves the threshold inside the test, through \code{edgeR::glmTreat}, \code{limma::treat} or the \code{lfcThreshold} of \code{DESeq2::results}, which is stricter and better calibrated than filtering the output afterwards. Default: \code{0}.
#' @param FDR Numeric value with the adjusted p-value cut-off used to fill the \code{diff.status} column. Default: \code{0.05}.
#' @param log2FC Numeric value with the absolute log2 fold change cut-off used to fill the \code{diff.status} column. Default: \code{0}.
#' @param adjustMethod String with the multiple testing correction, passed to \code{stats::p.adjust}. Default: \code{"BH"}.
#' @param regionSets Character vector with the names of the region sets to keep in the output. Default: \code{NULL}, all of them.
#' @param signalBy String with a column of the \code{colData}: every level of it gets an \code{average.signal.<level>} column in the results, the average signal over the samples of that level alone. \code{FALSE} adds none. Default: \code{NULL}, the column the contrast compares two levels of, which is known for the three-element form and for most coefficients and expressions, and none when it is not.
#' @param extraColumns Annotation carried by the regions that must be appended to the result, at the end of the table. Either \code{TRUE} for every column of the \code{rowData} beyond the ones the package writes itself, \code{FALSE} for none, or a character vector naming the ones wanted. Default: \code{TRUE}.
#' @param carryCounts Logical value to indicate whether the counts must travel inside the result, so that \code{\link{plotRegion}} and \code{\link{plotTopHeatmap}} can draw the values without being handed the counts object again. Several contrasts run on one fit share the same copy in memory. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.results} object. Its table, read with \code{\link{resultsTable}}, holds one row per region with the coordinates, then:
#' \itemize{
#'   \item \code{log2FC}: the log2 fold change, first group of the contrast over the second.
#'   \item \code{average.signal}: the average abundance of the region over all the samples, on the scale of the engine: the log2 counts per million of \code{edgeR::aveLogCPM} for edgeR, the average of the log2 values the linear model was fitted on for voom, limma and dream, \code{log2(baseMean + 1)} for DESeq2.
#'   \item \code{average.signal.<level>}: the same quantity computed on the samples of one level of \code{signalBy} only, one column per level.
#'   \item \code{stat}: the test statistic of the engine: the quasi-likelihood F of edgeR, or its likelihood ratio when the dispersion was fixed, the moderated t of voom, limma and dream, the Wald statistic of DESeq2.
#'   \item \code{stat.distribution}: the distribution \code{stat} follows under the null hypothesis, one among \code{"f"}, \code{"chisq"}, \code{"t"} and \code{"norm"}. \code{NA} for the threshold tests run when \code{lfcThreshold > 0}, whose null is not centred on zero.
#'   \item \code{df1}, \code{df2}: the degrees of freedom of that distribution. For the F of edgeR, \code{df1} is the numerator, the number of coefficients tested, and \code{df2} the denominator, the residual degrees of freedom plus the prior ones. For the likelihood ratio \code{df1} is the number of coefficients tested. For the moderated t \code{df1} is the total degrees of freedom, residual plus prior. The normal distribution of DESeq2 has none. \code{NA} where the distribution does not use them.
#'   \item \code{p.value}, \code{FDR}: the p-value and its adjustment over all the rows of the contrast.
#'   \item \code{diff.status}: \code{"up"}, \code{"down"} or \code{"null"}, from \code{FDR} and \code{log2FC}.
#' }
#' On a tiled object the statistics, the degrees of freedom and the averages come from the tile carrying the p-value of the region, followed by the columns the combination adds.
#'
#' @details The multiple testing correction is applied over all the rows of the object, across the region sets, and \code{regionSets} subsets the output afterwards. Correcting inside each set separately would make the FDR of a set depend on how many other sets were loaded, which is not a property anyone wants in a result.
#'
#' Two things follow from the combination step. With the default \code{combineMethod} the p-value of a tiled region is a Simes combination, so it answers "does any part of this region change" rather than "does the whole region change", and a long domain that moves over one tile out of forty will come out with a small p-value and a small overall fold change. The \code{log2FC} reported for a combined region is the fold change of the most significant tile, not an average, which is the quantity that matches the p-value. The tile level table stays available in the \code{tiles} slot, and \code{\link{plotRegion}} draws it.
#'
#' A design written as \code{~ condition} spends one coefficient per level except the first, so a level can be a coefficient in the design or the reference the others are measured against, depending on how the factor was ordered. Naming a coefficient that turns out to be the reference is the usual source of confusion, and it is what \code{c("column", "groupA", "groupB")} avoids: that form averages the design rows of each group and takes the difference, which gives the same contrast whatever the reference is and whether the design was written as \code{~ condition} or \code{~ 0 + condition}. With other covariates in the design the averaging picks up their imbalance between the two groups, so it describes what it says only when the design is reasonably balanced.
#'
#' Whatever the regions were loaded with travels through to the result. A gene name, a peak score or any other column attached to the \code{rowData} comes out at the end of the table, which is what makes \code{topRegions()} readable and lets \code{plotVolcano(labelColumn = )} label the points with something other than an identifier. On a tiled object the value is read off the tile the combination reported, the same one the fold change comes from, so a row describes one place rather than an average over several.
#'
#' The statistic, the distribution it follows and its degrees of freedom are in the table so that the test can be taken further, into a power or sample size analysis for instance. \code{\link{contrastInfo}} gathers the rest of what such an analysis needs, the engine and the number of samples in each group of the contrast, which is also stored in the \code{n.samples} element of the \code{contrast.groups} slot.
#'
#' The \code{diff.status} column is a labelling convenience, not a claim. It is filled from \code{FDR} and \code{log2FC} and used by the plotting functions; the thresholds are stored in the object so that a figure can state them.
#'
#' @examples
#' \dontrun{
#' fit <- fitRegions(counts, design = ~ replicate + condition, engine = "edgeR")
#'
#' res <- testRegions(fit, contrast = "conditionCOMBO")
#'
#' # Two levels of a column, whichever of them the design took as reference
#' res <- testRegions(fit, contrast = c("condition", "COMBO", "DMSO"))
#'
#' # Difference between two coefficients of the design
#' res <- testRegions(fit, contrast = "conditionCOMBO - conditionEPZ")
#'
#' # Several contrasts on the same fit
#' resList <- testRegions(fit, contrast = list(combo = c("condition", "COMBO", "DMSO"),
#'                                             epz = c("condition", "EPZ", "DMSO")))
#' resList
#' topRegions(resList, contrast = "combo")
#'
#' # Threshold inside the test rather than on the output
#' resStrict <- testRegions(fit, contrast = "conditionCOMBO", lfcThreshold = 1)
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{fitRegions}}, \code{\link{topRegions}}, \code{\link{plotVolcano}}
#'
#' @importFrom SummarizedExperiment colData rowData rowRanges
#' @importFrom GenomicRanges GRanges split
#' @importFrom BiocGenerics width
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom stats p.adjust
#' @importFrom dplyr mutate filter select arrange left_join case_when
#' @importFrom rlang .data
#' @importFrom methods is new
#'
#' @export testRegions

testRegions <-
  function(fit,
           contrast,
           combine = TRUE,
           combineMethod = "simes",
           lfcThreshold = 0,
           FDR = 0.05,
           log2FC = 0,
           adjustMethod = "BH",
           regionSets = NULL,
           signalBy = NULL,
           extraColumns = TRUE,
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

      resultsList <-
        lapply(names(contrast),
               function(contrastName) {
                 if (isTRUE(verbose)) {
                   message(paste0("--- ", contrastName, " ---"))
                 }
                 return(testRegions(fit = fit, contrast = contrast[[contrastName]], combine = combine,
                                    combineMethod = combineMethod, lfcThreshold = lfcThreshold, FDR = FDR,
                                    log2FC = log2FC, adjustMethod = adjustMethod, regionSets = regionSets,
                                    signalBy = signalBy, extraColumns = extraColumns, carryCounts = carryCounts,
                                    verbose = verbose))
               })

      names(resultsList) <- names(contrast)

      return(new(Class = "RegionSetDE.resultsList",
                 results = resultsList,
                 contrasts = names(contrast)))
    }

    if (!(combineMethod %in% c("simes", "holm-min"))) {
      stop("The 'combineMethod' parameter must be either 'simes' or 'holm-min'.", call. = FALSE)
    }

    if (lfcThreshold < 0) {
      stop("The 'lfcThreshold' parameter cannot be negative.", call. = FALSE)
    }

    contrastObject <- .resolveContrast(contrast = contrast,
                                       design = fit@design,
                                       colData = SummarizedExperiment::colData(fit@counts))

    if (isTRUE(verbose)) {
      message(paste0("Testing '", contrastObject$label, "' on ", nrow(fit@counts), " ", fit@counting.level, "s."))

      # The p-values are conditional on a number that was assumed rather than measured
      if (isTRUE(fit@dispersion$no.replicates)) {
        message("The fit has no replicates, so the p-values rest entirely on the supplied dispersion. Read them next to checkNullCalibration().")
      }
    }

    #-------------------------------#
    # Engine specific test          #
    #-------------------------------#
    rawTable <- switch(fit@engine,
                       "edgeR" = .testEdgeR(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = lfcThreshold),
                       "voom" = .testVoom(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = lfcThreshold),
                       "limma" = .testVoom(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = lfcThreshold, trend = TRUE),
                       "dream" = .testDream(fit = fit, contrastObject = contrastObject, lfcThreshold = lfcThreshold, verbose = verbose),
                       "deseq2" = .testDESeq2(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = lfcThreshold))

    #-------------------------------#
    # Average signal per group      #
    #-------------------------------#
    # The column the contrast compares, unless another one was asked for or none at all
    signalColumn <- if (isFALSE(signalBy)) {NULL} else if (is.null(signalBy)) {contrastObject$column} else {signalBy}
    groupSignal <- .groupAverageSignal(fit = fit, column = signalColumn)
    groupSignalColumns <- colnames(groupSignal)

    if (length(groupSignalColumns) > 0) {
      rawTable <- cbind(rawTable, groupSignal)
    }

    # How many samples stand behind each side of the contrast, which a power analysis starts from
    if (!is.null(contrastObject$column)) {
      groupValues <- as.character(SummarizedExperiment::colData(fit@counts)[[contrastObject$column]])
      contrastObject$n.samples <- vapply(contrastObject$groups, function(group) {sum(groupValues == group, na.rm = TRUE)}, integer(1))
    }

    #-------------------------------#
    # Attach the row annotation     #
    #-------------------------------#
    rowTable <- as.data.frame(SummarizedExperiment::rowData(fit@counts))
    rowRangesObject <- SummarizedExperiment::rowRanges(fit@counts)

    rawTable <- dplyr::mutate(rawTable,
                              region.set = rowTable$region.set,
                              region.id = rowTable$region.id,
                              tile.id = rowTable$tile.id,
                              region.key = paste(rowTable$region.set, rowTable$region.id, sep = "|"),
                              seqnames = as.character(GenomeInfoDb::seqnames(rowRangesObject)),
                              start = BiocGenerics::start(rowRangesObject),
                              end = BiocGenerics::end(rowRangesObject),
                              width = BiocGenerics::width(rowRangesObject))

    #-------------------------------#
    # Annotation carried by the rows #
    #-------------------------------#
    extraTable <- .extraRowColumns(rowTable = rowTable,
                                   extraColumns = extraColumns,
                                   reserved = c(colnames(rawTable), .resultColumnNames()),
                                   verbose = verbose)

    if (ncol(extraTable) > 0) {
      rawTable <- cbind(rawTable, extraTable)
    }

    isTiled <- fit@counting.level == "tile"

    #-------------------------------#
    # Tiles to regions, or not      #
    #-------------------------------#
    if (isTiled & isTRUE(combine)) {
      combinedList <- .combineTiles(tileTable = rawTable,
                                    tileRanges = rowRangesObject,
                                    extraColumns = c("stat.distribution", "df1", "df2", groupSignalColumns, colnames(extraTable)),
                                    method = combineMethod,
                                    adjustMethod = adjustMethod,
                                    verbose = verbose)

      resultTable <- combinedList$results
      regionRanges <- combinedList$regions
      tileTable <- dplyr::mutate(rawTable, FDR = stats::p.adjust(rawTable$p.value, method = adjustMethod))
      combinationInfo <- list(applied = TRUE, method = combineMethod)

    } else {
      resultTable <- dplyr::mutate(rawTable, FDR = stats::p.adjust(rawTable$p.value, method = adjustMethod))
      regionRanges <- rowRangesObject
      S4Vectors::mcols(regionRanges) <- NULL
      tileTable <- if (isTiled) {resultTable} else {data.frame()}
      combinationInfo <- list(applied = FALSE, method = NA_character_)
    }

    #-------------------------------#
    # Label and tidy the output     #
    #-------------------------------#
    # 'FDR' and 'log2FC' also name two columns of the table, and dplyr reads the column before the argument
    FDRthreshold <- FDR
    log2FCthreshold <- log2FC

    resultTable <- dplyr::mutate(resultTable,
                                 diff.status = dplyr::case_when(.data$FDR < FDRthreshold & .data$log2FC > log2FCthreshold ~ "up",
                                                                .data$FDR < FDRthreshold & .data$log2FC < (-log2FCthreshold) ~ "down",
                                                                TRUE ~ "null"))
    resultTable$diff.status <- factor(resultTable$diff.status, levels = c("down", "null", "up"))

    if (!is.null(regionSets)) {
      absentSets <- setdiff(regionSets, unique(resultTable$region.set))
      if (length(absentSets) > 0) {
        stop(paste0("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), "."), call. = FALSE)
      }

      # The correction has already run over every row, subsetting here does not change the FDR of what is kept
      keptRows <- resultTable$region.set %in% regionSets
      resultTable <- resultTable[keptRows, , drop = FALSE]
      regionRanges <- regionRanges[keptRows]

      if (nrow(tileTable) > 0) {
        tileTable <- dplyr::filter(tileTable, .data$region.set %in% regionSets)
      }
    }

    # The statistics first, then what the combination added, then whatever the regions came with
    statisticColumns <- c("region.set", "region.id", "tile.id", "seqnames", "start", "end", "width",
                          "log2FC", "average.signal", groupSignalColumns,
                          "stat", "stat.distribution", "df1", "df2", "p.value", "FDR", "diff.status")
    combinationColumns <- c("n.tiles", "n.tiles.up", "n.tiles.down", "direction", "rep.tile.start", "rep.tile.end")

    annotationColumns <- setdiff(colnames(resultTable),
                                 c(statisticColumns, combinationColumns, "region.key"))

    resultTable <- resultTable[, c(intersect(statisticColumns, colnames(resultTable)),
                                   intersect(combinationColumns, colnames(resultTable)),
                                   annotationColumns),
                               drop = FALSE]
    rownames(resultTable) <- NULL

    #-------------------------------#
    # Counts travelling along        #
    #-------------------------------#
    carriedCounts <- new(Class = "RegionSetDE.counts")

    if (isTRUE(carryCounts)) {
      carriedCounts <- fit@counts

      # The rows stay at the level the model was fitted on, which is what plotRegion needs to draw a tiled profile
      if (!is.null(regionSets)) {
        carriedCounts <- carriedCounts[SummarizedExperiment::rowData(carriedCounts)$region.set %in% regionSets, ]
      }
    }

    #-------------------------------#
    # Assemble the object           #
    #-------------------------------#
    resultsObject <- new(Class = "RegionSetDE.results",
                         counts = carriedCounts,
                         results = resultTable,
                         tiles = tileTable,
                         regions = regionRanges,
                         contrast = contrastObject$label,
                         contrast.vector = contrastObject$vector,
                         contrast.groups = contrastObject[intersect(c("column", "groups", "n.samples"), names(contrastObject))],
                         engine = fit@engine,
                         counting.level = fit@counting.level,
                         combination = combinationInfo,
                         thresholds = list(FDR = FDR, log2FC = log2FC, lfcThreshold = lfcThreshold, adjust.method = adjustMethod),
                         blacklist = fit@blacklist,
                         whitelist = fit@whitelist,
                         genome.assembly = fit@genome.assembly,
                         seqlevels.style = fit@seqlevels.style,
                         filtering.log = fit@filtering.log,
                         parameters = c(fit@parameters,
                                        list(testRegions = list(contrast = contrastObject$label,
                                                                combine = combine,
                                                                combineMethod = combineMethod,
                                                                lfcThreshold = lfcThreshold,
                                                                FDR = FDR,
                                                                log2FC = log2FC,
                                                                adjustMethod = adjustMethod,
                                                                signalBy = signalColumn,
                                                                carryCounts = carryCounts))))

    if (isTRUE(verbose)) {
      statusTable <- table(resultTable$diff.status)
      message(paste0("Done. ", statusTable[["up"]], " up and ", statusTable[["down"]],
                     " down out of ", nrow(resultTable), " regions (FDR < ", FDR,
                     if (log2FC > 0) {paste0(", |log2FC| > ", log2FC)} else {""}, ")."))
    }

    return(resultsObject)
  } # END function




#' @title .resolveContrast
#'
#' @description Turns the \code{contrast} argument of \code{\link{testRegions}} into a numeric vector over the columns of the design.
#'
#' @param contrast String with a coefficient name or an expression over the coefficients, a character vector of length three naming a column and two of its levels, or a numeric vector.
#' @param design Design matrix.
#' @param colData \code{DataFrame} or data.frame with the sample metadata, needed by the three-element form. Default: \code{NULL}.
#'
#' @return A list with the \code{vector} of coefficients, a \code{label} describing the contrast, and the \code{column} of the \code{colData} and the two \code{groups} it separates when the contrast turns out to be a difference between two levels of one variable.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom limma makeContrasts
#'
#' @keywords internal

.resolveContrast <-
  function(contrast,
           design,
           colData = NULL) {

    coefficientNames <- colnames(design)

    #-------------------------------#
    # A vector is taken as it is    #
    #-------------------------------#
    if (is.numeric(contrast)) {
      if (length(contrast) != length(coefficientNames)) {
        stop(paste0("The contrast vector must have one value per design column (", length(coefficientNames), ")."), call. = FALSE)
      }
      contrastVector <- as.numeric(contrast)
      names(contrastVector) <- coefficientNames

      nonZero <- coefficientNames[contrastVector != 0]
      return(c(list(vector = contrastVector, label = paste(nonZero, collapse = " vs ")),
               .contrastGroups(contrastVector = contrastVector, design = design, colData = colData)))
    }

    if (!is.character(contrast)) {
      stop("The \'contrast\' parameter must be a string, a character vector of length three, or a numeric vector.", call. = FALSE)
    }

    #-------------------------------#
    # c("column", "groupA", "groupB")
    #-------------------------------#
    if (length(contrast) == 3) {
      if (is.null(colData)) {
        stop("The three-element form of \'contrast\' needs the sample metadata, which is not available here.", call. = FALSE)
      }

      colTable <- as.data.frame(colData)
      columnName <- contrast[1]

      if (!(columnName %in% colnames(colTable))) {
        stop(paste0("The column \'", columnName, "\' is absent from the colData. Available: ",
                    paste(colnames(colTable), collapse = ", "), "."), call. = FALSE)
      }

      columnValues <- as.character(colTable[[columnName]])
      absentGroups <- setdiff(contrast[2:3], unique(columnValues))
      if (length(absentGroups) > 0) {
        stop(paste0("The following levels are absent from \'", columnName, "\': ", paste(absentGroups, collapse = ", "),
                    ". Available: ", paste(unique(columnValues), collapse = ", "), "."), call. = FALSE)
      }

      firstRows <- which(columnValues == contrast[2])
      secondRows <- which(columnValues == contrast[3])

      # Averaging the design rows of each group gives the same contrast whatever the reference level is,
      # and works identically for ~ condition and ~ 0 + condition
      contrastVector <- colMeans(design[firstRows, , drop = FALSE]) - colMeans(design[secondRows, , drop = FALSE])
      contrastVector[abs(contrastVector) < 1e-10] <- 0
      names(contrastVector) <- coefficientNames

      if (all(contrastVector == 0)) {
        stop(paste0("The design does not separate \'", contrast[2], "\' from \'", contrast[3],
                    "\', the two groups share the same coefficients."), call. = FALSE)
      }

      return(list(vector = contrastVector,
                  label = paste0(columnName, ": ", contrast[2], " vs ", contrast[3]),
                  column = columnName,
                  groups = c(contrast[2], contrast[3])))
    }

    if (length(contrast) != 1) {
      stop("The \'contrast\' parameter must hold one string, or three when naming a column and two of its levels.", call. = FALSE)
    }

    #-------------------------------#
    # A plain coefficient name      #
    #-------------------------------#
    if (contrast %in% coefficientNames) {
      contrastVector <- as.numeric(coefficientNames == contrast)
      names(contrastVector) <- coefficientNames
      return(c(list(vector = contrastVector, label = contrast),
               .contrastGroups(contrastVector = contrastVector, design = design, colData = colData)))
    }

    #-------------------------------#
    # An expression on the names    #
    #-------------------------------#
    # makeContrasts needs syntactic names, the design columns are renamed and put back afterwards
    safeNames <- make.names(coefficientNames)
    safeContrast <- contrast
    for (i in order(nchar(coefficientNames), decreasing = TRUE)) {
      safeContrast <- gsub(pattern = coefficientNames[i], replacement = safeNames[i], x = safeContrast, fixed = TRUE)
    }

    contrastMatrix <- try(limma::makeContrasts(contrasts = safeContrast, levels = safeNames), silent = TRUE)

    if (inherits(contrastMatrix, "try-error")) {
      stop(paste0("The contrast \'", contrast, "\' could not be read. Available coefficients: ",
                  paste(coefficientNames, collapse = ", "), ".",
                  .contrastSuggestion(contrast = contrast, coefficientNames = coefficientNames, colData = colData)), call. = FALSE)
    }

    contrastVector <- as.numeric(contrastMatrix[, 1])
    names(contrastVector) <- coefficientNames

    return(c(list(vector = contrastVector, label = contrast),
             .contrastGroups(contrastVector = contrastVector, design = design, colData = colData)))
  } # END function




#' @title .contrastGroups
#'
#' @description Works out which variable of the sample metadata a contrast separates, and which two of its levels, by comparing the contrast against the difference between the design rows of every pair of levels. The variables of the design are tried first. A column holding a different value for every sample, such as the sample names, is never tried unless it is in the design, since any two samples from two groups reproduce a contrast between those groups.
#'
#' @param contrastVector Numeric vector with the contrast.
#' @param design Design matrix.
#' @param colData Sample metadata, or \code{NULL}.
#' @param maxLevels Numeric value with the number of levels above which a column is not considered a grouping variable. Default: \code{20}.
#'
#' @return A list with the \code{column} and the two \code{groups}, the first one being the level the contrast is positive for. An empty list when the contrast is not a difference between two levels of one variable.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.contrastGroups <-
  function(contrastVector,
           design,
           colData = NULL,
           maxLevels = 20) {

    if (is.null(colData)) {
      return(list())
    }

    colTable <- as.data.frame(colData)

    # A variable of the design names its coefficients after itself and its levels, e.g. conditionSHR
    isDesignColumn <- vapply(colnames(colTable),
                             function(columnName) {
                               any(paste0(columnName, unique(as.character(colTable[[columnName]]))) %in% colnames(design))
                             },
                             logical(1))

    # A column with one value per sample only names the samples, and any two of them taken from two groups match the contrast
    isSampleLabel <- vapply(colnames(colTable),
                            function(columnName) {
                              columnValues <- as.character(colTable[[columnName]])
                              anyDuplicated(columnValues[!is.na(columnValues)]) == 0
                            },
                            logical(1))

    candidateColumns <- c(colnames(colTable)[isDesignColumn],
                          colnames(colTable)[!isDesignColumn & !isSampleLabel])

    for (columnName in candidateColumns) {
      columnValues <- as.character(colTable[[columnName]])
      columnLevels <- unique(columnValues[!is.na(columnValues)])

      if (length(columnLevels) < 2 | length(columnLevels) > maxLevels) {
        next
      }

      # The contrast built by the three-element form is exactly this difference, whatever the reference level was
      levelMeans <- lapply(columnLevels, function(x) {colMeans(design[columnValues == x, , drop = FALSE])})
      names(levelMeans) <- columnLevels

      for (firstLevel in columnLevels) {
        for (secondLevel in setdiff(columnLevels, firstLevel)) {
          candidateVector <- levelMeans[[firstLevel]] - levelMeans[[secondLevel]]

          if (max(abs(candidateVector - contrastVector)) < 1e-8) {
            return(list(column = columnName, groups = c(firstLevel, secondLevel)))
          }
        }
      }
    }

    return(list())
  } # END function




#' @title .contrastSuggestion
#'
#' @description Builds the second half of the error message raised when a contrast cannot be read, pointing at the reference level when that is what went wrong.
#'
#' @param contrast String with the contrast the user wrote.
#' @param coefficientNames Character vector with the columns of the design.
#' @param colData Sample metadata, or \code{NULL}.
#'
#' @return A string, empty when nothing useful can be said.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.contrastSuggestion <-
  function(contrast,
           coefficientNames,
           colData = NULL) {

    if (is.null(colData)) {
      return("")
    }

    colTable <- as.data.frame(colData)

    #-------------------------------#
    # Look for a reference level    #
    #-------------------------------#
    # A level that the design took as reference has no coefficient, so writing its name never resolves
    for (columnName in colnames(colTable)) {
      columnValues <- unique(as.character(colTable[[columnName]]))
      if (length(columnValues) < 2) {
        next
      }

      writtenLevels <- columnValues[vapply(columnValues,
                                           function(x) {grepl(paste0(columnName, x), contrast, fixed = TRUE)},
                                           logical(1))]
      missingLevels <- writtenLevels[!(paste0(columnName, writtenLevels) %in% coefficientNames)]

      if (length(missingLevels) > 0) {
        otherLevels <- setdiff(writtenLevels, missingLevels)
        secondLevel <- if (length(otherLevels) > 0) {otherLevels[1]} else {setdiff(columnValues, missingLevels)[1]}

        return(paste0("\n  \'", missingLevels[1], "\' is the reference level of \'", columnName,
                      "\' and has no coefficient of its own. Write contrast = c(\"", columnName, "\", \"",
                      missingLevels[1], "\", \"", secondLevel, "\") instead."))
      }
    }

    return("")
  } # END function




#' @title .testEdgeR
#'
#' @description Runs the quasi-likelihood F test, or the threshold test, on an \code{edgeR} fit.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrastVector Numeric vector with the contrast.
#' @param lfcThreshold Numeric value with the log2 fold change of the null hypothesis.
#'
#' @return A data.frame with the \code{log2FC}, \code{average.signal}, \code{stat}, \code{stat.distribution}, \code{df1}, \code{df2} and \code{p.value} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom edgeR glmQLFTest glmLRT glmTreat topTags
#'
#' @keywords internal

.testEdgeR <-
  function(fit,
           contrastVector,
           lfcThreshold = 0) {

    # A fixed dispersion leaves nothing for the quasi-likelihood F to account for, so the test is a likelihood ratio
    isLikelihoodRatio <- identical(fit@fit$test, "lrt")

    if (lfcThreshold > 0) {
      testObject <- edgeR::glmTreat(glmfit = fit@fit$object, contrast = contrastVector, lfc = lfcThreshold)
    } else if (isTRUE(isLikelihoodRatio)) {
      testObject <- edgeR::glmLRT(glmfit = fit@fit$object, contrast = contrastVector)
    } else {
      testObject <- edgeR::glmQLFTest(glmfit = fit@fit$object, contrast = contrastVector)
    }

    # sort.by = "none" keeps the rows aligned with the object, the annotation is bound by position later
    topTable <- edgeR::topTags(object = testObject, n = Inf, sort.by = "none", adjust.method = "none")$table

    statColumn <- if ("F" %in% colnames(topTable)) {
      topTable$F
    } else if ("LR" %in% colnames(topTable)) {
      topTable$LR
    } else {
      rep(NA_real_, nrow(topTable))
    }

    # The F of the quasi-likelihood test has the coefficients tested on top and the residual plus prior
    # degrees of freedom below, the likelihood ratio only the former. glmTreat has no statistic to describe.
    rowNumber <- nrow(topTable)
    statDistribution <- if ("F" %in% colnames(topTable)) {"f"} else if ("LR" %in% colnames(topTable)) {"chisq"} else {NA_character_}
    firstDf <- if (is.null(testObject$df.test)) {rep(NA_real_, rowNumber)} else {rep_len(as.numeric(testObject$df.test), rowNumber)}
    secondDf <- if (identical(statDistribution, "f") & !is.null(testObject$df.total)) {rep_len(as.numeric(testObject$df.total), rowNumber)} else {rep(NA_real_, rowNumber)}

    return(data.frame(log2FC = topTable$logFC,
                      average.signal = topTable$logCPM,
                      stat = statColumn,
                      stat.distribution = rep(statDistribution, rowNumber),
                      df1 = firstDf,
                      df2 = secondDf,
                      p.value = topTable$PValue,
                      stringsAsFactors = FALSE))
  } # END function




#' @title .testVoom
#'
#' @description Runs the moderated t test on a \code{limma} fit.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrastVector Numeric vector with the contrast.
#' @param lfcThreshold Numeric value with the log2 fold change of the null hypothesis.
#'
#' @return A data.frame with the \code{log2FC}, \code{average.signal}, \code{stat}, \code{stat.distribution}, \code{df1}, \code{df2} and \code{p.value} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom limma contrasts.fit eBayes treat topTable topTreat
#'
#' @keywords internal

.testVoom <-
  function(fit,
           contrastVector,
           lfcThreshold = 0,
           trend = FALSE) {

    contrastFit <- limma::contrasts.fit(fit = fit@fit$object, contrasts = contrastVector)

    # voom carries its mean-variance relationship in the weights, limma-trend carries it here instead
    if (lfcThreshold > 0) {
      contrastFit <- limma::treat(fit = contrastFit, lfc = lfcThreshold, robust = isTRUE(fit@fit$robust), trend = trend)
      topTable <- limma::topTreat(fit = contrastFit, number = Inf, sort.by = "none", adjust.method = "none")
    } else {
      contrastFit <- limma::eBayes(fit = contrastFit, robust = isTRUE(fit@fit$robust), trend = trend)
      topTable <- limma::topTable(fit = contrastFit, number = Inf, sort.by = "none", adjust.method = "none")
    }

    # The moderated t follows a t with the residual and the prior degrees of freedom together, TREAT does not
    rowNumber <- nrow(topTable)

    return(data.frame(log2FC = topTable$logFC,
                      average.signal = topTable$AveExpr,
                      stat = topTable$t,
                      stat.distribution = rep(if (lfcThreshold > 0) {NA_character_} else {"t"}, rowNumber),
                      df1 = rep_len(as.numeric(contrastFit$df.total), rowNumber),
                      df2 = rep(NA_real_, rowNumber),
                      p.value = topTable$P.Value,
                      stringsAsFactors = FALSE))
  } # END function




#' @title .testDream
#'
#' @description Runs the test on a \code{dream} fit. A contrast that is not a single coefficient of the design needs the mixed model to be fitted again, since \code{dream} builds the contrast at fit time.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrastObject List returned by \code{.resolveContrast}.
#' @param lfcThreshold Numeric value with the log2 fold change of the null hypothesis.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A data.frame with the \code{log2FC}, \code{average.signal}, \code{stat}, \code{stat.distribution}, \code{df1}, \code{df2} and \code{p.value} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.testDream <-
  function(fit,
           contrastObject,
           lfcThreshold = 0,
           verbose = TRUE) {

    if (lfcThreshold > 0) {
      stop("The 'lfcThreshold' parameter is not available for the 'dream' engine, filter on 'log2FC' instead.", call. = FALSE)
    }

    contrastVector <- contrastObject$vector
    isSingleCoefficient <- sum(contrastVector != 0) == 1 & all(contrastVector %in% c(0, 1))

    if (isSingleCoefficient) {
      mixedFit <- fit@fit$object
      coefficientName <- names(contrastVector)[contrastVector == 1]

    } else {
      # dream needs the contrast before fitting, so an arbitrary one costs a second pass over the rows
      if (isTRUE(verbose)) {
        message("The contrast is not a single coefficient, the mixed model is being fitted again.")
      }

      contrastMatrix <- variancePartition::makeContrastsDream(formula = fit@fit$formula,
                                                              data = fit@fit$data,
                                                              contrasts = stats::setNames(contrastObject$label, "contrast"))

      mixedFit <- variancePartition::dream(exprObj = fit@fit$voom,
                                           formula = fit@fit$formula,
                                           data = fit@fit$data,
                                           L = contrastMatrix,
                                           BPPARAM = fit@fit$BPPARAM,
                                           quiet = TRUE)
      mixedFit <- variancePartition::eBayes(mixedFit)
      coefficientName <- "contrast"
    }

    topTable <- variancePartition::topTable(fit = mixedFit, coef = coefficientName,
                                            number = Inf, sort.by = "none", adjust.method = "none")

    # dream gives every region its own Satterthwaite degrees of freedom, one column per coefficient
    totalDf <- if (!is.null(mixedFit$df.total)) {mixedFit$df.total} else {mixedFit$df.residual}
    if (is.matrix(totalDf)) {
      totalDf <- totalDf[, if (coefficientName %in% colnames(totalDf)) {coefficientName} else {1}]
    }

    rowNumber <- nrow(topTable)

    return(data.frame(log2FC = topTable$logFC,
                      average.signal = topTable$AveExpr,
                      stat = topTable$t,
                      stat.distribution = rep("t", rowNumber),
                      df1 = if (is.null(totalDf)) {rep(NA_real_, rowNumber)} else {rep_len(as.numeric(totalDf), rowNumber)},
                      df2 = rep(NA_real_, rowNumber),
                      p.value = topTable$P.Value,
                      stringsAsFactors = FALSE))
  } # END function




#' @title .testDESeq2
#'
#' @description Runs the Wald test on a \code{DESeq2} fit.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrastVector Numeric vector with the contrast.
#' @param lfcThreshold Numeric value with the log2 fold change of the null hypothesis.
#'
#' @return A data.frame with the \code{log2FC}, \code{average.signal}, \code{stat}, \code{stat.distribution}, \code{df1}, \code{df2} and \code{p.value} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.testDESeq2 <-
  function(fit,
           contrastVector,
           lfcThreshold = 0) {

    # independentFiltering removes rows from the correction, which would break the alignment with the region annotation
    resultTable <- DESeq2::results(object = fit@fit$object,
                                   contrast = as.numeric(contrastVector),
                                   lfcThreshold = lfcThreshold,
                                   independentFiltering = FALSE,
                                   cooksCutoff = FALSE,
                                   pAdjustMethod = "none")

    # The Wald statistic is a standard normal under the null, which has no degrees of freedom
    rowNumber <- nrow(resultTable)

    return(data.frame(log2FC = resultTable$log2FoldChange,
                      average.signal = log2(resultTable$baseMean + 1),
                      stat = resultTable$stat,
                      stat.distribution = rep(if (lfcThreshold > 0) {NA_character_} else {"norm"}, rowNumber),
                      df1 = rep(NA_real_, rowNumber),
                      df2 = rep(NA_real_, rowNumber),
                      p.value = resultTable$pvalue,
                      stringsAsFactors = FALSE))
  } # END function




#' @title .combineTiles
#'
#' @description Combines the tile level statistics into one row per region, through \code{csaw::combineTests}.
#'
#' @param tileTable Data.frame with one row per tile, as returned by the engine specific test.
#' @param tileRanges \code{GRanges} with the coordinates of the tiles.
#' @param extraColumns Character vector with the annotation columns carried over from the tiles.
#' @param method String with the combination method.
#' @param directionFDR Numeric value with the false discovery rate, within each region, below which a tile counts as moving up or down in \code{n.tiles.up} and \code{n.tiles.down}. It is the \code{fc.threshold} of csaw, which despite its name is not a fold change. Default: \code{0.05}.
#' @param adjustMethod String with the multiple testing correction.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A list with the \code{results} data.frame and the \code{regions} \code{GRanges}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom csaw combineTests minimalTests
#' @importFrom GenomicRanges split
#' @importFrom BiocGenerics unlist
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom stats p.adjust
#' @importFrom dplyr mutate arrange
#' @importFrom rlang .data
#'
#' @keywords internal

.combineTiles <-
  function(tileTable,
           tileRanges,
           extraColumns = character(0),
           method = "simes",
           directionFDR = 0.05,
           adjustMethod = "BH",
           verbose = TRUE) {

    # combineTests reads the grouping as a factor, its level order is what the output rows follow
    tileGroups <- factor(tileTable$region.key, levels = unique(tileTable$region.key))

    tileStatistics <- data.frame(logFC = tileTable$log2FC,
                                 PValue = tileTable$p.value,
                                 stringsAsFactors = FALSE)

    # csaw has one function per procedure rather than a method argument, and both return the same
    # columns, so everything below reads them the same way whichever was used
    combinedTable <-
      if (method == "holm-min") {
        csaw::minimalTests(ids = tileGroups,
                           tab = tileStatistics,
                           pval.col = "PValue",
                           fc.col = "logFC",
                           fc.threshold = directionFDR)
      } else {
        csaw::combineTests(ids = tileGroups,
                           tab = tileStatistics,
                           pval.col = "PValue",
                           fc.col = "logFC",
                           fc.threshold = directionFDR)
      }

    combinedTable <- as.data.frame(combinedTable)
    regionKeys <- rownames(combinedTable)

    #-------------------------------#
    # Region span from the tiles    #
    #-------------------------------#
    rangeList <- GenomicRanges::split(x = tileRanges, f = tileGroups)
    regionRanges <- unlist(range(rangeList), use.names = FALSE)
    regionRanges <- regionRanges[match(regionKeys, names(rangeList))]
    S4Vectors::mcols(regionRanges) <- NULL

    #-------------------------------#
    # Representative tile           #
    #-------------------------------#
    # The p-value comes from the best tile, so the fold change reported next to it has to come from the same tile
    representativeIndex <- combinedTable$rep.test
    keySplit <- strsplit(regionKeys, split = "|", fixed = TRUE)

    resultTable <- data.frame(region.set = vapply(keySplit, function(x) {x[1]}, character(1)),
                              region.id = vapply(keySplit, function(x) {paste(x[-1], collapse = "|")}, character(1)),
                              region.key = regionKeys,
                              seqnames = as.character(GenomeInfoDb::seqnames(regionRanges)),
                              start = BiocGenerics::start(regionRanges),
                              end = BiocGenerics::end(regionRanges),
                              width = BiocGenerics::width(regionRanges),
                              log2FC = tileTable$log2FC[representativeIndex],
                              average.signal = tileTable$average.signal[representativeIndex],
                              stat = tileTable$stat[representativeIndex],
                              p.value = combinedTable$PValue,
                              n.tiles = combinedTable$num.tests,
                              n.tiles.up = combinedTable$num.up.logFC,
                              n.tiles.down = combinedTable$num.down.logFC,
                              direction = combinedTable$direction,
                              rep.tile.start = BiocGenerics::start(tileRanges)[representativeIndex],
                              rep.tile.end = BiocGenerics::end(tileRanges)[representativeIndex],
                              stringsAsFactors = FALSE)

    # The statistics come from the representative tile, so its annotation is the one that describes the row
    if (length(extraColumns) > 0) {
      resultTable <- cbind(resultTable, tileTable[representativeIndex, extraColumns, drop = FALSE])
    }

    resultTable <- dplyr::mutate(resultTable, FDR = stats::p.adjust(.data$p.value, method = adjustMethod))

    if (isTRUE(verbose)) {
      message(paste0(nrow(tileTable), " tiles combined into ", nrow(resultTable), " regions by ", method, "."))
    }

    return(list(results = resultTable, regions = regionRanges))
  } # END function




#' @title .resultColumnNames
#'
#' @description Lists the column names \code{\link{testRegions}} writes itself, so that a column carried by the regions can be spotted before it overwrites one of them.
#'
#' @return A character vector.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resultColumnNames <-
  function() {

    return(c("log2FC", "average.signal", "stat", "stat.distribution", "df1", "df2", "p.value", "FDR", "diff.status",
             "n.tiles", "n.tiles.up", "n.tiles.down", "direction",
             "rep.tile.start", "rep.tile.end"))
  } # END function




#' @title .groupAverageSignal
#'
#' @description Computes the average signal of every region over the samples of each level of a column of the \code{colData}, with the same function the engine uses for \code{average.signal} over all the samples, so that the columns can be read side by side: \code{edgeR::aveLogCPM} for edgeR, the mean of the log2 values of the linear model for voom, limma and dream, and \code{log2} of the mean normalised count plus one for DESeq2.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param column String with the name of the column, or \code{NULL}.
#'
#' @return A data.frame with one \code{average.signal.<level>} column per level and one row per row of the fit, or a data.frame with no column when \code{column} is \code{NULL}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom edgeR aveLogCPM
#'
#' @keywords internal

.groupAverageSignal <-
  function(fit,
           column) {

    rowNumber <- nrow(fit@counts)

    if (is.null(column)) {
      return(data.frame(row.names = seq_len(rowNumber))[, 0, drop = FALSE])
    }

    sampleTable <- as.data.frame(SummarizedExperiment::colData(fit@counts), optional = TRUE)

    if (!is.character(column) | length(column) != 1 || !(column %in% colnames(sampleTable))) {
      stop("The 'signalBy' parameter must be the name of a column of the colData, or FALSE.", call. = FALSE)
    }

    columnValues <- sampleTable[[column]]
    columnLevels <- if (is.factor(columnValues)) {levels(droplevels(columnValues))} else {unique(as.character(columnValues[!is.na(columnValues)]))}
    columnValues <- as.character(columnValues)

    # The same quantity as average.signal, computed on a subset of the samples
    signalList <-
      lapply(columnLevels,
             function(level) {
               sampleIndex <- which(columnValues == level)

               switch(fit@engine,
                      "edgeR" = as.numeric(edgeR::aveLogCPM(fit@fit$dge[, sampleIndex])),
                      "voom" = ,
                      "limma" = ,
                      "dream" = as.numeric(rowMeans(fit@fit$voom$E[, sampleIndex, drop = FALSE], na.rm = TRUE)),
                      "deseq2" = as.numeric(log2(rowMeans(DESeq2::counts(fit@fit$object, normalized = TRUE)[, sampleIndex, drop = FALSE]) + 1)))
             })

    # The whole name is made syntactic rather than the level, so that a level such as 4h does not become X4h
    signalTable <- as.data.frame(signalList, col.names = make.names(paste0("average.signal.", columnLevels)), optional = TRUE)

    return(signalTable)
  } # END function




#' @title .extraRowColumns
#'
#' @description Picks the annotation columns of the \code{rowData} that must travel into a result, leaving out the ones the package writes itself and renaming any that would collide with a statistic.
#'
#' @param rowTable Data.frame with the \code{rowData} of the counts.
#' @param extraColumns \code{TRUE} for every annotation column, \code{FALSE} for none, or a character vector naming the ones wanted.
#' @param reserved Character vector with the names already spoken for.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A data.frame with one row per row of the counts, possibly with no column at all.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.extraRowColumns <-
  function(rowTable,
           extraColumns = TRUE,
           reserved = character(0),
           verbose = TRUE) {

    emptyTable <- data.frame(row.names = seq_len(nrow(rowTable)))

    if (isFALSE(extraColumns)) {
      return(emptyTable)
    }

    #-------------------------------#
    # Which columns                 #
    #-------------------------------#
    # These three are written by the counting and are already in the result under their own names
    candidateColumns <- setdiff(colnames(rowTable), c("region.set", "region.id", "tile.id"))

    if (is.character(extraColumns)) {
      absentColumns <- setdiff(extraColumns, colnames(rowTable))
      if (length(absentColumns) > 0) {
        stop("The following columns are absent from the region annotation: ",
             paste(absentColumns, collapse = ", "), ".", call. = FALSE)
      }
      candidateColumns <- intersect(extraColumns, candidateColumns)
    }

    if (length(candidateColumns) == 0) {
      return(emptyTable)
    }

    #-------------------------------#
    # What can be bound             #
    #-------------------------------#
    # A list column has no place in a flat table and would break the binding rather than the row
    isAtomic <- vapply(rowTable[candidateColumns], is.atomic, logical(1))

    if (any(!isAtomic) & isTRUE(verbose)) {
      message("The following region columns are not atomic and have been left out: ",
              paste(candidateColumns[!isAtomic], collapse = ", "), ".")
    }

    candidateColumns <- candidateColumns[isAtomic]

    if (length(candidateColumns) == 0) {
      return(emptyTable)
    }

    extraTable <- as.data.frame(rowTable[, candidateColumns, drop = FALSE], stringsAsFactors = FALSE)

    #-------------------------------#
    # Names already spoken for      #
    #-------------------------------#
    collidingColumns <- candidateColumns %in% reserved

    if (any(collidingColumns)) {
      colnames(extraTable)[collidingColumns] <- paste0(candidateColumns[collidingColumns], ".region")

      if (isTRUE(verbose)) {
        message("The following region columns share a name with another column of the table and carry the suffix '.region': ",
                paste(candidateColumns[collidingColumns], collapse = ", "), ".")
      }
    }

    return(extraTable)
  } # END function
