#' @title testRegions
#'
#' @description Tests a contrast on a \code{RegionSetDE.fit} object and returns one row per region. When the counts were tiled, every tile is tested on its own and the p-values are then combined back to the region, so that the region stays the unit of inference even though the signal was measured at a finer scale.
#'
#' @param fit \code{RegionSetDE.fit} object.
#' @param contrast Contrast to test, given in one of four ways. A character vector of length three, \code{c("column", "groupA", "groupB")}, naming a column of the \code{colData} and two of its levels, which is the form to reach for when the design uses a reference level. A string with the name of a design column, e.g. \code{"conditionCOMBO"}. A string written as an expression over the design columns, e.g. \code{"conditionCOMBO - conditionEPZ"}. Or a numeric vector with one coefficient per column of the design. A named list of any of these runs every contrast on the same fit and returns a \code{RegionSetDE.resultsList}.
#' @param combine Logical value to indicate whether the tile level p-values must be combined into one value per region. Ignored when the counts were not tiled. Default: \code{TRUE}.
#' @param combineMethod String with the method used to combine the tiles, among those accepted by \code{csaw::combineTests}: \code{"simes"}, \code{"holm-min"}, \code{"wilcoxon"} and \code{"stouffer"}. Default: \code{"simes"}.
#' @param lfcThreshold Numeric value with the log2 fold change against which the null hypothesis is tested. A value above zero moves the threshold inside the test, through \code{edgeR::glmTreat}, \code{limma::treat} or the \code{lfcThreshold} of \code{DESeq2::results}, which is stricter and better calibrated than filtering the output afterwards. Default: \code{0}.
#' @param FDR Numeric value with the adjusted p-value cut-off used to fill the \code{diff.status} column. Default: \code{0.05}.
#' @param log2FC Numeric value with the absolute log2 fold change cut-off used to fill the \code{diff.status} column. Default: \code{0}.
#' @param adjustMethod String with the multiple testing correction, passed to \code{stats::p.adjust}. Default: \code{"BH"}.
#' @param regionSets Character vector with the names of the region sets to keep in the output. Default: \code{NULL}, all of them.
#' @param extraColumns Annotation carried by the regions that must be appended to the result, at the end of the table. Either \code{TRUE} for every column of the \code{rowData} beyond the ones the package writes itself, \code{FALSE} for none, or a character vector naming the ones wanted. Default: \code{TRUE}.
#' @param carryCounts Logical value to indicate whether the counts must travel inside the result, so that \code{\link{plotRegion}} and \code{\link{plotTopHeatmap}} can draw the values without being handed the counts object again. Several contrasts run on one fit share the same copy in memory. Default: \code{TRUE}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.results} object.
#'
#' @details The multiple testing correction is applied over all the rows of the object, across the region sets, and \code{regionSets} subsets the output afterwards. Correcting inside each set separately would make the FDR of a set depend on how many other sets were loaded, which is not a property anyone wants in a result.
#'
#' Two things follow from the combination step. The p-value of a tiled region is a Simes combination, so it answers "does any part of this region change" rather than "does the whole region change", and a long domain that moves over one tile out of forty will come out with a small p-value and a small overall fold change. The \code{log2FC} reported for a combined region is the fold change of the most significant tile, not an average, which is the quantity that matches the p-value. The tile level table stays available in the \code{tiles} slot, and \code{\link{plotRegion}} draws it.
#'
#' A design written as \code{~ condition} spends one coefficient per level except the first, so a level can be a coefficient in the design or the reference the others are measured against, depending on how the factor was ordered. Naming a coefficient that turns out to be the reference is the usual source of confusion, and it is what \code{c("column", "groupA", "groupB")} avoids: that form averages the design rows of each group and takes the difference, which gives the same contrast whatever the reference is and whether the design was written as \code{~ condition} or \code{~ 0 + condition}. With other covariates in the design the averaging picks up their imbalance between the two groups, so it describes what it says only when the design is reasonably balanced.
#'
#' Whatever the regions were loaded with travels through to the result. A gene name, a peak score or any other column attached to the \code{rowData} comes out at the end of the table, which is what makes \code{topRegions()} readable and lets \code{plotVolcano(labelColumn = )} label the points with something other than an identifier. On a tiled object the value is read off the tile the combination reported, the same one the fold change comes from, so a row describes one place rather than an average over several.
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
                                    extraColumns = extraColumns, carryCounts = carryCounts, verbose = verbose))
               })

      names(resultsList) <- names(contrast)

      return(new(Class = "RegionSetDE.resultsList",
                 results = resultsList,
                 contrasts = names(contrast)))
    }

    if (!(combineMethod %in% c("simes", "holm-min", "wilcoxon", "stouffer"))) {
      stop("The 'combineMethod' parameter must be one of 'simes', 'holm-min', 'wilcoxon', 'stouffer'.", call. = FALSE)
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
                       "dream" = .testDream(fit = fit, contrastObject = contrastObject, lfcThreshold = lfcThreshold, verbose = verbose),
                       "deseq2" = .testDESeq2(fit = fit, contrastVector = contrastObject$vector, lfcThreshold = lfcThreshold))

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
                                    extraColumns = colnames(extraTable),
                                    method = combineMethod,
                                    lfcThreshold = lfcThreshold,
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
                          "log2FC", "average.signal", "stat", "p.value", "FDR", "diff.status")
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
#' @description Works out which variable of the sample metadata a contrast separates, and which two of its levels, by comparing the contrast against the difference between the design rows of every pair of levels.
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

    for (columnName in colnames(colTable)) {
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
#' @return A data.frame with the \code{log2FC}, \code{average.signal}, \code{stat} and \code{p.value} columns.
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

    return(data.frame(log2FC = topTable$logFC,
                      average.signal = topTable$logCPM,
                      stat = statColumn,
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
#' @return A data.frame with the \code{log2FC}, \code{average.signal}, \code{stat} and \code{p.value} columns.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom limma contrasts.fit eBayes treat topTable topTreat
#'
#' @keywords internal

.testVoom <-
  function(fit,
           contrastVector,
           lfcThreshold = 0) {

    contrastFit <- limma::contrasts.fit(fit = fit@fit$object, contrasts = contrastVector)

    if (lfcThreshold > 0) {
      contrastFit <- limma::treat(fit = contrastFit, lfc = lfcThreshold, robust = isTRUE(fit@fit$robust))
      topTable <- limma::topTreat(fit = contrastFit, number = Inf, sort.by = "none", adjust.method = "none")
    } else {
      contrastFit <- limma::eBayes(fit = contrastFit, robust = isTRUE(fit@fit$robust))
      topTable <- limma::topTable(fit = contrastFit, number = Inf, sort.by = "none", adjust.method = "none")
    }

    return(data.frame(log2FC = topTable$logFC,
                      average.signal = topTable$AveExpr,
                      stat = topTable$t,
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
#' @return A data.frame with the \code{log2FC}, \code{average.signal}, \code{stat} and \code{p.value} columns.
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

    return(data.frame(log2FC = topTable$logFC,
                      average.signal = topTable$AveExpr,
                      stat = topTable$t,
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
#' @return A data.frame with the \code{log2FC}, \code{average.signal}, \code{stat} and \code{p.value} columns.
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

    return(data.frame(log2FC = resultTable$log2FoldChange,
                      average.signal = log2(resultTable$baseMean + 1),
                      stat = resultTable$stat,
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
#' @param lfcThreshold Numeric value used to count the tiles moving in each direction.
#' @param adjustMethod String with the multiple testing correction.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A list with the \code{results} data.frame and the \code{regions} \code{GRanges}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom csaw combineTests
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
           lfcThreshold = 0,
           adjustMethod = "BH",
           verbose = TRUE) {

    # combineTests reads the grouping as a factor, its level order is what the output rows follow
    tileGroups <- factor(tileTable$region.key, levels = unique(tileTable$region.key))

    combinedTable <- csaw::combineTests(ids = tileGroups,
                                        tab = data.frame(logFC = tileTable$log2FC,
                                                         PValue = tileTable$p.value,
                                                         stringsAsFactors = FALSE),
                                        pval.col = "PValue",
                                        fc.col = "logFC",
                                        fc.threshold = lfcThreshold,
                                        method = method)

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

    return(c("log2FC", "average.signal", "stat", "p.value", "FDR", "diff.status",
             "n.tiles", "n.tiles.up", "n.tiles.down", "direction",
             "rep.tile.start", "rep.tile.end"))
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
        message("The following region columns share a name with a statistic and carry the suffix '.region': ",
                paste(candidateColumns[collidingColumns], collapse = ", "), ".")
      }
    }

    return(extraTable)
  } # END function
