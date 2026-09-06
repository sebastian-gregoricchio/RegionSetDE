#' @title fitRegions
#'
#' @description Fits a linear model on the counts of a \code{RegionSetDE.counts} object, one model per row. The row is the region, or the tile when the counts were tiled, and the model is the same one that the set level tests will read later, so the two levels never disagree on the design, on the offsets or on the dispersion.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param design Formula evaluated on the \code{colData}, e.g. \code{~ replicate + condition}, the same formula written as a string, e.g. \code{"~ replicate + condition"}, or a design matrix with one row per sample. For the \code{"dream"} engine the formula must be given as a formula or a string and may contain random terms, e.g. \code{~ condition + (1|donor)}.
#' @param engine String with the model to fit, one of \code{"edgeR"} (quasi-likelihood negative binomial), \code{"voom"} (limma on log-CPM with precision weights), \code{"dream"} (limma with random effects), \code{"deseq2"} (negative binomial Wald test) and \code{"limma"} (limma-trend on the log2 signal, for values that are not counts). Default: \code{"edgeR"}.
#' @param samples Character vector with the names of the samples to keep, or conditions already applied with \code{\link{selectSamples}}. Default: \code{NULL}, all the samples.
#' @param block String with the name of a \code{colData} column holding a blocking variable whose effect is estimated as a correlation rather than as a coefficient. Only for \code{engine = "voom"}. Default: \code{NULL}.
#' @param random Formula with the random terms, e.g. \code{~ (1|donor)}, appended to \code{design}. Only for \code{engine = "dream"}. Default: \code{NULL}.
#' @param assay String with the name of the assay holding the values to model. Default: \code{"counts"}.
#' @param useOffsets Logical value to indicate whether the normalisation stored in the object must enter the model as offsets. Default: \code{TRUE}.
#' @param dispersion Dispersion to fit with instead of taking it from the residual variation, which is what a design with no replicates needs. A numeric value, the list returned by \code{\link{estimateNullDispersion}}, or one of the strings \code{"background"} and \code{"regionSet"}, which run that function here. Only for \code{engine = "edgeR"}. Default: \code{NULL}, from the residual variation when there is any and from \code{nullSource} when there is none.
#' @param assumeCountLike Logical value overriding what the object says about its own values, so that a count engine can be run on a signal declared as continuous. Default: \code{FALSE}.
#' @param nullSource String with where the null rows come from when a dispersion has to be estimated, either \code{"background"} or \code{"regionSet"}. Default: \code{"background"}.
#' @param nullRegionSets Character vector with the names of the sets used as null rows. Only for \code{nullSource = "regionSet"}. Default: \code{NULL}.
#' @param robust Logical value to indicate whether the dispersion, or the prior variance, must be estimated robustly against outlier regions. Default: \code{TRUE}.
#' @param universe What every region set will be compared against at the set level. Either the string \code{"matched"}, which builds a universe matched on \code{matchOn} through \code{\link{makeSetUniverse}}, the string \code{"all"}, which takes every other row as it is, \code{NULL} to build none, or a \code{RegionSetDE.universe} object. Default: \code{"matched"}.
#' @param matchOn Character vector with the covariates the comparison rows are matched on, among \code{"width"} and \code{"abundance"}. Default: \code{c("width", "abundance")}.
#' @param universeRatio Numeric value with the number of comparison rows drawn per region of the set. Default: \code{5}.
#' @param universeSets Character vector with the names of the sets the comparison rows are drawn from, which is what the competitive p-values will be relative to. Default: \code{NULL}, every set other than the one being tested.
#' @param BPPARAM \code{BiocParallelParam} object passed to the \code{"dream"} engine. Default: \code{NULL}, sequential.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.fit} object.
#'
#' @details Only one assay type belongs in a model. Marks, antibodies and assays differ in dispersion, in dynamic range and in what their scaling factors mean, so a fit that pools them borrows information between rows that describe unrelated experiments. Split the object with \code{\link{splitSamples}} and normalise each piece on its own before coming here.
#' Offsets are read from the object rather than recomputed. When \code{\link{normalizeCounts}} produced a matrix of log offsets, in the \code{offset} assay, that matrix is used; otherwise the per-sample \code{scaling.factor} is expanded into one. In both cases \code{edgeR::scaleOffset} puts the offsets back on the scale of the library sizes, which leaves the coefficients readable as log2 fold changes and keeps the fitted values on the scale of the raw counts. Running without offsets is possible and is almost never what you want: the library sizes alone assume that the depth is the only difference between the samples.
#' A design with one sample per condition has no residual degree of freedom, so nothing in the data says how much two libraries differ for reasons unrelated to the treatment. In that case the number is taken from rows assumed not to respond, through \code{\link{estimateNullDispersion}}, and the fit switches from the quasi-likelihood F test to a likelihood ratio test with the dispersion held fixed. It happens here rather than being asked for, since the alternative is a fit that cannot be produced at all, but it is announced when it does: the result is conditional on that number, and \code{\link{checkNullCalibration}} is what turns the assumption into something checkable. \code{edgeR}'s own guidance, when no null rows exist either, is to pick a BCV by experience, around 0.4 for human samples, 0.1 for genetically identical model organisms and 0.01 for technical replicates, and to read the output as descriptive.
#' The engines answer slightly different questions. \code{"edgeR"} is the default and the safest with few replicates, since the quasi-likelihood F test carries the uncertainty of the dispersion estimate into the p-value. \code{"voom"} is faster on large objects and more flexible on the design, and it is the one that can absorb a repeated-measures structure through \code{block} without spending a coefficient on it. \code{"dream"} extends the same machinery to explicit random effects, which is what a design with several samples per donor asks for. \code{"deseq2"} is included for comparison and for the shrunken fold changes; on a few thousand regions it agrees with edgeR almost everywhere.
#' Four of the five are count models and the fifth is not, which is the distinction that decides the choice on bigWig input. \code{edgeR}, \code{DESeq2}, \code{voom} and \code{dream} all describe the number of fragments falling in an interval: the first two through a negative binomial likelihood, the other two through a mean-variance trend estimated on the count scale. Coverage read out of a bigWig is not that number, and rounding it to an integer does not make it one. \code{"limma"} models the log2 signal with an abundance trend on the residual variance and asks nothing of the values beyond their being continuous, which is what an already normalised track can supply. An object counted from BAM files may use any of the five; one built by \code{\link{countBigwig}} is refused by the count engines unless it was declared \code{countLike}, or unless \code{assumeCountLike} overrides it here.
#' The comparison universe of the set level tests is built here rather than there, because it depends on the rows and not on the contrast: one fit, one universe, however many contrasts are run on it afterwards. An object holding a single region set has nothing to compare that set against, and in that case the universe comes out empty with a message; the per-region analysis is unaffected.
#' Low-count rows are not removed here. Filter them with \code{\link{filterRegions}} first, or the dispersion trend is fitted on rows that carry no information.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#' counts <- normalizeCounts(counts, method = "background", verbose = FALSE)
#' counts <- filterRegions(counts, verbose = FALSE)
#'
#' fit <- fitRegions(counts, design = ~ condition, engine = "edgeR", verbose = FALSE)
#' fit
#'
#' \donttest{
#' # limma-voom on the same design
#' voomFit <- fitRegions(counts, design = ~ condition, engine = "voom", verbose = FALSE)
#' }
#'
#' \dontrun{
#' # Random effects need the dream engine and the formula given as such
#' mixedFit <- fitRegions(counts, design = ~ condition + (1|donor), engine = "dream")
#' }
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testRegions}}, \code{\link{testRegionSets}}, \code{\link{makeSetUniverse}}, \code{\link{selectSamples}}
#'
#' @importFrom SummarizedExperiment assay assayNames colData rowData
#' @importFrom S4Vectors metadata
#' @importFrom stats model.matrix as.formula
#' @importFrom methods is
#'
#' @export fitRegions

fitRegions <-
  function(counts,
           design,
           engine = "edgeR",
           samples = NULL,
           block = NULL,
           random = NULL,
           assay = "counts",
           useOffsets = TRUE,
           dispersion = NULL,
           assumeCountLike = FALSE,
           nullSource = "background",
           nullRegionSets = NULL,
           robust = TRUE,
           universe = "matched",
           matchOn = c("width", "abundance"),
           universeRatio = 5,
           universeSets = NULL,
           BPPARAM = NULL,
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    engine <- tolower(as.character(engine[1]))
    if (!(engine %in% c("edger", "voom", "dream", "deseq2", "limma"))) {
      stop("The 'engine' parameter must be one of 'edgeR', 'voom', 'dream', 'deseq2', 'limma'.", call. = FALSE)
    }
    engine <- c("edger" = "edgeR", "voom" = "voom", "dream" = "dream", "deseq2" = "deseq2", "limma" = "limma")[[engine]]

    enginePackage <- c("edgeR" = "edgeR", "voom" = "limma", "dream" = "variancePartition",
                       "deseq2" = "DESeq2", "limma" = "limma")[[engine]]
    if (!requireNamespace(enginePackage, quietly = TRUE)) {
      stop("The '", enginePackage, "' package is needed for the '", engine, "' engine.", call. = FALSE)
    }

    if (!is.null(block) & !(engine %in% c("voom", "limma"))) {
      stop("The 'block' parameter applies to the 'voom' and 'limma' engines only, add the variable to the design for the other engines.", call. = FALSE)
    }

    if (!is.null(random) & engine != "dream") {
      stop("The 'random' parameter applies to the 'dream' engine only.", call. = FALSE)
    }

    if (!(assay %in% SummarizedExperiment::assayNames(counts))) {
      stop("The assay '", assay, "' is absent from the object.", call. = FALSE)
    }

    #-------------------------------#
    # Restrict the samples          #
    #-------------------------------#
    if (!is.null(samples)) {
      counts <- selectSamples(counts = counts, samples = samples, dropNormalization = FALSE, verbose = FALSE)
    }

    if (ncol(counts) < 2) {
      stop("At least two samples are needed to fit a model.", call. = FALSE)
    }

    colTable <- as.data.frame(SummarizedExperiment::colData(counts))
    rownames(colTable) <- colnames(counts)

    #-------------------------------#
    # Build the design              #
    #-------------------------------#
    designObject <- .buildDesign(design = design, random = random, colTable = colTable, engine = engine)
    designMatrix <- designObject$matrix
    designFormula <- designObject$formula

    #-------------------------------#
    # Dispersion from outside       #
    #-------------------------------#
    residualDegrees <- nrow(designMatrix) - ncol(designMatrix)

    dispersionObject <- .resolveDispersion(counts = counts,
                                           dispersion = dispersion,
                                           residualDegrees = residualDegrees,
                                           engine = engine,
                                           nullSource = nullSource,
                                           nullRegionSets = nullRegionSets,
                                           verbose = verbose)

    dispersion <- dispersionObject$dispersion

    if (residualDegrees < 1 & engine != "dream" & isTRUE(verbose)) {
      message("No residual degree of freedom: fitting with a fixed dispersion of ", signif(dispersion, 3),
              " (BCV ", signif(sqrt(dispersion), 3), ") and a likelihood ratio test.")
    }

    if (isTRUE(verbose)) {
      message("Fitting ", nrow(counts), " ", counts@counting.level, "s over ", ncol(counts),
              " samples with '", engine, "' (", paste(colnames(designMatrix), collapse = ", "), ").")
    }

    #-------------------------------#
    # Counts, sizes and offsets     #
    #-------------------------------#
    countMatrix <- as.matrix(SummarizedExperiment::assay(counts, assay))

    librarySizes <- SummarizedExperiment::colData(counts)$library.size
    if (is.null(librarySizes) | any(is.na(librarySizes))) {
      librarySizes <- colSums(countMatrix)
    }

    offsetMatrix <- .fitOffsets(counts = counts, useOffsets = useOffsets, verbose = verbose)

    #-------------------------------#
    # Counts, or something else     #
    #-------------------------------#
    countEngines <- c("edgeR", "voom", "dream", "deseq2")
    signalType <- S4Vectors::metadata(counts)$signal.type
    countLikeFlag <- S4Vectors::metadata(counts)$count.like

    # Objects built before the flag was recorded: alignments are counts, anything else has to say so
    if (is.null(countLikeFlag)) {
      countLikeFlag <- is.null(signalType) | identical(signalType, "reads")
    }

    countLike <- isTRUE(countLikeFlag) | isTRUE(assumeCountLike)

    # Rounding coverage produces integers, which is not the same thing as producing counts, and every
    # count engine here reads its input as a number of fragments whether or not that is what it holds
    if (engine %in% countEngines & !is.null(signalType) & !countLike) {
      stop("The values in this object were not declared as counts (signal.type '", signalType, "'), and the '",
           engine, "' engine models a fragment count. Fit them with engine = 'limma', which models the log2 ",
           "signal directly, or set 'assumeCountLike = TRUE' to state that the values are raw coverage and ",
           "that the count model is an approximation you are willing to make.", call. = FALSE)
    }

    # A negative binomial likelihood is defined on integers, a normalised assay silently breaks it
    if (engine %in% countEngines & any(abs(countMatrix - round(countMatrix)) > 1e-8)) {
      if (engine == "deseq2") {
        stop("The 'deseq2' engine needs integer counts, use the raw 'counts' assay and let the offsets carry the normalisation.", call. = FALSE)
      }
      warning("The assay does not contain integers, the count model is being applied to non-integer values.", call. = FALSE)
    }

    #-------------------------------#
    # Fit                           #
    #-------------------------------#
    fitList <- switch(engine,
                      "edgeR" = .fitEdgeR(countMatrix = countMatrix, designMatrix = designMatrix, librarySizes = librarySizes,
                                          offsetMatrix = offsetMatrix, colTable = colTable, robust = robust,
                                          dispersion = dispersion),
                      "voom" = .fitVoom(countMatrix = countMatrix, designMatrix = designMatrix, librarySizes = librarySizes,
                                        offsetMatrix = offsetMatrix, colTable = colTable, block = block, robust = robust, verbose = verbose),
                      "dream" = .fitDream(countMatrix = countMatrix, designMatrix = designMatrix, designFormula = designFormula,
                                          librarySizes = librarySizes, offsetMatrix = offsetMatrix, colTable = colTable, BPPARAM = BPPARAM),
                      "deseq2" = .fitDESeq2(countMatrix = countMatrix, designMatrix = designMatrix, offsetMatrix = offsetMatrix,
                                            librarySizes = librarySizes, colTable = colTable, verbose = verbose),
                      "limma" = .fitLimma(countMatrix = countMatrix, designMatrix = designMatrix, librarySizes = librarySizes,
                                          offsetMatrix = offsetMatrix, colTable = colTable, block = block,
                                          robust = robust, verbose = verbose))

    fitList$dispersion$source <- dispersionObject$source
    fitList$dispersion$holdout.index <- dispersionObject$holdout.index
    fitList$dispersion$holdout.type <- dispersionObject$holdout.type

    #-------------------------------#
    # Universe of the set level     #
    #-------------------------------#
    # It reads the rows and the abundances, neither of which the contrast touches, so once is enough
    universeObject <- .resolveUniverse(object = counts,
                                       universe = universe,
                                       matchOn = matchOn,
                                       universeRatio = universeRatio,
                                       universeSets = universeSets,
                                       soft = TRUE,
                                       verbose = verbose)

    #-------------------------------#
    # Assemble the object           #
    #-------------------------------#
    fitObject <- new(Class = "RegionSetDE.fit",
                     fit = fitList$fit,
                     engine = engine,
                     design = designMatrix,
                     design.formula = if (is.null(designFormula)) {list()} else {list(formula = designFormula)},
                     blocking = fitList$blocking,
                     dispersion = fitList$dispersion,
                     counts = counts,
                     universe = universeObject,
                     samples = colnames(counts),
                     counting.level = counts@counting.level,
                     blacklist = counts@blacklist,
                     whitelist = counts@whitelist,
                     genome.assembly = counts@genome.assembly,
                     seqlevels.style = counts@seqlevels.style,
                     filtering.log = counts@filtering.log,
                     parameters = c(counts@parameters,
                                    list(fitRegions = list(engine = engine,
                                                           design = designMatrix,
                                                           block = block,
                                                           random = random,
                                                           assay = assay,
                                                           useOffsets = useOffsets,
                                                           assumeCountLike = assumeCountLike,
                                                           dispersion = dispersion,
                                                           dispersion.source = dispersionObject$source,
                                                           robust = robust,
                                                           universe = universe,
                                                           matchOn = matchOn,
                                                           universeRatio = universeRatio))))

    if (isTRUE(verbose)) {
      if (!is.null(fitList$dispersion$common)) {
        message("Done. Common dispersion: ", signif(fitList$dispersion$common, 3),
                " (BCV ", signif(sqrt(fitList$dispersion$common), 3), ").")
      } else {
        message("Done.")
      }
    }

    return(fitObject)
  } # END function




#' @title .buildDesign
#'
#' @description Turns the \code{design} argument of \code{\link{fitRegions}} into a model matrix, keeping the formula aside when the engine needs it.
#'
#' @param design Formula, string holding a formula, or design matrix.
#' @param random Formula, or string holding a formula, with the random terms, or \code{NULL}.
#' @param colTable Data.frame with the sample metadata.
#' @param engine String with the engine being used.
#'
#' @return A list with the \code{matrix} of the fixed effects and the \code{formula}, \code{NULL} when a matrix was supplied.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats model.matrix as.formula terms
#'
#' @keywords internal

.buildDesign <-
  function(design,
           random = NULL,
           colTable,
           engine) {

    #-------------------------------#
    # A matrix is taken as it is    #
    #-------------------------------#
    if (is.matrix(design)) {
      if (engine == "dream") {
        stop("The 'dream' engine needs the design as a formula, so that the random terms can be read.", call. = FALSE)
      }
      if (nrow(design) != nrow(colTable)) {
        stop("The design matrix must have one row per sample.", call. = FALSE)
      }
      if (is.null(colnames(design))) {
        stop("The columns of the design matrix must be named.", call. = FALSE)
      }
      return(list(matrix = design, formula = NULL))
    }

    design <- .asFormula(x = design, parameterName = "design")

    if (!inherits(design, "formula")) {
      stop("The 'design' parameter must be a formula, a string holding a formula, or a design matrix.", call. = FALSE)
    }

    #-------------------------------#
    # Merge the random terms        #
    #-------------------------------#
    designFormula <- design
    if (!is.null(random)) {
      random <- .asFormula(x = random, parameterName = "random")
      if (!inherits(random, "formula")) {
        stop("The 'random' parameter must be a formula, e.g. ~ (1|donor).", call. = FALSE)
      }
      designFormula <- stats::as.formula(paste("~", paste(c(attr(stats::terms(design), "term.labels"),
                                                            attr(stats::terms(random), "term.labels")), collapse = " + ")))
    }

    hasRandomTerm <- any(grepl("\\|", attr(stats::terms(designFormula), "term.labels")))

    if (hasRandomTerm & engine != "dream") {
      stop("Random terms are supported by the 'dream' engine only.", call. = FALSE)
    }

    if (engine == "dream" & !hasRandomTerm) {
      warning("The formula carries no random term, the 'voom' engine gives the same answer faster.", call. = FALSE)
    }

    #-------------------------------#
    # Fixed effect model matrix     #
    #-------------------------------#
    fixedFormula <- designFormula
    if (hasRandomTerm) {
      if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("The 'lme4' package is needed to separate the random terms from the design.", call. = FALSE)
      }
      fixedFormula <- lme4::nobars(designFormula)
    }

    designVariables <- all.vars(fixedFormula)
    absentVariables <- setdiff(designVariables, colnames(colTable))
    if (length(absentVariables) > 0) {
      stop("The following design variables are absent from the colData: ", paste(absentVariables, collapse = ", "), ".", call. = FALSE)
    }

    designMatrix <- stats::model.matrix(object = fixedFormula, data = colTable)

    # A rank deficient design gives coefficients that cannot be told apart, better to catch it here than inside the engine
    if (qr(designMatrix)$rank < ncol(designMatrix)) {
      stop("The design matrix is rank deficient, some coefficients are confounded.", call. = FALSE)
    }

    return(list(matrix = designMatrix, formula = designFormula))
  } # END function




#' @title .fitOffsets
#'
#' @description Extracts the normalisation stored in a counts object as a matrix of log offsets.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param useOffsets Logical value indicating whether the normalisation must be used.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A numeric matrix of log offsets, or \code{NULL}.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom SummarizedExperiment assay assayNames colData
#'
#' @keywords internal

.fitOffsets <-
  function(counts,
           useOffsets = TRUE,
           verbose = TRUE) {

    if (isFALSE(useOffsets)) {
      return(NULL)
    }

    # The loess normalisation writes one value per row and per sample, nothing else can reproduce it
    if ("offset" %in% SummarizedExperiment::assayNames(counts)) {
      return(as.matrix(SummarizedExperiment::assay(counts, "offset")))
    }

    scalingFactors <- SummarizedExperiment::colData(counts)$scaling.factor

    if (is.null(scalingFactors) | all(is.na(scalingFactors))) {
      if (isTRUE(verbose)) {
        warning("No normalisation is stored in the object, the model uses the library sizes alone.", call. = FALSE)
      }
      return(NULL)
    }

    if (any(is.na(scalingFactors)) | any(scalingFactors <= 0)) {
      stop("The scaling factors must be finite and strictly positive.", call. = FALSE)
    }

    # One factor per sample, expanded so that every engine receives the same object shape
    offsetMatrix <- matrix(data = rep(log(scalingFactors), each = nrow(counts)),
                           nrow = nrow(counts), ncol = ncol(counts),
                           dimnames = list(rownames(counts), colnames(counts)))

    return(offsetMatrix)
  } # END function




#' @title .fitEdgeR
#'
#' @description Fits the quasi-likelihood negative binomial model of \code{edgeR}.
#'
#' @param countMatrix Numeric matrix of counts.
#' @param designMatrix Design matrix.
#' @param librarySizes Numeric vector with the library sizes.
#' @param offsetMatrix Matrix of log offsets, or \code{NULL}.
#' @param colTable Data.frame with the sample metadata.
#' @param robust Logical value passed to the dispersion estimation.
#' @param dispersion Numeric value with a dispersion held fixed, or \code{NULL} to estimate one.
#'
#' @return A list with the \code{fit}, the \code{blocking} and the \code{dispersion} elements.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom edgeR DGEList scaleOffset estimateDisp glmQLFit glmFit
#'
#' @keywords internal

.fitEdgeR <-
  function(countMatrix,
           designMatrix,
           librarySizes,
           offsetMatrix,
           colTable,
           robust = TRUE,
           dispersion = NULL) {

    dgeList <- edgeR::DGEList(counts = countMatrix, samples = colTable, lib.size = librarySizes)

    # scaleOffset puts the offsets back on the scale of the library sizes, so the intercept stays interpretable
    if (!is.null(offsetMatrix)) {
      dgeList <- edgeR::scaleOffset(y = dgeList, offset = offsetMatrix)
    }

    #-------------------------------#
    # Dispersion held fixed         #
    #-------------------------------#
    # glmQLFit carries the uncertainty of an estimated dispersion, which there is none of when the value is given
    if (!is.null(dispersion)) {
      dgeList$common.dispersion <- dispersion
      modelFit <- edgeR::glmFit(y = dgeList, design = designMatrix, dispersion = dispersion)

      return(list(fit = list(object = modelFit, dge = dgeList, test = "lrt"),
                  blocking = list(),
                  dispersion = list(common = dispersion,
                                    fixed = TRUE,
                                    no.replicates = (nrow(designMatrix) - ncol(designMatrix)) < 1)))
    }

    dgeList <- edgeR::estimateDisp(y = dgeList, design = designMatrix, robust = robust)
    quasiFit <- edgeR::glmQLFit(y = dgeList, design = designMatrix, robust = robust)

    return(list(fit = list(object = quasiFit, dge = dgeList, test = "ql"),
                blocking = list(),
                dispersion = list(common = dgeList$common.dispersion,
                                  trended = dgeList$trended.dispersion,
                                  fixed = FALSE,
                                  no.replicates = FALSE,
                                  prior.df = quasiFit$df.prior)))
  } # END function




#' @title .fitVoom
#'
#' @description Fits the \code{limma} model on log-CPM values with the precision weights of \code{voom}.
#'
#' @param countMatrix Numeric matrix of counts.
#' @param designMatrix Design matrix.
#' @param librarySizes Numeric vector with the library sizes.
#' @param offsetMatrix Matrix of log offsets, or \code{NULL}.
#' @param colTable Data.frame with the sample metadata.
#' @param block String with the name of the blocking column, or \code{NULL}.
#' @param robust Logical value passed to the empirical Bayes step.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A list with the \code{fit}, the \code{blocking} and the \code{dispersion} elements.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom edgeR DGEList scaleOffset
#' @importFrom limma voom lmFit duplicateCorrelation
#'
#' @keywords internal

.fitVoom <-
  function(countMatrix,
           designMatrix,
           librarySizes,
           offsetMatrix,
           colTable,
           block = NULL,
           robust = TRUE,
           verbose = TRUE) {

    dgeList <- edgeR::DGEList(counts = countMatrix, samples = colTable, lib.size = librarySizes)
    if (!is.null(offsetMatrix)) {
      dgeList <- edgeR::scaleOffset(y = dgeList, offset = offsetMatrix)
    }

    voomObject <- limma::voom(counts = dgeList, design = designMatrix, plot = FALSE)

    #-------------------------------#
    # Blocked, or not               #
    #-------------------------------#
    if (is.null(block)) {
      linearFit <- limma::lmFit(object = voomObject, design = designMatrix)
      blockingInfo <- list()

    } else {
      if (!(block %in% colnames(colTable))) {
        stop("The column '", block, "' is absent from the colData.", call. = FALSE)
      }
      blockVector <- colTable[[block]]

      # The first correlation is estimated on unweighted values, the second pass feeds it back into voom
      firstCorrelation <- limma::duplicateCorrelation(object = voomObject, design = designMatrix, block = blockVector)$consensus.correlation
      voomObject <- limma::voom(counts = dgeList, design = designMatrix, plot = FALSE,
                                block = blockVector, correlation = firstCorrelation)
      consensusCorrelation <- limma::duplicateCorrelation(object = voomObject, design = designMatrix, block = blockVector)$consensus.correlation

      if (isTRUE(verbose)) {
        message("Consensus correlation within '", block, "': ", signif(consensusCorrelation, 3), ".")
      }

      linearFit <- limma::lmFit(object = voomObject, design = designMatrix,
                                block = blockVector, correlation = consensusCorrelation)
      blockingInfo <- list(block = block, correlation = consensusCorrelation)
    }

    return(list(fit = list(object = linearFit, voom = voomObject, robust = robust),
                blocking = blockingInfo,
                dispersion = list(sigma.median = stats::median(linearFit$sigma, na.rm = TRUE))))
  } # END function




#' @title .logSignalMatrix
#'
#' @description Puts a matrix of values on a log2 scale comparable across samples, dividing out the offsets when there are any and the library sizes otherwise.
#'
#' @param countMatrix Numeric matrix with the values.
#' @param librarySizes Numeric vector with the library size of each sample.
#' @param offsetMatrix Numeric matrix of log offsets, or \code{NULL}.
#' @param priorCount Numeric value added before the logarithm, damping the variance of the low rows. Default: \code{2}.
#'
#' @return A numeric matrix of log2 values on a per-million scale.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.logSignalMatrix <-
  function(countMatrix,
           librarySizes,
           offsetMatrix = NULL,
           priorCount = 2) {

    if (is.null(offsetMatrix)) {
      offsetMatrix <- matrix(data = rep(log(librarySizes), each = nrow(countMatrix)),
                             nrow = nrow(countMatrix), ncol = ncol(countMatrix))
    } else {
      # Centring on the row puts the offsets back on the scale of the library sizes, so the values stay readable
      offsetMatrix <- offsetMatrix - rowMeans(offsetMatrix) + mean(log(librarySizes))
    }

    logMatrix <- log2(countMatrix + priorCount) - offsetMatrix / log(2) + log2(1e6)
    dimnames(logMatrix) <- dimnames(countMatrix)

    return(logMatrix)
  } # END function




#' @title .fitLimma
#'
#' @description Fits the \code{limma} model on the log2 signal, with the residual variance trended on abundance, for values that are not counts.
#'
#' @param countMatrix Numeric matrix with the values.
#' @param designMatrix Design matrix.
#' @param librarySizes Numeric vector with the library size of each sample.
#' @param offsetMatrix Numeric matrix of log offsets, or \code{NULL}.
#' @param colTable Data.frame with the sample annotation.
#' @param block String naming a \code{colData} column holding a blocking variable, or \code{NULL}. Default: \code{NULL}.
#' @param robust Logical value to indicate whether the prior variance must be estimated robustly. Default: \code{TRUE}.
#' @param priorCount Numeric value added before the logarithm. Default: \code{2}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A list with the fit, the blocking information and the median residual standard deviation.
#'
#' @details Where \code{voom} estimates a weight per observation from a mean-variance trend fitted on the count scale, this fits the trend on the residual variance itself and applies no weights. Nothing about the values then has to be a fragment count, which is what makes it the engine for coverage, for an externally normalised matrix, or for anything else that arrives already on a continuous scale.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom limma lmFit duplicateCorrelation
#' @importFrom stats median
#'
#' @keywords internal

.fitLimma <-
  function(countMatrix,
           designMatrix,
           librarySizes,
           offsetMatrix,
           colTable,
           block = NULL,
           robust = TRUE,
           priorCount = 2,
           verbose = TRUE) {

    logMatrix <- .logSignalMatrix(countMatrix = countMatrix,
                                  librarySizes = librarySizes,
                                  offsetMatrix = offsetMatrix,
                                  priorCount = priorCount)

    #-------------------------------#
    # Blocked, or not               #
    #-------------------------------#
    if (is.null(block)) {
      linearFit <- limma::lmFit(object = logMatrix, design = designMatrix)
      blockingInfo <- list()

    } else {
      if (!(block %in% colnames(colTable))) {
        stop("The column '", block, "' is absent from the colData.", call. = FALSE)
      }
      blockVector <- colTable[[block]]

      consensusCorrelation <- limma::duplicateCorrelation(object = logMatrix, design = designMatrix, block = blockVector)$consensus.correlation

      if (isTRUE(verbose)) {
        message("Consensus correlation within '", block, "': ", signif(consensusCorrelation, 3), ".")
      }

      linearFit <- limma::lmFit(object = logMatrix, design = designMatrix,
                                block = blockVector, correlation = consensusCorrelation)
      blockingInfo <- list(block = block, correlation = consensusCorrelation)
    }

    # The trend on the residual variance replaces the observation weights voom would have supplied
    return(list(fit = list(object = linearFit, voom = list(E = logMatrix), robust = robust, trend = TRUE),
                blocking = blockingInfo,
                dispersion = list(sigma.median = stats::median(linearFit$sigma, na.rm = TRUE))))
  } # END function




#' @title .fitDream
#'
#' @description Fits the mixed model of \code{variancePartition}, with precision weights estimated under the same formula.
#'
#' @param countMatrix Numeric matrix of counts.
#' @param designMatrix Design matrix of the fixed effects, kept for the contrast resolution.
#' @param designFormula Formula including the random terms.
#' @param librarySizes Numeric vector with the library sizes.
#' @param offsetMatrix Matrix of log offsets, or \code{NULL}.
#' @param colTable Data.frame with the sample metadata.
#' @param BPPARAM \code{BiocParallelParam} object, or \code{NULL}.
#'
#' @return A list with the \code{fit}, the \code{blocking} and the \code{dispersion} elements.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom edgeR DGEList scaleOffset
#'
#' @keywords internal

.fitDream <-
  function(countMatrix,
           designMatrix,
           designFormula,
           librarySizes,
           offsetMatrix,
           colTable,
           BPPARAM = NULL) {

    dgeList <- edgeR::DGEList(counts = countMatrix, samples = colTable, lib.size = librarySizes)
    if (!is.null(offsetMatrix)) {
      dgeList <- edgeR::scaleOffset(y = dgeList, offset = offsetMatrix)
    }

    if (is.null(BPPARAM)) {
      BPPARAM <- BiocParallel::SerialParam()
    }

    # The weights have to come from the same formula as the fit, a voom run on the fixed part alone underestimates them
    voomObject <- variancePartition::voomWithDreamWeights(counts = dgeList,
                                                          formula = designFormula,
                                                          data = colTable,
                                                          BPPARAM = BPPARAM,
                                                          quiet = TRUE)

    mixedFit <- variancePartition::dream(exprObj = voomObject,
                                         formula = designFormula,
                                         data = colTable,
                                         BPPARAM = BPPARAM,
                                         quiet = TRUE)
    mixedFit <- variancePartition::eBayes(mixedFit)

    return(list(fit = list(object = mixedFit, voom = voomObject, data = colTable, formula = designFormula, BPPARAM = BPPARAM),
                blocking = list(random = designFormula),
                dispersion = list(sigma.median = stats::median(mixedFit$sigma, na.rm = TRUE))))
  } # END function




#' @title .fitDESeq2
#'
#' @description Fits the negative binomial model of \code{DESeq2} on a user supplied model matrix.
#'
#' @param countMatrix Integer matrix of counts.
#' @param designMatrix Design matrix.
#' @param offsetMatrix Matrix of log offsets, or \code{NULL}.
#' @param librarySizes Numeric vector with the library sizes.
#' @param colTable Data.frame with the sample metadata.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A list with the \code{fit}, the \code{blocking} and the \code{dispersion} elements.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.fitDESeq2 <-
  function(countMatrix,
           designMatrix,
           offsetMatrix,
           librarySizes,
           colTable,
           verbose = TRUE) {

    storage.mode(countMatrix) <- "integer"

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = countMatrix,
                                          colData = colTable,
                                          design = ~ 1)

    #-------------------------------#
    # Offsets, in the DESeq2 idiom  #
    #-------------------------------#
    if (is.null(offsetMatrix)) {
      DESeq2::sizeFactors(dds) <- librarySizes / exp(mean(log(librarySizes)))
    } else {
      # DESeq2 asks the factors to have geometric mean one on each row, which is where the row centring comes from
      DESeq2::normalizationFactors(dds) <- exp(offsetMatrix - rowMeans(offsetMatrix))
    }

    # The model matrix is passed as 'full', which keeps the coefficient names identical to the other engines
    dds <- DESeq2::DESeq(object = dds, full = designMatrix, betaPrior = FALSE, quiet = !verbose)

    return(list(fit = list(object = dds),
                blocking = list(),
                dispersion = list(common = mean(DESeq2::dispersions(dds), na.rm = TRUE))))
  } # END function




#' @title .asFormula
#'
#' @description Parses a design written as a string, leaving anything else untouched.
#'
#' @param x Object passed by the user.
#' @param parameterName String with the name of the parameter, used in the error message.
#'
#' @return A formula when the input was a string, the input itself otherwise.
#'
#' @author Sebastian Gregoricchio
#'
#' @importFrom stats as.formula
#'
#' @keywords internal

.asFormula <-
  function(x,
           parameterName = "design") {

    if (!is.character(x) | length(x) != 1) {
      return(x)
    }

    # A formula stored in a table or read from a configuration file arrives as a string, parsing it here costs nothing
    parsedFormula <- try(stats::as.formula(x), silent = TRUE)

    if (inherits(parsedFormula, "try-error")) {
      stop("The '", parameterName, "' parameter is a string that cannot be read as a formula.", call. = FALSE)
    }

    return(parsedFormula)
  } # END function




#' @title .resolveDispersion
#'
#' @description Works out which dispersion \code{\link{fitRegions}} should hold fixed, estimating one from rows assumed not to respond when the design leaves no residual to take it from.
#'
#' @param counts \code{RegionSetDE.counts} object.
#' @param dispersion Numeric value, list from \code{\link{estimateNullDispersion}}, keyword, or \code{NULL}.
#' @param residualDegrees Numeric value with the residual degrees of freedom of the design.
#' @param engine String with the engine being used.
#' @param nullSource String with where the null rows come from.
#' @param nullRegionSets Character vector with the sets used as null rows, or \code{NULL}.
#' @param verbose Logical value to indicate whether the messages must be printed.
#'
#' @return A list with the \code{dispersion}, \code{NULL} when it must come from the residual variation, the \code{source} it came from, and the rows held out of the estimate.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.resolveDispersion <-
  function(counts,
           dispersion = NULL,
           residualDegrees,
           engine = "edgeR",
           nullSource = "background",
           nullRegionSets = NULL,
           verbose = TRUE) {

    #-------------------------------#
    # Whatever was handed in        #
    #-------------------------------#
    if (is.list(dispersion)) {
      if (is.null(dispersion$dispersion)) {
        stop("The \'dispersion\' list must carry a \'dispersion\' element, as returned by estimateNullDispersion.", call. = FALSE)
      }
      dispersionList <- .checkDispersion(dispersion = dispersion$dispersion,
                                         source = if (is.null(dispersion$source)) {"supplied"} else {dispersion$source},
                                         engine = engine)
      dispersionList$holdout.index <- dispersion$holdout.index
      dispersionList$holdout.type <- if (is.null(dispersion$holdout.type)) {"none"} else {dispersion$holdout.type}
      return(dispersionList)
    }

    if (is.character(dispersion)) {
      if (!(dispersion %in% c("background", "regionSet"))) {
        stop("The \'dispersion\' parameter must be a value, a list from estimateNullDispersion, \'background\' or \'regionSet\'.", call. = FALSE)
      }
      nullSource <- dispersion
      dispersion <- NULL

    } else if (!is.null(dispersion)) {
      return(.checkDispersion(dispersion = dispersion, source = "supplied", engine = engine))
    }

    #-------------------------------#
    # Nothing to estimate it from   #
    #-------------------------------#
    if (residualDegrees >= 1 | engine %in% c("dream", "limma")) {
      return(list(dispersion = NULL, source = "residual", holdout.index = integer(0), holdout.type = "none"))
    }

    if (engine != "edgeR") {
      stop("The design uses one sample per coefficient, which leaves no residual to estimate the dispersion from. ",
           "Only the \'edgeR\' engine can be fitted with a dispersion brought in from outside.", call. = FALSE)
    }

    # A fit that cannot be produced at all is worse than one that announces where its dispersion came from
    if (isTRUE(verbose)) {
      message("No residual degree of freedom, estimating the dispersion from rows assumed not to respond.")
    }

    nullDispersion <- try(estimateNullDispersion(counts = counts,
                                                 source = nullSource,
                                                 regionSets = nullRegionSets,
                                                 verbose = verbose),
                          silent = TRUE)

    if (inherits(nullDispersion, "try-error")) {
      stop("The dispersion could not be estimated: ", sub("^Error[^:]*: ", "", nullDispersion[1]),
           "\n  Supply one through \'dispersion\', as a BCV squared, and check it with `checkNullCalibration()`.", call. = FALSE)
    }

    return(list(dispersion = nullDispersion$dispersion,
                source = nullSource,
                holdout.index = nullDispersion$holdout.index,
                holdout.type = nullDispersion$holdout.type))
  } # END function




#' @title .checkDispersion
#'
#' @description Validates a dispersion given by the user.
#'
#' @param dispersion Numeric value.
#' @param source String with where it came from.
#' @param engine String with the engine being used.
#'
#' @return A list with the \code{dispersion}, its \code{source} and an empty holdout.
#'
#' @author Sebastian Gregoricchio
#'
#' @keywords internal

.checkDispersion <-
  function(dispersion,
           source = "supplied",
           engine = "edgeR") {

    dispersion <- as.numeric(dispersion[1])

    if (is.na(dispersion) | dispersion < 0) {
      stop("The \'dispersion\' parameter must be a single non-negative value.", call. = FALSE)
    }

    if (engine != "edgeR") {
      stop("A dispersion brought in from outside is only available for the \'edgeR\' engine.", call. = FALSE)
    }

    return(list(dispersion = dispersion, source = source, holdout.index = integer(0)))
  } # END function
