#' @title scoreRegionSets
#'
#' @description Computes one signal score per region set per library, without a contrast, so that the sets can be compared to each other in an experiment holding a single condition. The score is the summarised signal over the regions of a set divided by a reference measured in the same library, and the sets are then compared library by library, which puts the replication in the libraries rather than in the regions. It answers whether one set carries more signal than another. It does not answer whether that difference comes from the factor or from what the sets are made of.
#'
#' @param counts \code{RegionSetDE.counts} object returned by \code{\link{countReads}}.
#' @param regionSets Character vector with the names of the sets to score. Default: \code{NULL}, all of them.
#' @param comparisons List of character vectors of length two, naming the sets compared to each other. Default: \code{NULL}, every pair.
#' @param reference String with what the signal of a set is divided by, one among \code{"background"} (the bins of \code{\link{countBackground}}), \code{"regions"} (every region loaded in the object) and \code{"none"} (no division). Default: \code{"background"}.
#' @param assay String with the name of the assay holding the signal, ignored when \code{reference} is \code{"background"}. Default: \code{NULL}, the normalised assay when present and the raw counts otherwise.
#' @param perBasepair Logical value indicating whether the signal must be divided by the width of the region before it is summarised. Default: \code{TRUE}.
#' @param summary String with how the regions of a set are summarised, either \code{"mean"} or \code{"median"}. Default: \code{"mean"}.
#' @param pseudoCount Numeric value added to the signal before the ratio is taken. Default: \code{0.5}.
#' @param minRegions Numeric value with the number of regions a set needs to be scored. Default: \code{10}.
#' @param confLevel Numeric value with the level of the confidence interval on the difference between two sets. Default: \code{0.95}.
#' @param adjustMethod String with the multiple testing correction applied across the comparisons. Default: \code{"BH"}.
#' @param verbose Logical value to indicate whether the messages must be printed. Default: \code{TRUE}.
#'
#' @return A \code{RegionSetDE.setScores} object. \code{\link{resultsTable}} returns the comparisons between the sets, \code{\link{scoreTable}} the per-library scores they were computed from, and \code{\link{plotSetSignal}} draws both.
#'
#' @details Every test in the package runs on a contrast, so a design holding one condition and three replicates leaves \code{\link{testSetContrast}} with nothing to compare. This function covers that case, and the price of covering it is worth stating before the output is read.
#'
#' The score of a set is a ratio taken inside one library, so the sequencing depth and any scaling factor shared by the numerator and the denominator cancel before the score exists. The background bins carry raw counts, so with \code{reference = "background"} the regions are read from the \code{counts} assay as well and the normalisation never enters, which also means that changing it cannot move the result. With \code{reference = "none"} that protection is gone: the score becomes the signal itself, the libraries sit on their own scales again, and the normalised assay is the one to pass.
#'
#' The reference cancels a second time when two sets are compared, since both scores of a library were divided by the same number. What the comparison tests is the log ratio between the summarised signal of the two sets, taken library by library. It depends on neither the reference nor the normalisation. It does depend on the libraries, and there are usually three of them, so the interval rests on two degrees of freedom and comes out wide. That width is what three replicates support.
#'
#' The regions are not the replication here and are deliberately not treated as it. Thirty thousand promoters collapse to one number per library, and the variance reaching the test is the variance between libraries. Summarising over the regions instead of testing across them is why the p-value is not the vanishing one a per-region test over the same data would return.
#'
#' What no arrangement of these numbers separates is the factor from the regions. Two sets differing in width, mappability, GC content or accessibility differ in coverage in a library where nothing is bound, and that difference repeats across replicates in the same way real binding does. \code{perBasepair} takes away the crudest part of it and leaves the rest standing. A difference between two sets is therefore a difference in signal, and calling it a difference in factor is a separate argument, made with the composition of the sets shown next to the result. An input, an IgG or a spike-in turns the score into an enrichment and removes the need to make that argument at all; where one exists, \code{\link{testSetContrast}} on the IP against it is the better instrument.
#'
#' @examples
#' counts <- loadExampleData("counts", verbose = FALSE)
#'
#' # One score per library per set, plus every pairwise comparison
#' setScores <- scoreRegionSets(counts, verbose = FALSE)
#' setScores
#'
#' resultsTable(setScores)
#' head(scoreTable(setScores))
#'
#' @author Sebastian Gregoricchio
#'
#' @seealso \code{\link{testSetContrast}}, \code{\link{countBackground}}, \code{\link{plotSetSignal}}
#'
#' @importFrom SummarizedExperiment assay assayNames colData rowData rowRanges
#' @importFrom S4Vectors metadata
#' @importFrom BiocGenerics width
#' @importFrom dplyr arrange count filter group_by left_join mutate n n_distinct summarise
#' @importFrom rlang .data
#' @importFrom methods is new validObject
#' @importFrom stats median p.adjust pt qt sd
#' @importFrom utils combn
#'
#' @export scoreRegionSets

scoreRegionSets <-
  function(counts,
           regionSets = NULL,
           comparisons = NULL,
           reference = "background",
           assay = NULL,
           perBasepair = TRUE,
           summary = "mean",
           pseudoCount = 0.5,
           minRegions = 10,
           confLevel = 0.95,
           adjustMethod = "BH",
           verbose = TRUE) {

    #------------------------#
    # Check of the arguments #
    #------------------------#
    if (!methods::is(counts, "RegionSetDE.counts")) {
      stop("The 'counts' parameter must be a RegionSetDE.counts object.", call. = FALSE)
    }

    if (!(summary %in% c("mean", "median"))) {
      stop("The 'summary' parameter must be either 'mean' or 'median'.", call. = FALSE)
    }

    if (!(reference %in% c("background", "regions", "none"))) {
      stop("The 'reference' parameter must be one among 'background', 'regions' or 'none'.", call. = FALSE)
    }

    if (ncol(counts) < 2) {
      stop("At least two libraries are required, the difference between two sets is measured across them.", call. = FALSE)
    }

    if (!is.numeric(confLevel) | confLevel[1] <= 0 | confLevel[1] >= 1) {
      stop("The 'confLevel' parameter must be a number strictly between 0 and 1.", call. = FALSE)
    }

    if (reference == "background") {
      backgroundCounts <- S4Vectors::metadata(counts)$background

      if (is.null(backgroundCounts)) {
        stop("No background bins are stored in the object, run countBackground() first or set 'reference' to 'regions'.", call. = FALSE)
      }

      # The bins hold raw counts, so a ratio against them has to be taken on raw counts too
      if (!is.null(assay) && assay != "counts") {
        stop("The background bins hold raw counts, the 'assay' parameter cannot be set together with reference = 'background'.", call. = FALSE)
      }

      assay <- "counts"

    } else {
      assay <- .defaultAssay(counts = counts, assay = assay)
    }

    #-----------------------------#
    # Table of the scored regions #
    #-----------------------------#
    regionTable <- data.frame(row.index = seq_len(nrow(counts)),
                              region.set = as.character(SummarizedExperiment::rowData(counts)$region.set),
                              width = as.numeric(BiocGenerics::width(SummarizedExperiment::rowRanges(counts))),
                              stringsAsFactors = FALSE)

    if (!is.null(regionSets)) {
      absentSets <- setdiff(regionSets, unique(regionTable$region.set))

      if (length(absentSets) > 0) {
        stop("The following region sets are absent from the object: ", paste(absentSets, collapse = ", "), ".", call. = FALSE)
      }

      regionTable <- dplyr::filter(regionTable, .data$region.set %in% regionSets)
    }

    # A summary over a handful of regions is one that a single locus can move on its own
    setSizes <- dplyr::count(regionTable, .data$region.set, name = "n.regions")
    smallSets <- dplyr::filter(setSizes, .data$n.regions < minRegions)

    if (nrow(smallSets) > 0) {
      warning("Dropped for holding fewer than ", minRegions, " regions: ", paste(smallSets$region.set, collapse = ", "), ".", call. = FALSE)
      regionTable <- dplyr::filter(regionTable, !(.data$region.set %in% smallSets$region.set))
    }

    if (dplyr::n_distinct(regionTable$region.set) < 2) {
      stop("At least two region sets must survive the filters, there is nothing to compare otherwise.", call. = FALSE)
    }

    #------------------------------#
    # Summarise the signal per set #
    #------------------------------#
    signalMatrix <- as.matrix(SummarizedExperiment::assay(counts, assay))[regionTable$row.index, , drop = FALSE] + pseudoCount

    # Wide regions collect more reads before any factor is involved, this removes the crudest part of that
    if (isTRUE(perBasepair)) {
      signalMatrix <- signalMatrix / regionTable$width
    }

    summaryFunction <- if (summary == "mean") {mean} else {stats::median}

    signalTable <- data.frame(sample = rep(colnames(counts), each = nrow(signalMatrix)),
                              region.set = rep(regionTable$region.set, times = ncol(signalMatrix)),
                              signal = as.numeric(signalMatrix),
                              stringsAsFactors = FALSE)

    setSignal <-
      dplyr::summarise(dplyr::group_by(signalTable, .data$sample, .data$region.set),
                       n.regions = dplyr::n(),
                       set.signal = summaryFunction(.data$signal),
                       .groups = "drop")

    #-------------------------------------------#
    # Reference measured in the same library    #
    #-------------------------------------------#
    referenceSignal <-
      switch(reference,
             "background" = {
               backgroundMatrix <- as.matrix(SummarizedExperiment::assay(backgroundCounts, "counts")) + pseudoCount

               if (isTRUE(perBasepair)) {
                 backgroundMatrix <- backgroundMatrix / as.numeric(BiocGenerics::width(SummarizedExperiment::rowRanges(backgroundCounts)))
               }

               data.frame(sample = colnames(backgroundCounts),
                          reference.signal = as.numeric(apply(backgroundMatrix, 2, summaryFunction)),
                          stringsAsFactors = FALSE)
             },
             "regions" = {
               as.data.frame(dplyr::summarise(dplyr::group_by(signalTable, .data$sample),
                                              reference.signal = summaryFunction(.data$signal),
                                              .groups = "drop"))
             },
             "none" = {
               data.frame(sample = colnames(counts),
                          reference.signal = 1,
                          stringsAsFactors = FALSE)
             })

    if (!all(colnames(counts) %in% referenceSignal$sample)) {
      stop("The reference is missing for some libraries, the regions and the bins were not counted from the same samples.", call. = FALSE)
    }

    setScores <-
      as.data.frame(dplyr::arrange(dplyr::mutate(dplyr::left_join(setSignal, referenceSignal, by = "sample"),
                                                 score = log2(.data$set.signal / .data$reference.signal)),
                                   .data$region.set, .data$sample))

    #------------------------------------------#
    # Compare the sets, one library at a time  #
    #------------------------------------------#
    scoredSets <- sort(unique(setScores$region.set))

    if (is.null(comparisons)) {
      comparisons <- utils::combn(scoredSets, m = 2, simplify = FALSE)
    }

    if (!is.list(comparisons) | any(lengths(comparisons) != 2)) {
      stop("The 'comparisons' parameter must be a list of character vectors of length two.", call. = FALSE)
    }

    unknownSets <- setdiff(unlist(comparisons), scoredSets)

    if (length(unknownSets) > 0) {
      stop("The following sets are named in 'comparisons' but were not scored: ", paste(unknownSets, collapse = ", "), ".", call. = FALSE)
    }

    if (isTRUE(verbose)) {
      message("Comparing ", length(comparisons), " pair(s) of region sets across ", ncol(counts), " libraries...")
    }

    comparisonTable <-
      do.call(rbind,
              lapply(comparisons,
                     function(pair) {
                       firstScores <- dplyr::filter(setScores, .data$region.set == pair[1])
                       secondScores <- dplyr::filter(setScores, .data$region.set == pair[2])

                       # Both scores of a library were divided by the same reference, so it drops out of the difference
                       sharedSamples <- intersect(firstScores$sample, secondScores$sample)
                       pairedDifference <- firstScores$score[match(sharedSamples, firstScores$sample)] -
                         secondScores$score[match(sharedSamples, secondScores$sample)]

                       nLibraries <- length(pairedDifference)
                       meanDifference <- mean(pairedDifference)
                       standardError <- stats::sd(pairedDifference) / sqrt(nLibraries)

                       # Differences identical to the last digit leave nothing to divide by, which is degenerate and not decisive
                       if (is.na(standardError) || standardError == 0) {
                         tStatistic <- NA_real_
                         pValue <- NA_real_
                         marginOfError <- NA_real_
                       } else {
                         tStatistic <- meanDifference / standardError
                         pValue <- 2 * stats::pt(-abs(tStatistic), df = nLibraries - 1)
                         marginOfError <- stats::qt(1 - (1 - confLevel) / 2, df = nLibraries - 1) * standardError
                       }

                       data.frame(set.1 = pair[1],
                                  set.2 = pair[2],
                                  n.libraries = nLibraries,
                                  mean.delta.score = meanDifference,
                                  CI.lower = meanDifference - marginOfError,
                                  CI.upper = meanDifference + marginOfError,
                                  t.statistic = tStatistic,
                                  df = nLibraries - 1,
                                  p.value = pValue,
                                  stringsAsFactors = FALSE)
                     }))

    comparisonTable$FDR <- stats::p.adjust(comparisonTable$p.value, method = adjustMethod)
    comparisonTable <- dplyr::arrange(comparisonTable, .data$p.value)

    #---------------------#
    # Assemble the object #
    #---------------------#
    # The colData travels with the scores, so the libraries can be grouped later without the counts object
    sampleMetadata <- data.frame(sample = colnames(counts),
                                 as.data.frame(SummarizedExperiment::colData(counts)),
                                 row.names = NULL,
                                 stringsAsFactors = FALSE)

    scoresObject <-
      methods::new(Class = "RegionSetDE.setScores",
                   scores = setScores,
                   comparisons = comparisonTable,
                   sample.metadata = sampleMetadata,
                   reference = reference,
                   assay = assay,
                   summary = summary,
                   per.basepair = perBasepair,
                   blacklist = counts@blacklist,
                   whitelist = counts@whitelist,
                   genome.assembly = counts@genome.assembly,
                   seqlevels.style = counts@seqlevels.style,
                   filtering.log = counts@filtering.log,
                   parameters = c(counts@parameters,
                                  list(scoreRegionSets = list(reference = reference,
                                                              assay = assay,
                                                              perBasepair = perBasepair,
                                                              summary = summary,
                                                              pseudoCount = pseudoCount,
                                                              minRegions = minRegions,
                                                              confLevel = confLevel,
                                                              adjustMethod = adjustMethod))))

    methods::validObject(scoresObject)

    if (isTRUE(verbose)) {
      message("Done. What separates two sets here is signal, the composition of the sets is not controlled for.")
    }

    return(scoresObject)
  } # END function
